# /// script
# requires-python = ">=3.10"
# dependencies = ["gerbonara", "cairosvg", "numpy", "pillow"]
# ///
"""把 EasyEDA 匯出的 Gerber zip 轉成 PCB 正反兩面的 PBR 貼圖。

用法:
    uv run scripts/make_textures.py <Gerber.zip> <輸出資料夾> [--size 4096] [--fit-outline]

--fit-outline:依板框(BoardOutlineLayer)自動決定貼圖範圍,給板子不在 Gerber 原點、
或不是 Ø40 的板子用。板框中太短的線段(例如零件外形的小框)不計入。

輸出:
    top_basecolor.png / bottom_basecolor.png      板面顏色(防焊、焊墊、絲印、孔)
    top_mr.png / bottom_mr.png                    glTF metallicRoughness (G=粗糙度, B=金屬度)
    texture.json                                  貼圖範圍 center_mm / half_mm(build_glb.py 會讀)

貼圖範圍預設為 Gerber 座標 [-HALF, HALF] mm(板子中心在原點),
影像第 0 列是 +Y(上方),對應 glTF UV 的 v=0。
"""
import argparse
import io
import json
import re
import warnings
import zipfile
from pathlib import Path

import cairosvg
import numpy as np
from gerbonara import ExcellonFile, GerberFile
from gerbonara.utils import MM
from PIL import Image

HALF = 20.32  # 貼圖涵蓋 ±20.32 mm(板半徑 20 mm,留一點邊)

# sRGB 顏色
MASK = np.array([22, 42, 110])          # 防焊(底下是基材)
MASK_ON_CU = np.array([38, 72, 150])    # 防焊(底下有銅,看得出走線)
GOLD = np.array([214, 170, 70])         # 沉金焊墊
FR4 = np.array([150, 140, 95])          # 開窗但沒銅的基材
SILK = np.array([242, 242, 238])        # 絲印
HOLE = np.array([18, 18, 18])           # 鑽孔

LAYERS = {
    'top': dict(cu='Gerber_TopLayer.GTL', mask='Gerber_TopSolderMaskLayer.GTS',
                silk='Gerber_TopSilkscreenLayer.GTO'),
    'bottom': dict(cu='Gerber_BottomLayer.GBL', mask='Gerber_BottomSolderMaskLayer.GBS',
                   silk='Gerber_BottomSilkscreenLayer.GBO'),
}
DRILLS = ['Drill_PTH_Through.DRL', 'Drill_PTH_Through_Via.DRL', 'Drill_NPTH_Through.DRL']


SLOT = re.compile(r'^X(\d+)Y(\d+)G85X(\d+)Y(\d+)$')


def expand_slots(text, step=0.02):
    """gerbonara 不認得 G85 槽孔(X..Y..G85X..Y..),展開成沿槽中心線密排的圓孔。
    座標格式依 EasyEDA 的 METRIC,LZ,000.000(3 位整數 + 3 位小數)。"""
    out = []
    for line in text.splitlines():
        m = SLOT.match(line.strip())
        if not m:
            out.append(line)
            continue
        x1, y1, x2, y2 = (int(v) / 1000 for v in m.groups())
        n = max(2, int(np.hypot(x2 - x1, y2 - y1) / step) + 1)
        out += [f'X{round(x * 1000):06d}Y{round(y * 1000):06d}'
                for x, y in zip(np.linspace(x1, x2, n), np.linspace(y1, y2, n))]
    return '\n'.join(out)


def outline_window(outline, margin=0.32):
    """由板框線段求貼圖中心與半寬(mm)。短於最長線段 10% 的線段不計入。"""
    segs = [(o.x1, o.y1, o.x2, o.y2) for o in outline.objects if hasattr(o, 'x2')]
    length = lambda s: np.hypot(s[2] - s[0], s[3] - s[1])
    lmax = max(map(length, segs))
    pts = np.array([p for s in segs if length(s) >= 0.1 * lmax for p in (s[:2], s[2:])])
    lo, hi = pts.min(0), pts.max(0)
    return (lo + hi) / 2, float((hi - lo).max() / 2 + margin)


def rasterize(cam, size, center=(0.0, 0.0), half=HALF):
    """把一個 Gerber/Excellon 物件畫成 0~1 的灰階遮罩(1 = 有圖形)。"""
    cx, cy = center
    bounds = ((cx - half, cy - half), (cx + half, cy + half))
    svg = str(cam.to_svg(force_bounds=bounds, arg_unit=MM, fg='white', bg='black'))
    png = cairosvg.svg2png(bytestring=svg.encode(), output_width=size, output_height=size,
                           background_color='black')
    img = Image.open(io.BytesIO(png)).convert('L')
    return np.asarray(img, dtype=np.float32) / 255.0


def lerp(base, color, a):
    return base * (1 - a[..., None]) + color[None, None, :] * a[..., None]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('gerber_zip')
    ap.add_argument('out_dir')
    ap.add_argument('--size', type=int, default=4096)
    ap.add_argument('--fit-outline', action='store_true', help='依板框決定貼圖範圍')
    args = ap.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    zf = zipfile.ZipFile(args.gerber_zip)
    read = lambda name: zf.read(name).decode('utf-8', errors='replace')
    names_in_zip = set(zf.namelist())

    center, half = np.zeros(2), HALF
    if args.fit_outline:
        ko = 'Gerber_BoardOutlineLayer.GKO'
        center, half = outline_window(GerberFile.from_string(read(ko), filename=ko))
    (out / 'texture.json').write_text(json.dumps(dict(center_mm=list(map(float, center)), half_mm=half)))
    print(f'貼圖範圍:中心 ({center[0]:.3f}, {center[1]:.3f}) mm,半寬 {half:.3f} mm,'
          f'{2 * half / args.size * 1000:.1f} µm/px')
    ras = lambda cam: rasterize(cam, args.size, center, half)
    blank = np.zeros((args.size, args.size), np.float32)

    holes = blank.copy()
    for d in DRILLS:
        if d in names_in_zip:
            with warnings.catch_warnings():  # EasyEDA 鑽孔檔的格式小問題(重複定義刀具等)
                warnings.simplefilter('ignore', SyntaxWarning)
                exc = ExcellonFile.from_string(expand_slots(read(d)), filename=d)
            holes = np.maximum(holes, ras(exc))

    for side, names in LAYERS.items():
        m = {}
        for k, v in names.items():
            if v in names_in_zip:
                m[k] = ras(GerberFile.from_string(read(v), filename=v))
            else:
                print(f'  {side}: 沒有 {v},當成空白')
                m[k] = blank
        cu, opening, silk = m['cu'], m['mask'], m['silk']

        # 板面顏色:防焊 -> 有銅處較亮 -> 開窗(銅=金、無銅=基材) -> 絲印(不印在開窗上) -> 孔
        rgb = np.broadcast_to(MASK.astype(np.float32), (args.size, args.size, 3)).copy()
        rgb = lerp(rgb, MASK_ON_CU, cu)
        exposed = lerp(np.broadcast_to(FR4.astype(np.float32), rgb.shape).copy(), GOLD, cu)
        rgb = rgb * (1 - opening[..., None]) + exposed * opening[..., None]
        silk_vis = silk * (1 - opening)
        rgb = lerp(rgb, SILK, silk_vis)
        rgb = lerp(rgb, HOLE, holes)
        Image.fromarray(np.clip(rgb, 0, 255).astype(np.uint8)).save(out / f'{side}_basecolor.png', optimize=True)

        # metallicRoughness:R 不用,G=粗糙度,B=金屬度
        metal_mask = opening * cu * (1 - holes)
        rough = 0.35 + (0.7 - 0.35) * opening          # 防焊亮面 0.35,開窗基材 0.7
        rough = rough * (1 - metal_mask) + 0.28 * metal_mask
        rough = rough * (1 - silk_vis) + 0.8 * silk_vis
        mr = np.zeros((args.size, args.size, 3), np.float32)
        mr[..., 1] = rough
        mr[..., 2] = metal_mask
        mr_img = Image.fromarray((mr * 255).astype(np.uint8))
        mr_img.resize((args.size // 2, args.size // 2), Image.LANCZOS).save(out / f'{side}_mr.png', optimize=True)
        print(side, 'done')


if __name__ == '__main__':
    main()
