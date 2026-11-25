# /// script
# requires-python = ">=3.10"
# dependencies = ["gerbonara", "cairosvg", "numpy", "pillow"]
# ///
"""把 EasyEDA 匯出的 Gerber zip 轉成 PCB 正反兩面的 PBR 貼圖。

用法:
    uv run scripts/make_textures.py <Gerber.zip> <輸出資料夾> [--size 4096]

輸出:
    top_basecolor.png / bottom_basecolor.png      板面顏色(防焊、焊墊、絲印、孔)
    top_mr.png / bottom_mr.png                    glTF metallicRoughness (G=粗糙度, B=金屬度)

貼圖範圍固定為 Gerber 座標 [-HALF, HALF] mm(板子中心在原點),
影像第 0 列是 +Y(上方),對應 glTF UV 的 v=0。
"""
import argparse
import io
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


def rasterize(cam, size):
    """把一個 Gerber/Excellon 物件畫成 0~1 的灰階遮罩(1 = 有圖形)。"""
    bounds = ((-HALF, -HALF), (HALF, HALF))
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
    args = ap.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    zf = zipfile.ZipFile(args.gerber_zip)
    read = lambda name: zf.read(name).decode('utf-8', errors='replace')

    holes = np.zeros((args.size, args.size), np.float32)
    for d in DRILLS:
        holes = np.maximum(holes, rasterize(ExcellonFile.from_string(read(d), filename=d), args.size))

    for side, names in LAYERS.items():
        m = {k: rasterize(GerberFile.from_string(read(v), filename=v), args.size)
             for k, v in names.items()}
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
