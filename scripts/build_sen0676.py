# /// script
# requires-python = ">=3.10,<3.13"
# dependencies = ["numpy", "trimesh", "scipy", "networkx", "manifold3d", "pillow", "rtree", "cadquery"]
# ///
"""DFRobot SEN0676(80 GHz 雷達液位感測器)模組:外殼用官方 STEP,電路板依尺寸圖重建。

用法:
    uv run scripts/build_sen0676.py <SEN0676_stpfile.zip 或 .stp> <輸出.glb> [--step 輸出.step]

官方 STEP(Pro/E 匯出,料號 LDRR04C-01)只有外殼一個實體,沒有電路板;電路板依
dimension_V1.1.pdf 第 1 頁重建:
- 35 × 35 mm(四角 R2:外殼方槽內角約 R1.7,方角放不進去),4 個 Ø2.6 安裝孔、孔距 31 × 31
- 左側 5P 排針 TX / RX / DE / GND / VCC(2.54 mm,VCC 離下緣 9 mm),右側 4P、3P 焊孔
- 雷達晶片(AiP,天線在封裝上)與下方黑色晶片、幾顆被動元件、Ø30 絲印圓
  (位置依圖面比例量,約 ±0.3 mm;元件高度是估計值)

組裝:電路板的元件面(雷達晶片)朝鏡頭,放在外殼內 4 根凸柱(頂面有 Ø2.8 × 0.78 沉孔、Ø2 × 5.5 導孔)上,
4 顆 M2.5 × 6 盤頭自攻螺絲從背面(外殼開口端)鎖入。5P 排針裝在背面,朝外殼開口。
板子轉角:試 0 / 90 / 180 / 270°,用 OpenCASCADE 布林運算算每個角度下零件與外殼的重疊體積,取沒有干涉的角度
(凸柱那一層有一側的支撐邊,排針腳若壓在上面就會干涉)。

座標(mm,Z 朝上):鏡頭朝下,鏡頭頂點 z=0,軸心在原點。
"""
import argparse
import sys
import tempfile
import zipfile
from pathlib import Path

import cadquery as cq
import numpy as np
import trimesh

sys.path.insert(0, str(Path(__file__).parent))
from build_glb import GLB, SMOOTH_ANGLE, gltf_arrays, srgb_to_linear  # noqa: E402
from build_gimbal import gltf_matrix, make_extras  # noqa: E402
from hardware import pan_screw, shade  # noqa: E402

# ---------------- 電路板(板子座標:元件面 z=0、背面 z=-PCB_T;x 右、y 上,同尺寸圖第 1 頁) ----------------
PCB_W, PCB_T = 35.0, 1.6
PCB_CORNER_R = 2.0               # 外殼 35.28 方槽的內角約 R1.7,板子四角要倒圓角才放得進去
HOLE_PITCH, HOLE_D = 31.0, 2.6
PX_MM = 14.0                         # 尺寸圖第 1 頁:35 mm = 490 px
PX_C = (550.0, 470.0)                # 板子中心的 px 座標


def px(x, y):
    """尺寸圖 px -> 板子座標 mm。"""
    return (x - PX_C[0]) / PX_MM, (PX_C[1] - y) / PX_MM


def px_box(x0, x1, y0, y1):
    (a, b), (c, d) = px(x0, y1), px(x1, y0)
    return (a + c) / 2, (b + d) / 2, c - a, d - b


PITCH = 2.54
HDR5_X, HDR5_Y0 = -16.0, -17.5 + 9.0          # 5P:VCC 離下緣 9 mm,往上 GND / DE / RX / TX
HDR5 = [(HDR5_X, HDR5_Y0 + k * PITCH) for k in range(5)]          # VCC, GND, DE, RX, TX(TX 為方形焊墊)
HDR4 = [(16.0, 17.5 - 9.0 + 2.54 - k * PITCH) for k in range(4)]  # 右側 4P(第 2 腳離上緣 9 mm)
HDR3 = [(17.5 - 4.0, 17.5 - 9.0 - k * PITCH) for k in range(3)]   # 3P(離右緣 4 mm)
CHIPS = {   # 名稱: (中心 x, 中心 y, 寬, 高, 厚, 材質)
    'radar_chip': (*px_box(507, 597, 310, 395), 1.0, 'chip'),
    'ic': (*px_box(483, 615, 397, 508), 0.9, 'chip_dark'),
}
PASSIVES = [px_box(*b) for b in (                                    # 0805 / 0603 被動元件(依圖面方塊)
    (441, 471, 305, 319), (441, 471, 320, 334), (472, 501, 290, 318), (504, 531, 272, 289), (504, 531, 290, 303),
    (708, 737, 497, 512), (708, 737, 540, 554), (708, 737, 560, 575))]
PASSIVE_H = 0.6
SILK_R, SILK_W = 15.0, 0.25          # Ø30 絲印圓
HDR_BASE_H, PIN_ABOVE, PIN_TAIL = 2.5, 6.0, 1.5   # 5P 排針:塑膠座高、座上針長、元件面露出的焊腳

SCREW = (2.5, 6.0)                   # M2.5 × 6 盤頭自攻(凸柱導孔 Ø2)

PARTS = {
    'housing': ('外殼', '外購件', 1, 'SEN0676 官方 STEP(LDRR04C-01),G1½ 螺紋、六角法蘭 62'),
    'pcb': ('感測器電路板', '電路板', 1, '35 × 35 mm,依尺寸圖重建;雷達晶片朝鏡頭'),
    'header': ('5P 排針', '電子', 1, 'TX / RX / DE / GND / VCC,2.54 mm,裝在背面'),
    'pcb_screws': ('電路板螺絲', '五金', 1, 'M2.5 × 6 盤頭自攻,鎖進外殼凸柱 Ø2 導孔'),
}
EXPLODE = {'housing': (0, 0, 0), 'pcb': (0, 0, 45), 'header': (0, 0, 55), 'pcb_screws': (0, 0, 75)}
MATS = {   # sRGB, 金屬度, 粗糙度
    'housing': ((0.06, 0.06, 0.07), 0.0, 0.55),       # 實品為黑色(STEP 原色是淺灰)
    'pcb': ((0.05, 0.32, 0.16), 0.0, 0.45),
    'silk': ((0.95, 0.95, 0.93), 0.0, 0.8),
    'gold': ((0.83, 0.66, 0.28), 1.0, 0.3),
    'chip': ((0.20, 0.20, 0.22), 0.2, 0.4),
    'chip_dark': ((0.06, 0.06, 0.07), 0.0, 0.5),
    'smd': ((0.30, 0.25, 0.20), 0.0, 0.6),
    'plastic': ((0.05, 0.05, 0.05), 0.0, 0.6),
    'tin': ((0.80, 0.80, 0.82), 1.0, 0.3),
    'screw': ((0.80, 0.80, 0.82), 1.0, 0.25),
}


def cbox(cx, cy, w, h, z0, z1):
    return cq.Workplane('XY').box(w, h, z1 - z0, centered=(True, True, False)).translate((cx, cy, z0)).val()


def ccyl(cx, cy, r, z0, z1):
    return cq.Workplane('XY').circle(r).extrude(z1 - z0).translate((cx, cy, z0)).val()


def fuse(shapes):
    out = shapes[0]
    for s in shapes[1:]:
        out = out.fuse(s)
    return out.clean()


def make_board():
    """回傳 {(零件 id, 材質): cadquery Shape}(板子座標)。"""
    tht = HDR5 + HDR4 + HDR3
    holes = [(s * HOLE_PITCH / 2, t * HOLE_PITCH / 2) for s in (-1, 1) for t in (-1, 1)]
    board = (cq.Workplane('XY').rect(PCB_W, PCB_W).extrude(PCB_T).edges('|Z').fillet(PCB_CORNER_R)
             .translate((0, 0, -PCB_T)).val())
    board = board.cut(fuse([ccyl(x, y, HOLE_D / 2, -PCB_T - 1, 1) for x, y in holes] +
                           [ccyl(x, y, 0.5, -PCB_T - 1, 1) for x, y in tht]))
    pads = fuse([ccyl(x, y, 0.85, z0, z1).cut(ccyl(x, y, 0.5, z0 - 1, z1 + 1))
                 for x, y in tht for z0, z1 in ((0, 0.035), (-PCB_T - 0.035, -PCB_T))])
    silk = ccyl(0, 0, SILK_R + SILK_W / 2, 0, 0.02).cut(ccyl(0, 0, SILK_R - SILK_W / 2, -1, 1))
    parts = {('pcb', 'pcb'): board, ('pcb', 'gold'): pads, ('pcb', 'silk'): silk}
    for name, (cx, cy, w, h, t, mat) in CHIPS.items():
        parts[('pcb', mat)] = parts[('pcb', mat)].fuse(cbox(cx, cy, w, h, 0, t)) if ('pcb', mat) in parts \
            else cbox(cx, cy, w, h, 0, t)
    bodies, caps = [], []
    for cx, cy, w, h in PASSIVES:
        horiz = w >= h
        L, W = (w, h) if horiz else (h, w)
        cap = 0.3
        bodies.append(cbox(cx, cy, w - (2 * cap if horiz else 0), h - (0 if horiz else 2 * cap), 0, PASSIVE_H))
        for s in (-1, 1):
            if horiz:
                caps.append(cbox(cx + s * (w / 2 - cap / 2), cy, cap, h, 0, PASSIVE_H))
            else:
                caps.append(cbox(cx, cy + s * (h / 2 - cap / 2), w, cap, 0, PASSIVE_H))
    parts[('pcb', 'smd')] = fuse(bodies)
    parts[('pcb', 'tin')] = fuse(caps)
    # 5P 排針(背面):塑膠座 + 方針,焊腳從元件面露出 PIN_TAIL
    ys = [y for _, y in HDR5]
    base = cbox(HDR5_X, (min(ys) + max(ys)) / 2, 2.5, PITCH * 5, -PCB_T - HDR_BASE_H, -PCB_T)
    pins = fuse([cbox(x, y, 0.64, 0.64, -PCB_T - HDR_BASE_H - PIN_ABOVE, PIN_TAIL) for x, y in HDR5])
    parts[('header', 'plastic')] = base
    parts[('header', 'tin')] = pins
    return parts, holes


def read_step(path):
    path = Path(path)
    if path.suffix.lower() == '.zip':
        with zipfile.ZipFile(path) as z:
            name = next(n for n in z.namelist() if n.lower().endswith(('.stp', '.step')))
            tmp = Path(tempfile.mkdtemp()) / Path(name).name
            tmp.write_bytes(z.read(name))
            path = tmp
    return cq.importers.importStep(str(path)).val()


def to_trimesh(shape, tol=0.05, ang=0.3):
    v, t = shape.tessellate(tol, ang)
    m = trimesh.Trimesh(np.array([p.toTuple() for p in v]), np.array(t), process=True)
    trimesh.repair.fix_normals(m, multibody=True)
    return m


def place(shape, M):
    """4x4 齊次矩陣 -> cadquery Location 套用。"""
    from OCP.gp import gp_Trsf
    T = gp_Trsf()
    T.SetValues(*M[0, :4], *M[1, :4], *M[2, :4])
    return cq.Shape.cast(shape.wrapped.Moved(cq.Location(T).wrapped))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('step')
    ap.add_argument('out')
    ap.add_argument('--step-out', help='整個模組另存 STEP')
    args = ap.parse_args()

    housing = read_step(args.step)
    hm = to_trimesh(housing, 0.05, 0.3)
    z_tip = float(hm.bounds[0][2])
    lift = np.eye(4)
    lift[2, 3] = -z_tip                                   # 鏡頭頂點 -> z=0
    housing = place(housing, lift)
    hm.apply_transform(lift)
    print(f'外殼:{np.ptp(hm.bounds[:, 0]):.2f} × {np.ptp(hm.bounds[:, 1]):.2f} × {np.ptp(hm.bounds[:, 2]):.2f} mm')

    parts, holes = make_board()
    # 凸柱頂面:在孔旁 2.2 mm(實心處)往下打射線
    tops = []
    for x, y in holes:
        h, _, _ = hm.ray.intersects_location([[x + 2.2, y, 200]], [[0, 0, -1]])
        tops.append(float(h[:, 2].max()))
    assert np.ptp(tops) < 0.05, f'凸柱頂面不一致 {tops}'
    z_face = tops[0]
    print(f'凸柱頂面 z={z_face:.2f}(元件面貼這裡),板子背面 z={z_face + PCB_T:.2f}')

    flip = np.diag([1.0, -1.0, -1.0])                    # 元件面朝下(朝鏡頭)

    def M_for(yaw):
        a = np.radians(yaw)
        Rz = np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])
        M = np.eye(4)
        M[:3, :3] = Rz @ flip
        M[2, 3] = z_face
        return M

    # 試 4 個角度,算零件與外殼重疊體積
    scores = {}
    for yaw in (0, 90, 180, 270):
        M = M_for(yaw)
        vol = 0.0
        for shape in parts.values():
            vol += housing.intersect(place(shape, M)).Volume()
        scores[yaw] = vol
        print(f'  板子轉 {yaw:3d}°:零件與外殼重疊 {vol:.3f} mm³')
    yaw = min(scores, key=lambda k: (round(scores[k], 3), k))
    assert scores[yaw] < 0.01, '每個角度都會干涉'
    M = M_for(yaw)
    print(f'採用 {yaw}°')

    # 螺絲:座面在板子背面,頭朝 +Z(外殼開口),螺桿往下鎖進凸柱
    screw_l = SCREW[1]
    engage = screw_l - PCB_T
    hole_bottom = []
    for x, y in holes:
        wx, wy = (M[:3, :3] @ [x, y, 0])[:2]
        h, _, _ = hm.ray.intersects_location([[wx, wy, z_face - 0.01]], [[0, 0, -1]])
        hole_bottom.append(float(h[:, 2].max()))
    depth = z_face - max(hole_bottom)
    print(f'螺絲 M{SCREW[0]} × {screw_l:.0f}:鎖入凸柱 {engage:.1f} mm(導孔深 {depth:.2f} mm)')
    assert engage < depth - 0.3

    # ---- GLB ----
    g = GLB()
    nodes = g.j['nodes']
    nodes[0]['name'] = 'SEN0676'
    mats = {}

    def material(key):
        if key not in mats:
            srgb, metal, rough = MATS[key]
            mats[key] = g.material(dict(name=key, pbrMetallicRoughness=dict(
                baseColorFactor=[*srgb_to_linear(srgb).tolist(), 1.0], metallicFactor=metal, roughnessFactor=rough)))
        return mats[key]

    def add(name, tm_mm, key, part, matrix=None, mesh=None):
        if mesh is None:
            tm = tm_mm.copy()
            tm.apply_scale(0.001)
            pos, nrm, idx = gltf_arrays(shade(tm, key, SMOOTH_ANGLE))
            g.mesh(name, pos, nrm, idx, material(key))
        else:
            nodes.append(dict(name=name, mesh=mesh))
            nodes[0]['children'].append(len(nodes) - 1)
        nodes[-1]['extras'] = make_extras(PARTS, EXPLODE, part)
        if matrix is not None:
            nodes[-1]['matrix'] = matrix
        return nodes[-1]['mesh']

    add('housing', hm, 'housing', 'housing')
    for (part, key), shape in parts.items():
        add(f'{part}_{key}', to_trimesh(place(shape, M), 0.02, 0.2), key, part)
    mg = pan_screw(*SCREW).to_mesh()
    screw = trimesh.Trimesh(np.asarray(mg.vert_properties)[:, :3], np.asarray(mg.tri_verts), process=True)
    mid = None
    for i, (x, y) in enumerate(holes):
        wx, wy = (M[:3, :3] @ [x, y, 0])[:2]
        mid = add(f'pcb_screw_{i}', screw, 'screw', 'pcb_screws', gltf_matrix(t_mm=[wx, wy, z_face + PCB_T]), mid)
    g.save(args.out)
    print('輸出', args.out, f'{Path(args.out).stat().st_size / 1e6:.1f} MB')

    if args.step_out:
        asm = cq.Assembly(name='SEN0676')
        asm.add(housing, name='housing', color=cq.Color(*MATS['housing'][0]))
        for (part, key), shape in parts.items():
            asm.add(place(shape, M), name=f'{part}_{key}', color=cq.Color(*MATS[key][0]))
        asm.export(args.step_out)
        print('輸出', args.step_out)


if __name__ == '__main__':
    main()
