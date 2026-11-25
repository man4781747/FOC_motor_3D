# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "trimesh", "scipy", "networkx", "manifold3d", "pillow", "rtree", "lxml", "shapely"]
# ///
"""水情盒:把底板(3MF)放進 RS PRO 1862264 防水盒底部,用 4 顆螺絲鎖在盒內的凸柱上;
可選擇再把 RS485 擴展板用 12 mm 塑膠六角柱裝在底板上。

用法:
    uv run scripts/build_water_box.py <防水盒 3MF 或 zip> <底板 3MF> <輸出.glb> [--rs485 output/rs485/rs485_v1.glb]

盒子底部內側有 4 根凸柱(Ø7.9,中心 Ø2.4 自攻導孔,頂面比內底高 5.9 mm),位置 (±62.75, ±61.4);
底板的 4 個 Ø4.0 孔剛好對準(誤差 < 0.1 mm),四角與兩側的缺口避開盒蓋螺絲柱。
螺絲用 ST3.0 × 6 盤頭自攻螺絲:穿過 2 mm 底板,鎖進凸柱約 4 mm(導孔深 5.9 mm,不會頂到底)。

RS485 擴展板(--rs485):底板上有兩組 95 × 95 mm 的 Ø3.0 孔,用離底板中心較近的那組。
板子四角 Ø3.0 孔間距 95.1 × 95.1 mm,裝在 M3 × 12 尼龍母母六角柱上;下方用 M3 × 6 盤頭螺絲從底板背面
鎖入(要在底板放進盒子前先鎖好),上方用 M3 × 30 + 6 金屬公母六角柱取代螺絲:公頭穿過板子鎖進尼龍柱,
母牙端朝上。

水情盒主板(--main):裝在金屬六角柱頂端(最上層),用四角 Ø3.0 孔(同樣 95.1 × 95.1 mm),
M3 × 6 盤頭螺絲從主板正面鎖進金屬柱母牙;轉角與 RS485 板相同。板子轉 270°,4 組端子朝 -X(盒子中央方向)。

SEN0676 雷達(--sen0676):在盒身 -X 短邊的牆上開 Ø48.5 孔(G1½ 螺紋外徑 47.8 要穿過),孔中心在 (-X, +Y)
角落盒蓋螺絲柱的邊緣往 -Y 70 mm、盒身上緣往下 27 mm。雷達從盒內裝,鏡頭(尖頭)朝外:
- 螺紋根部的凹槽(Ø44.6 × 2.5)先套 O 型環(GB/T 3452.1 43.7 × 2.65),法蘭把 O 型環壓在牆內面(壓縮 25%)
- 螺紋穿出牆外,用 KSS AGJL-112 尼龍螺帽(PF 1½",對邊 58.5、高 9.8)從外面鎖緊
- 六角法蘭與螺帽的一組對邊平行於盒底
目前只列出干涉,不擋。

座標(mm,Z 朝上):盒子長邊沿 X、外底面 z=0、中心在原點。
"""
import argparse
import io
import itertools
import json
import sys
import zipfile
from pathlib import Path

import numpy as np
import trimesh
from manifold3d import Manifold, Mesh
from shapely.geometry import Polygon
from shapely.ops import unary_union

sys.path.insert(0, str(Path(__file__).parent))
from build_assembly import load_glb  # noqa: E402
from build_glb import GLB, SMOOTH_ANGLE, gltf_arrays, srgb_to_linear  # noqa: E402
from build_gimbal import gltf_matrix, make_extras  # noqa: E402
from hardware import pan_screw, shade, standoff  # noqa: E402

# 防水盒 3MF 的網格名稱(同 build_car.py);其餘 6 個是盒蓋螺牙嵌件
BOX_MESHES = {'1862264_1': ('box_body', 'box_gray'), '1862264_3': ('box_lid', 'box_gray'),
              '1862264_10': ('box_gasket', 'box_gasket')}
# 盒子 3MF 座標(Y 朝上、長邊沿 Z)-> 本檔座標:(x, y, z)_盒 -> (z, x, y)
BOX_ROT = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], float)
# 底板 3MF 座標(長邊沿 Y)-> 本檔座標:長邊沿 X;(x, y, z)_板 -> (y, -x, z)
PLATE_ROT = np.array([[0, 1, 0], [-1, 0, 0], [0, 0, 1]], float)

PLATE_HOLE_D = 4.0           # 底板上的固定孔(用來找孔位)
SCREW_D, SCREW_L = 3.0, 6.0  # ST3.0 × 6 盤頭自攻
SCREW_HEAD_D, SCREW_HEAD_H = 5.6, 2.2

PCB_HOLE_D = 3.0             # RS485 板四角孔、底板上裝電路板的孔
PCB_HOLE_PITCH = 95.1        # RS485 板孔距(底板孔組要和它差 0.3 mm 內)
SPACER_L, SPACER_FLATS = 12.0, 5.5   # M3 × 12 尼龍母母六角柱,對邊 5.5
SPACER_THREAD_DEPTH = 5.5    # 尼龍六角柱通常整支有牙,模型只做到 5.5 mm
PCB_SCREW_L = 6.0            # M3 × 6 盤頭(底板背面)
TOP_POST_L, TOP_POST_STUD = 30.0, 6.0   # 板子正面:M3 × 30 + 6 金屬公母六角柱(公頭取代螺絲)
TOP_POST_FLATS = 5.5
TOP_POST_THREAD = 4.0        # 金屬柱母牙深度(hardware.standoff 預設)
PCB_YAW = 270.0              # 板子轉角:端子(板子 -Y 邊)朝 -X(盒子中央方向)
MAIN_YAW = PCB_YAW           # 主板轉角,與 RS485 板相同
SEN_HOLE_D = 48.5            # 雷達孔徑(G1½ 螺紋外徑 47.8 穿過)
SEN_ORING = (43.7, 2.65)     # 螺紋根部凹槽的 O 型環(GB/T 3452.1 內徑 × 線徑)
SEN_ORING_SQUEEZE = 0.25     # O 型環被法蘭壓在牆內面的壓縮率
SEN_NUT = (58.5, 9.8)        # KSS AGJL-112 尼龍螺帽:對邊 W、高 H(PF 1½")
SEN_HOLE_FROM_POST = 70.0    # 孔中心:(-X, +Y) 角落盒蓋螺絲柱邊緣往 -Y
SEN_HOLE_FROM_TOP = 27.0     # 孔中心:盒身上緣往下

PARTS = {   # id -> (名稱, 分類, 每個節點件數, 說明)
    'plate': ('底板', '列印件', 1, '底板_3.3mf · 155 × 205 × 2 mm'),
    'plate_screws': ('底板螺絲', '五金', 1, 'ST3.0 × 6 盤頭自攻,鎖進盒底凸柱的 Ø2.4 導孔'),
    'box_body': ('防水盒 盒身', '外購件', 1, 'RS PRO 1862264,225 × 175 × 81 mm;盒底內有 4 根底板凸柱'),
    'box_lid': ('防水盒 盒蓋', '外購件', 1, 'RS PRO 1862264 附'),
    'box_gasket': ('防水盒 密封條', '外購件', 1, 'RS PRO 1862264 附'),
    'box_inserts': ('防水盒 螺牙嵌件', '外購件', 6, 'RS PRO 1862264 附,盒蓋螺絲鎖入處'),
    'rs485': ('RS485 擴展板', '電路板', 1, 'RS485_V1,100 × 100 mm;轉 270°,端子朝盒子中央(-X)'),
    'rs485_spacers': ('RS485 六角柱', '五金', 1, 'M3 × 12 尼龍母母六角柱,對邊 5.5'),
    'main_pcb': ('水情盒主板', '電路板', 1, 'Main_PCB V2.4,100 × 100 mm;裝在金屬六角柱頂端'),
    'main_screws': ('主板螺絲', '五金', 1, 'M3 × 6 盤頭,從主板正面鎖進金屬六角柱'),
    'sen_oring': ('SEN0676 O 型環', '外購件', 1, 'GB/T 3452.1 43.7 × 2.65(NBR),套在螺紋根部凹槽,法蘭壓在牆內面'),
    'sen_nut': ('SEN0676 鎖緊螺帽', '外購件', 1, 'KSS AGJL-112 尼龍螺帽 PF 1½"(對邊 58.5、高 9.8,黑),盒外鎖緊'),
    'rs485_posts_top': ('RS485 上方金屬六角柱', '五金', 1, 'M3 × 30 + 6 公母六角柱(對邊 5.5),公頭穿過板子鎖進尼龍柱'),
    'rs485_screws_bottom': ('RS485 底板螺絲', '五金', 1, 'M3 × 6 盤頭,從底板背面鎖進六角柱(底板裝進盒子前先鎖)'),
}
EXPLODE = {   # mm,Z 朝上:照組裝順序往上拆
    'box_body': (0, 0, 0), 'plate': (0, 0, 45), 'plate_screws': (0, 0, 75),
    'rs485_screws_bottom': (0, 0, 25), 'rs485_spacers': (0, 0, 70), 'rs485': (0, 0, 100),
    'rs485_posts_top': (0, 0, 125),
    'main_pcb': (0, 0, 160), 'main_screws': (0, 0, 185), 'sen0676': (70, 0, 0), 'sen_oring': (40, 0, 0), 'sen_nut': (-45, 0, 0),
    'box_inserts': (0, 0, 205), 'box_gasket': (0, 0, 220), 'box_lid': (0, 0, 250),
}
MATS = {   # sRGB, 金屬度, 粗糙度
    'box_gray': ((0.74, 0.75, 0.76), 0.0, 0.55),
    'box_gasket': ((0.12, 0.12, 0.13), 0.0, 0.85),
    'brass': ((0.80, 0.62, 0.32), 1.0, 0.3),
    'screw': ((0.80, 0.80, 0.82), 1.0, 0.25),
    'print_white': ((0.93, 0.93, 0.91), 0.0, 0.6),
    'nylon_black': ((0.05, 0.05, 0.05), 0.0, 0.55),
    'rubber': ((0.04, 0.04, 0.04), 0.0, 0.9),
    'nylon': ((0.96, 0.95, 0.90), 0.0, 0.45),
}


def read_3mf(path):
    path = Path(path)
    data = path.read_bytes()
    if path.suffix == '.zip':
        with zipfile.ZipFile(path) as z:
            data = z.read(next(n for n in z.namelist() if n.endswith('.3mf')))
    return trimesh.load(io.BytesIO(data), file_type='3mf', force='scene')


def affine(R, t=(0, 0, 0)):
    M = np.eye(4)
    M[:3, :3], M[:3, 3] = R, t
    return M


def load_box(path):
    """{(零件 id, 材質): [trimesh]},外底面 z=0、中心在原點。"""
    scene = read_3mf(path)
    out = {}
    for n in scene.graph.nodes_geometry:
        T, g = scene.graph[n]
        m = scene.geometry[g].copy()
        m.apply_transform(affine(BOX_ROT) @ T)
        out.setdefault(BOX_MESHES.get(n, ('box_inserts', 'brass')), []).append(m)
    allm = [m for ms in out.values() for m in ms]
    lo = np.min([m.bounds[0] for m in allm], 0)
    hi = np.max([m.bounds[1] for m in allm], 0)
    shift = [-(lo[0] + hi[0]) / 2, -(lo[1] + hi[1]) / 2, -lo[2]]
    for m in allm:
        m.apply_translation(shift)
    return out


def loops(path3d):
    """平面截面 -> [(中心 xy, 直徑)](只取小的封閉圈,也就是孔或凸柱)。"""
    res = []
    for e in path3d.discrete:
        d = np.ptp(e[:, :2], 0)
        if d.max() < 12:
            res.append(((e[:, :2].min(0) + e[:, :2].max(0)) / 2, float(d.mean())))
    return res


def find_bosses(body, plate_holes):
    """在每個底板孔的正下方找凸柱:回傳 [(x, y, 頂面 z, 導孔直徑)]。"""
    bosses = []
    for (x, y), _ in plate_holes:
        # 往下打一條射線偏離導孔中心 2.2 mm(凸柱實心處),第一個交點是凸柱頂面
        hits, _, _ = body.ray.intersects_location([[x + 2.2, y, 200]], [[0, 0, -1]])
        top = float(hits[:, 2].max())
        sec = body.section([0, 0, 1], [0, 0, top - 1.0])
        near = [(c, d) for c, d in loops(sec) if np.hypot(*(c - [x, y])) < 1.0]
        assert near, f'孔 ({x:.1f}, {y:.1f}) 下方找不到凸柱'
        c, pilot = min(near, key=lambda cd: cd[1])
        bosses.append((float(c[0]), float(c[1]), top, pilot))
    return bosses


def check_clearance(body, plate, z0, z1):
    """底板厚度範圍內,盒身截面與底板截面不能重疊(凸柱頂面剛好貼底板,所以從 z0 + 0.05 開始查)。"""
    worst = 0.0
    for z in np.linspace(z0 + 0.05, z1 - 0.05, 5):
        ps = plate.section([0, 0, 1], [0, 0, z])
        bs = body.section([0, 0, 1], [0, 0, z])
        if ps is None or bs is None:
            continue
        p_poly = unary_union([Polygon(e[:, :2]).buffer(0) for e in ps.discrete if len(e) > 3])
        # 盒身截面:外牆(兩圈之間)與各凸柱。用每一圈的多邊形做 symmetric difference 得到實體區域
        solid = Polygon()
        for e in bs.discrete:
            if len(e) > 3:
                solid = solid.symmetric_difference(Polygon(e[:, :2]).buffer(0))
        worst = max(worst, p_poly.intersection(solid).area)
    return worst


def square_sets(holes, pitch, tol=0.3):
    """從孔清單找出所有「邊長 ≈ pitch 的正方形」四孔組,回傳 [(中心, 4 個孔中心)]。"""
    out = []
    for combo in itertools.combinations([c for c, _ in holes], 4):
        pts = np.array(combo)
        c = pts.mean(0)
        ok = all(np.abs(np.abs(p - c) - pitch / 2).max() < tol for p in pts)
        if ok:
            out.append((c, pts))
    return out


def pcb_holes(pcb_glb):
    """PCB GLB -> 板子座標(mm,Z 朝上,背面 z=0)的四角孔中心(95.1 × 95.1 那組)與板厚。"""
    sc = trimesh.load(pcb_glb, force='scene')
    edge = sc.geometry[next(n for n in sc.geometry if n.startswith('PCB_edge'))].copy()
    edge.apply_transform(np.array([[1, 0, 0, 0], [0, 0, -1, 0], [0, 1, 0, 0], [0, 0, 0, 1]], float))  # Y-up -> Z-up
    edge.apply_scale(1000)
    thick = float(np.ptp(edge.bounds[:, 2]))
    hs = [(c, d) for c, d in loops(edge.section([0, 0, 1], [0, 0, thick / 2])) if abs(d - PCB_HOLE_D) < 0.2]
    sets = square_sets(hs, PCB_HOLE_PITCH, tol=0.1)
    assert len(sets) == 1, f'{pcb_glb}:應有一組 95.1 × 95.1 的 Ø{PCB_HOLE_D} 孔,找到 {len(sets)} 組'
    return sets[0][1], thick


def merge_glb(g, path):
    """把另一個 GLB(build_glb.py 的輸出)併進 g,回傳它根節點底下的節點編號(未掛到場景上)。"""
    other = load_glb(path)
    j, o = g.j, other.j
    while len(g.bin) % 4:
        g.bin.append(0)
    off = {k: len(j[k]) for k in ('bufferViews', 'accessors', 'images', 'samplers', 'textures', 'materials',
                                  'meshes', 'nodes')}
    for bv in o['bufferViews']:
        j['bufferViews'].append(dict(bv, byteOffset=bv.get('byteOffset', 0) + len(g.bin)))
    g.bin += other.bin
    for a in o['accessors']:
        j['accessors'].append(dict(a, bufferView=a['bufferView'] + off['bufferViews']))
    for im in o.get('images', []):
        j['images'].append(dict(im, bufferView=im['bufferView'] + off['bufferViews']))
    j.setdefault('samplers', []).extend(o.get('samplers', []))
    for t in o.get('textures', []):
        j['textures'].append(dict(t, source=t['source'] + off['images'], sampler=t.get('sampler', 0) + off['samplers']))
    for m in o['materials']:
        m = json.loads(json.dumps(m))
        pbr = m.get('pbrMetallicRoughness', {})
        for k in ('baseColorTexture', 'metallicRoughnessTexture'):
            if k in pbr:
                pbr[k]['index'] += off['textures']
        j['materials'].append(m)
    for me in o['meshes']:
        prims = [dict(p, attributes={k: v + off['accessors'] for k, v in p['attributes'].items()},
                      indices=p['indices'] + off['accessors'], material=p['material'] + off['materials'])
                 for p in me['primitives']]
        j['meshes'].append(dict(me, primitives=prims))
    for n in o['nodes']:
        n = dict(n)
        if 'mesh' in n:
            n['mesh'] += off['meshes']
        if 'children' in n:
            n['children'] = [c + off['nodes'] for c in n['children']]
        j['nodes'].append(n)
    return [c + off['nodes'] for c in o['nodes'][0]['children']]


def rz(deg):
    a = np.radians(deg)
    return np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])


def sensor_profile(sen_glb):
    """SEN0676 GLB(鏡頭頂點 z=0、軸 = Z)-> 外殼各高度的最大半徑 [(z, r)]。"""
    sc = trimesh.load(sen_glb, force='scene')
    T, gname = sc.graph['housing']
    h = sc.geometry[gname].copy()
    h.apply_transform(np.array([[1, 0, 0, 0], [0, 0, -1, 0], [0, 1, 0, 0], [0, 0, 0, 1]], float) @ T)
    h.apply_scale(1000)
    prof = []
    for z in np.arange(0.2, h.bounds[1][2], 0.1):
        sec = h.section([0, 0, 1], [0, 0, z])
        if sec is not None:
            prof.append((float(z), max(float(np.hypot(*e[:, :2].T).max()) for e in sec.discrete)))
    return prof


def sen_nut_mesh(W, H, bore, collar_h=2.8):
    """AGJL 尼龍螺帽(雷達座標:軸 = Z,z ∈ [0, H],z=H 端是貼牆的圓形墊圈面);六角平邊垂直於 x。"""
    from manifold3d import CrossSection
    R = W / np.sqrt(3)                                     # 外接圓半徑
    hexa = CrossSection([[(R * np.cos(np.radians(a)), R * np.sin(np.radians(a))) for a in range(30, 390, 60)]])
    body = hexa.extrude(H - collar_h)
    body = body ^ Manifold.cylinder(H - collar_h, R - 0.4, R + 1.5, 96)          # 六角外緣倒角
    collar = Manifold.cylinder(collar_h, W / 2 + 1.5, W / 2 + 1.5, 96).translate([0, 0, H - collar_h])
    hole = Manifold.cylinder(H + 2, bore / 2, bore / 2, 96).translate([0, 0, -1])
    out = ((body + collar) - hole).to_mesh()
    return trimesh.Trimesh(np.asarray(out.vert_properties)[:, :3], np.asarray(out.tri_verts), process=True)


def cut_hole_x(mesh, y, z, r, x0, x1):
    """在網格上沿 X 打圓孔(manifold3d 布林)。"""
    mm = mesh.copy()
    mm.merge_vertices(merge_tex=True, merge_norm=True)
    m = Manifold(Mesh(vert_properties=np.asarray(mm.vertices, np.float32), tri_verts=np.asarray(mm.faces, np.uint32)))
    cyl = Manifold.cylinder(x1 - x0, r, r, 96).rotate([0, 90, 0]).translate([x0, y, z])
    out = (m - cyl).to_mesh()
    return trimesh.Trimesh(np.asarray(out.vert_properties)[:, :3], np.asarray(out.tri_verts), process=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('box')
    ap.add_argument('plate')
    ap.add_argument('out')
    ap.add_argument('--rs485', help='RS485 擴展板 GLB(build_glb.py 的輸出),裝在底板上')
    ap.add_argument('--main', help='水情盒主板 GLB,裝在 RS485 板上方金屬六角柱頂端(需要 --rs485)')
    ap.add_argument('--sen0676', help='SEN0676 雷達 GLB(build_sen0676.py 的輸出),裝在 -Y 牆的 Ø40 孔')
    args = ap.parse_args()

    box = load_box(args.box)
    body = trimesh.util.concatenate(box[('box_body', 'box_gray')])

    # ---- 底板:轉成長邊沿 X,用 4 個 Ø4 孔對準盒底凸柱 ----
    plate = trimesh.load(args.plate, file_type='3mf', force='mesh')
    plate.apply_transform(affine(PLATE_ROT))
    plate.apply_translation(-(plate.bounds[0] + plate.bounds[1]) / 2 * [1, 1, 0] - [0, 0, plate.bounds[0][2]])
    thick = float(np.ptp(plate.bounds[:, 2]))
    holes = [(c, d) for c, d in loops(plate.section([0, 0, 1], [0, 0, thick / 2]))
             if abs(d - PLATE_HOLE_D) < 0.2]
    assert len(holes) == 4, f'底板應有 4 個 Ø{PLATE_HOLE_D} 固定孔,找到 {len(holes)} 個'
    bosses = find_bosses(body, holes)
    tops = [b[2] for b in bosses]
    assert np.ptp(tops) < 0.05, f'凸柱頂面高度不一致 {tops}'
    # 平移:孔中心平均對到凸柱中心平均
    dxy = np.mean([b[:2] for b in bosses], 0) - np.mean([c for c, _ in holes], 0)
    plate.apply_translation([dxy[0], dxy[1], tops[0]])
    err = max(np.hypot(*(np.array(b[:2]) - (c + dxy))) for b, (c, _) in zip(bosses, holes))
    z_top = tops[0] + thick
    floor_hits, _, _ = body.ray.intersects_location([[bosses[0][0], bosses[0][1], 200]], [[0, 0, -1]])
    pilot_bottom = float(floor_hits[:, 2].max())          # 導孔底(= 盒內底面)
    print(f'底板 {np.ptp(plate.bounds[:, 0]):.1f} × {np.ptp(plate.bounds[:, 1]):.1f} × {thick:.1f} mm,'
          f'放在凸柱頂 z={tops[0]:.2f},平移 ({dxy[0]:+.2f}, {dxy[1]:+.2f}),孔位最大誤差 {err:.3f} mm')
    for x, y, top, pilot in bosses:
        print(f'  凸柱 ({x:+7.2f}, {y:+7.2f}) 頂 z={top:.2f},導孔 Ø{pilot:.2f} 深 {top - pilot_bottom:.1f} mm')

    engage = SCREW_L - thick
    depth = tops[0] - pilot_bottom
    assert engage < depth - 0.5, f'螺絲太長:鎖入 {engage} mm,導孔只有 {depth:.1f} mm'
    assert SCREW_D > bosses[0][3], '螺絲直徑要大於導孔才咬得住'
    overlap = check_clearance(body, plate, tops[0], z_top)
    print(f'鎖入凸柱 {engage:.1f} mm(導孔深 {depth:.1f});底板與盒身重疊面積 {overlap:.3f} mm²')
    assert overlap < 0.5, '底板和盒身干涉'
    lid_z = min(m.bounds[0][2] for m in box[('box_lid', 'box_gray')])
    print(f'底板頂面 z={z_top:.1f},盒身開口 z={lid_z:.1f},可用高度約 {lid_z - z_top:.0f} mm(到盒緣)')

    # ---- RS485 擴展板:12 mm 六角柱,裝在底板較置中的那組 Ø3 孔 ----
    pcb = None
    if args.rs485:
        loops_d3 = [(c, d) for c, d in loops(plate.section([0, 0, 1], [0, 0, z_top - thick / 2]))
                    if abs(d - PCB_HOLE_D) < 0.2]
        sets = square_sets(loops_d3, PCB_HOLE_PITCH)
        assert sets, '底板上找不到 95 × 95 的 Ø3 孔組'
        plate_c = (plate.bounds[0][:2] + plate.bounds[1][:2]) / 2
        sets.sort(key=lambda s: np.hypot(*(s[0] - plate_c)))
        for c, _ in sets:
            print(f'  底板 Ø3 孔組:中心 ({c[0]:+.2f}, {c[1]:+.2f}),離底板中心 {np.hypot(*(c - plate_c)):.1f} mm')
        set_c, set_pts = sets[0]
        local, pcb_t = pcb_holes(args.rs485)
        R = rz(PCB_YAW)
        z_pcb = z_top + SPACER_L
        world = (R[:2, :2] @ (local - local.mean(0)).T).T + set_c
        err = max(np.min(np.hypot(*(set_pts - w).T)) for w in world)
        t = np.array([*(set_c - R[:2, :2] @ local.mean(0)), z_pcb])
        pcb = dict(R=R, t=t, holes=set_pts, thick=pcb_t)
        # 底板背面的螺絲頭不能碰到盒底
        head_bot = tops[0] - SCREW_HEAD_H
        for x, y in set_pts:
            hits, _, _ = body.ray.intersects_location([[x, y, tops[0] - 0.01]], [[0, 0, -1]])
            floor = float(hits[:, 2].max())
            assert floor < head_bot - 0.5, f'({x:.1f}, {y:.1f}) 底板背面螺絲頭碰到盒底'
        gap = head_bot - max(float(body.ray.intersects_location([[x, y, tops[0] - 0.01]], [[0, 0, -1]])[0][:, 2].max())
                             for x, y in set_pts)
        stud_in = TOP_POST_STUD - pcb_t
        assert 0 < stud_in < SPACER_THREAD_DEPTH and PCB_SCREW_L - thick <= SPACER_THREAD_DEPTH
        # 上方金屬柱頂端要在盒蓋下方(盒蓋內面:從柱頂往上打射線)
        post_top = z_pcb + pcb_t + TOP_POST_L
        lid = trimesh.util.concatenate(box[('box_lid', 'box_gray')])
        lid_gap = min(float(lid.ray.intersects_location([[x, y, post_top]], [[0, 0, 1]])[0][:, 2].min()) - post_top
                      for x, y in set_pts)
        rim_gap = lid_z - post_top
        assert lid_gap > 1.0, f'上方六角柱頂到盒蓋只剩 {lid_gap:.1f} mm'
        print(f'RS485 板:裝在中心 ({set_c[0]:+.2f}, {set_c[1]:+.2f}) 那組孔,轉 {PCB_YAW:.0f}°,'
              f'孔位最大誤差 {err:.3f} mm;板子背面 z={z_pcb:.1f},正面 z={z_pcb + pcb_t:.1f}')
        print(f'  尼龍六角柱 M3 × {SPACER_L:.0f};下螺絲鎖入 {PCB_SCREW_L - thick:.1f} mm;底板背面螺絲頭離盒底 {gap:.1f} mm')
        print(f'  上方金屬公母六角柱 M3 × {TOP_POST_L:.0f} + {TOP_POST_STUD:.0f}:公頭鎖入尼龍柱 {stud_in:.2f} mm,'
              f'柱頂 z={post_top:.1f},離盒緣 {rim_gap:.1f} mm、離盒蓋內面 {lid_gap:.1f} mm')

    # ---- 水情盒主板:金屬六角柱頂端 ----
    main = None
    if args.main:
        assert pcb, '--main 需要同時給 --rs485(主板裝在 RS485 板上方的金屬六角柱上)'
        mlocal, main_t = pcb_holes(args.main)
        Rm = rz(MAIN_YAW)
        world = (Rm[:2, :2] @ (mlocal - mlocal.mean(0)).T).T + set_c
        merr = max(np.min(np.hypot(*(set_pts - w).T)) for w in world)
        assert merr < 0.3, f'主板孔位對不上六角柱({merr:.2f} mm)'
        tm = np.array([*(set_c - Rm[:2, :2] @ mlocal.mean(0)), post_top])
        main = dict(R=Rm, t=tm, thick=main_t)
        assert PCB_SCREW_L - main_t < TOP_POST_THREAD
        # 主板與零件(世界座標)不能碰到盒蓋 / 盒身
        msc = trimesh.load(args.main, force='scene')
        Zup = np.array([[1, 0, 0, 0], [0, 0, -1, 0], [0, 1, 0, 0], [0, 0, 0, 1]], float)
        M = affine(Rm, tm)
        mm = []
        for n in msc.graph.nodes_geometry:
            T, gname = msc.graph[n]
            m = msc.geometry[gname].copy()
            m.apply_transform(Zup @ T)
            m.apply_scale(1000)
            m.apply_transform(M)
            mm.append(m)
        mmesh = trimesh.util.concatenate(mm)
        top_z = float(mmesh.bounds[1][2])
        lid_mesh = trimesh.util.concatenate(box[('box_lid', 'box_gray')])
        # 各零件分開切(合併後切,不同零件的截面線會被串成錯誤的多邊形)
        enclosure = [lid_mesh, body] + box[('box_gasket', 'box_gasket')] + box[('box_inserts', 'brass')]
        foot = Polygon(np.array([[*mmesh.bounds[0][:2]], [mmesh.bounds[1][0], mmesh.bounds[0][1]],
                                 [*mmesh.bounds[1][:2]], [mmesh.bounds[0][0], mmesh.bounds[1][1]]]))
        hit = 0.0
        for z in np.linspace(float(mmesh.bounds[0][2]) + 0.05, top_z, 12):
            for part_mesh in enclosure:
                sec = part_mesh.section([0, 0, 1], [0, 0, z])
                if sec is None:
                    continue
                for e in sec.discrete:
                    if len(e) > 3:
                        hit = max(hit, Polygon(e[:, :2]).buffer(0).boundary.intersection(foot).length)
        ceil = min(float(lid_mesh.ray.intersects_location([[x, y, top_z]], [[0, 0, 1]])[0][:, 2].min())
                   for x, y in [mmesh.bounds.mean(0)[:2]] + list(set_pts))
        print(f'主板:裝在金屬柱頂 z={post_top:.1f},轉 {MAIN_YAW:.0f}°,孔位誤差 {merr:.3f} mm;'
              f'板子正面 z={post_top + main_t:.1f},最高零件 z={top_z:.1f},離盒蓋內面 {ceil - top_z:.1f} mm')
        print(f'  外形範圍 X {mmesh.bounds[0][0]:.1f}~{mmesh.bounds[1][0]:.1f},Y {mmesh.bounds[0][1]:.1f}~'
              f'{mmesh.bounds[1][1]:.1f};盒壁穿過主板外框的長度 {hit:.2f} mm;螺絲鎖入金屬柱 {PCB_SCREW_L - main_t:.2f} mm')
        assert hit < 0.01 and ceil - top_z > 0.5, '主板碰到盒蓋或盒壁'

    # ---- SEN0676 雷達:-X 牆開 Ø40 孔,鏡頭朝外,六角一組對邊平行盒底 ----
    sen = None
    if args.sen0676:
        top = float(body.bounds[1][2])
        hz = top - SEN_HOLE_FROM_TOP
        hits, _, _ = body.ray.intersects_location([[0, 0, hz]], [[-1, 0, 0]])
        xs = sorted(hits[:, 0], reverse=True)
        inner_x, outer_x = xs[0], xs[-1]
        # (-X, +Y) 角落盒蓋螺絲柱:孔高度的截面上,貼著 -X 牆內面、y > 40 的凸出部分,取它朝 -Y 的邊
        pts = np.vstack(body.section([0, 0, 1], [0, 0, hz]).discrete)
        bump = pts[(pts[:, 0] > inner_x + 0.5) & (pts[:, 0] < inner_x + 15) & (pts[:, 1] > 40)]
        post_y = float(bump[:, 1].min())
        hy = post_y - SEN_HOLE_FROM_POST
        prof = sensor_profile(args.sen0676)
        # 螺紋(最大半徑)、螺紋根部凹槽、法蘭面
        r_thread = max(r for z, r in prof if z < 24.5)
        z_thread_end = max(z for z, r in prof if z < 26 and r > r_thread - 0.05)
        groove = [(z, r) for z, r in prof if z > z_thread_end and r < r_thread - 0.5]
        z_flange = next(z for z, r in prof if z > z_thread_end and r > r_thread + 1.0)
        r_groove = min(r for z, r in groove if z < z_flange)
        assert 2 * r_thread < SEN_HOLE_D, f'螺紋 Ø{2 * r_thread:.1f} 穿不過 Ø{SEN_HOLE_D} 孔'
        o_id, o_cs = SEN_ORING
        o_cs_on = o_cs * np.sqrt((o_id + o_cs) / (2 * r_groove + o_cs))   # 拉伸後線徑
        gap = o_cs_on * (1 - SEN_ORING_SQUEEZE)                           # 法蘭面到牆內面
        wall_t = inner_x - outer_x
        nut_w, nut_h = SEN_NUT
        z_wall_out = z_flange - gap - wall_t                              # 牆外面(雷達座標)
        z_nut_end = z_wall_out - nut_h
        full_thread = min(z for z, r in prof if r > r_thread - 0.05)
        # 雷達座標:軸 +Z(鏡頭 -> 背面),六角對邊垂直於 x(x = ±31 是平邊)。
        # 軸 -> 盒內 +X;x -> 上(+Z),平邊就平行盒底;y -> -Y(右手座標)
        Rs = np.array([[0, 0, 1], [0, -1, 0], [1, 0, 0]], float)
        ts = np.array([inner_x + gap - z_flange, hy, hz])           # 法蘭面在牆內面往盒內 gap
        sen = dict(R=Rs, t=ts)
        sc = trimesh.load(args.sen0676, force='scene')
        T, gname = sc.graph['housing']
        hm = sc.geometry[gname].copy()
        hm.apply_transform(np.array([[1, 0, 0, 0], [0, 0, -1, 0], [0, 1, 0, 0], [0, 0, 0, 1]], float) @ T)
        hm.apply_scale(1000)
        hm.apply_transform(affine(Rs, ts))
        lo, hi = hm.bounds
        print(f'SEN0676:-X 牆內面 x={inner_x:.2f}、外面 x={outer_x:.2f}(厚 {inner_x - outer_x:.2f}),'
              f'(-X, +Y) 角落螺絲柱邊緣 y={post_y:.2f}')
        print(f'  Ø{SEN_HOLE_D:.1f} 孔中心 y={hy:.2f}、z={hz:.2f}(盒身上緣 z={top:.2f} 往下 {SEN_HOLE_FROM_TOP:.0f}),'
              f'孔範圍 z {hz - SEN_HOLE_D / 2:.1f}~{hz + SEN_HOLE_D / 2:.1f}')
        print(f'  螺紋 Ø{2 * r_thread:.1f}(離頂點 {full_thread:.1f}~{z_thread_end:.1f}),根部凹槽 Ø{2 * r_groove:.1f}'
              f',法蘭面離頂點 {z_flange:.1f}')
        print(f'  O 型環 {o_id} × {o_cs}:套上拉伸 {(2 * r_groove / o_id - 1) * 100:.1f}%、線徑 {o_cs_on:.2f},'
              f'外緣 Ø{2 * r_groove + 2 * o_cs_on:.1f}(壓在孔邊外 {r_groove + o_cs_on - SEN_HOLE_D / 2:.2f} mm 寬),'
              f'壓縮 {SEN_ORING_SQUEEZE * 100:.0f}% -> 法蘭面離牆內面 {gap:.2f} mm')
        print(f'  螺帽 AGJL-112(對邊 {nut_w}、高 {nut_h}):鎖在螺紋 {z_nut_end:.1f}~{z_wall_out:.1f}'
              f'(完整螺紋從 {full_thread:.1f} 起,{"OK" if z_nut_end >= full_thread else "螺紋不夠長"});'
              f'雷達凸出牆外 {z_wall_out:.1f} mm(螺帽外 {z_nut_end:.1f} mm)')
        # 牆內面有拔模角:O 型環在孔上緣 / 下緣的實際壓縮率
        squeeze = []
        for zz in (hz + r_groove + o_cs_on / 2, hz - r_groove - o_cs_on / 2):
            h_, _, _ = body.ray.intersects_location([[0, hy, zz]], [[-1, 0, 0]])
            g_ = gap + (inner_x - float(h_[:, 0].max()))
            squeeze.append(1 - g_ / o_cs_on)
        print(f'  牆內面拔模:O 型環壓縮率 上緣 {squeeze[0] * 100:.0f}%、下緣 {squeeze[1] * 100:.0f}%')
        hole_top = hz + SEN_HOLE_D / 2
        gk = trimesh.util.concatenate(box[('box_gasket', 'box_gasket')])
        print(f'  孔頂 z={hole_top:.2f},密封條底 z={gk.bounds[0][2]:.2f}(相距 {gk.bounds[0][2] - hole_top:.2f} mm)')
        # 盒外:螺帽是否碰到盒蓋 / 盒身外側
        nut_top = hz + nut_w / 2
        lid_m = trimesh.util.concatenate(box[('box_lid', 'box_gray')])
        worst = 0.0
        for zz in np.linspace(hz - nut_w / 2, nut_top, 25):
            for yy in np.linspace(hy - nut_w / 2, hy + nut_w / 2, 9):
                for mesh_ in (lid_m, body):
                    h_, _, _ = mesh_.ray.intersects_location([[outer_x - 30, yy, zz]], [[1, 0, 0]])
                    if len(h_):
                        worst = max(worst, outer_x - float(h_[:, 0].min()))
        # 盒身外牆有拔模角(越往上越凸出),螺帽平面會先頂到最凸的地方 -> 螺帽往外移 worst
        print(f'  螺帽範圍內外牆最多凸出 {worst:.2f} mm(拔模角),螺帽貼在最凸處、往外移 {worst:.2f} mm;'
              f'螺帽鎖在螺紋 {z_nut_end - worst:.1f}~{z_wall_out - worst:.1f}'
              f'({"OK" if z_nut_end - worst >= full_thread else "螺紋不夠長"})')
        # 干涉(只列出,不擋)
        floor_z = float(body.ray.intersects_location([[inner_x + 25, hy, 100]], [[0, 0, -1]])[0][:, 2].max())
        print(f'  雷達外殼在盒內佔 x {inner_x:.1f}~{hi[0]:.1f}、y {lo[1]:.1f}~{hi[1]:.1f}、z {lo[2]:.1f}~{hi[2]:.1f}')
        print(f'  [{"干涉" if lo[2] < floor_z else "OK"}] 最低點 z={lo[2]:.2f},盒內底 z={floor_z:.2f}'
              f'({"壓進" if lo[2] < floor_z else "間隙"} {abs(floor_z - lo[2]):.2f} mm)')
        if plate.bounds[0][0] < hi[0] and lo[2] < z_top:
            print(f'  [干涉] 底板 -X 邊緣 x={plate.bounds[0][0]:.1f}(z {tops[0]:.1f}~{z_top:.1f})在雷達範圍內')
        if pcb and min(x for x, _ in set_pts) - SPACER_FLATS / 2 < hi[0]:
            print('  [干涉] 電路板六角柱在雷達範圍內')
        elif pcb:
            print(f'  [OK] 離電路板六角柱 {min(x for x, _ in set_pts) - SPACER_FLATS / 2 - hi[0]:.1f} mm')
        lid = trimesh.util.concatenate(box[('box_lid', 'box_gray')])
        up = lid.ray.intersects_location([[float(hi[0]) - 1, hy, float(hi[2])]], [[0, 0, 1]])[0]
        if len(up):
            print(f'  雷達頂 z={hi[2]:.1f},盒蓋內面 z={float(up[:, 2].min()):.1f}')
        box[('box_body', 'box_gray')] = [cut_hole_x(body, hy, hz, SEN_HOLE_D / 2, outer_x - 1, inner_x + 1)]
        tor = trimesh.creation.torus(r_groove + o_cs_on / 2, o_cs_on / 2, major_sections=96, minor_sections=16)
        tor.apply_scale([1, 1, 1 - SEN_ORING_SQUEEZE])              # 軸向被壓扁
        tor.apply_translation([0, 0, z_flange - gap / 2])
        sen['oring'] = tor
        nut = sen_nut_mesh(nut_w, nut_h, 2 * r_thread)
        nut.apply_translation([0, 0, z_nut_end - worst])
        sen['nut'] = nut

    # ---- 輸出 GLB ----
    if pcb:
        g = load_glb(args.rs485)               # 保留板子的貼圖與材質,其他零件加在後面
        nodes = g.j['nodes']
        pcb_children = nodes[0]['children']
        nodes.append(dict(name='rs485', children=pcb_children, matrix=gltf_matrix(pcb['R'], pcb['t']),
                          extras=make_extras(PARTS, EXPLODE, 'rs485')))
        nodes[0]['children'] = [len(nodes) - 1]
    else:
        g = GLB()
        nodes = g.j['nodes']
    nodes[0]['name'] = 'Water_box'
    if main:
        kids = merge_glb(g, args.main)
        nodes.append(dict(name='main_pcb', children=kids, matrix=gltf_matrix(main['R'], main['t']),
                          extras=make_extras(PARTS, EXPLODE, 'main_pcb')))
        nodes[0]['children'].append(len(nodes) - 1)
    if sen:
        n0 = len(nodes)
        kids = merge_glb(g, args.sen0676)
        for n in nodes[n0:]:                        # 零件清單:加上「SEN0676」前綴,避免和盒子零件混在一起
            part = n.get('extras', {}).get('part')
            if part:
                part['id'] = 'sen_' + part['id']
                part['name'] = 'SEN0676 ' + part['name']
        nodes.append(dict(name='sen0676', children=kids, matrix=gltf_matrix(sen['R'], sen['t']),
                          extras=make_extras(PARTS, EXPLODE, explode_key='sen0676')))
        nodes[0]['children'].append(len(nodes) - 1)
    mat_cache = {}

    def material(key):
        if key not in mat_cache:
            srgb, metal, rough = MATS[key]
            mat_cache[key] = g.material(dict(name=key, pbrMetallicRoughness=dict(
                baseColorFactor=[*srgb_to_linear(srgb).tolist(), 1.0], metallicFactor=metal, roughnessFactor=rough)))
        return mat_cache[key]

    def add_mesh(name, tm_mm, key, part, matrix=None, mesh=None):
        if mesh is None:
            tm = tm_mm.copy()
            tm.apply_scale(0.001)
            pos, nrm, idx = gltf_arrays(shade(tm, key, SMOOTH_ANGLE))
            g.mesh(name, pos, nrm, idx, material(key))
            node = nodes[-1]
        else:
            nodes.append(dict(name=name, mesh=mesh))
            nodes[0]['children'].append(len(nodes) - 1)
            node = nodes[-1]
        node['extras'] = make_extras(PARTS, EXPLODE, part)
        if matrix is not None:
            node['matrix'] = matrix
        return node['mesh']

    for (part, key), ms in box.items():
        add_mesh(part, trimesh.util.concatenate(ms), key, part)
    add_mesh('plate', plate, 'print_white', 'plate')
    if sen:
        add_mesh('sen_oring', sen['oring'], 'rubber', 'sen_oring', gltf_matrix(sen['R'], sen['t']))
        add_mesh('sen_nut', sen['nut'], 'nylon_black', 'sen_nut', gltf_matrix(sen['R'], sen['t']))

    mg = pan_screw(SCREW_D, SCREW_L, SCREW_HEAD_D, SCREW_HEAD_H).to_mesh()
    screw = trimesh.Trimesh(np.asarray(mg.vert_properties)[:, :3], np.asarray(mg.tri_verts), process=True)
    mesh_id = None
    for i, (x, y, _, _) in enumerate(bosses):   # 4 個節點共用一個網格,件數各記 1
        mesh_id = add_mesh(f'plate_screw_{i}', screw, 'screw', 'plate_screws',
                           gltf_matrix(t_mm=[x, y, z_top]), mesh_id)

    if pcb:
        def to_tm(m):
            mg = m.to_mesh()
            return trimesh.Trimesh(np.asarray(mg.vert_properties)[:, :3], np.asarray(mg.tri_verts), process=True)
        spacer = to_tm(standoff(3.0, SPACER_L, SPACER_FLATS, hole_depth=SPACER_THREAD_DEPTH))
        m3 = to_tm(pan_screw(3.0, PCB_SCREW_L))
        post = to_tm(standoff(3.0, TOP_POST_L, TOP_POST_FLATS, stud=TOP_POST_STUD))
        flip = np.diag([1.0, -1.0, -1.0])                  # 螺絲頭朝下(從底板背面鎖)
        ids = dict(s=None, top=None, bot=None)
        z_pcb = pcb['t'][2]
        for i, (x, y) in enumerate(pcb['holes']):
            ids['s'] = add_mesh(f'rs485_spacer_{i}', spacer, 'nylon', 'rs485_spacers',
                                gltf_matrix(t_mm=[x, y, z_pcb]), ids['s'])
            ids['top'] = add_mesh(f'rs485_post_top_{i}', post, 'brass', 'rs485_posts_top',
                                  gltf_matrix(t_mm=[x, y, z_pcb + pcb['thick'] + TOP_POST_L]), ids['top'])
            ids['bot'] = add_mesh(f'rs485_screw_bottom_{i}', m3, 'screw', 'rs485_screws_bottom',
                                  gltf_matrix(flip, [x, y, tops[0]]), ids['bot'])
            if main:
                ids['main'] = add_mesh(f'main_screw_{i}', m3, 'screw', 'main_screws',
                                       gltf_matrix(t_mm=[x, y, main['t'][2] + main['thick']]), ids.get('main'))

    g.save(args.out)
    print('輸出', args.out, f'{Path(args.out).stat().st_size / 1e6:.1f} MB')


if __name__ == '__main__':
    main()
