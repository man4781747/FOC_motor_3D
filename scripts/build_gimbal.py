# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "trimesh", "scipy", "networkx", "manifold3d", "pillow"]
# ///
"""用兩組「PCB + 2804 馬達」模組組成雙軸(水平 Pan + 俯仰 Tilt)雲台,
並輸出三個 3D 列印結構件的 STL(已擺好列印方向)。

用法:
    uv run scripts/build_gimbal.py <pcb_motor_assembly.glb> <輸出.glb> [--stl-dir output/print]

世界座標(mm,Z 朝上),水平軸 = Z 軸。

安裝方式(每個接合都確認過鎖得到):
1. 底座:八角底板 + 4 根凸柱。水平模組 PCB「正面朝下」放在凸柱上,
   M2.5x12 螺絲從底板下方沉頭孔穿過凸柱、PCB,鎖進模組原本的銅柱
   (取代模組正面原本的 4 顆短螺絲)。馬達朝上、轉子面朝上;
   定子與 PCB 固定在底座,馬達線不會被扭到。
2. 框架(平台 + 立板一體):平台用 4 顆 M2.5 從上方鎖進水平馬達轉子面螺孔(Ø16)。
3. 俯仰模組:PCB 正面朝立板,以 4 支 M2.5 公母銅柱 8 mm 連接 ——
   公頭穿過 PCB 鎖進模組銅柱,螺絲從立板外側鎖進母頭。馬達朝 +X 水平伸出。
4. 相機支架:4 顆 M2.5 從支架外側鎖進俯仰馬達轉子面(Ø16);相機用 4 顆 M2 鎖在前板。

3D 列印考量:每件都有大平面貼床、懸空 <= 45°(內角圓角、補強肋、孔頂皆 45°)、
水平大孔用水滴形、壁厚 >= 2.4 mm、M2.5 孔 Ø2.9 / M3 孔 Ø3.4 / M2 孔 Ø2.3。

節點階層(three.js / Blender 可轉動):
    Gimbal_2axis
      ├ base_fixed     底座、水平模組(轉子以外)
      └ pan_group      繞 Z 軸(glTF Y 軸):水平馬達轉子、框架、俯仰模組(轉子以外)
          └ tilt_group 以俯仰軸為樞紐繞 X 軸:俯仰馬達轉子、相機支架、相機
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import trimesh
from manifold3d import CrossSection, JoinType, Manifold

sys.path.insert(0, str(Path(__file__).parent))
from build_assembly import load_glb  # noqa: E402
from build_glb import SMOOTH_ANGLE, gltf_arrays, srgb_to_linear  # noqa: E402
from build_motor import box, cyl, polar, revolve, union  # noqa: E402

# ---------------- 模組尺寸(模組座標:PCB 背面 z=0、正面 z=2.54,馬達朝 -Z) ----------------
MODULE_PCB_TOP = 2.54
MODULE_ROTOR_Z = -21.5                                    # 轉子面
MOUNT = [polar(9.5, a) for a in (0, 90, 180, 270)]        # PCB 安裝孔 = 模組銅柱位置
ROTOR_HOLES = [polar(8.0, a) for a in (0, 90, 180, 270)]  # 轉子面 M2.5 孔(Ø16)

# ---------------- 列印件與組裝尺寸(mm) ----------------
M25, M3, M2 = 2.9, 3.4, 2.3     # 螺絲孔徑(含列印公差)
FOOT_T = 4.0                    # 底板厚
BOSS_H = 8.0                    # 凸柱高(PCB 正面零件最高約 3.5 mm)
PCB_FACE_Z = FOOT_T + BOSS_H    # 水平模組 PCB 正面(朝下)的高度
R_PAN = np.diag([1.0, -1.0, -1.0])                        # 繞 X 翻轉 180°:馬達朝上
PAN_T = np.array([0, 0, PCB_FACE_Z + MODULE_PCB_TOP])
ROTOR_FACE_Z = PAN_T[2] - MODULE_ROTOR_Z                  # 水平馬達轉子面 = 36.04
PLATFORM_Z0 = ROTOR_FACE_Z
PLATFORM_T = 4.0
PLATFORM_TOP = PLATFORM_Z0 + PLATFORM_T
TILT_Z = PLATFORM_TOP + 36.0                              # 俯仰軸高度
SPACER = 8.0                                              # 俯仰模組 PCB 到立板的公母銅柱
ROTOR_X = -20.0                                           # 俯仰馬達轉子面的世界 X
TILT_T = np.array([ROTOR_X + MODULE_ROTOR_Z, 0.0, TILT_Z])
PLATE_X1 = TILT_T[0] - (MODULE_PCB_TOP + SPACER)          # 立板內側面
PLATE_X0 = PLATE_X1 - 5.0                                 # 立板厚 5
R_TILT = np.array([[0, 0, -1], [-1, 0, 0], [0, 1, 0]], float)  # 模組 -Z(馬達)-> 世界 +X
TILT_PIVOT = np.array([0.0, 0.0, TILT_Z])
CAM_X = 3.0                                               # 相機中心 X(避開支架內角圓角)

C = np.array([[1, 0, 0], [0, 0, 1], [0, -1, 0]], float)   # Z-up -> glTF Y-up

MATS = {   # sRGB, 金屬度, 粗糙度
    'print_white': ((0.92, 0.92, 0.90), 0.0, 0.6),      # 3D 列印 PLA
    'brass': ((0.80, 0.62, 0.32), 1.0, 0.3),
    'screw': ((0.80, 0.80, 0.82), 1.0, 0.25),
    'cam_pcb': ((0.08, 0.09, 0.10), 0.0, 0.45),
    'cam_holder': ((0.05, 0.05, 0.05), 0.0, 0.5),
    'lens_barrel': ((0.10, 0.10, 0.11), 0.5, 0.35),
    'lens_glass': ((0.05, 0.07, 0.14), 0.0, 0.05),
}


# ---------------- 幾何小工具 ----------------
def along_x(m):
    """沿 +Z 建立的形狀 -> 沿 +X(x' = z)。"""
    return m.rotate([0, 90, 0])


def along_y(m):
    """沿 +Z 建立的形狀 -> 沿 +Y(y' = z)。"""
    return m.rotate([-90, 0, 0])


def hexagon(flats, h):
    r = flats / np.sqrt(3)
    return Manifold.cylinder(h, r, circular_segments=6).rotate([0, 0, 30])


def rounded_rect(w, h, t, r):
    cs = CrossSection.square([w - 2 * r, h - 2 * r], center=True).offset(r, JoinType.Round, 2.0, 32)
    return cs.extrude(t)


def teardrop_x(r, x0, length, y, z):
    """沿 X 的水滴形孔(尖端朝 +Z),水平列印時孔頂 45° 可自我支撐。"""
    cs = CrossSection.batch_hull([CrossSection.circle(r, 48),
                                  CrossSection.square([0.01, 0.01]).translate([0, r * np.sqrt(2)])])
    m = cs.extrude(length).rotate([90, 0, 0]).rotate([0, 0, 90])   # (u, v, w) -> (w, u, v)
    return m.translate([x0, y, z])


def ccw(pts):
    """CrossSection 需要逆時針頂點;順時針會變成無效實體,合併時把結果整個清空。"""
    p = np.asarray(pts, float)
    area = np.sum(p[:, 0] * np.roll(p[:, 1], -1) - np.roll(p[:, 0], -1) * p[:, 1])
    return [tuple(q) for q in (p if area > 0 else p[::-1])]


def prism_y(tri_xz, y0, y1):
    """XZ 平面三角形沿 Y 擠出。"""
    m = CrossSection([ccw(tri_xz)]).extrude(y1 - y0)     # (x, z2d, w)
    return m.rotate([90, 0, 0]).translate([0, y1, 0])    # -> (x, -w, z2d),y ∈ [y0, y1]


def prism_z(tri_xy, z0, z1):
    return CrossSection([ccw(tri_xy)]).extrude(z1 - z0).translate([0, 0, z0])


def pan_head(d, h):
    """十字盤頭螺絲頭(頭頂朝 +Z,座面 z=0)。"""
    r = d / 2
    head = revolve([(0, 0), (r, 0), (r, h * 0.55), (r * 0.8, h * 0.9), (r * 0.5, h), (0, h)])
    s = 0.5
    cross = union([box(-r * 0.55, r * 0.55, -s / 2, s / 2, h * 0.45, h + 1),
                   box(-s / 2, s / 2, -r * 0.55, r * 0.55, h * 0.45, h + 1)])
    return head - cross


def tilt_world(p_module):
    return R_TILT @ np.asarray(p_module, float) + TILT_T


def pan_world(p_module):
    return R_PAN @ np.asarray(p_module, float) + PAN_T


# ---------------- 列印件 1:底座 ----------------
def make_base():
    """八角底板 + 4 根凸柱(頂住水平模組 PCB 正面)。列印方向:底板朝下。"""
    foot = Manifold.cylinder(0.6, 33.4, 34, circular_segments=8)            # 底緣倒角(防象腳)
    foot += Manifold.cylinder(FOOT_T - 0.6, 34, circular_segments=8).translate([0, 0, 0.6])
    foot = foot.rotate([0, 0, 22.5])
    bosses, cuts = [], [cyl(6, 10, -1)]                                     # 中央走線孔
    cuts += [cyl(M3 / 2, 10, -1, *polar(27, a), seg=32) for a in (22.5, 112.5, 202.5, 292.5)]
    for m in MOUNT:
        x, y, _ = pan_world([*m, MODULE_PCB_TOP])
        # 頂端 Ø4.6(與原本螺絲頭一樣大,避開 LQFP / 晶振),往下 45° 擴大到 Ø7
        bosses.append(revolve([(0, FOOT_T - 0.01), (3.5, FOOT_T - 0.01), (3.5, PCB_FACE_Z - 3.2),
                               (2.3, PCB_FACE_Z - 2.0), (2.3, PCB_FACE_Z), (0, PCB_FACE_Z)]).translate([x, y, 0]))
        cuts.append(cyl(M25 / 2, 20, -1, x, y, 32))
        cuts.append(cyl(5.6 / 2, 7.0, -1, x, y, 32))                       # 底面沉頭孔(深 6 mm)
    return (foot + union(bosses)) - union(cuts)


def base_screws():
    """M2.5x12 螺絲頭在底座沉頭孔內(頭朝下)。"""
    heads = []
    for m in MOUNT:
        x, y, _ = pan_world([*m, 0])
        heads.append(pan_head(4.5, 1.7).mirror([0, 0, 1]).translate([x, y, 6.0]))
    return union(heads)


# ---------------- 列印件 2:框架(平台 + 立板) ----------------
def make_frame():
    """平台鎖在水平馬達轉子上,立板支撐俯仰模組。列印方向:平台朝下。"""
    zb = PLATFORM_TOP
    px0, px1 = PLATE_X0, 32.0
    plat = rounded_rect(px1 - px0, 64, PLATFORM_T, 4).translate([(px0 + px1) / 2, 0, PLATFORM_Z0])
    cuts = [cyl(4.5, 20, PLATFORM_Z0 - 1)]                                  # 中空軸 / 走線孔
    for x, y in ROTOR_HOLES:
        cuts.append(cyl(M25 / 2, 20, PLATFORM_Z0 - 1, x, y, 32))
        cuts.append(cyl(5.6 / 2, 5, PLATFORM_TOP - 2.0, x, y, 32))         # 頂面沉頭
    # 立板(YZ 平面),上緣圓角
    top = TILT_Z + 26
    plate = rounded_rect(top - zb + 10, 52, PLATE_X1 - PLATE_X0, 6)
    plate = plate.rotate([0, 90, 0]).translate([PLATE_X0, 0, (top + zb) / 2 - 5])
    plate = plate.trim_by_plane([0, 0, 1], zb - 0.01)
    cuts.append(teardrop_x(6.0, PLATE_X0 - 1, 8, 0, TILT_Z))               # 走線孔(水滴形)
    for m in MOUNT:
        _, wy, wz = tilt_world([*m, 0])
        cuts.append(teardrop_x(M25 / 2, PLATE_X0 - 1, 8, wy, wz))
        cuts.append(teardrop_x(5.6 / 2, PLATE_X0 - 1, 2.5, wy, wz))         # 外側沉頭
    # 45° 內角圓角 + 兩側 45° 補強肋
    fillet = prism_y([(PLATE_X1 - 0.01, zb - 0.01), (PLATE_X1 + 6, zb - 0.01), (PLATE_X1 - 0.01, zb + 6)], -26, 26)
    ribs = [prism_y([(PLATE_X1 - 0.01, zb - 0.01), (PLATE_X1 + 18, zb - 0.01), (PLATE_X1 - 0.01, zb + 18)],
                    y0, y0 + 3) for y0 in (-26, 23)]
    return (plat + plate + fillet + union(ribs)) - union(cuts)


def frame_hardware():
    heads = [pan_head(4.5, 1.5).translate([x, y, PLATFORM_TOP - 2.0]) for x, y in ROTOR_HOLES]
    standoffs = []
    for m in MOUNT:
        wx, wy, wz = tilt_world([*m, MODULE_PCB_TOP])
        standoffs.append(along_x(hexagon(4.0, SPACER)).mirror([1, 0, 0]).translate([wx, wy, wz]))  # 往 -X 到立板
        heads.append(along_x(pan_head(4.5, 1.5)).mirror([1, 0, 0]).translate([PLATE_X0 + 1.5, wy, wz]))
    return {'screw': union(heads), 'brass': union(standoffs)}


# ---------------- 列印件 3:相機支架 ----------------
def make_camera_bracket():
    """L 形支架:法蘭鎖轉子面,前板裝相機。列印方向:法蘭外側面(-X)朝下。"""
    zc = TILT_Z
    x0, x1 = ROTOR_X, ROTOR_X + 4.0           # 法蘭(厚 4)
    yf0, yf1 = 18.0, 21.2                     # 前板(厚 3.2)
    flange = box(x0, x1, -18, yf1, zc - 18, zc + 18)
    front = box(x0, CAM_X + 19, yf0, yf1, zc - 18, zc + 18)
    fillet = prism_z([(x1 - 0.01, yf0 + 0.01), (x1 + 3, yf0 + 0.01), (x1 - 0.01, yf0 - 3)], zc - 18, zc + 18)
    gussets = [prism_z([(x1 - 0.01, yf0 + 0.01), (x1 + 10, yf0 + 0.01), (x1 - 0.01, yf0 - 10)], z0, z0 + 2.4)
               for z0 in (zc - 18, zc + 15.6)]
    cuts = [along_x(cyl(4.0, 10)).translate([x0 - 2, 0, zc])]               # 中空軸孔
    for v in ROTOR_HOLES:
        _, wy, wz = tilt_world([*v, MODULE_ROTOR_Z])
        cuts.append(along_x(cyl(M25 / 2, 10, seg=32)).translate([x0 - 2, wy, wz]))
    # 鏡頭座方孔:列印時朝上的一側(+X)做 45° 尖頂,不需支撐
    hole = CrossSection([[(CAM_X - 7.5, -7.5), (CAM_X + 7.5, -7.5), (CAM_X + 15.0, 0.0),
                          (CAM_X + 7.5, 7.5), (CAM_X - 7.5, 7.5)]])
    cuts.append(hole.extrude(10).rotate([90, 0, 0]).translate([0, yf1 + 3, zc]))   # (u,v,w)->(u,-w,v)
    for dx in (-12, 12):
        for dz in (-12, 12):
            cuts.append(along_y(cyl(M2 / 2, 10, seg=24)).translate([CAM_X + dx, yf0 - 2, zc + dz]))
    return (flange + front + fillet + union(gussets)) - union(cuts)


def camera_parts():
    zc, p = TILT_Z, {}
    p['cam_pcb'] = box(CAM_X - 15, CAM_X + 15, 16.4, 18.0, zc - 15, zc + 15)
    p['cam_holder'] = box(CAM_X - 7, CAM_X + 7, 18.0, 25.0, zc - 7, zc + 7)
    p['lens_barrel'] = along_y(cyl(6.0, 9.0)).translate([CAM_X, 25.0, zc]) - \
        along_y(cyl(4.5, 2.0)).translate([CAM_X, 33.2, zc])
    p['lens_glass'] = along_y(cyl(4.5, 0.6)).translate([CAM_X, 33.0, zc])
    heads = []
    for dx in (-12, 12):
        for dz in (-12, 12):
            heads.append(along_y(pan_head(3.8, 1.3)).translate([CAM_X + dx, 21.2, zc + dz]))
    for v in ROTOR_HOLES:
        _, wy, wz = tilt_world([*v, MODULE_ROTOR_Z])
        heads.append(along_x(pan_head(4.5, 1.5)).translate([ROTOR_X + 4.0, wy, wz]))
    p['screw'] = union(heads)
    return p


# ---------------- 列印方向 / 檢查 ----------------
PRINT_ORIENT = {   # 世界座標 -> 列印方向(旋轉後自動貼床、置中)
    'base': lambda m: m,
    'frame': lambda m: m,
    'camera_bracket': lambda m: m.rotate([0, -90, 0]),   # (x,y,z)->(-z,y,x):法蘭外側面朝下
}


def export_stl(name, man, out_dir):
    m = PRINT_ORIENT[name](man)
    bb = np.asarray(m.bounding_box(), float)
    m = m.translate([-(bb[0] + bb[3]) / 2, -(bb[1] + bb[4]) / 2, -bb[2]])
    mg = m.to_mesh()
    tm = trimesh.Trimesh(np.asarray(mg.vert_properties)[:, :3], np.asarray(mg.tri_verts), process=True)
    tm.export(out_dir / f'{name}.stl')
    # 懸空檢查:朝下且比 45° 更陡、又不在床面上的面積
    n, c, a = tm.face_normals, tm.triangles_center, tm.area_faces
    over = (n[:, 2] < -np.cos(np.radians(45)) - 1e-3) & (c[:, 2] > 0.3)
    size = tm.extents
    print(f'  {name:15} {size[0]:5.1f} x {size[1]:5.1f} x {size[2]:5.1f} mm  體積 {tm.volume / 1000:5.1f} cm³  '
          f'水密 {tm.is_watertight}  懸空 {a[over].sum():5.1f} mm²(最高 z {c[over, 2].max() if over.any() else 0:.1f})')


# ---------------- glTF 輔助 ----------------
def gltf_matrix(R=np.eye(3), t_mm=np.zeros(3)):
    M = np.eye(4)
    M[:3, :3] = C @ R @ C.T
    M[:3, 3] = C @ (np.asarray(t_mm) / 1000.0)
    return M.T.ravel().tolist()


def to_trimesh(m, offset_mm=np.zeros(3)):
    mg = m.to_mesh()
    v = (np.asarray(mg.vert_properties)[:, :3] - offset_mm) / 1000.0
    return trimesh.Trimesh(v, np.asarray(mg.tri_verts), process=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('src')
    ap.add_argument('out')
    ap.add_argument('--stl-dir', default='output/print')
    args = ap.parse_args()

    g = load_glb(args.src)
    nodes = g.j['nodes']
    module_nodes = list(nodes[0]['children'])
    by_name = {nodes[i]['name']: i for i in module_nodes}
    rotor_names = {'motor_rotor', 'shaft'}
    skip = {'screws'}           # 模組正面原本的短螺絲 -> 改由底座長螺絲 / 公母銅柱取代
    nodes[0]['children'] = []
    mat_cache = {}

    def material(key):
        if key not in mat_cache:
            srgb, metal, rough = MATS[key]
            mat_cache[key] = g.material(dict(name=key, pbrMetallicRoughness=dict(
                baseColorFactor=[*srgb_to_linear(srgb).tolist(), 1.0],
                metallicFactor=metal, roughnessFactor=rough)))
        return mat_cache[key]

    def add_node(node, parent):
        nodes.append(node)
        nodes[parent].setdefault('children', []).append(len(nodes) - 1)
        return len(nodes) - 1

    def add_mesh(name, man, key, parent, offset_mm=np.zeros(3)):
        mesh = to_trimesh(man, offset_mm)
        pos, nrm, idx = gltf_arrays(trimesh.graph.smooth_shade(mesh, angle=SMOOTH_ANGLE))
        g.mesh(name, pos, nrm, idx, material(key))
        nodes[0]['children'].pop()
        nodes[parent].setdefault('children', []).append(len(nodes) - 1)

    def module_copy(prefix, names, parent, matrix):
        holder = add_node(dict(name=prefix, matrix=matrix), parent)
        for n in names:
            add_node(dict(name=f'{prefix}_{n}', mesh=nodes[by_name[n]]['mesh']), holder)

    nodes[0]['name'] = 'Gimbal_2axis'
    nodes[0]['rotation'] = [0, 1, 0, 0]          # 繞垂直軸 180°,讓鏡頭朝向檢視器預設視角
    base = add_node(dict(name='base_fixed'), 0)
    pan = add_node(dict(name='pan_group'), 0)
    tilt = add_node(dict(name='tilt_group', translation=(C @ TILT_PIVOT / 1000).tolist()), pan)

    rest = [n for n in by_name if n not in rotor_names | skip]
    rotor = [n for n in by_name if n in rotor_names]
    module_copy('pan_module', rest, base, gltf_matrix(R_PAN, PAN_T))               # 定子 + PCB 固定
    module_copy('pan_module_rotor', rotor, pan, gltf_matrix(R_PAN, PAN_T))         # 轉子帶動上方
    module_copy('tilt_module', rest, pan, gltf_matrix(R_TILT, TILT_T))
    module_copy('tilt_module_rotor', rotor, tilt, gltf_matrix(R_TILT, TILT_T - TILT_PIVOT))

    printed = {'base': make_base(), 'frame': make_frame(), 'camera_bracket': make_camera_bracket()}
    add_mesh('print_base', printed['base'], 'print_white', base)
    add_mesh('base_screws', base_screws(), 'screw', base)
    add_mesh('print_frame', printed['frame'], 'print_white', pan)
    for key, man in frame_hardware().items():
        add_mesh(f'frame_{key}', man, key, pan)
    add_mesh('print_camera_bracket', printed['camera_bracket'], 'print_white', tilt, TILT_PIVOT)
    for key, man in camera_parts().items():
        add_mesh(f'camera_{key}', man, key, tilt, TILT_PIVOT)

    g.save(args.out)
    print('輸出', args.out, f'{Path(args.out).stat().st_size / 1e6:.1f} MB')
    print(f'凸柱頂 z={PCB_FACE_Z},水平轉子面 z={ROTOR_FACE_Z:.2f},平台頂 z={PLATFORM_TOP:.2f},'
          f'俯仰軸 z={TILT_Z:.2f},立板 x={PLATE_X0:.2f}~{PLATE_X1:.2f}')

    stl_dir = Path(args.stl_dir)
    stl_dir.mkdir(parents=True, exist_ok=True)
    print('列印件 STL(已擺好列印方向):')
    for name, man in printed.items():
        export_stl(name, man, stl_dir)


if __name__ == '__main__':
    main()
