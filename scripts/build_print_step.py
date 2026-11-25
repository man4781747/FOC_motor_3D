# /// script
# requires-python = ">=3.10,<3.13"
# dependencies = ["numpy", "trimesh", "scipy", "networkx", "manifold3d", "pillow", "rtree", "cadquery"]
# ///
"""把雲台的三個 3D 列印件輸出成 STEP(真正的 B-rep 實體),
可在 FreeCAD / SolidWorks / Fusion 等 CAD 軟體裡直接修改。

用法:
    uv run scripts/build_print_step.py [--out-dir output/step] [--no-check]

- 尺寸常數全部從 build_gimbal.py 匯入,改尺寸只要改那邊,STL 與 STEP 會一致
- 幾何與 build_gimbal.py 的 make_base / make_frame / make_camera_bracket 逐一對應,
  但圓孔、圓角、水滴孔都是真正的曲面(manifold3d 版本是多邊形近似)
- STEP 為「組裝位置」的世界座標(mm,Z 朝上),方便在 CAD 裡與其他零件對位;
  要列印用的方向請用 build_gimbal.py 輸出的 STL
- 預設會把 STEP 細分成網格,與 manifold3d 版本做布林差集,檢查兩者形狀一致
"""
import argparse
import sys
from pathlib import Path

import cadquery as cq
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from build_gimbal import (CAM_BOSS_R, CAM_BOSS_Y1, CAM_HOLES, CAM_X, FOOT_T, M3, M25, MODULE_PCB_TOP,  # noqa: E402
                          MODULE_ROTOR_Z, MOUNT,
                          PCB_FACE_Z, PLATE_X0, PLATE_X1, PLATFORM_T, PLATFORM_TOP, PLATFORM_Z0,
                          ROTOR_HOLES, ROTOR_X, TILT_Z, pan_world, polar, tilt_world)

V = cq.Vector


# ---------------- 幾何小工具(對應 build_gimbal / build_motor 的同名函式) ----------------
def cyl(r, h, z0=0.0, x=0.0, y=0.0):
    return cq.Solid.makeCylinder(r, h, V(x, y, z0))


def cyl_x(r, h, x0, y, z):
    """沿 +X 的圓柱,x ∈ [x0, x0 + h]。"""
    return cq.Solid.makeCylinder(r, h, V(x0, y, z), V(1, 0, 0))


def cyl_y(r, h, x, y0, z):
    return cq.Solid.makeCylinder(r, h, V(x, y0, z), V(0, 1, 0))


def box(x0, x1, y0, y1, z0, z1):
    return cq.Solid.makeBox(x1 - x0, y1 - y0, z1 - z0, V(x0, y0, z0))


def revolve(profile, x=0.0, y=0.0):
    """profile: [(r, z), ...] 繞 Z 軸旋轉一圈。"""
    wp = cq.Workplane('XZ', origin=(x, 0, 0)).polyline(profile).close()
    return wp.revolve(360, (0, 0, 0), (0, 1, 0)).val().translate(V(0, y, 0))


def prism_y(tri_xz, y0, y1):
    """XZ 平面多邊形沿 Y 擠出,y ∈ [y0, y1](XZ 工作平面法向量是 -Y,所以從 y1 往回擠)。"""
    return cq.Workplane('XZ', origin=(0, y1, 0)).polyline(tri_xz).close().extrude(y1 - y0).val()


def prism_z(tri_xy, z0, z1):
    return cq.Workplane('XY', origin=(0, 0, z0)).polyline(tri_xy).close().extrude(z1 - z0).val()


def teardrop_x(r, x0, length, y, z):
    """沿 X 的水滴形孔(尖端朝 +Z),與 build_gimbal.teardrop_x 相同:圓 + 往上 r·√2 的尖頂。"""
    s = r / np.sqrt(2)
    wp = (cq.Workplane('YZ', origin=(x0, y, z))        # 區域座標 (u, v) = 世界 (Y, Z)
          .moveTo(s, s).threePointArc((0, -r), (-s, s)).lineTo(0, r * np.sqrt(2)).close())
    return wp.extrude(length).val()


def teardrop_y(r, y0, length, x, z, tip):
    """沿 Y 的水滴形,y ∈ [y0, y0 + length];tip = +1 尖端朝 +X,-1 朝 -X(同 build_gimbal.teardrop_y)。"""
    s = r / np.sqrt(2)
    wp = (cq.Workplane('XZ', origin=(0, y0 + length, 0))       # 區域座標 (u, v) = 世界 (X, Z),往 -Y 擠出
          .moveTo(x + tip * s, z + s).threePointArc((x - tip * r, z), (x + tip * s, z - s))
          .lineTo(x + tip * r * np.sqrt(2), z).close())
    return wp.extrude(length).val()


def fuse(*shapes):
    out = shapes[0]
    for s in shapes[1:]:
        out = out.fuse(s)
    return out


def finish(body, cuts):
    return body.cut(*cuts).clean()


# ---------------- 列印件 1:底座 ----------------
def make_base():
    def octagon(r, z):
        pts = [V(*polar(r, 22.5 + 45 * k), z) for k in range(8)]
        return cq.Wire.makePolygon(pts, close=True)

    top = cq.Face.makeFromWires(octagon(34.0, 0.6))
    # 底緣倒角(防象腳):從 z=0.6 往下帶拔模角擠出,外接圓 34 -> 33.4。
    # 不用 loft,因為 loft 會把斜面存成 B-spline,CAD 裡不好選、FeatureWorks 也辨識不了
    taper = np.degrees(np.arctan((34.0 - 33.4) * np.cos(np.radians(22.5)) / 0.6))
    foot = cq.Solid.extrudeLinear(top, V(0, 0, -0.6), taper=taper)
    foot = foot.fuse(cq.Solid.extrudeLinear(top, V(0, 0, FOOT_T - 0.6)))
    bosses, cuts = [], [cyl(6, 10, -1)]                                        # 中央走線孔
    cuts += [cyl(M3 / 2, 10, -1, *polar(27, a)) for a in (22.5, 112.5, 202.5, 292.5)]
    for m in MOUNT:
        x, y, _ = pan_world([*m, MODULE_PCB_TOP])
        # 頂端 Ø4.6(避開 LQFP / 晶振),往下 45° 擴大到 Ø7
        bosses.append(revolve([(0, FOOT_T - 0.01), (3.5, FOOT_T - 0.01), (3.5, PCB_FACE_Z - 3.2),
                               (2.3, PCB_FACE_Z - 2.0), (2.3, PCB_FACE_Z), (0, PCB_FACE_Z)], x, y))
        cuts.append(cyl(M25 / 2, 20, -1, x, y))
        cuts.append(cyl(5.6 / 2, 7.0, -1, x, y))                               # 底面沉頭孔(深 6 mm)
    return finish(fuse(foot, *bosses), cuts)


# ---------------- 列印件 2:框架(平台 + 立板) ----------------
def make_frame():
    zb = PLATFORM_TOP
    px0, px1 = PLATE_X0, 32.0
    plat = (cq.Workplane('XY', origin=((px0 + px1) / 2, 0, PLATFORM_Z0))
            .rect(px1 - px0, 64).extrude(PLATFORM_T).edges('|Z').fillet(4).val())
    cuts = [cyl(4.5, 20, PLATFORM_Z0 - 1)]                                     # 中空軸 / 走線孔
    for x, y in ROTOR_HOLES:
        cuts.append(cyl(M25 / 2, 20, PLATFORM_Z0 - 1, x, y))
        cuts.append(cyl(5.6 / 2, 5, PLATFORM_TOP - 2.0, x, y))                 # 頂面沉頭
    # 立板(YZ 平面),上緣兩角 R6
    top = TILT_Z + 26
    plate = (cq.Workplane().add(box(PLATE_X0, PLATE_X1, -26, 26, zb - 0.01, top))
             .edges('|X').edges('>Z').fillet(6).val())
    cuts.append(teardrop_x(6.0, PLATE_X0 - 1, 8, 0, TILT_Z))                  # 走線孔(水滴形)
    for m in MOUNT:
        _, wy, wz = tilt_world([*m, 0])
        cuts.append(teardrop_x(M25 / 2, PLATE_X0 - 1, 8, wy, wz))
        cuts.append(teardrop_x(5.6 / 2, PLATE_X0 - 1, 2.5, wy, wz))            # 外側沉頭
    # 45° 內角圓角 + 兩側 45° 補強肋
    fillet = prism_y([(PLATE_X1 - 0.01, zb - 0.01), (PLATE_X1 + 6, zb - 0.01), (PLATE_X1 - 0.01, zb + 6)], -26, 26)
    ribs = [prism_y([(PLATE_X1 - 0.01, zb - 0.01), (PLATE_X1 + 18, zb - 0.01), (PLATE_X1 - 0.01, zb + 18)],
                    y0, y0 + 3) for y0 in (-26, 23)]
    return finish(fuse(plat, plate, fillet, *ribs), cuts)


# ---------------- 列印件 3:相機支架 ----------------
def make_camera_bracket():
    zc = TILT_Z
    x0, x1 = ROTOR_X, ROTOR_X + 4.0           # 法蘭(厚 4)
    yf0, yf1 = 18.0, 21.2                     # 前板(厚 3.2)
    flange = box(x0, x1, -18, yf1, zc - 18, zc + 18)
    front = box(x0, CAM_X + 19, yf0, yf1, zc - 18, zc + 18)
    fillet = prism_z([(x1 - 0.01, yf0 + 0.01), (x1 + 3, yf0 + 0.01), (x1 - 0.01, yf0 - 3)], zc - 18, zc + 18)
    gussets = [prism_z([(x1 - 0.01, yf0 + 0.01), (x1 + 10, yf0 + 0.01), (x1 - 0.01, yf0 - 10)], z0, z0 + 2.4)
               for z0 in (zc - 18, zc + 15.6)]
    cuts = [cyl_x(4.0, 10, x0 - 2, 0, zc)]                                    # 中空軸孔
    for v in ROTOR_HOLES:
        _, wy, wz = tilt_world([*v, MODULE_ROTOR_Z])
        cuts.append(cyl_x(M25 / 2, 10, x0 - 2, wy, wz))
    # 鏡頭座方孔:列印時朝上的一側(+X)做 45° 尖頂;XZ 平面的五邊形沿 -Y 擠出,y ∈ [yf1 - 7, yf1 + 3]
    hole = [(CAM_X - 7.5, -7.5), (CAM_X + 7.5, -7.5), (CAM_X + 15.0, 0.0), (CAM_X + 7.5, 7.5), (CAM_X - 7.5, 7.5)]
    cuts.append(cq.Workplane('XZ', origin=(0, yf1 + 3, zc)).polyline(hole).close().extrude(10).val())
    # 前板正面 4 個水滴形螺絲柱,背面開 M3 熱熔螺母孔
    from hardware import INSERT_M3
    bosses = []
    for x, z in CAM_HOLES:
        bosses.append(teardrop_y(CAM_BOSS_R, yf1 - 0.01, CAM_BOSS_Y1 - yf1 + 0.01, x, z, -1))
        cuts.append(teardrop_y(INSERT_M3['hole_d'] / 2, yf0 - 1, INSERT_M3['hole_depth'] + 1, x, z, +1))
    return finish(fuse(flange, front, fillet, *gussets, *bosses), cuts)


PARTS = {'base': make_base, 'frame': make_frame, 'camera_bracket': make_camera_bracket}


# ---------------- 與 manifold3d 版本比對 ----------------
def compare(name, shape):
    """STEP 實體細分後與 build_gimbal 的 manifold3d 版本做布林差集,回報體積差。"""
    import build_gimbal
    from manifold3d import Manifold, Mesh

    ref = {'base': build_gimbal.make_base, 'frame': build_gimbal.make_frame,
           'camera_bracket': build_gimbal.make_camera_bracket}[name]()
    import trimesh

    verts, tris = shape.tessellate(0.01, 0.1)
    tm = trimesh.Trimesh(np.array([v.toTuple() for v in verts]), np.array(tris))   # 每個面各自有頂點 -> 合併
    tm.merge_vertices(digits_vertex=4)
    m = Manifold(Mesh(np.asarray(tm.vertices, np.float32), np.asarray(tm.faces, np.uint32)))
    # 無效網格做布林會靜靜回傳空結果,差集變成 0,所以一定要先確認
    assert not m.is_empty(), f'{name}:細分後的網格不是封閉實體({m.status()})'
    diff = (m - ref).volume() + (ref - m).volume()
    # 最大表面偏差:兩邊表面互相取樣量距離。多邊形近似圓的誤差約 0.01 mm,
    # 少一個孔或尺寸不同會是 mm 等級(只看體積的話,少一個 Ø2.3 小孔才差 0.15%,容易漏掉)
    rg = ref.to_mesh()
    tr = trimesh.Trimesh(np.asarray(rg.vert_properties)[:, :3], np.asarray(rg.tri_verts))
    dev = 0.0
    for a, b in ((tm, tr), (tr, tm)):
        pts = np.vstack([a.vertices, trimesh.sample.sample_surface(a, 20000, seed=0)[0]])
        dev = max(dev, max(trimesh.proximity.closest_point(b, pts[i:i + 5000])[1].max()
                           for i in range(0, len(pts), 5000)))
    print(f'  {name:15} STEP 體積 {shape.Volume() / 1000:6.2f} cm³  STL 版 {ref.volume() / 1000:6.2f} cm³  '
          f'差集 {diff:5.1f} mm³  最大表面偏差 {dev:.3f} mm')
    return dev


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out-dir', default='output/step')
    ap.add_argument('--no-check', action='store_true')
    args = ap.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    worst = 0.0     # 最大表面偏差(mm)
    for name, make in PARTS.items():
        shape = make()
        assert shape.isValid(), f'{name} 不是有效實體'
        assert len(shape.Solids()) == 1, f'{name} 有 {len(shape.Solids())} 個實體'
        path = out / f'{name}.step'
        asm = cq.Assembly(name='gimbal_print').add(shape, name=name)   # 讓 CAD 裡顯示零件名稱
        asm.export(str(path))
        print(f'輸出 {path}  ({path.stat().st_size / 1e3:.0f} kB, {len(shape.Faces())} 個面)')
        if not args.no_check:
            worst = max(worst, compare(name, shape))
    if not args.no_check and worst > 0.05:
        sys.exit(f'STEP 與 manifold3d 版本的表面偏差 {worst:.2f} mm 超過 0.05 mm,請檢查兩邊幾何是否同步')


if __name__ == '__main__':
    main()
