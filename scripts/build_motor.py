# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "trimesh", "scipy", "networkx", "manifold3d"]
# ///
"""建立 2804 雲台無刷馬達本體的 3D 模型(GLB)。

用法:
    uv run scripts/build_motor.py <輸出.glb>

尺寸依賣家圖片(單位 mm,原點 = 馬達底面中心,Z 朝上):
    外徑 Ø34.5、高 15(上方轉子 + 下方定子座,中間一圈細縫)
    中空軸孔 Ø6.5(內有軸承環)
    頂部 4 x M2.5 螺孔,孔位圓 Ø16;另有 2 個小孔
    底部 4 x M2.5 螺孔,孔位圓 Ø19
    出線 3 條(黑/紅/白)約 30 mm,末端 1.0 mm 間距 3P 接頭
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import trimesh
from manifold3d import CrossSection, Manifold, OpType

sys.path.insert(0, str(Path(__file__).parent))
from build_glb import GLB, SMOOTH_ANGLE, gltf_arrays, srgb_to_linear  # noqa: E402

SEG = 128                 # 圓的分段數

# ---------------- 尺寸(mm) ----------------
MOTOR_R = 34.5 / 2
HEIGHT = 15.0
BASE_H = 2.5              # 定子座高度(細縫以下)
GAP = 0.25                # 轉子與定子座之間的細縫
CHAMFER = 0.8             # 轉子頂緣倒角
SHAFT_R = 6.5 / 2
BOLT_TOP_R = 16 / 2
BOLT_BOTTOM_R = 19 / 2
M25_R = 2.5 / 2
DIAG = [45, 135, 225, 315]


def cyl(r, h, z0=0.0, x=0.0, y=0.0, seg=SEG):
    return Manifold.cylinder(h, r, circular_segments=seg).translate([x, y, z0])


def box(x0, x1, y0, y1, z0, z1):
    return Manifold.cube([x1 - x0, y1 - y0, z1 - z0]).translate([x0, y0, z0])


def revolve(profile):
    """profile: [(r, z), ...] 逆時針多邊形,繞 Z 軸旋轉。"""
    return CrossSection([profile]).revolve(SEG)


def polar(r, deg):
    a = np.radians(deg)
    return r * np.cos(a), r * np.sin(a)


def union(parts):
    return Manifold.batch_boolean(parts, OpType.Add)


def capsule_path(points, r):
    """沿折線產生圓管(每段是兩顆球的凸包)。"""
    balls = [Manifold.sphere(r, 16).translate(p) for p in points]
    return union([Manifold.batch_hull([a, b]) for a, b in zip(balls[:-1], balls[1:])])


def motor_parts():
    p = {}
    zg = BASE_H                     # 細縫高度
    # 定子座:底緣小倒角,細縫處內縮 0.3 mm
    stator = revolve([(SHAFT_R, 0), (MOTOR_R - 0.5, 0), (MOTOR_R, 0.5), (MOTOR_R, zg),
                      (MOTOR_R - 0.3, zg), (MOTOR_R - 0.3, zg + GAP), (SHAFT_R, zg + GAP)])
    # 底部:4 x M2.5 螺孔(Ø19)+ 中心一圈淺凹槽
    bottom_holes = [cyl(M25_R, 2.0, -0.5, *polar(BOLT_BOTTOM_R, a)) for a in DIAG]
    bottom_holes.append(revolve([(SHAFT_R, -1), (5.5, -1), (5.5, 0.4), (SHAFT_R, 0.4)]))
    stator = stator - union(bottom_holes)
    # 出線口護套(-Y 方向)
    p['motor_stator'] = stator + box(-2.2, 2.2, -MOTOR_R - 1.2, -MOTOR_R + 1.0, 0.4, 2.0)

    # 轉子:頂緣倒角
    zt = HEIGHT
    rotor = revolve([(SHAFT_R, zg + GAP), (MOTOR_R, zg + GAP), (MOTOR_R, zt - CHAMFER),
                     (MOTOR_R - CHAMFER, zt), (SHAFT_R, zt)])
    holes = [cyl(M25_R, 3.0, zt - 2.5, *polar(BOLT_TOP_R, a)) for a in DIAG]
    holes += [cyl(0.75, 2.0, zt - 1.5, *polar(12.0, a)) for a in (15, 195)]
    holes.append(revolve([(10.5, zt - 0.25), (11.0, zt - 0.25), (11.0, zt + 1), (10.5, zt + 1)]))
    p['motor_rotor'] = rotor - union(holes)

    # 中空軸 / 軸承內環(頂面略凹 0.2 mm,底部凸出 0.4 mm)
    p['shaft'] = revolve([(2.4, -0.4), (SHAFT_R, -0.4), (SHAFT_R, zt - 0.2), (2.4, zt - 0.2)])
    return p


def wire_parts():
    """3 條出線從定子側面(-Y)拉出,自然下垂到 3P 接頭。"""
    p = {}
    z = 1.2
    for name, x in (('wire_black', -1.1), ('wire_red', 0.0), ('wire_white', 1.1)):
        pts = [(x, -MOTOR_R + 0.5, z), (x, -21.5, z), (x, -25.0, z - 2.0),
               (x * 0.9, -28.0, z - 8.0), (x * 0.9, -30.0, z - 15.0)]
        p[name] = capsule_path(pts, 0.45)
    zc = z - 15.0
    body = box(-2.5, 2.5, -31.5, -28.5, zc - 4.2, zc + 0.2)
    p['connector'] = body - box(-1.9, 1.9, -30.9, -29.1, zc - 4.6, zc - 1.0)
    return p


MATERIALS = {   # sRGB 顏色, 金屬度, 粗糙度
    'motor_rotor': ((0.05, 0.05, 0.055), 0.6, 0.42),   # 黑色陽極處理鋁
    'motor_stator': ((0.06, 0.06, 0.065), 0.6, 0.45),
    'shaft': ((0.82, 0.82, 0.84), 1.0, 0.22),
    'wire_black': ((0.04, 0.04, 0.04), 0.0, 0.5),
    'wire_red': ((0.75, 0.06, 0.05), 0.0, 0.5),
    'wire_white': ((0.92, 0.92, 0.90), 0.0, 0.5),
    'connector': ((0.93, 0.91, 0.86), 0.0, 0.6),
}


def to_trimesh(m):
    mg = m.to_mesh()
    v = np.asarray(mg.vert_properties)[:, :3] / 1000.0   # mm -> m
    return trimesh.Trimesh(v, np.asarray(mg.tri_verts), process=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    out = Path(ap.parse_args().out)

    g = GLB()
    g.j['nodes'][0]['name'] = 'Motor_2804'
    parts = {**motor_parts(), **wire_parts()}
    for name, man in parts.items():
        mesh = to_trimesh(man)
        srgb, metal, rough = MATERIALS[name]
        mat = g.material(dict(name=name, pbrMetallicRoughness=dict(
            baseColorFactor=[*srgb_to_linear(srgb).tolist(), 1.0],
            metallicFactor=metal, roughnessFactor=rough)))
        pos, nrm, idx = gltf_arrays(trimesh.graph.smooth_shade(mesh, angle=SMOOTH_ANGLE))
        g.mesh(name, pos, nrm, idx, mat)
        print(f'{name:13} 三角形 {len(mesh.faces):6}')

    g.save(out)
    print('輸出', out, f'{out.stat().st_size / 1e6:.2f} MB')


if __name__ == '__main__':
    main()
