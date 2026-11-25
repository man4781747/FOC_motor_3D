# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "trimesh", "scipy", "networkx", "manifold3d", "pillow"]
# ///
"""把 2804 馬達用 6.5 mm 銅柱裝到 PCB 背面,出線焊到 U/V/W 焊墊,輸出成一個 GLB。

用法:
    uv run scripts/build_assembly.py <pcb.glb> <貼圖資料夾> <輸出.glb>

座標(mm,Z 朝上):原點 = PCB 中心、PCB 背面高度;正面在 z = 2.54。
- PCB 的 4 個 NPTH 安裝孔 (0, ±9.5)、(±9.5, 0) 正好是 M2.5、孔位圓 Ø19,
  與馬達底部螺孔相同 -> 馬達翻轉朝下並轉 45° 對孔
- 銅柱 z = -6.5 ~ 0,正面以 M2.5 盤頭螺絲鎖住
- U/V/W 焊墊位置由背面貼圖(Gerber 產生)自動找出
"""
import json
import struct
import sys
from pathlib import Path

import numpy as np
import trimesh
from manifold3d import Manifold
from PIL import Image
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).parent))
from build_glb import GLB, SMOOTH_ANGLE, TEX_HALF_MM, gltf_arrays, srgb_to_linear  # noqa: E402
from build_motor import MATERIALS, capsule_path, cyl, motor_parts, polar, union  # noqa: E402

BOARD_T = 2.54            # PCB 厚度(與 OBJ 相同)
STANDOFF = 6.5
MOUNT_R = 9.5             # 安裝孔孔位圓半徑
MOTOR_ROT = 45            # 馬達繞 Z 旋轉角度(翻轉後),讓螺孔對到 0/90/180/270°
EXIT_ANGLE = 90 + MOTOR_ROT   # 翻轉後出線口方向
WIRE_R = 0.45


def load_glb(path):
    """讀入既有 GLB,回傳可繼續加網格的 GLB 物件。"""
    data = Path(path).read_bytes()
    _, _, total = struct.unpack_from('<4sII', data, 0)
    off, g = 12, GLB()
    while off < total:
        length, kind = struct.unpack_from('<I4s', data, off)
        chunk = data[off + 8: off + 8 + length]
        if kind == b'JSON':
            g.j = json.loads(chunk)
        elif kind == b'BIN\0':
            g.bin = bytearray(chunk)
        off += 8 + length
    for k in ('textures', 'images', 'samplers'):
        g.j.setdefault(k, [])
    return g


def find_uvw_pads(tex_path):
    """從背面貼圖找出左上角 U/V/W 三個大焊墊的中心(mm)。"""
    img = np.asarray(Image.open(tex_path).convert('RGB')).astype(int)
    size = img.shape[0]
    gold = (np.abs(img - np.array([214, 170, 70])).sum(-1) < 60)
    # 金色焊墊上有 via 孔,先填補再分群
    gold = ndimage.binary_closing(gold, iterations=6)
    lab, n = ndimage.label(gold)
    mm = 2 * TEX_HALF_MM / size
    pads = []
    for i in range(1, n + 1):
        ys, xs = np.nonzero(lab == i)
        area = len(xs) * mm * mm
        x, y = xs.mean() * mm - TEX_HALF_MM, TEX_HALF_MM - ys.mean() * mm
        if area > 3.0 and x < -5 and y > 8:          # 大焊墊、位於左上角
            pads.append((x, y, area))
    assert len(pads) == 3, f'找到 {len(pads)} 個 U/V/W 焊墊,預期 3 個'
    return pads


def place_motor(m):
    """馬達座標 -> 板子座標:繞 X 翻轉 180°、繞 Z 轉 45°、往下移 6.5 mm。"""
    return m.rotate([180, 0, 0]).rotate([0, 0, MOTOR_ROT]).translate([0, 0, -STANDOFF])


def pan_head_top(x, y, d=4.5, h=1.5):
    from build_motor import revolve, box
    r = d / 2
    head = revolve([(0, 0), (r, 0), (r, h * 0.55), (r * 0.8, h * 0.9), (r * 0.5, h), (0, h)])
    s = 0.5
    cross = union([box(-r * 0.55, r * 0.55, -s / 2, s / 2, h * 0.45, h + 1),
                   box(-s / 2, s / 2, -r * 0.55, r * 0.55, h * 0.45, h + 1)])
    return (head - cross).translate([x, y, BOARD_T])


def hexagon(flats, h, z0, x, y):
    r = flats / np.sqrt(3)
    return Manifold.cylinder(h, r, circular_segments=6).rotate([0, 0, 30]).translate([x, y, z0])


def wires_to_pads(pads):
    """三條線從出線口往外、往上拉到板子背面,沿焊墊徑向鋪平後上錫。"""
    parts = {}
    z_exit = -STANDOFF - 1.2
    # 依角度排序,出線順序與焊墊順序一致才不會交叉
    pads = sorted(pads, key=lambda p: np.degrees(np.arctan2(p[1], p[0])))
    names = ['wire_white', 'wire_red', 'wire_black']   # 角度由小到大
    offsets = [-1.1, 0.0, 1.1]                          # 切向偏移(逆時針為正)
    a = np.radians(EXIT_ANGLE)
    radial, tangent = np.array([np.cos(a), np.sin(a)]), np.array([-np.sin(a), np.cos(a)])
    cores, solder = [], []
    for (px, py, _), name, t in zip(pads, names, offsets):
        phi = np.degrees(np.arctan2(py, px))
        rp = np.hypot(px, py)
        p0 = radial * 16.9 + tangent * t
        p1 = radial * 18.4 + tangent * t
        pts = [(*p0, z_exit), (*p1, z_exit),
               (*polar(19.0, phi), -3.5), (*polar(19.1, phi), -1.2),
               (*polar(rp + 1.0, phi), -WIRE_R - 0.02)]
        parts[name] = capsule_path(pts, WIRE_R)
        # 剝線後的銅芯 + 焊錫
        cores.append(capsule_path([(*polar(rp + 1.0, phi), -0.32), (*polar(rp - 1.0, phi), -0.32)], 0.3))
        blob = Manifold.sphere(1.0, 32).scale([1.5, 0.95, 0.45]).rotate([0, 0, phi])
        solder.append(blob.translate([*polar(rp - 0.2, phi), -0.1]))
    parts['wire_core'] = union(cores)
    parts['solder'] = union(solder)
    return parts


EXTRA_MATERIALS = {
    'brass': ((0.80, 0.62, 0.32), 1.0, 0.3),
    'screws': ((0.80, 0.80, 0.82), 1.0, 0.25),
    'wire_core': ((0.85, 0.50, 0.30), 1.0, 0.3),
    'solder': ((0.78, 0.78, 0.80), 1.0, 0.18),
}


def to_trimesh(m):
    mg = m.to_mesh()
    v = np.asarray(mg.vert_properties)[:, :3] / 1000.0
    return trimesh.Trimesh(v, np.asarray(mg.tri_verts), process=True)


def main():
    pcb_glb, tex_dir, out = sys.argv[1:4]
    g = load_glb(pcb_glb)
    g.j['nodes'][0]['name'] = 'PCB_with_Motor_2804'

    pads = find_uvw_pads(Path(tex_dir) / 'bottom_basecolor.png')
    for x, y, a in pads:
        print(f'U/V/W 焊墊 ({x:6.2f}, {y:6.2f}) mm,面積 {a:.1f} mm²,角度 {np.degrees(np.arctan2(y, x)):.1f}°')

    holes = [polar(MOUNT_R, a) for a in (0, 90, 180, 270)]
    parts = {k: place_motor(v) for k, v in motor_parts().items()}
    parts['brass'] = union([hexagon(4.0, STANDOFF, -STANDOFF, x, y) for x, y in holes])
    parts['screws'] = union([pan_head_top(x, y) for x, y in holes])
    parts.update(wires_to_pads(pads))

    mats = {**MATERIALS, **EXTRA_MATERIALS}
    for name, man in parts.items():
        mesh = to_trimesh(man)
        srgb, metal, rough = mats[name]
        mat = g.material(dict(name=name, pbrMetallicRoughness=dict(
            baseColorFactor=[*srgb_to_linear(srgb).tolist(), 1.0],
            metallicFactor=metal, roughnessFactor=rough)))
        pos, nrm, idx = gltf_arrays(trimesh.graph.smooth_shade(mesh, angle=SMOOTH_ANGLE))
        g.mesh(name, pos, nrm, idx, mat)
        print(f'{name:13} 三角形 {len(mesh.faces):6}')

    g.save(out)
    print('輸出', out, f'{Path(out).stat().st_size / 1e6:.1f} MB')


if __name__ == '__main__':
    main()
