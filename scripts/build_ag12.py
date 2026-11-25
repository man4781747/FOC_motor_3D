# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "trimesh", "scipy", "networkx", "manifold3d", "pillow"]
# ///
"""KSS AG-12 尼龍電纜固定頭(M12 × 1.25)。

用法:
    uv run scripts/build_ag12.py <輸出.glb> [--panel 3.0] [--stl-dir output/ag12]

尺寸來源:
- 官網 / 型錄 0604:螺紋 M12 × 1.25、螺紋長 L1 = 7.9、開孔 Ø12 +0.2、適用線徑 3~6.5
- 型錄 0611:螺帽 AGL-12 對邊 17.0、高 5.2;墊片 AGW-12 外徑 17.0(PVC)
- 實測:本體六角對邊 15.2、迫緊頭外徑 15.5、總高 30(螺紋末端到迫緊頭頂,不含螺帽)
- 其餘依產品照片比例估計:本體凸緣 Ø17 × 1、六角高 5、迫緊頭頂端圓角與 Ø8.5 出線孔、
  夾緊爪 / 夾緊圈(藏在迫緊頭內)、墊片厚 1.0、內孔 Ø7

零件(型錄爆炸圖,由下往上):螺帽 AGL、墊片 AGW、本體 AGB(與夾緊爪一體)、夾緊圈 AGS(TPE)、迫緊頭 AGN。
座標(mm,Z 朝上):面板外側面 z=0,固定頭在 +Z,螺紋穿過面板往 -Z;螺帽鎖在面板內側(厚度 --panel)。
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import trimesh
from manifold3d import CrossSection, Manifold

sys.path.insert(0, str(Path(__file__).parent))
from build_glb import GLB, SMOOTH_ANGLE, gltf_arrays, srgb_to_linear  # noqa: E402
from build_gimbal import make_extras  # noqa: E402
from hardware import revolve, shade, thread  # noqa: E402

SEG = 96
THREAD_D, PITCH, L1 = 12.0, 1.25, 7.9
TOTAL_H = 30.0                      # 螺紋末端到迫緊頭頂
BORE = 7.0                          # 內孔(線徑 3~6.5)
WASHER_OD, WASHER_ID, WASHER_T = 17.0, 12.2, 1.0   # 厚度未公開,估 1.0
FLANGE_D, FLANGE_H = 17.0, 1.0
BODY_HEX, BODY_HEX_H = 15.2, 5.0
CAP_D, CAP_TOP_HOLE = 15.5, 8.5
THREAD_END = WASHER_T - L1                           # 螺紋從本體凸緣底面(墊片上緣)往下 L1
CAP_Z0 = WASHER_T + FLANGE_H + BODY_HEX_H + 0.8      # 迫緊頭下緣(與本體六角留 0.8 縫)
CAP_TOP = THREAD_END + TOTAL_H                       # 總高從螺紋末端量
CLAW_D = 11.5                                        # 夾緊爪(迫緊頭內的外螺紋段)
SEAL_OD, SEAL_ID, SEAL_H = 9.5, 6.5, 7.0             # 夾緊圈 TPE
NUT_AF, NUT_H = 17.0, 5.2                            # AGL-12
FLUTES = 30

PARTS = {
    'body': ('本體 AGB', '固定頭', 1, 'NYLON 66 黑,M12 × 1.25、L1 7.9、六角 15.2;與夾緊爪一體'),
    'cap': ('迫緊頭 AGN', '固定頭', 1, 'NYLON 66 黑,外徑 15.5'),
    'seal': ('夾緊圈 AGS', '固定頭', 1, 'TPE,線徑 3~6.5'),
    'washer': ('墊片 AGW-12', '固定頭', 1, 'PVC 黑,外徑 17'),
    'nut': ('螺帽 AGL-12', '固定頭', 1, 'NYLON 66 黑,對邊 17、高 5.2,鎖在面板內側'),
}
EXPLODE = {'body': (0, 0, 0), 'washer': (0, 0, -6), 'nut': (0, 0, -22), 'seal': (0, 0, 22), 'cap': (0, 0, 34)}
MATS = {   # sRGB, 金屬度, 粗糙度
    'nylon': ((0.06, 0.06, 0.065), 0.0, 0.45),
    'tpe': ((0.14, 0.14, 0.15), 0.0, 0.85),
    'pvc': ((0.03, 0.03, 0.03), 0.0, 0.7),
}


def cyl(r, z0, z1, r1=None):
    return Manifold.cylinder(z1 - z0, r, r if r1 is None else r1, SEG).translate([0, 0, z0])


def hex_prism(af, z0, z1, chamfer=0.6):
    R = af / np.sqrt(3)
    h = CrossSection([[(R * np.cos(np.radians(a)), R * np.sin(np.radians(a))) for a in range(0, 360, 60)]])
    m = h.extrude(z1 - z0).translate([0, 0, z0])
    # 兩端 45° 倒角:與兩端削角的圓柱取交集(六角角落被削掉)
    return m ^ revolve([(0, z0), (R - chamfer, z0), (R + 0.01, z0 + chamfer), (R + 0.01, z1 - chamfer),
                        (R - chamfer, z1), (0, z1)])


def make_body():
    z_hex0 = WASHER_T + FLANGE_H
    z_hex1 = z_hex0 + BODY_HEX_H
    thr = thread(THREAD_D, L1, PITCH).mirror([0, 0, 1]).translate([0, 0, WASHER_T])   # 螺紋到墊片下緣
    flange = cyl(FLANGE_D / 2, WASHER_T, z_hex0)
    hexa = hex_prism(BODY_HEX, z_hex0, z_hex1)
    claws = cyl(CLAW_D / 2, z_hex1, CAP_TOP - 2.5, CLAW_D / 2 - 1.2)          # 夾緊爪(往上收)
    slots = [Manifold.cube([0.9, CLAW_D, 6.0]).translate([-0.45, -CLAW_D / 2, CAP_TOP - 8.5])
             .rotate([0, 0, a]) for a in range(0, 180, 30)]                    # 爪與爪之間的縫
    body = thr + flange + hexa + claws
    for s in slots:
        body = body - s
    seal_z0 = CAP_TOP - 3.5 - SEAL_H
    return body - cyl(BORE / 2, -L1 - 5, CAP_TOP + 1) - cyl(SEAL_OD / 2 + 0.05, seal_z0, CAP_TOP + 1)   # 夾緊圈裝在爪內


def make_cap():
    r, z0, z1 = CAP_D / 2, CAP_Z0, CAP_TOP
    fillet = 3.0
    prof = [(CLAW_D / 2 + 0.05, z0), (r, z0), (r, z1 - fillet)]
    prof += [(r - fillet + fillet * np.cos(t), z1 - fillet + fillet * np.sin(t)) for t in np.linspace(0, np.pi / 2, 9)[1:]]
    prof += [(CAP_TOP_HOLE / 2, z1), (CAP_TOP_HOLE / 2, z1 - 2.0), (CLAW_D / 2 - 1.2 + 0.05, z1 - 3.5),
             (CLAW_D / 2 + 0.05, z1 - 6.0)]
    cap = revolve(prof)
    for k in range(FLUTES):                                                     # 外圈直紋
        a = 360 * k / FLUTES
        cap = cap - Manifold.cube([0.6, 1.2, z1 - fillet - z0 - 1.0]).translate([r - 0.45, -0.3, z0 + 0.5]) \
            .rotate([0, 0, a])
    return cap


def make_seal():
    z0 = CAP_TOP - 3.5 - SEAL_H
    return cyl(SEAL_OD / 2, z0, z0 + SEAL_H) - cyl(SEAL_ID / 2, z0 - 1, z0 + SEAL_H + 1)


def make_washer():
    return cyl(WASHER_OD / 2, 0, WASHER_T) - cyl(WASHER_ID / 2, -1, WASHER_T + 1)


def make_nut(panel):
    z1 = -panel
    nut = hex_prism(NUT_AF, z1 - NUT_H, z1 - 1.0) + cyl(NUT_AF / 2 + 0.5, z1 - 1.0, z1)   # 貼面板那側有圓墊面
    return nut - cyl(THREAD_D / 2 - 0.1, z1 - NUT_H - 1, z1 + 1)


def to_tm(m):
    mg = m.to_mesh()
    return trimesh.Trimesh(np.asarray(mg.vert_properties)[:, :3], np.asarray(mg.tri_verts), process=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('--panel', type=float, default=3.0, help='面板厚度(螺帽位置)')
    ap.add_argument('--stl-dir', help='另存各零件 STL')
    args = ap.parse_args()

    parts = {'body': (make_body(), 'nylon'), 'cap': (make_cap(), 'nylon'), 'seal': (make_seal(), 'tpe'),
             'washer': (make_washer(), 'pvc'), 'nut': (make_nut(args.panel), 'nylon')}
    for k, (m, _) in parts.items():
        assert m.status().name == 'NoError' and m.volume() > 0, f'{k} 不是有效實體'
    all_m = sum((m for m, _ in parts.values()), Manifold())
    bb = to_tm(all_m).bounds
    thread_end = THREAD_END
    print(f'AG-12:螺紋 M{THREAD_D:.0f} × {PITCH}、L1 {L1}(z {thread_end:.1f}~{WASHER_T:.1f});'
          f'迫緊頭頂 z={CAP_TOP:.1f} -> 總高(螺紋末端到頂){CAP_TOP - thread_end:.1f} mm')
    print(f'  本體六角 {BODY_HEX}、迫緊頭 Ø{CAP_D}、整體外形 {np.ptp(bb[:, 0]):.1f} × {np.ptp(bb[:, 1]):.1f} × {np.ptp(bb[:, 2]):.1f}')
    engage = -args.panel - THREAD_END                     # 面板內側以下的螺紋長
    print(f'  面板 {args.panel} mm:面板內側以下剩 {engage:.1f} mm 螺紋,螺帽高 {NUT_H}'
          f' -> {"全部咬合" if engage >= NUT_H else f"咬合 {engage:.1f} mm,螺帽超出螺紋末端 {NUT_H - engage:.1f} mm"}')

    g = GLB()
    g.j['nodes'][0]['name'] = 'AG-12'
    mats = {}
    for k, (m, key) in parts.items():
        if key not in mats:
            srgb, metal, rough = MATS[key]
            mats[key] = g.material(dict(name=key, pbrMetallicRoughness=dict(
                baseColorFactor=[*srgb_to_linear(srgb).tolist(), 1.0], metallicFactor=metal, roughnessFactor=rough)))
        tm = to_tm(m)
        if args.stl_dir:
            Path(args.stl_dir).mkdir(parents=True, exist_ok=True)
            tm.export(Path(args.stl_dir) / f'ag12_{k}.stl')
        tm.apply_scale(0.001)
        pos, nrm, idx = gltf_arrays(shade(tm, key if k != 'body' else 'screw', SMOOTH_ANGLE))
        g.mesh(k, pos, nrm, idx, mats[key])
        g.j['nodes'][-1]['extras'] = make_extras(PARTS, EXPLODE, k)
    g.save(args.out)
    print('輸出', args.out, f'{Path(args.out).stat().st_size / 1e6:.2f} MB')


if __name__ == '__main__':
    main()
