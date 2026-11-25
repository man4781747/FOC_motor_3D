# /// script
# requires-python = ">=3.10,<3.13"
# dependencies = ["numpy", "trimesh", "scipy", "networkx", "manifold3d", "pillow", "rtree", "lxml", "cadquery"]
# ///
"""用 8 組「PCB + 2804 馬達」模組做 4 輪舵輪(swerve)平台車:
每個角 1 顆驅動(輪轂直驅)+ 1 顆轉向。輸出組合模型 GLB,以及列印件的 STL(已擺好列印方向)與 STEP。

用法:
    uv run scripts/build_car.py <pcb_motor_assembly.glb> <輸出.glb> [--out-dir output/car]

世界座標(mm,Z 朝上,地面 z=0),車頭朝 +X。每個角的「舵輪座標」:原點在轉向軸與地面交點,
輪軸沿 +X(輪轂朝外),再依左右側轉 ±90° 放到車身四角。

結構(由下往上,每個角):
1. 車輪:白色輪轂(杯狀,包住會轉的馬達轉子)+ 黑色 TPU 輪胎,Ø64 × 14。
   輪轂用 4 顆 M2.5 鎖在驅動馬達轉子面(Ø16),轉向軸通過輪子中心面,原地轉向時輪胎不拖地。
2. 驅動模組:PCB 正面朝轉向叉立板,以 4 支 M2.5 公母銅柱 8 mm 連接(同雲台俯仰模組)。
3. 轉向叉(L 形,列印件):上方法蘭用 4 顆 M2.5 從下方鎖進轉向馬達轉子面。
4. 轉向模組:轉子面朝下穿過底板 Ø40 孔;PCB 正面朝上,由「橋架」(列印件)的 4 根凸柱頂住,
   M2.5 螺絲從橋架上方鎖進模組銅柱(取代模組正面原本的短螺絲)。
5. 車身:黑色底板 + 白色上蓋。橋架腳與上蓋螺絲柱都埋 M3 熱熔螺母,M3 × 8 沉頭螺絲從底板下方鎖入;
   上蓋四角各有 4 個孔,對準橋架螺絲,不用拆蓋就能鎖。上蓋頂部中央有 4 個 M3 熱熔螺母,
   可直接鎖雙軸雲台底座;電池放車尾。

五金規則:鎖進馬達螺孔(本身有牙)與 PCB 模組銅柱的用 M2.5;其他鎖進列印件的地方一律用
M3 埋頭型熱熔螺母(M3 × 5、外徑 4.5,列印件開 Ø4.0 × 深 6)。

節點階層(three.js / Blender 可轉動):
    Car_swerve
      └ body                    底板、上蓋、電池
          └ corner_FL/FR/RL/RR  橋架、轉向模組(轉子以外)
              └ steer_*         繞垂直軸:轉向轉子、轉向叉、驅動模組(轉子以外)
                  └ wheel_*     繞輪軸:驅動轉子、輪轂、輪胎
    extras:steer_*.spin_steer = 原地旋轉時的轉向角;wheel_*.spin_sign = 前進時車輪正轉為 +1 / 反轉為 -1;
            part = 零件清單資料(PARTS);explode = 爆炸圖位移(EXPLODE,glTF 座標、公尺)
"""
import argparse
import sys
from pathlib import Path

import cadquery as cq
import numpy as np
import trimesh

sys.path.insert(0, str(Path(__file__).parent))
import build_print_step as q  # noqa: E402  CadQuery 小工具
from build_assembly import load_glb  # noqa: E402
from build_gimbal import (M3, M25, MODULE_PCB_TOP, MODULE_ROTOR_Z, MOUNT, R_TILT, ROTOR_HOLES,  # noqa: E402
                          MODULE_EXPLODE, MODULE_ROTOR, along_x, gltf_matrix, make_extras, module_groups,
                          rotor_explode)
from build_glb import SMOOTH_ANGLE, gltf_arrays, srgb_to_linear  # noqa: E402
from build_motor import polar, union  # noqa: E402
from hardware import INSERT_M3, flat_screw, heat_insert, pan_screw, shade, standoff  # noqa: E402

# ---------------- 舵輪尺寸(舵輪座標,mm) ----------------
WHEEL_R = 32.0                       # 輪胎外半徑
AXLE_Z = WHEEL_R                     # 輪軸高度
RIM_R0, RIM_R1, LIP_R = 23.5, 26.5, 28.0   # 輪轂杯內 / 外半徑、擋邊半徑
XR = 4.0                             # 驅動轉子面的 X;輪子 x ∈ [-7, 7],中心面在轉向軸上
R_D = R_TILT                         # 驅動模組:模組 -Z(馬達)-> +X
T_D = np.array([XR, 0, AXLE_Z]) - R_D @ np.array([0, 0, MODULE_ROTOR_Z])
SPACER = 8.0                         # 驅動模組 PCB 到叉立板的公母銅柱
FORK_X1 = T_D[0] - (MODULE_PCB_TOP + SPACER)   # 叉立板內側面
FORK_X0 = FORK_X1 - 4.0              # 叉立板厚 4
FORK_HALF_W = 22.0
FLANGE_X1 = 22.0
FLANGE_T = 4.0
FLANGE_Z1 = AXLE_Z + WHEEL_R + 3.0 + FLANGE_T  # 轉向轉子面 = 法蘭頂(輪頂留 3 mm)
FLANGE_Z0 = FLANGE_Z1 - FLANGE_T
STEER_T = np.array([0, 0, FLANGE_Z1 - MODULE_ROTOR_Z])   # 轉向模組:不旋轉,馬達朝下
STEER_PCB_FACE = STEER_T[2] + MODULE_PCB_TOP

# ---------------- 車身尺寸(世界座標,mm) ----------------
BODY_L, BODY_W, BODY_R = 250.0, 200.0, 20.0
PLATE_Z0 = FLANGE_Z1 + 2.0           # 底板底面(與轉向叉法蘭留 2 mm)
PLATE_T = 4.0
PLATE_Z1 = PLATE_Z0 + PLATE_T
BOSS_GAP = 6.0                       # 轉向模組 PCB 正面到橋架底面(正面零件最高約 3.5 mm)
BRIDGE_Z0 = STEER_PCB_FACE + BOSS_GAP
BRIDGE_T = 4.0
BRIDGE_Z1 = BRIDGE_Z0 + BRIDGE_T
BRIDGE_R = 32.0                       # 腳外緣 30,留 2 mm(不可剛好相切,細分會破面)
LEG_R, LEG_POS_R = 4.0, 26.0         # 橋架腳(避開 Ø40 PCB)
COVER_TOP_T, WALL = 3.0, 2.4
COVER_Z0 = PLATE_Z1
COVER_Z1 = BRIDGE_Z1 + 1.0 + COVER_TOP_T
CORNER_X, CORNER_Y = 82.0, 57.0      # 轉向軸位置
SWEEP_R = float(np.hypot(FORK_X0, FORK_HALF_W))   # 舵輪旋轉半徑(叉立板外角)
COVER_SCREWS = [(0, 88), (0, -88), (113, 0), (-113, 0)]
BATTERY = (75.0, 35.0, 25.0)         # 3S 1300 mAh 佔位
BATTERY_X = -70.0                    # 電池中心:放車尾,讓出頂部中央雲台孔位的螺絲柱
GIMBAL_HOLES = [polar(27, a) for a in (22.5, 112.5, 202.5, 292.5)]   # 雙軸雲台底座的 4 個 M3 孔

CORNERS = {   # 名稱: (x, y, 偏擺角, 前進時車輪轉向)
    'FL': (CORNER_X, CORNER_Y, 90, 1), 'RL': (-CORNER_X, CORNER_Y, 90, 1),
    'FR': (CORNER_X, -CORNER_Y, -90, -1), 'RR': (-CORNER_X, -CORNER_Y, -90, -1),
}

# 零件清單(檢視頁的「零件拆分」):id -> (名稱, 分類, 每個節點的件數, 說明)
# 列印件的說明是 STL 檔名;分類順序就是清單顯示順序
PARTS = {
    'top_cover': ('上蓋', '列印件', 1, 'top_cover.stl · PLA 白'),
    'bottom_plate': ('底板', '列印件', 1, 'bottom_plate.stl · PLA 黑'),
    'bridge': ('橋架', '列印件', 1, 'bridge.stl · PLA 白'),
    'fork': ('轉向叉', '列印件', 1, 'fork.stl · PLA 白'),
    'hub': ('輪轂', '列印件', 1, 'hub.stl · PLA 白'),
    'tire': ('輪胎', '列印件', 1, 'tire.stl · TPU 黑'),
    'steer_motor': ('轉向馬達', '馬達', 1, '2804 雲台馬達(轉子與定子一體,不可拆);轉子帶動轉向叉'),
    'drive_motor': ('驅動馬達', '馬達', 1, '2804 雲台馬達(轉子與定子一體,不可拆);轉子帶動輪轂'),
    'steer_pcb': ('轉向 PCB', '電路板', 1, 'Ø40 驅動板,正面朝上由橋架凸柱頂住'),
    'drive_pcb': ('驅動 PCB', '電路板', 1, 'Ø40 驅動板,正面朝轉向叉立板'),
    'battery': ('電池', '電子', 1, '3S 1300 mAh 尺寸佔位'),
    'steer_standoffs': ('轉向模組銅柱', '五金', 4, 'M2.5 公母 6.5 mm,公頭鎖進馬達底部'),
    'drive_standoffs': ('驅動模組銅柱', '五金', 4, 'M2.5 公母 6.5 mm,公頭鎖進馬達底部'),
    'bridge_screws': ('橋架螺絲', '五金', 4, 'M2.5 × 14 盤頭,穿過凸柱、PCB 鎖進模組銅柱'),
    'leg_screws': ('橋架腳螺絲', '五金', 4, 'M3 × 8 沉頭,從底板下方鎖進橋架腳的熱熔螺母'),
    'cover_screws': ('上蓋螺絲', '五金', 1, 'M3 × 8 沉頭,從底板下方鎖進上蓋螺絲柱的熱熔螺母'),
    'leg_inserts': ('橋架腳熱熔螺母', '五金', 4, 'M3 × 5 外徑 4.5 埋頭型,腳底開 Ø4.0 × 深 6'),
    'cover_inserts': ('上蓋熱熔螺母', '五金', 1, 'M3 × 5 外徑 4.5 埋頭型,螺絲柱開 Ø4.0 × 深 6'),
    'box_body': ('防水盒 盒身', '外購件', 1, 'RS PRO 1862264,225 × 175 × 81 mm;直接放在上蓋上,未固定'),
    'box_lid': ('防水盒 盒蓋', '外購件', 1, 'RS PRO 1862264 附'),
    'box_gasket': ('防水盒 密封條', '外購件', 1, 'RS PRO 1862264 附'),
    'box_inserts': ('防水盒 螺牙嵌件', '外購件', 6, 'RS PRO 1862264 附,盒蓋螺絲鎖入處'),
    'gimbal_inserts': ('雲台孔位熱熔螺母', '五金', 4, 'M3 × 5 外徑 4.5 埋頭型,頂面開 Ø4.0 × 深 6;鎖雙軸雲台底座'),
    'flange_screws': ('法蘭螺絲', '五金', 4, 'M2.5 × 5 盤頭,從下方鎖進轉向馬達轉子面'),
    'plate_screws': ('立板螺絲', '五金', 4, 'M2.5 × 6 盤頭,鎖進間隔銅柱母牙'),
    'spacers': ('間隔銅柱', '五金', 4, 'M2.5 公母 8 mm,公頭 6 mm 穿過 PCB 鎖進模組銅柱'),
    'hub_screws': ('輪轂螺絲', '五金', 4, 'M2.5 × 5 盤頭,鎖進驅動馬達轉子面'),
}
# 爆炸圖位移(mm,Z 朝上,在父節點座標中):照組裝順序往外拆
EXPLODE = {
    'top_cover': (0, 0, 110), 'battery': (0, 0, 50), 'bottom_plate': (0, 0, 0),
    'box_body': (0, 0, 150), 'box_inserts': (0, 0, 170), 'box_gasket': (0, 0, 190), 'box_lid': (0, 0, 215),
    'cover_screws': (0, 0, -25), 'cover_inserts': (0, 0, 100), 'gimbal_inserts': (0, 0, 125),
    'bridge': (0, 0, 70), 'bridge_screws': (0, 0, 95), 'leg_screws': (0, 0, -25), 'leg_inserts': (0, 0, 62),
    'steer_module': (0, 0, 30),
    'steer': (0, 0, -30),                       # 整組舵輪往下
    'fork': (0, 0, 0),
    'flange_screws': (0, 0, -16), 'plate_screws': (-18, 0, 0), 'spacers': (6, 0, 0),
    'drive_module': (25, 0, 0), 'wheel': (50, 0, 0),   # 沿輪軸(舵輪座標 +X)往外
    'hub': (10, 0, 0), 'tire': (22, 0, 0), 'hub_screws': (34, 0, 0),
}

MATS = {   # sRGB, 金屬度, 粗糙度
    'print_white': ((0.93, 0.93, 0.91), 0.0, 0.6),
    'print_black': ((0.07, 0.07, 0.08), 0.0, 0.7),
    'tire': ((0.04, 0.04, 0.04), 0.0, 0.9),
    'brass': ((0.80, 0.62, 0.32), 1.0, 0.3),
    'screw': ((0.80, 0.80, 0.82), 1.0, 0.25),
    'battery': ((0.12, 0.22, 0.45), 0.0, 0.5),
    'box_gray': ((0.74, 0.75, 0.76), 0.0, 0.55),        # RS PRO 淺灰
    'box_gasket': ((0.12, 0.12, 0.13), 0.0, 0.85),
}


def drive_world(p_module):
    return R_D @ np.asarray(p_module, float) + T_D


def rz(deg):
    a = np.radians(deg)
    return np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])


def corner_points(local_xy):
    """舵輪座標的 XY 點 -> 四個角的世界座標(橋架螺絲、腳位等)。"""
    out = []
    for cx, cy, yaw, _ in CORNERS.values():
        for p in local_xy:
            x, y, _ = rz(yaw) @ np.array([*p, 0.0])
            out.append((cx + x, cy + y))
    return out


def rounded_box(length, width, r, z0, h, cx=0.0, cy=0.0):
    return cq.Workplane('XY', origin=(cx, cy, z0)).rect(length, width).extrude(h).edges('|Z').fillet(r)


def teardrop_x(r, x0, length, y, z, tip=1):
    """沿 X 的水滴形孔;tip=+1 尖端朝 +Z,-1 朝 -Z(配合列印時的「上方」)。"""
    m = q.teardrop_x(r, 0, length, 0, 0)
    if tip < 0:
        m = m.mirror('XY')
    return m.translate(q.V(x0, y, z))


def countersink(x, y, z_bottom, d=3.4, dk=6.4):
    """M3 沉頭孔(從下方鎖入,90° 錐)。"""
    cone = cq.Solid.makeCone(dk / 2, d / 2, (dk - d) / 2, q.V(x, y, z_bottom - 0.01))
    return cone.fuse(q.cyl(d / 2, 30, z_bottom - 1, x, y))


# ---------------- 列印件(CadQuery) ----------------
def make_tire():
    """黑色 TPU 輪胎(驅動模組座標:輪軸 = Z,朝外 = -Z)。內側倒角配合輪轂擋邊。"""
    z0, z1 = MODULE_ROTOR_Z - 1.5, MODULE_ROTOR_Z + 9.5      # -23 ~ -12
    prof = [(RIM_R1, z0), (WHEEL_R - 1, z0), (WHEEL_R, z0 + 1), (WHEEL_R, z1 - 1), (WHEEL_R - 1, z1),
            (LIP_R, z1), (RIM_R1, z1 - (LIP_R - RIM_R1))]
    tire = q.revolve(prof)
    grooves = [q.box(-0.6, 0.6, WHEEL_R - 0.8, WHEEL_R + 1, z0 - 1, z1 + 1).rotate(q.V(0, 0, 0), q.V(0, 0, 1), k * 15)
               for k in range(24)]
    return q.finish(tire, grooves)


def make_hub():
    """白色輪轂(驅動模組座標)。外側面朝下列印;內側擋邊下緣 45°,不需支撐。"""
    zf = MODULE_ROTOR_Z                                           # 轉子面 -21.5
    plate = q.cyl(RIM_R1, 3.0, zf - 3.0)                           # 輪轂盤 -24.5 ~ -21.5
    rim = q.cyl(RIM_R1, 14.0, zf - 3.0).cut(q.cyl(RIM_R0, 20, zf))  # 杯緣到 -10.5(離定子座 1.25)
    outer_lip = q.cyl(LIP_R, 1.5, zf - 3.0)
    inner_lip = q.revolve([(RIM_R0 + 0.1, zf + 7.9), (RIM_R1 - 0.1, zf + 7.9), (LIP_R, zf + 9.5),
                           (LIP_R, zf + 11.0), (RIM_R0 + 0.1, zf + 11.0)])
    cuts = [q.cyl(4.5, 10, zf - 5)]                                # 中空軸
    for x, y in ROTOR_HOLES:
        cuts += [q.cyl(M25 / 2, 10, zf - 5, x, y), q.cyl(5.6 / 2, 2.5, zf - 4.0, x, y)]   # 外側沉頭 1.5
    for k in range(6):                                             # 六角鏤空(輪轂外觀)
        x, y = polar(16.5, 30 + 60 * k)
        cuts.append(cq.Workplane('XY', origin=(x, y, zf - 4)).polygon(6, 9.0).extrude(5).val())
    return q.finish(q.fuse(plate, rim, outer_lip, inner_lip), cuts)


def make_fork():
    """L 形轉向叉(舵輪座標)。列印方向:法蘭頂面朝下,立板孔用尖端朝 -Z 的水滴形。"""
    flange = rounded_box(FLANGE_X1 - FORK_X0, 2 * FORK_HALF_W, 6, FLANGE_Z0, FLANGE_T,
                         (FLANGE_X1 + FORK_X0) / 2).val()
    plate = q.box(FORK_X0, FORK_X1, -FORK_HALF_W, FORK_HALF_W, AXLE_Z, FLANGE_Z1).fuse(
        q.cyl_x(FORK_HALF_W, FORK_X1 - FORK_X0, FORK_X0, 0, AXLE_Z))            # 下緣半圓
    zb = FLANGE_Z0
    fillet = q.prism_y([(FORK_X1 - 0.01, zb + 0.01), (FORK_X1 + 6, zb + 0.01), (FORK_X1 - 0.01, zb - 6)],
                       -FORK_HALF_W, FORK_HALF_W)
    ribs = [q.prism_y([(FORK_X1 - 0.01, zb + 0.01), (FORK_X1 + 18, zb + 0.01), (FORK_X1 - 0.01, zb - 18)],
                      y0, y0 + 3) for y0 in (-FORK_HALF_W, FORK_HALF_W - 3)]
    cuts = [q.cyl(4.5, 10, zb - 1)]                                              # 中空軸 / 走線
    for x, y in ROTOR_HOLES:
        cuts += [q.cyl(M25 / 2, 10, zb - 1, x, y), q.cyl(5.6 / 2, 3.0, zb - 1, x, y)]   # 下方沉頭 2
    cuts.append(teardrop_x(5.0, FORK_X0 - 1, 6, 0, AXLE_Z, tip=-1))              # 驅動模組走線孔
    for m in MOUNT:
        _, wy, wz = drive_world([*m, 0])
        cuts.append(teardrop_x(M25 / 2, FORK_X0 - 1, 6, wy, wz, tip=-1))
        cuts.append(teardrop_x(5.6 / 2, FORK_X0 - 1, 2.5, wy, wz, tip=-1))       # 外側沉頭 1.5
    return q.finish(q.fuse(flange, plate, fillet, *ribs), cuts)


def make_bridge():
    """轉向模組橋架(舵輪座標):4 根凸柱頂住 PCB 正面、4 隻腳鎖底板。列印方向:頂板朝下。"""
    pf = STEER_PCB_FACE
    top = q.cyl(BRIDGE_R, BRIDGE_T, BRIDGE_Z0)
    legs = [q.cyl(LEG_R, BRIDGE_Z0 - PLATE_Z1 + 0.01, PLATE_Z1, *polar(LEG_POS_R, 45 + 90 * k)) for k in range(4)]
    bosses, cuts = [], [q.cyl(5.0, 10, BRIDGE_Z0 - 1)]                            # 走線孔
    for x, y in MOUNT:
        # 與雲台底座凸柱相同:碰 PCB 的一端 Ø4.6(避開 LQFP / 晶振),往上 45° 擴大到 Ø7
        bosses.append(q.revolve([(0, pf), (2.3, pf), (2.3, pf + 2.0), (3.5, pf + 3.2),
                                 (3.5, BRIDGE_Z0 + 0.01), (0, BRIDGE_Z0 + 0.01)], x, y))
        cuts += [q.cyl(M25 / 2, 20, pf - 1, x, y), q.cyl(5.6 / 2, 3.0, BRIDGE_Z1 - 2.0, x, y)]   # 頂面沉頭 2
    for k in range(4):
        cuts.append(q.cyl(INSERT_M3['hole_d'] / 2, INSERT_M3['hole_depth'] + 1, PLATE_Z1 - 1,
                          *polar(LEG_POS_R, 45 + 90 * k)))                          # M3 熱熔螺母孔
    return q.finish(q.fuse(top, *legs, *bosses), cuts)


def make_bottom_plate():
    plate = rounded_box(BODY_L, BODY_W, BODY_R, PLATE_Z0, PLATE_T).edges('<Z').chamfer(0.6).val()
    cuts = [q.cyl(20.0, 10, PLATE_Z0 - 1, cx, cy) for cx, cy, _, _ in CORNERS.values()]   # 轉向轉子穿過
    legs = [polar(LEG_POS_R, 45 + 90 * k) for k in range(4)]
    cuts += [countersink(x, y, PLATE_Z0) for x, y in corner_points(legs) + COVER_SCREWS]
    for sx in (BATTERY_X - 20, BATTERY_X + 20):                                # 電池綁帶孔
        for sy in (-1, 1):
            cuts.append(rounded_box(22, 3, 1.4, PLATE_Z0 - 1, 10, sx, sy * (BATTERY[1] / 2 + 4)).val())
    return q.finish(plate, cuts)


def make_cover():
    """上蓋。列印方向:頂面朝下(頂緣用 45° 倒角,不用圓角,第一層才不會懸空)。"""
    h = COVER_Z1 - COVER_Z0
    shell = rounded_box(BODY_L, BODY_W, BODY_R, COVER_Z0, h).edges('>Z').chamfer(3.0).val()
    shell = shell.cut(rounded_box(BODY_L - 2 * WALL, BODY_W - 2 * WALL, BODY_R - WALL,
                                  COVER_Z0 - 1, h - COVER_TOP_T + 1).val())
    hd, hdep = INSERT_M3['hole_d'] / 2, INSERT_M3['hole_depth']
    under = COVER_Z1 - COVER_TOP_T                                              # 頂板底面
    bosses = [q.cyl(4.0, under - COVER_Z0 + 0.01, COVER_Z0, x, y) for x, y in COVER_SCREWS]
    # 雲台孔位:頂板只有 3 mm,往下加螺絲柱讓 6 mm 深的熱熔螺母孔留 1 mm 底
    bosses += [q.cyl(4.0, under - (COVER_Z1 - hdep - 1) + 0.01, COVER_Z1 - hdep - 1, x, y) for x, y in GIMBAL_HOLES]
    cuts = [q.cyl(hd, hdep + 1, COVER_Z0 - 1, x, y) for x, y in COVER_SCREWS]  # 熱熔螺母孔(從下方)
    cuts += [q.cyl(hd, hdep + 1, COVER_Z1 - hdep, x, y) for x, y in GIMBAL_HOLES]   # 熱熔螺母孔(從頂面)
    cuts += [q.cyl(5.5 / 2, 10, COVER_Z1 - 5, x, y) for x, y in corner_points(MOUNT)]   # 橋架螺絲維修孔
    cuts.append(rounded_box(8, 24, 3.9, COVER_Z1 - 5, 10, -100, 0).val())         # 開關 / 充電口
    return q.finish(q.fuse(shell, *bosses), cuts)


# 名稱: (建立函式, 材質, 列印方向, 數量)
PRINTED = {
    'top_cover': (make_cover, 'print_white', 'flip', 1),
    'bottom_plate': (make_bottom_plate, 'print_black', 'up', 1),
    'bridge': (make_bridge, 'print_white', 'flip', 4),
    'fork': (make_fork, 'print_white', 'flip', 4),
    'hub': (make_hub, 'print_white', 'up', 4),
    'tire': (make_tire, 'tire', 'up', 4),
}


# ---------------- 五金(manifold3d) ----------------
def bridge_screws():
    """M2.5 × 14:頭在橋架頂面沉頭孔,穿過凸柱、PCB,鎖進模組銅柱約 3.5 mm。"""
    return union([pan_screw(2.5, 14, 4.5, 1.5).translate([x, y, BRIDGE_Z1 - 2.0]) for x, y in MOUNT])


def leg_screws():
    """M3 × 8 沉頭:頭頂與底板底面齊平,往上鎖進橋架腳的熱熔螺母(咬合約 4 mm)。"""
    return union([flat_screw(3.0, 8, 6.0).mirror([0, 0, 1]).translate([*polar(LEG_POS_R, 45 + 90 * k), PLATE_Z0])
                  for k in range(4)])


def leg_inserts():
    """熱熔螺母頂面與橋架腳底面齊平,往上埋入。"""
    return union([heat_insert().mirror([0, 0, 1]).translate([*polar(LEG_POS_R, 45 + 90 * k), PLATE_Z1])
                  for k in range(4)])


def cover_screw():
    return flat_screw(3.0, 8, 6.0).mirror([0, 0, 1]).translate([0, 0, PLATE_Z0])


def cover_insert():
    return heat_insert().mirror([0, 0, 1]).translate([0, 0, COVER_Z0])


def gimbal_inserts():
    """上蓋頂面的雲台孔位:熱熔螺母頂面與上蓋頂面齊平。"""
    return union([heat_insert().translate([x, y, COVER_Z1]) for x, y in GIMBAL_HOLES])


def fork_hardware():
    flange = [pan_screw(2.5, 5, 4.5, 1.5).mirror([0, 0, 1]).translate([x, y, FLANGE_Z0 + 2.0])
              for x, y in ROTOR_HOLES]                                     # 頭朝下,螺桿往上鎖進轉子面
    standoffs, plate = [], []
    for m in MOUNT:
        wx, wy, wz = drive_world([*m, MODULE_PCB_TOP])
        # 公母銅柱:母端在立板,本體到 PCB 正面,公頭穿過 PCB 鎖進模組銅柱
        standoffs.append(along_x(standoff(2.5, SPACER, stud=6.0)).mirror([1, 0, 0]).translate([wx - SPACER, wy, wz]))
        plate.append(along_x(pan_screw(2.5, 6, 4.5, 1.5)).mirror([1, 0, 0]).translate([FORK_X0 + 1.5, wy, wz]))
    return {'flange_screws': union(flange), 'plate_screws': union(plate), 'spacers': union(standoffs)}


def hub_screws():
    """驅動模組座標:頭在輪轂外側沉頭孔內(朝 -Z),螺桿鎖進轉子面。"""
    return union([pan_screw(2.5, 5, 4.5, 1.5).mirror([0, 0, 1]).translate([x, y, MODULE_ROTOR_Z - 1.5])
                  for x, y in ROTOR_HOLES])


def battery():
    from build_motor import box
    L, W, H = BATTERY
    return box(BATTERY_X - L / 2, BATTERY_X + L / 2, -W / 2, W / 2, PLATE_Z1, PLATE_Z1 + H)


# ---------------- 輸出 ----------------
# 防水盒(外購件,TraceParts 的 3MF)。模型座標:Y 朝上、長邊沿 Z、底面 y = -1.5
BOX_MESHES = {'1862264_1': ('box_body', 'box_gray'), '1862264_3': ('box_lid', 'box_gray'),
              '1862264_10': ('box_gasket', 'box_gasket')}          # 其餘 6 個是螺牙嵌件


def load_box(path):
    """讀 3MF(可在 zip 裡),轉到車身座標:盒子長邊沿 X(車頭方向)、底面貼上蓋頂面。
    回傳 {零件 id: [trimesh(mm), ...]}。"""
    import io
    import zipfile
    path = Path(path)
    if path.suffix == '.zip':
        with zipfile.ZipFile(path) as z:
            name = next(n for n in z.namelist() if n.endswith('.3mf'))
            data = z.read(name)
    else:
        data = path.read_bytes()
    scene = trimesh.load(io.BytesIO(data), file_type='3mf', force='scene')
    M = np.eye(4)
    M[:3, :3] = [[0, 0, 1], [1, 0, 0], [0, 1, 0]]          # (x, y, z)_盒 -> (z, x, y):長邊 -> 車身 X,Y 朝上 -> Z
    out = {}
    for n in scene.graph.nodes_geometry:
        T, g = scene.graph[n]
        m = scene.geometry[g].copy()
        m.apply_transform(M @ T)
        part = BOX_MESHES.get(n, ('box_inserts', 'brass'))
        out.setdefault(part, []).append(m)
    z0 = min(m.bounds[0][2] for ms in out.values() for m in ms)
    for ms in out.values():
        for m in ms:
            m.apply_translation([0, 0, COVER_Z1 - z0])
    return out


def man_to_trimesh(m):
    """manifold3d -> trimesh(mm)。"""
    mg = m.to_mesh()
    return trimesh.Trimesh(np.asarray(mg.vert_properties)[:, :3], np.asarray(mg.tri_verts), process=True)


def cq_to_trimesh(shape, tol=0.02, ang=0.15):
    v, t = shape.tessellate(tol, ang)
    tm = trimesh.Trimesh(np.array([p.toTuple() for p in v]), np.array(t))
    tm.merge_vertices(digits_vertex=4)
    return tm


def export_print(name, shape, orient, out_dir):
    s = shape.rotate(q.V(0, 0, 0), q.V(1, 0, 0), 180) if orient == 'flip' else shape
    bb = s.BoundingBox()
    s = s.translate(q.V(-bb.center.x, -bb.center.y, -bb.zmin))
    tm = cq_to_trimesh(s, 0.01, 0.1)
    tm.export(out_dir / 'print' / f'{name}.stl')
    cq.Assembly(name='car_print').add(shape, name=name).export(str(out_dir / 'step' / f'{name}.step'))
    n, c, a = tm.face_normals, tm.triangles_center, tm.area_faces
    over = (n[:, 2] < -np.cos(np.radians(44))) & (c[:, 2] > 0.3)   # 比 46° 更平才算懸空(容許 45° 倒角的細分誤差)
    e = tm.extents
    print(f'  {name:13} {e[0]:5.1f} x {e[1]:5.1f} x {e[2]:5.1f} mm  體積 {tm.volume / 1000:6.1f} cm³  '
          f'水密 {tm.is_watertight}  懸空 {a[over].sum():5.1f} mm²')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('src')
    ap.add_argument('out')
    ap.add_argument('--out-dir', default='output/car')
    ap.add_argument('--box', help='防水盒 3MF(或含 3MF 的 zip),放在上蓋頂面')
    args = ap.parse_args()

    assert CORNER_X + SWEEP_R + 3 <= BODY_L / 2 and CORNER_Y + SWEEP_R + 3 <= BODY_W / 2, '舵輪轉動範圍超出車身'

    g = load_glb(args.src)
    nodes = g.j['nodes']
    module_nodes = list(nodes[0]['children'])
    by_name = {nodes[i]['name']: i for i in module_nodes}
    rotor_names = MODULE_ROTOR
    skip = {'screws'}           # 模組正面原本的短螺絲 -> 改由橋架長螺絲 / 公母銅柱取代
    nodes[0]['children'] = []
    mat_cache, mesh_cache = {}, {}

    def material(key):
        if key not in mat_cache:
            srgb, metal, rough = MATS[key]
            mat_cache[key] = g.material(dict(name=key, pbrMetallicRoughness=dict(
                baseColorFactor=[*srgb_to_linear(srgb).tolist(), 1.0], metallicFactor=metal, roughnessFactor=rough)))
        return mat_cache[key]

    def add_node(node, parent):
        nodes.append(node)
        nodes[parent].setdefault('children', []).append(len(nodes) - 1)
        return len(nodes) - 1

    def mesh_index(name, make, key):
        """同一個零件只寫一次網格,四個角共用。make() 回傳 mm 單位的 trimesh。"""
        if name not in mesh_cache:
            tm = make()
            tm.apply_scale(0.001)
            pos, nrm, idx = gltf_arrays(shade(tm, key, SMOOTH_ANGLE))
            g.mesh(name, pos, nrm, idx, material(key))
            nodes.pop()
            nodes[0]['children'].pop()
            mesh_cache[name] = len(g.j['meshes']) - 1
        return mesh_cache[name]

    def extras(part=None, explode=None, **kw):
        return make_extras(PARTS, EXPLODE, part, explode, **kw)

    def place(name, make, key, parent, matrix=None, part=None):
        node = dict(name=name, mesh=mesh_index(name.split('@')[0], make, key))
        if part:
            node['extras'] = extras(part)
        if matrix is not None:
            node['matrix'] = matrix
        add_node(node, parent)

    def module_copy(prefix, names, parent, matrix, role):
        """模組(轉子以外)拆成 PCB / 銅柱 / 馬達三個零件;role = steer / drive。"""
        holder = add_node(dict(name=prefix, matrix=matrix, extras=extras(explode=f'{role}_module')), parent)
        for sub, members in module_groups(names).items():
            sub_node = add_node(dict(name=f'{prefix}_{sub}', extras=extras(
                f'{role}_{sub}', explode_vec=MODULE_EXPLODE[sub])), holder)
            for n in members:
                add_node(dict(name=f'{prefix}_{n}', mesh=nodes[by_name[n]]['mesh']), sub_node)

    def rotor_copy(prefix, names, parent, matrix, role, R, group_key):
        """轉子:與定子同一個零件(件數記 0),爆炸時補償群組位移,與定子停在一起。"""
        holder = add_node(dict(name=prefix, matrix=matrix, extras=extras(
            f'{role}_motor', pcs=0, explode_vec=rotor_explode(EXPLODE, f'{role}_module', R, group_key))), parent)
        for n in names:
            add_node(dict(name=f'{prefix}_{n}', mesh=nodes[by_name[n]]['mesh']), holder)

    rest = [n for n in by_name if n not in rotor_names | skip]
    rotor = [n for n in by_name if n in rotor_names]
    shapes = {name: make() for name, (make, *_rest) in PRINTED.items()}

    nodes[0]['name'] = 'Car_swerve'
    body = add_node(dict(name='body'), 0)
    place('print_top_cover', lambda: cq_to_trimesh(shapes['top_cover']), 'print_white', body, part='top_cover')
    place('print_bottom_plate', lambda: cq_to_trimesh(shapes['bottom_plate']), 'print_black', body,
          part='bottom_plate')
    place('battery', lambda: man_to_trimesh(battery()), 'battery', body, part='battery')
    for i, (x, y) in enumerate(COVER_SCREWS):
        place(f'cover_screw@{i}', lambda: man_to_trimesh(cover_screw()), 'screw', body,
              gltf_matrix(t_mm=[x, y, 0]), part='cover_screws')
        place(f'cover_insert@{i}', lambda: man_to_trimesh(cover_insert()), 'brass', body,
              gltf_matrix(t_mm=[x, y, 0]), part='cover_inserts')
    place('gimbal_inserts', lambda: man_to_trimesh(gimbal_inserts()), 'brass', body, part='gimbal_inserts')
    if args.box:
        box_parts = load_box(args.box)
        box_top = max(m.bounds[1][2] for ms in box_parts.values() for m in ms)
        for (part, mat), ms in box_parts.items():
            place(f'{part}', lambda ms=ms: trimesh.util.concatenate(ms), mat, body,
                  part=part if part != 'box_inserts' else None)
            if part == 'box_inserts':          # 6 個嵌件合成一個網格,件數記 6
                nodes[-1]['extras'] = extras('box_inserts')
        print(f'防水盒:放在上蓋頂面 z={COVER_Z1:.1f},頂部 z={box_top:.1f}(整車高 {box_top:.0f} mm)')
    fork_hw = fork_hardware()
    for cname, (cx, cy, yaw, sign) in CORNERS.items():
        corner = add_node(dict(name=f'corner_{cname}', matrix=gltf_matrix(rz(yaw), [cx, cy, 0])), body)
        place(f'print_bridge@{cname}', lambda: cq_to_trimesh(shapes['bridge']), 'print_white', corner, part='bridge')
        place(f'bridge_screws@{cname}', lambda: man_to_trimesh(bridge_screws()), 'screw', corner,
              part='bridge_screws')
        place(f'leg_screws@{cname}', lambda: man_to_trimesh(leg_screws()), 'screw', corner, part='leg_screws')
        place(f'leg_inserts@{cname}', lambda: man_to_trimesh(leg_inserts()), 'brass', corner, part='leg_inserts')
        module_copy(f'steer_module_{cname}', rest, corner, gltf_matrix(np.eye(3), STEER_T), 'steer')
        # 原地旋轉:滾動方向(舵輪座標 +Y)與以車身中心為圓心的圓相切,取 ±90° 以內
        spin = (np.degrees(np.arctan2(cy, cx)) - yaw + 90) % 180 - 90
        steer = add_node(dict(name=f'steer_{cname}', extras=extras(explode='steer', spin_steer=round(float(spin), 2))),
                         corner)
        rotor_copy(f'steer_module_rotor_{cname}', rotor, steer, gltf_matrix(np.eye(3), STEER_T), 'steer', np.eye(3),
                   'steer')
        place(f'print_fork@{cname}', lambda: cq_to_trimesh(shapes['fork']), 'print_white', steer, part='fork')
        for key, man in fork_hw.items():
            place(f'fork_{key}@{cname}', lambda man=man: man_to_trimesh(man),
                  'brass' if key == 'spacers' else 'screw', steer, part=key)
        module_copy(f'drive_module_{cname}', rest, steer, gltf_matrix(R_D, T_D), 'drive')
        axle = np.array([0, 0, AXLE_Z])
        wheel = add_node(dict(name=f'wheel_{cname}', translation=gltf_matrix(t_mm=axle)[12:15],
                              extras=extras(explode='wheel', spin_sign=sign)), steer)
        m = gltf_matrix(R_D, T_D - axle)
        rotor_copy(f'drive_module_rotor_{cname}', rotor, wheel, m, 'drive', R_D, 'wheel')
        place(f'print_hub@{cname}', lambda: cq_to_trimesh(shapes['hub']), 'print_white', wheel, m, part='hub')
        place(f'print_tire@{cname}', lambda: cq_to_trimesh(shapes['tire']), 'tire', wheel, m, part='tire')
        place(f'hub_screws@{cname}', lambda: man_to_trimesh(hub_screws()), 'screw', wheel, m, part='hub_screws')

    g.save(args.out)
    print('輸出', args.out, f'{Path(args.out).stat().st_size / 1e6:.1f} MB')
    print(f'車身 {BODY_L:.0f} x {BODY_W:.0f},總高 {COVER_Z1:.1f} mm(底板 z={PLATE_Z0:.1f}~{PLATE_Z1:.1f},'
          f'上蓋頂 {COVER_Z1:.1f});輪 Ø{2 * WHEEL_R:.0f},舵輪旋轉半徑 {SWEEP_R:.1f} mm')

    out = Path(args.out_dir)
    for sub in ('print', 'step'):
        (out / sub).mkdir(parents=True, exist_ok=True)
    print(f'列印件(STL 已擺好列印方向,STEP 為零件座標)-> {out}/')
    for name, (_, _, orient, qty) in PRINTED.items():
        shape = shapes[name]
        assert shape.isValid() and len(shape.Solids()) == 1, f'{name} 不是單一有效實體'
        export_print(name, shape, orient, out)
        print(f'{"":17}x {qty}')


if __name__ == '__main__':
    main()
