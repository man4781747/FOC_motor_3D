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

AG-12 電纜固定頭(--ag12):同一面 -X 牆、SEN0676 右側(從盒外看 = -Y 側)上下兩顆,Ø12.2 孔;
孔中心在盒身外圍 -Y 側往 +Y 34 mm,下孔離盒外最底部 19 mm、上孔 44 mm。固定頭在盒外、螺帽在盒內。

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
from shapely.geometry import Point, Polygon
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
SIM_SPACER_L = 14.0          # SIM7000 板:M3 × 14 尼龍母母六角柱(對邊 5.5),裝在主板中段那組孔
SIM_HOLES = (55.1, 61.1)     # SIM7000 板四角孔距(板子 x × y);主板中段孔組是 61.1 × 55.1
SIM_REL_YAW = 270.0          # SIM7000 板相對主板的轉角(孔組轉 90° 或 270° 才對得上;用 270)
SIM_CABLE_L = 100.0          # SIM7000 5P XH 排針座 -> 主板 5P XH 排針座的排線長(兩端膠殼後緣之間的線長)
XH_PITCH, XH_WIRE_D = 2.5, 1.3   # XH 2.5 mm 間距;AWG26 線外徑
XHP_LEN = 9.8                # XHP 膠殼長度(插入排針座 + 露出)
XH_WIRE_COLORS = [(0.75, 0.08, 0.06), (0.9, 0.75, 0.1), (0.1, 0.55, 0.2), (0.1, 0.3, 0.8), (0.05, 0.05, 0.05)]
# 插卡與 OLED(板子座標:mm,板子中心為原點,z=0 是板子背面;位置量自 OBJ / 絲印)
SD_SLOT = dict(x_front=-49.62, y_c=5.6, z0=3.52)      # 主板 microSD 卡座:開口在 -X 板邊,卡底面 z(接點彈片頂)
SD_PROTRUDE = 1.0            # 推入式卡座插好後,卡尾露出卡座前緣
SIM_SLOTS = [(-0.03, -22.98), (16.97, -22.98)]        # SIM7000 板兩個 Micro SIM 掀蓋卡座的中心
SIM_CARD_Z0 = 3.40           # 卡底面(接點彈片頂 3.45,被卡片壓下;卡頂 4.16 < 卡座蓋內面 4.2)
OLED_HEADER = dict(x_pin1=13.25, y=-14.99, z_top=11.04, pins=4)   # 主板 OLED 1×4 母座(pin1 = 3.3V,往 +X)
OLED_PCB = (35.4, 33.5, 1.2)     # 1.3 吋 I2C OLED 模組板(SH1106):寬 × 高 × 厚
OLED_GLASS = (34.5, 23.0, 1.4)   # 面板玻璃
OLED_ACTIVE = (29.42, 14.7)      # 顯示區
OLED_PIN_FROM_TOP = 1.5          # 排針中心離模組板上緣
OLED_HOLE = (2.6, 30.4, 28.5)    # 四角孔徑、孔距(估計)
DS3231_HEADER = dict(x=-24.85, y_pin1=-13.37, z_top=5.08, pins=5)   # 主板 DS3231 1×5 公排針:塑膠頂 z、pin1(+3.3V)在 +Y 端
DS3231_PCB = (13.64, 13.71, 1.0)   # 時鐘模組板(排針方向 × 垂直方向 × 厚),量自賣場照片(游標卡尺 13.64 / 13.71)
DS3231_FEMALE_H = 8.5              # 模組上的 1×5 母座高
# 跳線帽:主板 1×3 排針(+3.3V / IN / +5V),跳在 IN-+5V(選 5V);值為 (x, IN 腳 y, +5V 腳 y),排針塑膠頂 z
JUMPERS = {'RANGE': (-23.59, -33.01, -35.55), 'SPEED': (-6.07, -32.80, -35.33), 'D_SENSOR': (-17.36, 37.10, 34.55)}
JUMPER_Z0 = 5.08
# 主板 LORA:原 1×7 母座(2.54 mm)與 LORA-SUP 1×3 排針已在 fixes/main_pcb_v2.4.json 刪除;
# 改用 XH2.5 7P 90° 公座(S7B-XH-A 類),裝在主板底面,開口朝板子 +Y(LORA-SUP / 上緣那側)
LORA_XH = dict(x=-34.76, y=15.745, pins=7)   # 原 7 個焊墊中心(x 置中)、焊墊列 y
LORA_CABLE_L = 50.0          # 主板 LORA XH 7P(底面臥式)-> RS485 板 XH 7P(立式)排線長(兩端膠殼後緣之間)
RS485_XH = dict(x=-33.75, y=16.12, pins=7)   # RS485 板 XH 7P 立式公座中心(板子座標;量自 OBJ:外形 20.0 × 5.9 × 7.0)
XH_WIRE_COLORS7 = [(0.05, 0.05, 0.05), (0.75, 0.08, 0.06), (0.95, 0.45, 0.05), (0.9, 0.75, 0.1), (0.1, 0.55, 0.2),
                   (0.1, 0.3, 0.8), (0.5, 0.2, 0.7)]
# RS485 板 BAT XH 2P(立式公座)插 XHP-2 母頭,只接一條線(+ 腳,走線接到旁邊 IC 的 BAT 腳)到主板 DC 接頭 + 腳(中心針,
# 離開口最遠的那支)的底面焊點
RS485_BAT = dict(pins=((20.48, -19.76), (22.99, -19.76)), wire_pin=1, housing_w=7.5)
DC_PLUS = (39.2, -36.12, -1.46)   # 主板 DC 接頭 + 腳(板子座標;z = 針腳尾端,板子背面 z=0)
BAT_WIRE_K = 9.0                  # S 形線兩端垂直段的控制長度
# ESP32-S3 Wi-Fi 天線:2.4 GHz 銅管天線(全向 3 dBi),1 代 IPEX(U.FL)、線長 50 mm(賣場圖)
# 外形依賣場照片比例估計:粗段 Ø5.5 × 30、細段 Ø4.6 × 17、頂端鞭狀線 Ø0.8 × 28、同軸線 Ø1.13
ANT_IPEX = (10.48, 2.86, 3.35)    # 主板上 ESP32-S3 模組的 IPEX 座(板子座標;z = 模組板面)
ANT_BASE = (12.0, 29.0, 9.6)      # 天線管底端中心(板子座標):平躺在 SIM7000 板與主板之間,往板子 -X 伸出
ANT_DIR = (-1.0, 0.0, 0.0)
ANT_TAIL_YAW = -90.0              # IPEX 公頭線尾方向:-90° = 板子 +X(朝模組邊緣,不壓在模組屏蔽罩上)
ANT_CABLE_L = 50.0
AG12_HOLE_D = 12.2           # AG-12 開孔(型錄 Ø12 +0.2)
AG12_FROM_SIDE = 34.0        # 孔中心:盒身外圍 -Y 側(從盒外看 -X 牆的右側)往 +Y
AG12_FROM_BOTTOM = (19.0, 44.0)   # 孔中心:盒外最底部往上(下孔 19、上孔再往上 25)

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
    'sim7000': ('SIM7000 IoT 模組板', '電路板', 1, 'My-SIM7000 V1.1,59.7 × 65.7 mm;裝在主板中段孔組上方'),
    'sim_spacers': ('SIM7000 六角柱', '五金', 1, 'M3 × 14 尼龍母母六角柱,對邊 5.5'),
    'sim_screws': ('SIM7000 螺絲', '五金', 1, 'M3 × 6 盤頭,上:鎖 SIM7000 板;下:從主板背面鎖'),
    'sim_cable': ('SIM7000 排線', '線材', 1, f'XH 2.5 5P 雙頭排線 {SIM_CABLE_L:.0f} mm(AWG26),SIM7000 板 -> 主板'),
    'sim_cable_plugs': ('XHP-5 膠殼', '線材', 1, 'JST XH 2.5 mm 5P 母端膠殼(排線兩端,外形簡化)'),
    'sd_card': ('microSD 卡', '外購件', 1, '插在主板 microSD 卡座(推入式,卡尾露出約 1 mm)'),
    'sim_cards': ('Micro SIM 卡', '外購件', 1, '15 × 12 × 0.76 mm,裝在 SIM7000 板兩個掀蓋卡座,接點朝下'),
    'oled': ('1.3 吋 OLED 模組', '外購件', 1, 'SH1106 128 × 64、I2C 4P(VCC / GND / SCL / SDA),35.4 × 33.5 mm;'
             '插在主板 OLED 母座,平躺、螢幕朝上'),
    'ds3231': ('DS3231 時鐘模組', '外購件', 1, 'DS3231SN「For Pi」迷你模組 13.6 × 13.7 mm,1×5 母座(+ / D / C / NC / −),'
               'LIR1220 充電電池;倒插在主板 DS3231 排針上,晶片朝上、電池朝下'),
    'lora_xh': ('LORA XH 7P 90° 公座', '電子', 1, 'JST XH 2.5 mm 7P 臥式公座(S7B-XH-A 類),裝在主板底面、開口朝板子 +Y;'
                '取代原 LORA 1×7 母座'),
    'lora_cable': ('LORA 排線', '線材', 1, f'XH 2.5 7P 雙頭排線 {LORA_CABLE_L:.0f} mm(AWG26),主板 LORA -> RS485 板 7P'),
    'lora_cable_plugs': ('XHP-7 膠殼', '線材', 1, 'JST XH 2.5 mm 7P 母端膠殼(排線兩端,外形簡化)'),
    'bat_wire': ('BAT 電源線', '線材', 1, 'AWG26 單線(紅):RS485 板 BAT 2P 母頭(只接一條) -> 主板 DC 接頭 + 腳底面焊點'),
    'bat_plug': ('XHP-2 膠殼', '線材', 1, 'JST XH 2.5 mm 2P 母端膠殼,只壓一支端子'),
    'wifi_ant': ('Wi-Fi 天線', '外購件', 1, '2.4 GHz 銅管天線(全向 3 dBi)+ 1 代 IPEX 線 50 mm,插在 ESP32-S3 模組 IPEX 座;'
                 '平躺在 SIM7000 板下方'),
    'jumpers': ('跳線帽', '五金', 1, '2.54 mm 跳線帽(黑),主板 RANGE / SPEED / D_SENSOR(濁度電壓)選 +5V:跳在 IN-+5V'),
    'ag12_body': ('AG-12 本體', '外購件', 1, 'KSS AG-12 尼龍電纜固定頭 M12 × 1.25(黑),-X 牆 SEN0676 右側上下兩顆'),
    'ag12_cap': ('AG-12 迫緊頭', '外購件', 1, 'KSS AG-12 附,外徑 15.5'),
    'ag12_seal': ('AG-12 夾緊圈', '外購件', 1, 'KSS AG-12 附,TPE,線徑 3~6.5'),
    'ag12_washer': ('AG-12 墊片', '外購件', 1, 'KSS AGW-12,盒外貼牆'),
    'ag12_nut': ('AG-12 螺帽', '外購件', 1, 'KSS AGL-12,對邊 17,盒內鎖緊'),
    'rs485_posts_top': ('RS485 上方金屬六角柱', '五金', 1, 'M3 × 30 + 6 公母六角柱(對邊 5.5),公頭穿過板子鎖進尼龍柱'),
    'rs485_screws_bottom': ('RS485 底板螺絲', '五金', 1, 'M3 × 6 盤頭,從底板背面鎖進六角柱(底板裝進盒子前先鎖)'),
}
EXPLODE = {   # mm,Z 朝上:照組裝順序往上拆
    'box_body': (0, 0, 0), 'plate': (0, 0, 45), 'plate_screws': (0, 0, 75),
    'rs485_screws_bottom': (0, 0, 25), 'rs485_spacers': (0, 0, 70), 'rs485': (0, 0, 100),
    'rs485_posts_top': (0, 0, 125),
    'main_pcb': (0, 0, 160), 'main_screws': (0, 0, 185),
    'sim_spacers': (0, 0, 200), 'sim7000': (0, 0, 225), 'sim_screws': (0, 0, 240), 'sim_cable': (0, -60, 0), 'sim_cable_plugs': (0, -40, 0),
    'sd_card': (0, 30, 160), 'sim_cards': (0, 0, 245), 'oled': (0, 0, 200), 'ds3231': (0, 0, 205), 'jumpers': (0, 0, 190), 'lora_xh': (0, 0, 135), 'lora_cable': (40, 0, 60), 'lora_cable_plugs': (30, 0, 60), 'bat_wire': (0, -30, 50), 'bat_plug': (0, -20, 30), 'wifi_ant': (0, 0, 175),
    'sen0676': (70, 0, 0), 'sen_oring': (40, 0, 0), 'sen_nut': (-45, 0, 0),
    'ag12_body': (-40, 0, 0), 'ag12_cap': (-70, 0, 0), 'ag12_seal': (-55, 0, 0), 'ag12_washer': (-25, 0, 0),
    'ag12_nut': (30, 0, 0),
    'box_inserts': (0, 0, 255), 'box_gasket': (0, 0, 270), 'box_lid': (0, 0, 300),
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
    'xh_housing': ((0.93, 0.92, 0.86), 0.0, 0.5),
    'sd_black': ((0.06, 0.06, 0.07), 0.0, 0.45),
    'gold': ((0.85, 0.68, 0.30), 1.0, 0.3),
    'sim_white': ((0.92, 0.92, 0.90), 0.0, 0.5),
    'oled_pcb': ((0.05, 0.22, 0.62), 0.0, 0.5),
    'oled_glass': ((0.10, 0.10, 0.11), 0.2, 0.08),
    'oled_screen': ((0.01, 0.01, 0.015), 0.0, 0.05),
    'header_black': ((0.04, 0.04, 0.04), 0.0, 0.6),
    'fpc': ((0.80, 0.45, 0.10), 0.0, 0.5),
    'pcb_black': ((0.03, 0.03, 0.035), 0.0, 0.55),
    'chip_black': ((0.07, 0.07, 0.075), 0.0, 0.6),
    'battery_yellow': ((0.95, 0.78, 0.10), 0.0, 0.5),
    'nickel': ((0.75, 0.74, 0.72), 1.0, 0.3),
    'smd_tan': ((0.80, 0.68, 0.48), 0.0, 0.6),
    'wire_red': ((0.78, 0.06, 0.05), 0.0, 0.5),
    'solder': ((0.78, 0.78, 0.80), 1.0, 0.25),
    'ant_black': ((0.04, 0.04, 0.045), 0.0, 0.55),
    'coax': ((0.03, 0.03, 0.03), 0.0, 0.6),
    'whip': ((0.85, 0.85, 0.84), 0.6, 0.35),
    **{f'wire{i}': (c, 0.0, 0.5) for i, c in enumerate(XH_WIRE_COLORS)},
    **{f'wire7_{i}': (c, 0.0, 0.5) for i, c in enumerate(XH_WIRE_COLORS7)},
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


def pcb_edge_holes(pcb_glb):
    """PCB GLB -> 板子座標(mm,Z 朝上,背面 z=0)的所有 Ø3.0 孔 [(中心, 直徑)] 與板厚。"""
    sc = trimesh.load(pcb_glb, force='scene')
    edge = sc.geometry[next(n for n in sc.geometry if n.startswith('PCB_edge'))].copy()
    edge.apply_transform(np.array([[1, 0, 0, 0], [0, 0, -1, 0], [0, 1, 0, 0], [0, 0, 0, 1]], float))  # Y-up -> Z-up
    edge.apply_scale(1000)
    thick = float(np.ptp(edge.bounds[:, 2]))
    hs = [(c, d) for c, d in loops(edge.section([0, 0, 1], [0, 0, thick / 2])) if abs(d - PCB_HOLE_D) < 0.2]
    return hs, thick


def rect_holes(holes, px, py, tol=0.15):
    """從孔清單找出 x 間距 px、y 間距 py 的矩形四孔組(和軸平行)。"""
    out = []
    for combo in itertools.combinations([c for c, _ in holes], 4):
        pts = np.array(combo)
        c = pts.mean(0)
        if all(abs(abs(p[0] - c[0]) - px / 2) < tol and abs(abs(p[1] - c[1]) - py / 2) < tol for p in pts):
            out.append(pts)
    return out


def pcb_holes(pcb_glb):
    """PCB GLB -> 板子座標(mm,Z 朝上,背面 z=0)的四角孔中心(95.1 × 95.1 那組)與板厚。"""
    hs, thick = pcb_edge_holes(pcb_glb)
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


def cut_holes_x(mesh, holes, x0, x1):
    """在網格上沿 X 打圓孔(manifold3d 布林);holes = [(y, z, r)]。"""
    mm = mesh.copy()
    mm.merge_vertices(merge_tex=True, merge_norm=True)
    m = Manifold(Mesh(vert_properties=np.asarray(mm.vertices, np.float32), tri_verts=np.asarray(mm.faces, np.uint32)))
    for y, z, r in holes:
        m = m - Manifold.cylinder(x1 - x0, r, r, 96).rotate([0, 90, 0]).translate([x0, y, z])
    out = m.to_mesh()
    return trimesh.Trimesh(np.asarray(out.vert_properties)[:, :3], np.asarray(out.tri_verts), process=True)


def find_xh_header(meshes, n_pins=5):
    """在零件網格(世界座標)裡找 XH 臥式排針座:外形約 (n-1)*2.5+4.9 寬、11.5 深、5.75 高的相連區塊。
    回傳 [(lo, hi, 開口方向 (dx, dy), 腔深, 腔中心 z)]。"""
    w = (n_pins - 1) * XH_PITCH + 4.9
    out = []
    for m in meshes:
        for c in m.split(only_watertight=False):
            lo, hi = c.bounds
            d = hi - lo
            if not (abs(max(d[0], d[1]) - w) < 0.4 and abs(min(d[0], d[1]) - 11.5) < 1.0 and 4.5 < d[2] < 7.0):
                continue
            mid = (lo + hi) / 2
            best = None
            for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):        # 從四個側面往內打射線,腔最深的那面是開口
                face = np.array([hi[0] if dx > 0 else lo[0] if dx < 0 else mid[0],
                                 hi[1] if dy > 0 else lo[1] if dy < 0 else mid[1]])
                deps = []
                for zz in np.linspace(lo[2] + 0.3, hi[2] - 0.3, 25):
                    o = [face[0] + dx * 5, face[1] + dy * 5, zz]
                    h = c.ray.intersects_location([o], [[-dx, -dy, 0]])[0]
                    deps.append(float(np.min(np.abs((h[:, :2] - face) @ [dx, dy]))) if len(h) else 0.0)
                deps = np.array(deps)
                k = deps > 1.0
                if k.sum() and (best is None or deps.max() > best[1]):
                    zs = np.linspace(lo[2] + 0.3, hi[2] - 0.3, 25)[k]
                    best = ((dx, dy), float(np.median(deps[k])), float((zs.min() + zs.max()) / 2))
            if best:
                out.append((lo, hi, *best))
    return out


def tube(path, r, seg=12):
    """沿路徑(N × 3)掃圓截面的管子(兩端封口)。"""
    path = np.asarray(path, float)
    t = np.gradient(path, axis=0)
    t /= np.linalg.norm(t, axis=1)[:, None]
    ref = np.array([1.0, 0, 0])
    verts = []
    for p, tt in zip(path, t):
        n = np.cross(tt, ref)
        if np.linalg.norm(n) < 1e-6:
            n = np.cross(tt, [0, 0, 1])
        n /= np.linalg.norm(n)
        b = np.cross(tt, n)
        a = np.linspace(0, 2 * np.pi, seg, endpoint=False)
        verts.append(p + r * (np.cos(a)[:, None] * n + np.sin(a)[:, None] * b))
    V = np.vstack(verts + [path[0], path[-1]])
    F = []
    N = len(path)
    for i in range(N - 1):
        for j in range(seg):
            a0, a1 = i * seg + j, i * seg + (j + 1) % seg
            F += [[a0, a1, a1 + seg], [a0, a1 + seg, a0 + seg]]
    c0, c1 = N * seg, N * seg + 1
    for j in range(seg):
        F += [[c0, (j + 1) % seg, j], [c1, (N - 1) * seg + j, (N - 1) * seg + (j + 1) % seg]]
    m = trimesh.Trimesh(V, F, process=True)
    trimesh.repair.fix_normals(m)
    return m


def j_loop(a, b, length, r1=5.0, r2=4.0, r_u=4.5, step=0.5):
    """XH 排線路徑(y-z 平面,開口朝 -Y):從上方插頭 a 往 -Y 出線、彎下、垂直往下、U 形迴轉、
    往上、彎成水平往 +Y 進下方插頭 b。外圈腿在內圈腿外 2*r_u;用 length 解出最低點。
    a、b = 插頭後緣中心 (x, y, z);回傳 (路徑點, 最低點 z)。"""
    ya, za = a[1], a[2]
    yb, zb = b[1], b[2]
    y_in = min(yb - r2 - 3.0, ya - 8.0)               # 內圈垂直腿(往上接 b)
    y_out = y_in - 2 * r_u                             # 外圈垂直腿(從 a 下來)
    fixed = (ya - (y_out + r1)) + np.pi / 2 * r1 + np.pi * r_u + np.pi / 2 * r2 + (yb - (y_in + r2))
    # 兩條垂直腿:外腿 (za - r1) -> (zl + r_u),內腿 (zl + r_u) -> (zb - r2)
    zl = ((za - r1) + (zb - r2) - (length - fixed)) / 2 - r_u
    pts = []

    def line(p, q):
        n = max(2, int(np.ceil(np.linalg.norm(np.subtract(q, p)) / step)))
        pts.extend(np.linspace(p, q, n, endpoint=False))

    def arc(c, rr, a0, a1):
        n = max(4, int(np.ceil(abs(a1 - a0) * rr / step)))
        for t in np.linspace(a0, a1, n, endpoint=False):
            pts.append([c[0] + rr * np.cos(t), c[1] + rr * np.sin(t)])
    # 在 (y, z) 平面組路徑
    line((ya, za), (y_out + r1, za))
    arc((y_out + r1, za - r1), r1, np.pi / 2, np.pi)
    line((y_out, za - r1), (y_out, zl + r_u))
    arc((y_out + r_u, zl + r_u), r_u, np.pi, 2 * np.pi)
    line((y_in, zl + r_u), (y_in, zb - r2))
    arc((y_in + r2, zb - r2), r2, np.pi, np.pi / 2)
    line((y_in + r2, zb), (yb, zb))
    pts.append((yb, zb))
    yz = np.array(pts, float)
    seglen = np.r_[0, np.cumsum(np.linalg.norm(np.diff(yz, axis=0), axis=1))]
    x = a[0] + (b[0] - a[0]) * seglen / seglen[-1]
    return np.column_stack([x, yz]), zl, float(seglen[-1])


def xhp_plug(n_pins, cav_w, cav_h, insert, total):
    """XHP 母端膠殼(簡化):軸 = +Y(插入方向),y ∈ [-total + insert, insert] 中 y>0 那段插進排針座;
    x 置中、z 置中(z=0 = 腔中心)。後端有 n 個入線孔,頂面有卡勾凸起。"""
    from manifold3d import Manifold as M
    w, h = cav_w - 0.2, cav_h - 0.2
    body = M.cube([w, total, h]).translate([-w / 2, insert - total, -h / 2])
    latch = M.cube([3.0, 2.5, 0.6]).translate([-1.5, insert - total + 1.2, h / 2 - 0.01])
    body = body + latch
    for k in range(n_pins):
        x = (k - (n_pins - 1) / 2) * XH_PITCH
        body = body - M.cylinder(3.0, 0.75, 0.75, 24).rotate([90, 0, 0]).translate([x, insert - total + 2.9, 0])
        body = body - M.cube([1.2, 3.0, 1.2]).translate([x - 0.6, insert - 2.9, -0.6])   # 前端端子孔
    o = body.to_mesh()
    return trimesh.Trimesh(np.asarray(o.vert_properties)[:, :3], np.asarray(o.tri_verts), process=True)


def _m2t(m):
    o = m.to_mesh()
    return trimesh.Trimesh(np.asarray(o.vert_properties)[:, :3], np.asarray(o.tri_verts), process=True)


def _rbox(x0, x1, y0, y1, z0, z1, r=0.0, seg=6):
    """垂直邊倒圓角的方塊(manifold)。"""
    from manifold3d import CrossSection
    if r <= 0:
        return Manifold.cube([x1 - x0, y1 - y0, z1 - z0]).translate([x0, y0, z0])
    cs = CrossSection.square([x1 - x0 - 2 * r, y1 - y0 - 2 * r]).offset(r, circular_segments=4 * seg)
    return cs.extrude(z1 - z0).translate([x0 + r, y0 + r, z0])


def microsd_card():
    """microSD 卡(SD 協會規格 11 × 15 × 0.7/1.0,簡化):卡座標 x = 插入方向(+X),x ∈ [0, 15],
    x=0 是卡尾(握把,厚 1.0);y 置中;z=0 是卡底(接點那面)。回傳 {材質: trimesh}。"""
    from manifold3d import CrossSection
    # 外形:右側(-y)有缺口的段差
    poly = [(0, -5.5), (15, -5.5), (15, 5.5), (0.6, 5.5), (0, 4.9)]
    poly = [(0, -5.5), (6.4, -5.5), (6.4, -4.7), (7.2, -4.7), (8.0, -5.5), (15, -5.5), (15, 5.5), (0.6, 5.5), (0, 4.9)]
    body = CrossSection([poly]).extrude(0.7)
    grip = _rbox(0, 2.0, -5.5, 5.5, 0.6, 1.0, 0.3) ^ CrossSection([poly]).extrude(1.0)
    card = body + grip
    pads = Manifold()
    for k in range(8):
        y = -4.4 + k * 1.1
        pads = pads + Manifold.cube([3.0, 0.8, 0.02]).translate([11.0, y - 0.4, -0.02])
    return {'sd_black': _m2t(card), 'gold': _m2t(pads)}


def micro_sim_card():
    """Micro SIM 卡(3FF:15 × 12 × 0.76,一角 45° 切角):x 置中(寬 12)、y 置中(長 15),z=0 是接點面。"""
    from manifold3d import CrossSection
    c = 2.5
    poly = [(-6, -7.5), (6, -7.5), (6, 7.5 - c), (6 - c, 7.5), (-6, 7.5)]
    card = CrossSection([poly]).extrude(0.76)
    chip = _rbox(-4.0, 4.0, -3.5, 4.5, -0.02, 0.0, 1.0)
    return {'sim_white': _m2t(card), 'gold': _m2t(chip)}


def oled_module():
    """1.3 吋 I2C OLED 模組(SH1106,簡化外形):原點 = 排針列中心、模組板底面;x = 排針方向(pin1 在 -x)、
    模組往 -y 延伸、z 朝上(螢幕面)。排針從板底往下:黑色塑膠 2.5 mm + 插針。回傳 {材質: trimesh}。"""
    W, H, T = OLED_PCB
    gw, gh, gt = OLED_GLASS
    aw, ah = OLED_ACTIVE
    hd, hx, hy = OLED_HOLE
    y_top = OLED_PIN_FROM_TOP
    y_bot = y_top - H
    pcb = _rbox(-W / 2, W / 2, y_bot, y_top, 0, T, 0.5)
    yc = (y_top + y_bot) / 2
    for sx in (-1, 1):
        for sy in (-1, 1):
            pcb = pcb - Manifold.cylinder(T + 2, hd / 2, hd / 2, 32).translate([sx * hx / 2, yc + sy * hy / 2, -1])
    n = OLED_HEADER['pins']
    xs = [(k - (n - 1) / 2) * 2.54 for k in range(n)]
    for x in xs:
        pcb = pcb - Manifold.cylinder(T + 2, 0.5, 0.5, 16).translate([x, 0, -1])
    g_top = y_top - 4.2                                   # 面板玻璃上緣(排針區下方)
    glass = Manifold.cube([gw, gh, gt]).translate([-gw / 2, g_top - gh, T + 0.1])
    ay = g_top - gh / 2 + 1.0
    screen = Manifold.cube([aw, ah, 0.02]).translate([-aw / 2, ay - ah / 2, T + 0.1 + gt])
    fpc = Manifold.cube([12.0, g_top - gh - y_bot + 0.6, 0.12]).translate([-6.0, y_bot - 0.6, T])
    fpc = fpc + Manifold.cube([12.0, 0.12, T + 0.3]).translate([-6.0, y_bot - 0.6, -0.1])
    plastic = _rbox(-n * 1.27, n * 1.27, -1.25, 1.25, -2.5, 0, 0.2)
    pins = Manifold()
    for x in xs:
        pins = pins + Manifold.cube([0.64, 0.64, 6.0 + 2.5 + T + 1.0]).translate([x - 0.32, -0.32, -8.5])
    return {'oled_pcb': _m2t(pcb), 'oled_glass': _m2t(glass), 'oled_screen': _m2t(screen), 'fpc': _m2t(fpc),
            'header_black': _m2t(plastic), 'gold': _m2t(pins)}


def ds3231_module():
    """DS3231 「For Pi」迷你時鐘模組(簡化,依賣場照片)。模組座標 (a, v, z):a 沿排針,從「−」(a = -5.08)
    到「+」(a = +5.08);v = 從排針往模組本體;z 朝上,z=0 是模組板底面(電池那面朝下)。
    母座從板底往下 DS3231_FEMALE_H(開口朝下套在主板公排針上);板頂面是 DS3231SN(SOIC-16W)與電阻電容。
    回傳 {材質: trimesh}。"""
    W, D, T = DS3231_PCB
    n = 5
    pcb = _rbox(-W / 2, W / 2, -1.37, D - 1.37, 0, T, 0.3)
    for k in range(n):
        pcb = pcb - Manifold.cylinder(T + 2, 0.5, 0.5, 16).translate([(k - 2) * 2.54, 0, -1])
    # 母座(開口在底面)
    hd = DS3231_FEMALE_H
    fem = _rbox(-n * 1.27, n * 1.27, -1.27, 1.27, -hd, 0, 0.15)
    for k in range(n):
        fem = fem - Manifold.cube([1.0, 1.0, hd - 0.8]).translate([(k - 2) * 2.54 - 0.5, -0.5, -hd - 0.1])
    # 板頂:DS3231SN(本體 10.3 × 7.5 × 2.35,長邊沿 a,接腳在 ±v 兩側)、焊點、2 顆電阻 + 1 顆電容
    cx, cv = -1.4, 7.4
    chip = _rbox(cx - 5.15, cx + 5.15, cv - 3.75, cv + 3.75, T + 0.15, T + 2.35, 0.2)
    leads = Manifold()
    for k in range(8):
        x = cx + (k - 3.5) * 1.27
        for sgn in (-1, 1):
            leads = leads + Manifold.cube([0.42, 1.4, 0.25]).translate([x - 0.21, cv + sgn * 3.75 + (0 if sgn > 0 else -1.4), T])
            leads = leads + Manifold.cube([0.42, 0.25, 1.2]).translate([x - 0.21, cv + sgn * 3.75 + (-0.05 if sgn > 0 else -0.2), T + 0.2])
    joints = Manifold()
    for k in range(n):
        joints = joints + Manifold.sphere(0.8, 24).scale([1, 1, 0.6]).translate([(k - 2) * 2.54, 0, T])
    smd = Manifold()
    for vv in (3.6, 6.4):
        smd = smd + Manifold.cube([1.25, 2.0, 0.5]).translate([5.4 - 0.62, vv - 1.0, T])
    cap = Manifold.cube([1.25, 2.0, 0.9]).translate([5.4 - 0.62, 9.7 - 1.0, T])
    # 板底:LIR1220 電池(Ø12.5 × 2.0,外圈黃色絕緣套,外露面是鎳色),鎳片從電池中心拉到 + 端板邊焊點
    ba, bv, br, bt, gap = 0.2, 6.6, 6.25, 2.0, 0.4
    z_b1 = -gap
    z_b0 = z_b1 - bt
    yellow = Manifold.cylinder(bt, br, br, 64).translate([ba, bv, z_b0]) - \
        Manifold.cylinder(bt, br - 1.6, br - 1.6, 64).translate([ba, bv, z_b0 - 0.05])
    face = Manifold.cylinder(bt - 0.05, br - 1.6, br - 1.6, 64).translate([ba, bv, z_b0 + 0.05])
    strap = Manifold.cube([W / 2 - ba + 0.3, 5.0, 0.15]).translate([ba, bv - 2.5, z_b0 - 0.15])
    strap = strap + Manifold.cube([0.15, 5.0, bt + gap + 0.15]).translate([W / 2 + 0.3, bv - 2.5, z_b0 - 0.15])
    tabs = Manifold.cube([0.15, 1.6, 2.0]).translate([W / 2 + 0.3, 1.8, -1.2]) + \
        Manifold.cube([0.15, 1.6, 2.0]).translate([W / 2 + 0.3, 3.8, -1.2])
    return {'pcb_black': _m2t(pcb), 'header_black': _m2t(fem), 'chip_black': _m2t(chip), 'nickel': _m2t(leads + face + strap + tabs),
            'gold': _m2t(joints), 'smd_tan': _m2t(smd + cap), 'battery_yellow': _m2t(yellow)}


def jumper_cap():
    """2.54 mm 跳線帽(簡化):5.0 × 2.5 × 6.0,沿 y 套住兩支針(y = ±1.27),底面 z=0、兩個 1.0 方孔深 5。"""
    cap = _rbox(-1.25, 1.25, -2.5, 2.5, 0, 6.0, 0.2)
    cap = (cap - Manifold.cube([2.5, 5.0, 0.6]).translate([-1.25, -2.5, 5.4])) + \
        _rbox(-1.05, 1.05, -2.3, 2.3, 5.4, 6.0, 0.2)                                   # 頂面一圈小段差
    for y in (-1.27, 1.27):
        cap = cap - Manifold.cube([1.0, 1.0, 5.0]).translate([-0.5, y - 0.5, -0.1])
    return _m2t(cap)


def xh_right_angle(n):
    """JST XH 2.5 mm 臥式公座(S nB-XH-A 類,簡化)。座標 (u, v, h):u 沿針排(置中)、v = 開口方向、
    h = 離開板面的方向(h=0 是板面);針排(穿板腳)在 v=0,本體 v ∈ [-2.5, 9.0]、高 5.75,
    開口面在 v=9.0(腔 4.1 高、深 5.0)。穿板腳往 -h 伸出 4.5 mm。回傳 {材質: trimesh}。"""
    W = (n - 1) * 2.5 + 4.9
    body = _rbox(-W / 2, W / 2, -2.5, 9.0, 0, 5.75, 0.0)
    cav = Manifold.cube([W - 1.0, 5.1, 4.1]).translate([-(W - 1.0) / 2, 4.0, 0.85])
    body = body - cav
    for k in (-1, 1):                                            # 側面卡榫缺口
        body = body - Manifold.cube([0.6, 2.4, 1.6]).translate([k * (W / 2 - 0.3) - 0.3, 6.0, 5.75 - 1.0])
    body = body + Manifold.cube([1.6, 4.0, 0.6]).translate([-0.8, 5.0, 4.4])        # 腔頂的鎖扣凸起
    pins = Manifold()
    for k in range(n):
        u = (k - (n - 1) / 2) * 2.5
        pins = pins + Manifold.cube([0.64, 8.5, 0.64]).translate([u - 0.32, -0.32, 2.58])     # 水平段(進腔)
        pins = pins + Manifold.cube([0.64, 0.64, 2.9 + 4.5]).translate([u - 0.32, -0.32, -4.5])  # 穿板腳
    return {'xh_housing': _m2t(body), 'gold': _m2t(pins)}


def hook_path(A, dv, B, length, r=4.0, r2=4.0, step=0.5):
    """排線路徑:從 A 沿水平方向 dv 出線、往下 U 形迴轉(半徑 r)往回走、再彎 90°(半徑 r2)垂直往下到 B。
    B 在 A 的 -dv 側下方;多出來的長度用「出線距離 a」吸收(兩段水平各加 a)。回傳 (N × 3 路徑, a)。"""
    A, B, dv = np.asarray(A, float), np.asarray(B, float), np.asarray(dv, float)
    sA = float((A - B) @ dv)
    zA, zB = A[2], B[2]
    zl = zA - 2 * r
    assert zl - r2 >= zB - 1e-6, f'高度不夠彎:{zA:.1f} -> {zB:.1f}'
    fixed = np.pi * r + (sA - r2) + np.pi / 2 * r2 + (zl - r2 - zB)
    a = (length - fixed) / 2
    assert a >= 0, f'線長 {length} 不夠(至少 {fixed:.1f} mm)'
    pts = []

    def line(p, q):
        n = max(2, int(np.ceil(np.hypot(q[0] - p[0], q[1] - p[1]) / step)))
        pts.extend(np.linspace(p, q, n, endpoint=False))

    def arc(c, rr, a0, a1):
        n = max(4, int(np.ceil(abs(a1 - a0) * rr / step)))
        for t in np.linspace(a0, a1, n, endpoint=False):
            pts.append([c[0] + rr * np.cos(t), c[1] + rr * np.sin(t)])
    line((sA, zA), (sA + a, zA))
    arc((sA + a, zA - r), r, np.pi / 2, -np.pi / 2)
    line((sA + a, zl), (r2, zl))
    arc((r2, zl - r2), r2, np.pi / 2, np.pi)
    line((0, zl - r2), (0, zB))
    pts.append((0, zB))
    sz = np.array(pts, float)
    seg = np.r_[0, np.cumsum(np.linalg.norm(np.diff(sz, axis=0), axis=1))]
    side = (A - B) - ((A - B) @ dv) * dv
    side[2] = 0
    t = 1 - seg / seg[-1]                                   # A 端 t=1、B 端 t=0:側向偏移線性內插
    P = B[None, :2] + sz[:, :1] * dv[None, :2] + t[:, None] * side[None, :2]
    return np.column_stack([P, sz[:, 1]]), a, float(seg[-1])


def wifi_antenna():
    """銅管天線本體(管軸 +X,底端 x=0 中心在原點)+ 頂端鞭狀線。回傳 {材質: trimesh}。"""
    from manifold3d import Manifold as M
    def cylx(r0, r1, x0, x1):
        return M.cylinder(x1 - x0, r0, r1, 48).rotate([0, 90, 0]).translate([x0, 0, 0])
    tube = cylx(1.2, 2.75, 0, 2.5) + cylx(2.75, 2.75, 2.5, 30) + cylx(2.75, 2.3, 30, 31.5) + cylx(2.3, 2.3, 31.5, 46) \
        + cylx(2.3, 0.8, 46, 47)
    whip = cylx(0.4, 0.4, 46.5, 47 + 28)
    return {'ant_black': _m2t(tube), 'whip': _m2t(whip)}


def ipex_plug():
    """1 代 IPEX(U.FL)公頭(簡化):圓帽 Ø2.6 套在座上(z 0.4~2.5,z=0 是座底面),線尾往 +y。回傳 (網格, 線尾端點)。"""
    from manifold3d import Manifold as M
    cap = M.cylinder(2.1, 1.3, 1.3, 32).translate([0, 0, 0.4])
    tail = _rbox(-0.75, 0.75, 0, 3.0, 1.3, 2.4, 0.2)
    ferrule = M.cylinder(1.0, 0.75, 0.75, 24).rotate([-90, 0, 0]).translate([0, 2.6, 1.85])
    return _m2t(cap + tail + ferrule), np.array([0, 3.6, 1.85])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('box')
    ap.add_argument('plate')
    ap.add_argument('out')
    ap.add_argument('--rs485', help='RS485 擴展板 GLB(build_glb.py 的輸出),裝在底板上')
    ap.add_argument('--main', help='水情盒主板 GLB,裝在 RS485 板上方金屬六角柱頂端(需要 --rs485)')
    ap.add_argument('--sim7000', help='SIM7000 IoT 模組板 GLB,用 14 mm 尼龍六角柱裝在主板上(需要 --main)')
    ap.add_argument('--sen0676', help='SEN0676 雷達 GLB(build_sen0676.py 的輸出),裝在 -X 牆的 Ø48.5 孔')
    ap.add_argument('--ag12', action='store_true', help='在 -X 牆 SEN0676 右側加上下兩顆 AG-12 電纜固定頭')
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

    # ---- SIM7000 IoT 模組板:M3 × 14 尼龍六角柱,裝在主板中段那組 Ø3 孔 ----
    sim = None
    if args.sim7000:
        assert main, '--sim7000 需要同時給 --main'
        mh, _ = pcb_edge_holes(args.main)
        msets = rect_holes(mh, SIM_HOLES[1], SIM_HOLES[0])
        assert len(msets) == 1, f'主板上應有一組 {SIM_HOLES[1]} × {SIM_HOLES[0]} 的 Ø3 孔,找到 {len(msets)} 組'
        sh, sim_t = pcb_edge_holes(args.sim7000)
        ssets = rect_holes(sh, *SIM_HOLES)
        assert len(ssets) == 1, f'SIM7000 板應有一組 {SIM_HOLES[0]} × {SIM_HOLES[1]} 的 Ø3 孔'
        main_top = main['t'][2] + main['thick']
        sim_holes = (Rm[:2, :2] @ msets[0].T).T + main['t'][:2]              # 主板中段孔(世界座標)
        Rsim = rz(MAIN_YAW + SIM_REL_YAW)
        slocal = ssets[0]
        z_sim = main_top + SIM_SPACER_L
        tsim = np.array([*(sim_holes.mean(0) - Rsim[:2, :2] @ slocal.mean(0)), z_sim])
        world = (Rsim[:2, :2] @ slocal.T).T + tsim[:2]
        serr = max(np.min(np.hypot(*(sim_holes - w).T)) for w in world)
        assert serr < 0.3, f'SIM7000 板孔位對不上主板中段孔({serr:.2f} mm)'
        sim = dict(R=Rsim, t=tsim, thick=sim_t, holes=sim_holes)
        ssc = trimesh.load(args.sim7000, force='scene')
        sm = []
        for n in ssc.graph.nodes_geometry:
            T, gname = ssc.graph[n]
            m = ssc.geometry[gname].copy()
            m.apply_transform(Zup @ T)
            m.apply_scale(1000)
            m.apply_transform(affine(Rsim, tsim))
            sm.append(m)
        smesh = trimesh.util.concatenate(sm)
        slo, shi = smesh.bounds
        # SIM7000 板與主板之間的實際最小垂直間隙:在 SIM 板範圍打一格網射線,
        # SIM 這側取最低點(往上打),主板這側取最高點(往下打)
        main_top_mesh = trimesh.util.concatenate([m for m in mm if m.bounds[1][2] > main_top - 0.01])
        gap_min, gap_at = 99.0, None
        for x in np.arange(slo[0] + 0.25, shi[0], 0.5):
            for y in np.arange(slo[1] + 0.25, shi[1], 0.5):
                hs = smesh.ray.intersects_location([[x, y, main_top - 1]], [[0, 0, 1]])[0]
                if not len(hs):
                    continue
                hm_ = main_top_mesh.ray.intersects_location([[x, y, float(hs[:, 2].min()) - 0.001]], [[0, 0, -1]])[0]
                if len(hm_):
                    gg = float(hs[:, 2].min()) - float(hm_[:, 2].max())
                    if gg < gap_min:
                        gap_min, gap_at = gg, (x, y, float(hm_[:, 2].max()))
        under = gap_at[2]
        sim_bottom = under + gap_min
        # 盒蓋內面(平的):在 SIM7000 板範圍的四角與中心往上打射線取最低
        lid_low = min(float(lid_mesh.ray.intersects_location([[x, y, 0]], [[0, 0, 1]])[0][:, 2].min())
                      for x, y in [shi[:2], slo[:2], [shi[0], slo[1]], [slo[0], shi[1]], (slo[:2] + shi[:2]) / 2])
        lid_gap = lid_low - float(shi[2])
        print(f'SIM7000 板:裝在主板中段孔組(中心 ({sim_holes.mean(0)[0]:+.2f}, {sim_holes.mean(0)[1]:+.2f})),'
              f'轉 {(MAIN_YAW + SIM_REL_YAW) % 360:.0f}°,孔位誤差 {serr:.3f} mm;M3 × {SIM_SPACER_L:.0f} 尼龍六角柱')
        print(f'  板子背面 z={z_sim:.1f}、正面 z={z_sim + sim_t:.1f},最高零件 z={shi[2]:.1f};'
              f'範圍 X {slo[0]:.1f}~{shi[0]:.1f},Y {slo[1]:.1f}~{shi[1]:.1f}')
        print(f'  [{"OK" if gap_min > 0 else "干涉"}] SIM7000 板(含背面焊腳)與主板零件最小垂直間隙 {gap_min:.1f} mm'
              f'(在 ({gap_at[0]:.1f}, {gap_at[1]:.1f}),主板零件頂 z={under:.1f})')
        print(f'  [{"OK" if lid_gap > 0.5 else "干涉"}] 盒蓋內面 z={lid_low:.1f},SIM7000 最高零件'
              f'{"離盒蓋 " + format(lid_gap, ".1f") if lid_gap > 0 else "伸進盒蓋 " + format(-lid_gap, ".1f")} mm')
        if lid_gap <= 0.5:
            need = 0.5 - lid_gap
            room = gap_min - 0.3
            print(f'    要放得下(上方留 0.5 mm):SIM7000 板要降 {need:.1f} mm;柱子最多能縮短 {room:.1f} mm'
                  f'(下方留 0.3 mm)-> {"柱子改 " + format(SIM_SPACER_L - need, ".1f") + " mm 以下即可" if room >= need else "只改這組柱子不夠,整疊要再降 " + format(need - room, ".1f") + " mm"}')

    # ---- SIM7000 <-> 主板:XH 2.5 5P 雙頭排線 ----
    cable = None
    if sim:
        hs_ = find_xh_header(sm)
        hm_ = find_xh_header(mm)
        assert len(hs_) == 1, f'SIM7000 板上應有一個 5P XH 臥式排針座,找到 {len(hs_)} 個'
        top_h = hs_[0]
        bot_h = min(hm_, key=lambda h: np.hypot(*((h[0][:2] + h[1][:2]) / 2 - (top_h[0][:2] + top_h[1][:2]) / 2)))
        assert top_h[2] == bot_h[2] == (0, -1), f'目前只處理兩個開口都朝 -Y 的情況({top_h[2]}, {bot_h[2]})'
        plugs, ends = [], []
        for lo, hi, _, depth, zc in (top_h, bot_h):
            xc = float((lo[0] + hi[0]) / 2)
            face_y = float(lo[1])
            insert = depth - 0.2
            cav_w = float(hi[0] - lo[0]) - 1.0                 # 腔寬(兩側壁各約 0.5)
            cav_h = 4.1
            plug = xhp_plug(5, cav_w, cav_h, insert, XHP_LEN)
            plug.apply_translation([xc, face_y, zc])
            plugs.append(plug)
            ends.append(np.array([xc, face_y + insert - XHP_LEN, zc]))
        path, zl, plen = j_loop(ends[0], ends[1], SIM_CABLE_L)
        wires = []
        for k in range(5):
            dx = (k - 2) * XH_PITCH
            wires.append(tube(path + [dx, 0, 0], XH_WIRE_D / 2))
        cable = dict(plugs=plugs, wires=wires)
        allw = trimesh.util.concatenate(wires)
        wl, wh = allw.bounds
        # 盒壁:-Y 牆內面
        wall = body.ray.intersects_location([[float(path[:, 0].mean()), 0, float(zl) + 5]], [[0, -1, 0]])[0]
        wall_y = float(wall[:, 1].max()) if len(wall) else -99
        print(f'SIM7000 排線:XH 2.5 5P,線長 {plen:.1f} mm(插頭後緣之間),XHP 膠殼長 {XHP_LEN}')
        print(f'  上端插頭:SIM7000 排針座 x={ends[0][0]:.1f},腔中心 z={ends[0][2]:.1f};'
              f'下端插頭:主板排針座 x={ends[1][0]:.1f},z={ends[1][2]:.1f}(兩座前緣 y 差 {top_h[0][1] - bot_h[0][1]:.1f})')
        print(f'  走線:往 -Y 出線、彎下、U 形迴轉回到主板(J 形),最低點 z={wl[2]:.1f},'
              f'最外 y={wl[1]:.1f}(-Y 牆內面 y={wall_y:.1f},間隙 {wl[1] - wall_y:.1f} mm)')
        # 盒內其他東西:電路板(RS485 / 主板外框)y 範圍
        print(f'  排線範圍 X {wl[0]:.1f}~{wh[0]:.1f}、Y {wl[1]:.1f}~{wh[1]:.1f}、Z {wl[2]:.1f}~{wh[2]:.1f};'
              f'電路板外緣 y={min(mmesh.bounds[0][1], -50):.1f}')

    # ---- 插卡與 OLED:microSD(主板)、Micro SIM × 2(SIM7000 板)、1.3 吋 OLED(主板 OLED 母座) ----
    extras = []          # [(零件 id, {材質: trimesh(世界座標 mm)})]
    if main:
        Mmain = affine(main['R'], main['t'])
        sd = microsd_card()
        x0 = SD_SLOT['x_front'] - SD_PROTRUDE
        for m in sd.values():
            m.apply_translation([x0, SD_SLOT['y_c'], SD_SLOT['z0']])
            m.apply_transform(Mmain)
        extras.append(('sd_card', sd))
        sdw = trimesh.util.concatenate(list(sd.values())).bounds
        print(f'microSD 卡:插進主板卡座(板子座標 x {x0:.1f}~{x0 + 15:.1f}、卡底 z={SD_SLOT["z0"]:.2f}),'
              f'卡尾露出卡座 {SD_PROTRUDE:.1f} mm;世界座標 X {sdw[0][0]:.1f}~{sdw[1][0]:.1f}、Y {sdw[0][1]:.1f}~{sdw[1][1]:.1f}')
        # OLED:排針插進母座,黑色塑膠座在母座頂上
        h = OLED_HEADER
        xc = h['x_pin1'] + (h['pins'] - 1) * 2.54 / 2
        z_pcb = h['z_top'] + 2.5
        ol = oled_module()
        for m in ol.values():
            m.apply_translation([xc, h['y'], z_pcb])
            m.apply_transform(Mmain)
        extras.append(('oled', ol))
        olm = trimesh.util.concatenate(list(ol.values()))
        ob = olm.bounds
        # 模組板底面與主板零件(正下方)最小垂直間隙;與 SIM7000 板 / 盒蓋
        pcb_only = ol['oled_pcb']
        main_parts = trimesh.util.concatenate([m for m in mm if m.bounds[1][2] > main['t'][2] + main['thick'] + 0.01])
        gmin, gat = 99.0, None
        pb = pcb_only.bounds
        for x in np.arange(pb[0][0] + 0.25, pb[1][0], 0.5):
            for y in np.arange(pb[0][1] + 0.25, pb[1][1], 0.5):
                hh = main_parts.ray.intersects_location([[x, y, pb[0][2] - 0.01]], [[0, 0, -1]])[0]
                if len(hh):
                    gg = pb[0][2] - float(hh[:, 2].max())
                    if gg < gmin:
                        gmin, gat = gg, (x, y)
        lid_o = min(float(lid_mesh.ray.intersects_location([[x, y, 0]], [[0, 0, 1]])[0][:, 2].min())
                    for x, y in [ob[0][:2], ob[1][:2], (ob[0][:2] + ob[1][:2]) / 2])
        hit_sim = sim is not None and not (ob[1][0] < slo[0] or ob[0][0] > shi[0] or ob[1][1] < slo[1] or ob[0][1] > shi[1])
        print(f'OLED 1.3 吋模組:插在主板 OLED 母座(pin1 = 3.3V 對 VCC),平躺、螢幕朝上,模組往板子 -Y(板子座標)延伸;'
              f'模組板底 z={pb[0][2]:.1f}、玻璃頂 z={ob[1][2]:.1f}')
        print(f'  範圍 X {ob[0][0]:.1f}~{ob[1][0]:.1f}、Y {ob[0][1]:.1f}~{ob[1][1]:.1f};'
              f'[{"OK" if gmin > 0.3 else "干涉"}] 模組板底與主板零件最小間隙 {gmin:.1f} mm(在 ({gat[0]:.1f}, {gat[1]:.1f}));'
              f'離盒蓋內面 {lid_o - ob[1][2]:.1f} mm;{"[干涉] 與 SIM7000 板範圍重疊" if hit_sim else "不在 SIM7000 板下方"}')
        # DS3231:倒插在主板 DS3231 公排針(母座坐在公排針塑膠上),模組本體往 -X(板子座標)
        h = DS3231_HEADER
        yc = h['y_pin1'] - (h['pins'] - 1) * 2.54 / 2
        Rd = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], float)      # a -> +Y、v -> -X、z -> z
        z_mod = h['z_top'] + DS3231_FEMALE_H
        ds = ds3231_module()
        for m in ds.values():
            m.apply_transform(affine(Rd, [h['x'], yc, z_mod]))
            m.apply_transform(Mmain)
        extras.append(('ds3231', ds))
        dsm = trimesh.util.concatenate(list(ds.values()))
        db = dsm.bounds
        bat = ds['battery_yellow'].bounds
        g2, g2at = 99.0, None
        for x in np.arange(bat[0][0] + 0.25, bat[1][0], 0.5):
            for y in np.arange(bat[0][1] + 0.25, bat[1][1], 0.5):
                hh = main_parts.ray.intersects_location([[x, y, bat[0][2] - 0.2]], [[0, 0, -1]])[0]
                if len(hh):
                    gg = bat[0][2] - 0.15 - float(hh[:, 2].max())
                    if gg < g2:
                        g2, g2at = gg, (x, y)
        lid_d = min(float(lid_mesh.ray.intersects_location([[x, y, 0]], [[0, 0, 1]])[0][:, 2].min())
                    for x, y in [db[0][:2], db[1][:2], (db[0][:2] + db[1][:2]) / 2])
        hit_o = not (db[1][0] < ob[0][0] or db[0][0] > ob[1][0] or db[1][1] < ob[0][1] or db[0][1] > ob[1][1])
        print(f'DS3231 時鐘模組:倒插在主板 DS3231 排針(+ / D / C / NC / − 對 3.3V / SDA / SCL / NC / GND),'
              f'晶片朝上、電池朝下;模組板底 z={z_mod + main["t"][2]:.1f}、最高 z={db[1][2]:.1f}')
        print(f'  範圍 X {db[0][0]:.1f}~{db[1][0]:.1f}、Y {db[0][1]:.1f}~{db[1][1]:.1f};'
              f'[{"OK" if g2 > 0.3 else "干涉"}] 電池(含鎳片)離下方主板零件最少 {g2:.1f} mm;'
              f'離盒蓋內面 {lid_d - db[1][2]:.1f} mm;{"[干涉] 與 OLED 重疊" if hit_o else "不碰 OLED"}')
        # 跳線帽:RANGE / SPEED / D_SENSOR 跳在 IN-+5V
        jc = jumper_cap()
        for name, (x, y_in, y_5v) in JUMPERS.items():
            m = jc.copy()
            m.apply_translation([x, (y_in + y_5v) / 2, JUMPER_Z0])
            m.apply_transform(Mmain)
            extras.append(('jumpers', {'header_black': m}))
        print(f'跳線帽 × {len(JUMPERS)}:{" / ".join(JUMPERS)} 跳在 IN-+5V(套在排針塑膠頂 z={JUMPER_Z0},帽頂與針頂齊)')
        # LORA:XH2.5 7P 90° 公座裝在主板底面(h -> -Z),開口朝板子 +Y;u = v × h -> -X
        L = LORA_XH
        Rl = np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]], float)
        xh = xh_right_angle(L['pins'])
        for m in xh.values():
            m.apply_transform(affine(Rl, [L['x'], L['y'], 0.0]))
            m.apply_transform(Mmain)
        extras.append(('lora_xh', xh))
        xb = trimesh.util.concatenate(list(xh.values())).bounds
        hous = xh['xh_housing'].bounds
        rsc = trimesh.load(args.rs485, force='scene')
        rparts = []
        for n_ in rsc.graph.nodes_geometry:
            T_, gname_ = rsc.graph[n_]
            m_ = rsc.geometry[gname_].copy()
            m_.apply_transform(Zup @ T_)
            m_.apply_scale(1000)
            m_.apply_transform(affine(pcb['R'], pcb['t']))
            rparts.append(m_)
        below = trimesh.util.concatenate(rparts)
        rparts_all = below
        gx = []
        for x in np.linspace(hous[0][0] + 0.5, hous[1][0] - 0.5, 6):
            for y in np.linspace(hous[0][1] + 0.5, hous[1][1] - 0.5, 6):
                hh = below.ray.intersects_location([[x, y, hous[0][2] - 0.01]], [[0, 0, -1]])[0]
                if len(hh):
                    gx.append(hous[0][2] - float(hh[:, 2].max()))
        print(f'LORA XH 7P 90° 公座:主板底面,開口朝板子 +Y;原焊墊 2.54 mm 間距、XH 是 2.5 mm(7P 總長差 '
              f'{6 * 0.04:.2f} mm,PCB 要改封裝);本體底 z={hous[0][2]:.1f},離下方 RS485 板零件最少 {min(gx):.1f} mm')
        print(f'  範圍 X {xb[0][0]:.1f}~{xb[1][0]:.1f}、Y {xb[0][1]:.1f}~{xb[1][1]:.1f}、Z {xb[0][2]:.1f}~{xb[1][2]:.1f}')
        # LORA 排線:主板底面臥式 XH 7P -> RS485 板立式 XH 7P
        Mrs = affine(pcb['R'], pcb['t'])
        n7 = L['pins']
        W7 = (n7 - 1) * 2.5 + 4.9
        # 主板端:插頭沿 -v 插入,前端在腔底(v = 9.0 - 5.0 + 0.2),後緣 = 前端 + XHP_LEN;腔中心 h = 0.85 + 2.05
        Mh = Mmain @ affine(Rl, [L['x'], L['y'], 0.0])          # 公座座標 (u, v, h) -> 世界
        dv = Mh[:3, :3] @ [0, 1, 0]
        h_c = 0.85 + 2.05
        plugA = xhp_plug(n7, W7 - 1.0, 4.1, 4.8, XHP_LEN)       # 插頭座標:x 沿針、+y 插入方向、z 厚度
        Pa = affine(np.column_stack([[1, 0, 0], [0, -1, 0], [0, 0, -1]]), [0, 9.0, h_c])   # 插頭 -> 公座座標
        plugA.apply_transform(Mh @ Pa)
        A = (Mh @ [0, 9.0 + XHP_LEN - 4.8, h_c, 1])[:3]
        # RS485 端:立式公座,插頭沿 -Z 插入;腔底用射線量
        rx, ry = RS485_XH['x'], RS485_XH['y']
        c_top = (Mrs @ [rx, ry, 0, 1])[:3]
        floor_ = rparts_all.ray.intersects_location([[*((Mrs @ [rx + 1.25, ry, 30, 1])[:2]), 200]], [[0, 0, -1]])[0]
        top_z = float((Mrs @ [0, 0, 9.5, 1])[2])
        cav_floor = float(floor_[:, 2].max()) if len(floor_) else top_z - 5.0     # 開口朝上:第一個交點就是腔底
        ins = top_z - cav_floor - 0.2
        plugB = xhp_plug(n7, W7 - 1.0, 4.1, ins, XHP_LEN)
        ex = Mrs[:3, :3] @ [1, 0, 0]
        Rb = np.column_stack([ex, [0, 0, -1], np.cross(ex, [0, 0, -1])])
        plugB.apply_transform(affine(Rb, [c_top[0], c_top[1], top_z]))
        B = np.array([c_top[0], c_top[1], top_z + XHP_LEN - ins])
        cpath, a_out, clen = hook_path(A, dv, B, LORA_CABLE_L)
        side_dir = np.cross([0, 0, 1], dv)                        # 排線寬度方向(水平、垂直 dv)
        ua = Mh[:3, :3] @ [1, 0, 0]
        wires7 = []
        for k in range(n7):
            off = (k - (n7 - 1) / 2) * 2.5
            # 主板端針排方向 ua、RS485 端 ex;排線寬度方向都是水平,取它們在 side_dir 上的正負一致
            sgn_a = np.sign(ua @ side_dir) or 1
            wires7.append(tube(cpath + off * sgn_a * side_dir, XH_WIRE_D / 2))
        extras.append(('lora_cable_plugs', {'xh_housing': trimesh.util.concatenate([plugA, plugB])}))
        for k, w in enumerate(wires7):
            extras.append(('lora_cable', {f'wire7_{k}': w}))
        cw = trimesh.util.concatenate(wires7).bounds
        print(f'LORA 排線:XH 2.5 7P {LORA_CABLE_L:.0f} mm(路徑長 {clen:.1f});主板端插頭後緣 z={A[2]:.1f},'
              f'RS485 端插頭插入 {ins:.1f} mm、後緣 z={B[2]:.1f};兩端水平距離 {np.hypot(*(A - B)[:2]):.1f} mm')
        print(f'  走線:沿主板 +Y 出線 {a_out:.1f} mm、往下 U 形迴轉、往回再彎下插入 RS485;範圍 X {cw[0][0]:.1f}~{cw[1][0]:.1f}、'
              f'Y {cw[0][1]:.1f}~{cw[1][1]:.1f}、Z {cw[0][2]:.1f}~{cw[1][2]:.1f}')
        # 與兩塊板子零件的最近距離(排線中心線)
        allb = trimesh.util.concatenate([main_parts, rparts_all] + [m for m in mm if m.bounds[1][2] <= main['t'][2] + 0.01])
        d = np.abs(trimesh.proximity.signed_distance(allb, cpath[::2]))
        print(f'  [{"OK" if d.min() > XH_WIRE_D / 2 + 0.2 else "注意"}] 排線中心線離兩塊板子 / 零件最近 {d.min():.1f} mm(線半徑 {XH_WIRE_D / 2})')
        # BAT 線:RS485 BAT 2P 母頭(單線)-> 主板 DC 接頭 + 腳底面
        bp = np.array(RS485_BAT['pins'])
        bc = bp.mean(0)
        W2 = RS485_BAT['housing_w']
        flo = rparts_all.ray.intersects_location([[*((Mrs @ [bc[0], bc[1] + 1.2, 30, 1])[:2]), 200]], [[0, 0, -1]])[0]
        btop = float(rparts_all.ray.intersects_location([[*((Mrs @ [bc[0] - W2 / 2 + 0.3, bc[1], 30, 1])[:2]), 200]],
                                                       [[0, 0, -1]])[0][:, 2].max())
        bfloor = float(flo[:, 2].max())
        bins = btop - bfloor - 0.2
        plug2 = xhp_plug(2, W2 - 1.0, 4.1, bins, XHP_LEN)
        cb = (Mrs @ [bc[0], bc[1], 0, 1])[:3]
        plug2.apply_transform(affine(Rb, [cb[0], cb[1], btop]))
        pw = (Mrs @ [*bp[RS485_BAT['wire_pin']], 0, 1])[:3]
        A2 = np.array([pw[0], pw[1], btop + XHP_LEN - bins])
        dcp = (Mmain @ [*DC_PLUS, 1])[:3]
        B2 = dcp + [0, 0, -0.4]                                   # 線頭包住針腳尾端
        k = BAT_WIRE_K
        P0, P1, P2, P3 = A2, A2 + [0, 0, k], B2 - [0, 0, k], B2
        tt = np.linspace(0, 1, 160)[:, None]
        bez = (1 - tt) ** 3 * P0 + 3 * (1 - tt) ** 2 * tt * P1 + 3 * (1 - tt) * tt ** 2 * P2 + tt ** 3 * P3
        blen = float(np.linalg.norm(np.diff(bez, axis=0), axis=1).sum())
        bw = tube(bez, XH_WIRE_D / 2)
        blob = trimesh.creation.icosphere(2, 1.1)
        blob.apply_scale([1, 1, 1.4])
        blob.apply_translation(dcp + [0, 0, 0.3])
        extras.append(('bat_plug', {'xh_housing': plug2}))
        extras.append(('bat_wire', {'wire_red': bw, 'solder': blob}))
        # 離其他東西:兩塊板子零件、六角柱
        posts = []
        for x_, y_ in pcb['holes']:
            posts.append((x_, y_))
        dpost = min(float(np.min(np.hypot(bez[:, 0] - x_, bez[:, 1] - y_))) for x_, y_ in posts) - TOP_POST_FLATS / np.sqrt(3)
        mid = bez[10:-10]                                          # 兩端會碰到插頭 / 焊點,不算
        dmin = float(np.abs(trimesh.proximity.signed_distance(allb, mid[::2])).min())
        print(f'BAT 線:RS485 BAT 2P 母頭(只接 x={bp[RS485_BAT["wire_pin"]][0]:.2f} 那支,+)插入 {bins:.1f} mm -> 主板 DC 接頭 + 腳'
              f'(板子座標 ({DC_PLUS[0]}, {DC_PLUS[1]}))底面焊點;線長約 {blen:.0f} mm(S 形)')
        print(f'  兩端水平距離 {np.hypot(*(A2 - B2)[:2]):.1f}、高度差 {B2[2] - A2[2]:.1f};'
              f'[{"OK" if dmin > XH_WIRE_D / 2 + 0.2 else "注意"}] 線離兩塊板子 / 零件最近 {dmin:.1f} mm,離六角柱 {dpost:.1f} mm')
        # Wi-Fi 天線:IPEX 公頭插在 ESP32-S3 模組 IPEX 座,同軸線 S 形繞到天線底端,天線平躺在兩塊板子之間
        ip, tail_end = ipex_plug()
        Rt = rz(ANT_TAIL_YAW)
        ip.apply_transform(affine(Rt, ANT_IPEX))
        E = np.array(ANT_IPEX) + Rt @ tail_end                 # 板子座標
        dT = Rt @ [0, 1, 0]
        Bb = np.array(ANT_BASE)
        dA = np.array(ANT_DIR) / np.linalg.norm(ANT_DIR)
        def bez_len(k):
            P0, P1, P2, P3 = E, E + dT * k, Bb - dA * k, Bb
            t = np.linspace(0, 1, 200)[:, None]
            pts = (1 - t) ** 3 * P0 + 3 * (1 - t) ** 2 * t * P1 + 3 * (1 - t) * t ** 2 * P2 + t ** 3 * P3
            return pts, float(np.linalg.norm(np.diff(pts, axis=0), axis=1).sum())
        want = ANT_CABLE_L - 3.6                               # 線長從公頭算起(扣掉線尾)
        lo_k, hi_k = 0.5, 60.0
        for _ in range(50):
            mk = (lo_k + hi_k) / 2
            if bez_len(mk)[1] < want:
                lo_k = mk
            else:
                hi_k = mk
        cpts, clen2 = bez_len((lo_k + hi_k) / 2)
        coax = tube(cpts, 1.13 / 2)
        ant = wifi_antenna()
        Ra = np.column_stack([dA, np.cross([0, 0, 1], dA), [0, 0, 1]])
        for m in ant.values():
            m.apply_transform(affine(Ra, Bb))
        parts_ant = {'nickel': ip, 'coax': coax, **ant}
        for m in parts_ant.values():
            m.apply_transform(Mmain)
        extras.append(('wifi_ant', parts_ant))
        # 檢查:天線 / 線離其他東西(兩塊板子零件、SIM7000 板)
        others = trimesh.util.concatenate([allb] + sm)
        body_pts = trimesh.util.concatenate(list(ant.values())).sample(1500)
        dant = float(np.abs(trimesh.proximity.signed_distance(others, body_pts)).min())
        cw_ = (Mmain @ np.c_[cpts, np.ones(len(cpts))].T).T[:, :3]
        dco = float(np.abs(trimesh.proximity.signed_distance(others, cw_[15:-5:3])).min()) - 1.13 / 2
        ab = trimesh.util.concatenate(list(ant.values())).bounds
        print(f'Wi-Fi 天線:IPEX 公頭插在 ESP32-S3 模組 IPEX 座(板子座標 ({ANT_IPEX[0]}, {ANT_IPEX[1]})),同軸線 {ANT_CABLE_L:.0f} mm'
              f'(S 形,路徑 {clen2 + 3.6:.1f});天線平躺在主板與 SIM7000 板之間,中心高 z={ab.mean(0)[2]:.1f}')
        print(f'  天線範圍 X {ab[0][0]:.1f}~{ab[1][0]:.1f}、Y {ab[0][1]:.1f}~{ab[1][1]:.1f}、Z {ab[0][2]:.1f}~{ab[1][2]:.1f};'
              f'[{"OK" if dant > 0.3 else "干涉"}] 天線離其他零件最近 {dant:.1f} mm、'
              f'[{"OK" if dco > 0.2 else "注意"}] 同軸線最近 {dco:.1f} mm')
    if sim:
        Msim = affine(sim['R'], sim['t'])
        for i, (x, y) in enumerate(SIM_SLOTS):
            sc_ = micro_sim_card()
            for m in sc_.values():
                m.apply_translation([x, y, SIM_CARD_Z0])
                m.apply_transform(Msim)
            extras.append(('sim_cards', sc_))
        print(f'Micro SIM 卡 × {len(SIM_SLOTS)}:裝在 SIM7000 板兩個掀蓋卡座內(接點朝下,卡頂 z(板子座標)'
              f'{SIM_CARD_Z0 + 0.76:.2f},卡座蓋內面 4.2)')

    # ---- SEN0676 雷達:-X 牆開 Ø40 孔,鏡頭朝外,六角一組對邊平行盒底 ----
    sen = None
    cut = []                                      # -X 牆要打的孔 [(y, z, r)]
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
        cut.append((hy, hz, SEN_HOLE_D / 2))
        tor = trimesh.creation.torus(r_groove + o_cs_on / 2, o_cs_on / 2, major_sections=96, minor_sections=16)
        tor.apply_scale([1, 1, 1 - SEN_ORING_SQUEEZE])              # 軸向被壓扁
        tor.apply_translation([0, 0, z_flange - gap / 2])
        sen['oring'] = tor
        nut = sen_nut_mesh(nut_w, nut_h, 2 * r_thread)
        nut.apply_translation([0, 0, z_nut_end - worst])
        sen['nut'] = nut

    # ---- AG-12 電纜固定頭:-X 牆、SEN0676 右側(從盒外看)上下兩顆,固定頭在盒外、螺帽在盒內 ----
    ag = None
    if args.ag12:
        import build_ag12 as A
        gy = float(body.bounds[0][1]) + AG12_FROM_SIDE
        zb = float(body.bounds[0][2])
        ag = dict(pos=[])
        # 固定頭座標:z 朝外(-X),面板外側面 z=0;轉 30° 讓六角一組平邊平行盒底
        R0 = np.array([[0, 0, -1], [0, 1, 0], [1, 0, 0]], float) @ rz(30)
        print(f'AG-12:-X 牆,盒身外圍 -Y 側 y={body.bounds[0][1]:.2f} 往 +Y {AG12_FROM_SIDE:.0f} -> 孔中心 y={gy:.2f},'
              f'Ø{AG12_HOLE_D} 孔')
        panels = []
        for dz in AG12_FROM_BOTTOM:
            gz = zb + dz
            xs = sorted(body.ray.intersects_location([[0, gy, gz]], [[-1, 0, 0]])[0][:, 0], reverse=True)
            ix, ox = float(xs[0]), float(xs[-1])
            # 外牆有拔模角:墊片(Ø17)範圍內最凸的點 -> 墊片貼在那裡
            prot = max(ox - float(body.ray.intersects_location([[ox - 20, gy + dy, gz + dz2]], [[1, 0, 0]])[0][:, 0].min())
                       for dy in np.linspace(-8.5, 8.5, 5) for dz2 in np.linspace(-8.5, 8.5, 5))
            # 盒內:螺帽(外接 Ø19.6)範圍內牆內面最凹的點
            inner = min(float(body.ray.intersects_location([[0, gy + dy, gz + dz2]], [[-1, 0, 0]])[0][:, 0].max())
                        for dy in np.linspace(-9.5, 9.5, 5) for dz2 in np.linspace(-9.5, 9.5, 5))
            panel = inner - (ox - prot)
            panels.append(panel)
            ag['pos'].append((gy, gz, ox - prot))
            cut.append((gy, gz, AG12_HOLE_D / 2))
            engage = A.L1 - A.WASHER_T - panel
            print(f'  孔中心 z={gz:.1f}(盒外底往上 {dz:.0f}):牆外面 x={ox:.2f}(墊片處最凸 {prot:.2f})、'
                  f'內面 x={ix:.2f},牆厚 {panel:.2f};螺帽可咬螺紋 {engage:.1f} mm(螺帽高 {A.NUT_H})')
            # 盒內:螺帽範圍與底板 / 盒底 / 雷達
            nut_x1 = inner + A.NUT_H
            lo_z = gz - A.NUT_AF / np.sqrt(3)
            if plate.bounds[0][0] < nut_x1 and lo_z < z_top and gz + A.NUT_AF / np.sqrt(3) > tops[0]:
                print(f'    [干涉] 螺帽(x {inner:.1f}~{nut_x1:.1f}、z {lo_z:.1f}~{gz + A.NUT_AF / np.sqrt(3):.1f})'
                      f'碰到底板 -X 邊緣 x={plate.bounds[0][0]:.1f}(z {tops[0]:.1f}~{z_top:.1f})')
            fl = body.ray.intersects_location([[inner + 2, gy, gz]], [[0, 0, -1]])[0]
            print(f'    螺帽最低點 z={lo_z:.1f},該處盒內底 z={float(fl[:, 2].max()):.1f}')
        if sen:
            d = abs(gy - sen['t'][1]) - SEN_NUT[0] / 2 - 1.5 - A.WASHER_OD / 2
            print(f'  盒外:AG-12 墊片離 SEN0676 螺帽 {d:.1f} mm(y 方向)')
        panel = float(np.mean(panels))
        ag['R'] = R0
        ag['parts'] = {'ag12_body': (A.make_body(), 'nylon_black'), 'ag12_cap': (A.make_cap(), 'nylon_black'),
                       'ag12_seal': (A.make_seal(), 'rubber'), 'ag12_washer': (A.make_washer(), 'nylon_black'),
                       'ag12_nut': (A.make_nut(panel), 'nylon_black')}
    if cut:
        xs = sorted(body.ray.intersects_location([[0, cut[0][0], cut[0][1]]], [[-1, 0, 0]])[0][:, 0], reverse=True)
        box[('box_body', 'box_gray')] = [cut_holes_x(body, cut, xs[-1] - 3, xs[0] + 3)]

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
    if sim:
        kids = merge_glb(g, args.sim7000)
        nodes.append(dict(name='sim7000', children=kids, matrix=gltf_matrix(sim['R'], sim['t']),
                          extras=make_extras(PARTS, EXPLODE, 'sim7000')))
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

    def add_mesh(name, tm_mm, key, part, matrix=None, mesh=None, pcs=None):
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
        node['extras'] = make_extras(PARTS, EXPLODE, part, pcs=pcs)
        if matrix is not None:
            node['matrix'] = matrix
        return node['mesh']

    for (part, key), ms in box.items():
        add_mesh(part, trimesh.util.concatenate(ms), key, part)
    add_mesh('plate', plate, 'print_white', 'plate')
    if sen:
        add_mesh('sen_oring', sen['oring'], 'rubber', 'sen_oring', gltf_matrix(sen['R'], sen['t']))
        add_mesh('sen_nut', sen['nut'], 'nylon_black', 'sen_nut', gltf_matrix(sen['R'], sen['t']))
    if ag:
        for part, (m, key) in ag['parts'].items():
            mg = m.to_mesh()
            tm = trimesh.Trimesh(np.asarray(mg.vert_properties)[:, :3], np.asarray(mg.tri_verts), process=True)
            mid = None
            for i, (y, z, x) in enumerate(ag['pos']):
                mid = add_mesh(f'{part}_{i}', tm, key, part, gltf_matrix(ag['R'], [x, y, z]), mid)

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

    for part, ms in extras:
        first = True
        for key, m in ms.items():
            add_mesh(f'{part}_{key}', m, key, part, pcs=1 if first else 0)
            first = False

    if sim:
        def to_tm2(m):
            mg = m.to_mesh()
            return trimesh.Trimesh(np.asarray(mg.vert_properties)[:, :3], np.asarray(mg.tri_verts), process=True)
        if cable:
            pid = None
            for i, pl in enumerate(cable['plugs']):
                add_mesh(f'sim_cable_plug_{i}', pl, 'xh_housing', 'sim_cable_plugs')
            for k, w in enumerate(cable['wires']):
                add_mesh(f'sim_cable_wire_{k}', w, f'wire{k}', 'sim_cable', None, None)
        sp = to_tm2(standoff(3.0, SIM_SPACER_L, SPACER_FLATS, hole_depth=SPACER_THREAD_DEPTH))
        m3s = to_tm2(pan_screw(3.0, PCB_SCREW_L))
        ids = dict(s=None, top=None, bot=None)
        for i, (x, y) in enumerate(sim['holes']):
            ids['s'] = add_mesh(f'sim_spacer_{i}', sp, 'nylon', 'sim_spacers',
                                gltf_matrix(t_mm=[x, y, sim['t'][2]]), ids['s'])
            ids['top'] = add_mesh(f'sim_screw_top_{i}', m3s, 'screw', 'sim_screws',
                                  gltf_matrix(t_mm=[x, y, sim['t'][2] + sim['thick']]), ids['top'])
            ids['bot'] = add_mesh(f'sim_screw_bottom_{i}', m3s, 'screw', 'sim_screws',
                                  gltf_matrix(np.diag([1.0, -1.0, -1.0]), [x, y, main['t'][2]]), ids['bot'])

    g.save(args.out)
    print('輸出', args.out, f'{Path(args.out).stat().st_size / 1e6:.1f} MB')


if __name__ == '__main__':
    main()
