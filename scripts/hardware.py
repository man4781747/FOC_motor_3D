"""五金零件模型(manifold3d,mm):有螺紋的盤頭 / 沉頭螺絲、公母 / 母母銅柱、六角螺帽。

螺紋做法:截面的半徑沿圓周呈三角波(牙底 -> 牙頂 -> 牙底),每升高一個螺距扭轉 360°,
得到單頭 V 形螺紋。只是外觀用,牙型比例不是精確的 ISO 規格。

座標慣例:
- 盤頭螺絲:頭的座面在 z=0,頭往 +Z,螺桿往 -Z(長度 L 不含頭)
- 沉頭螺絲:頭頂面在 z=0,頭與螺桿往 -Z(長度 L 含頭,與市售規格相同)
- 銅柱:本體 z ∈ [-L, 0];z=0 端為母牙孔,z=-L 端為母牙孔或公頭(往 -Z 伸出)
"""
import numpy as np
from manifold3d import CrossSection, Manifold

from build_motor import box, union

SEG = 32                  # 五金的圓分段數(小零件不需要馬達用的 128 段)
HW_MATERIALS = {'screw', 'screws', 'brass'}


def shade(mesh, material, angle):
    """平滑著色。五金不經過 trimesh.graph.smooth_shade:它會把螺紋網格的頂點拆成約 5 倍,
    檔案大好幾倍;小零件直接用合併後的頂點法向量就夠了。"""
    import trimesh
    if material in HW_MATERIALS:
        return mesh
    return trimesh.graph.smooth_shade(mesh, angle=angle)


def revolve(profile):
    """繞 Z 軸旋轉;截面自動轉成逆時針(順時針會變成無效實體)。"""
    p = np.asarray(profile, float)
    area = np.sum(p[:, 0] * np.roll(p[:, 1], -1) - np.roll(p[:, 0], -1) * p[:, 1])
    return CrossSection([[tuple(q) for q in (p if area > 0 else p[::-1])]]).revolve(SEG)

PITCH = {2.0: 0.4, 2.5: 0.45, 3.0: 0.5}
DIV_PER_TURN = 6          # 每圈切幾段(越多越平滑,三角形越多)
PROFILE_PTS = 12          # 截面點數


def thread(d, length, pitch=None):
    """沿 +Z 的外螺紋,z ∈ [0, length],z=length 端倒角。"""
    pitch = pitch or PITCH[d]
    r_maj, r_min = d / 2, d / 2 - 0.6 * pitch
    th = np.linspace(0, 2 * np.pi, PROFILE_PTS, endpoint=False)
    tri = 1 - np.abs(2 * (th / (2 * np.pi)) - 1)                # 0 -> 1 -> 0 繞一圈
    r = r_min + (r_maj - r_min) * tri
    cs = CrossSection([list(zip(r * np.cos(th), r * np.sin(th)))])
    turns = length / pitch
    m = cs.extrude(length, int(np.ceil(turns * DIV_PER_TURN)), 360.0 * turns)
    tip = revolve([(0, -0.01), (r_maj + 0.1, -0.01), (r_maj + 0.1, length - pitch),
                   (r_min - 0.1, length + 0.01), (0, length + 0.01)])
    return m ^ tip


def _cross_recess(r, h):
    s = 0.5
    return union([box(-r * 0.55, r * 0.55, -s / 2, s / 2, h * 0.45, h + 1),
                  box(-s / 2, s / 2, -r * 0.55, r * 0.55, h * 0.45, h + 1)])


def pan_screw(d, length, head_d=None, head_h=None):
    """十字盤頭螺絲:座面 z=0,頭往 +Z,螺桿往 -Z。"""
    head_d = head_d or d * 1.8
    head_h = head_h or d * 0.6
    r = head_d / 2
    head = revolve([(0, 0), (r, 0), (r, head_h * 0.55), (r * 0.8, head_h * 0.9), (r * 0.5, head_h), (0, head_h)])
    head = head - _cross_recess(r, head_h)
    shank = thread(d, length).mirror([0, 0, 1])               # z ∈ [-length, 0],尖端在 -length
    return head + shank


def flat_screw(d, length, head_d=None):
    """十字沉頭螺絲(90°):頭頂面 z=0,往 -Z;length 含頭。"""
    head_d = head_d or d * 2.0
    hk = (head_d - d) / 2                                      # 90° 錐的高度
    head = revolve([(0, -hk), (d / 2, -hk), (head_d / 2, 0), (0, 0)])
    head = head - _cross_recess(head_d / 2, hk).mirror([0, 0, 1]).translate([0, 0, hk])   # 槽口在頭頂面
    shank = thread(d, length - hk + 0.01).mirror([0, 0, 1]).translate([0, 0, -hk + 0.01])
    return head + shank


def _hex(flats, h):
    r = flats / np.sqrt(3)
    return Manifold.cylinder(h, r, circular_segments=6).rotate([0, 0, 30])


def standoff(d, length, flats=4.0, stud=0.0, hole_depth=None):
    """六角銅柱,本體 z ∈ [-length, 0]。z=0 端母牙孔;stud > 0 時 z=-length 端為公頭(往 -Z 伸出),
    否則也是母牙孔(母母銅柱)。"""
    hole_depth = hole_depth or min(length / 2 - 0.2, 4.0)
    body = _hex(flats, length).translate([0, 0, -length])
    # 六角兩端倒角
    r = flats / np.sqrt(3)
    body = body ^ revolve([(0, -length - 0.01), (r - 0.3, -length - 0.01), (r + 0.01, -length + 0.3),
                           (r + 0.01, -0.3), (r - 0.3, 0.01), (0, 0.01)])
    r_hole = d / 2 - 0.3
    holes = [revolve([(0, -hole_depth), (r_hole, -hole_depth), (r_hole, -0.35), (r_hole + 0.35, 0.01), (0, 0.01)])]
    if stud > 0:
        body = body + thread(d, stud + 0.01).mirror([0, 0, 1]).translate([0, 0, -length + 0.01])
    else:
        holes.append(revolve([(0, -length - 0.01), (r_hole + 0.35, -length - 0.01), (r_hole, -length + 0.35),
                              (r_hole, -length + hole_depth), (0, -length + hole_depth)]))
    return body - union(holes)


def nut(d, flats=None, h=None):
    """六角螺帽,z ∈ [0, h]。"""
    flats = flats or d * 2.0
    h = h or d * 0.8
    r = flats / np.sqrt(3)
    body = _hex(flats, h) ^ revolve([(0, -0.01), (r - 0.25, -0.01), (r + 0.01, 0.3),
                                      (r + 0.01, h - 0.3), (r - 0.25, h + 0.01), (0, h + 0.01)])
    return body - Manifold.cylinder(h + 2, d / 2 - 0.25, circular_segments=SEG).translate([0, 0, -1])


INSERT_M3 = dict(d=3.0, od=4.5, length=5.0, hole_d=4.0, hole_depth=6.0)   # M3 × 5 熱熔螺母與列印件開孔


def heat_insert(d=3.0, od=4.5, length=5.0):
    """埋頭型(滾花)熱熔螺母:頂面 z=0(與列印件表面齊平),本體往 -Z。
    上下兩段斜向滾花、中間一圈凹槽,內孔畫成 M 牙的小徑。"""
    def knurl(z0, h, twist):
        th = np.linspace(0, 2 * np.pi, 40, endpoint=False)
        r = np.where(np.arange(40) % 2 == 0, od / 2, od / 2 - 0.25)
        cs = CrossSection([list(zip(r * np.cos(th), r * np.sin(th)))])
        return cs.extrude(h, 4, twist).translate([0, 0, z0])
    band = (length - 0.8) / 2
    body = knurl(-length, band, 25) + knurl(-band, band, -25)
    body = body + Manifold.cylinder(0.8 + 0.02, od / 2 - 0.45, circular_segments=SEG).translate([0, 0, -band - 0.81])
    bore = revolve([(0, -length - 0.1), (d / 2 - 0.25, -length - 0.1), (d / 2 - 0.25, -0.35), (d / 2 + 0.1, 0.1), (0, 0.1)])
    return body - bore
