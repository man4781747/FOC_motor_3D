# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "trimesh", "scipy", "networkx"]
# ///
"""把 EasyEDA 匯出的 OBJ zip + Gerber 貼圖組合成一個 GLB(glTF 2.0)。

用法:
    uv run scripts/build_glb.py <OBJ.zip> <貼圖資料夾> <輸出.glb> [修正檔.json]

修正檔(選用):EasyEDA 有些零件的 3D 模型本身就壞掉(例如面被拉成大三角形),
可以用 JSON 刪掉該零件的材質群組,換成簡單的方塊:
    {"drop": ["mtl17"],
     "add": [{"name": "L1_body", "color": [0.2, 0.2, 0.2], "metal": 0, "rough": 0.6,
              "box": [x0, x1, y0, y1, z0, z1], "round": 0.8}]}
座標是 mm,xy 以板子中心為原點(與貼圖一致,+Y 朝上),z=0 是板子背面(板子正面 z=2.54)。
round > 0 時四條垂直邊倒圓角。

處理內容:
- EasyEDA OBJ 單位是 10 mil(0.254 mm),換算成公尺(glTF 標準單位)
- 以板子(mtl1)中心為原點;OBJ 的 Z 軸朝上 -> glTF 的 Y 軸朝上
- 板子拆成 正面 / 背面 / 側邊,正背面依 XY 座標貼上 Gerber 貼圖
- 其餘零件依 OBJ 材質分成獨立網格,套用 PBR 材質
"""
import json
import struct
import sys
import zipfile
from collections import OrderedDict
from pathlib import Path

import numpy as np
import trimesh

UNIT_MM = 0.254          # EasyEDA OBJ 1 單位 = 10 mil
TEX_HALF_MM = 20.32      # 與 make_textures.py 的 HALF 一致(貼圖資料夾沒有 texture.json 時使用)
BOARD_MTL = 'mtl1'


def srgb_to_linear(c):
    c = np.asarray(c, dtype=np.float64)
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def parse_mtl(text):
    mats, cur = {}, None
    for line in text.splitlines():
        t = line.split()
        if not t:
            continue
        if t[0] == 'newmtl':
            cur = t[1]
            mats[cur] = {}
        elif t[0] in ('Kd', 'Ks') and cur:
            mats[cur][t[0]] = [float(x) for x in t[1:4]]
    return mats


def parse_obj(text):
    """回傳 頂點、法線、以及 {材質: [(v, n) * 3 的三角形, ...]}。"""
    vs, ns = [], []
    groups = OrderedDict()
    cur = None
    for line in text.splitlines():
        if not line:
            continue
        c = line[0]
        if c == 'v':
            t = line.split()
            (vs if t[0] == 'v' else ns).append([float(x) for x in t[1:4]])
        elif c == 'f':
            tri = []
            for p in line.split()[1:]:
                a = p.split('/')
                tri.append((int(a[0]) - 1, int(a[2]) - 1 if len(a) > 2 and a[2] else -1))
            for i in range(1, len(tri) - 1):  # 多邊形扇形三角化
                groups.setdefault(cur, []).append((tri[0], tri[i], tri[i + 1]))
        elif line.startswith('usemtl'):
            cur = line.split()[1]
    return np.array(vs), np.array(ns), groups


def pbr_for(name, kd, ks):
    """依 MTL 顏色推估 PBR 參數:中性灰且反光強的當成金屬(引腳、外殼)。"""
    kd = np.array(kd)
    ks = np.array(ks if ks is not None else [0, 0, 0])
    sat = kd.max() - kd.min()
    spec = ks.max() / max(kd.max(), 1e-6)
    metallic = 0.0
    rough = float(np.clip(0.85 - spec * 0.6, 0.3, 0.9))
    if sat < 0.05 and kd.max() >= 0.6 and ks.max() >= 0.3:
        metallic, rough = 0.85, 0.32
    return dict(name=name, doubleSided=True, pbrMetallicRoughness=dict(
        baseColorFactor=[*srgb_to_linear(kd).tolist(), 1.0],
        metallicFactor=metallic, roughnessFactor=rough))


class GLB:
    def __init__(self):
        self.bin = bytearray()
        self.j = dict(asset=dict(version='2.0', generator='easyEDA-obj-to-glb'),
                      scene=0, scenes=[dict(nodes=[0])],
                      nodes=[dict(name='PCB', children=[])], meshes=[], materials=[],
                      accessors=[], bufferViews=[], buffers=[], textures=[], images=[],
                      samplers=[dict(magFilter=9729, minFilter=9987, wrapS=33071, wrapT=33071)])

    def _view(self, data, target=None):
        while len(self.bin) % 4:
            self.bin.append(0)
        bv = dict(buffer=0, byteOffset=len(self.bin), byteLength=len(data))
        if target:
            bv['target'] = target
        self.bin += data
        self.j['bufferViews'].append(bv)
        return len(self.j['bufferViews']) - 1

    def accessor(self, arr, kind):
        if kind == 'index':
            arr = arr.astype(np.uint32).ravel()
            view = self._view(arr.tobytes(), 34963)
            acc = dict(bufferView=view, componentType=5125, count=arr.size, type='SCALAR')
        else:
            arr = arr.astype(np.float32)
            view = self._view(arr.tobytes(), 34962)
            acc = dict(bufferView=view, componentType=5126, count=len(arr),
                       type={2: 'VEC2', 3: 'VEC3'}[arr.shape[1]])
            if kind == 'pos':
                acc['min'] = arr.min(0).tolist()
                acc['max'] = arr.max(0).tolist()
        self.j['accessors'].append(acc)
        return len(self.j['accessors']) - 1

    def image(self, path):
        view = self._view(Path(path).read_bytes())
        self.j['images'].append(dict(bufferView=view, mimeType='image/png'))
        self.j['textures'].append(dict(source=len(self.j['images']) - 1, sampler=0))
        return len(self.j['textures']) - 1

    def material(self, m):
        self.j['materials'].append(m)
        return len(self.j['materials']) - 1

    def mesh(self, name, pos, nrm, idx, mat, uv=None):
        attrs = dict(POSITION=self.accessor(pos, 'pos'), NORMAL=self.accessor(nrm, 'nrm'))
        if uv is not None:
            attrs['TEXCOORD_0'] = self.accessor(uv, 'uv')
        self.j['meshes'].append(dict(name=name, primitives=[dict(
            attributes=attrs, indices=self.accessor(idx, 'index'), material=mat)]))
        self.j['nodes'].append(dict(name=name, mesh=len(self.j['meshes']) - 1))
        self.j['nodes'][0]['children'].append(len(self.j['nodes']) - 1)

    def save(self, path):
        while len(self.bin) % 4:
            self.bin.append(0)
        self.j['buffers'] = [dict(byteLength=len(self.bin))]
        if not self.j['textures']:  # glTF 規範不允許空陣列
            for k in ('textures', 'images', 'samplers'):
                self.j.pop(k, None)
        js = json.dumps(self.j, separators=(',', ':')).encode()
        js += b' ' * (-len(js) % 4)
        total = 12 + 8 + len(js) + 8 + len(self.bin)
        with open(path, 'wb') as f:
            f.write(struct.pack('<4sII', b'glTF', 2, total))
            f.write(struct.pack('<I4s', len(js), b'JSON') + js)
            f.write(struct.pack('<I4s', len(self.bin), b'BIN\0') + bytes(self.bin))


SMOOTH_ANGLE = np.radians(35)


def to_mesh(V, tris, center):
    """取出三角形,換算成公尺(仍是 OBJ 的 Z-up),合併重複頂點。"""
    vi = np.array(tris)[:, :, 0]
    uniq, inv = np.unique(vi, return_inverse=True)
    verts = (V[uniq] - center) * UNIT_MM / 1000.0
    return trimesh.Trimesh(verts, inv.reshape(-1, 3), process=True)


def clean_overlapping_planes(mesh, ratio=1.3):
    """EasyEDA 會把含孔洞(例如晶片印字)的平面三角化成互相重疊的三角形,
    渲染時會閃爍。偵測同一水平面上「三角形面積總和 > 外框凸包面積 * ratio」的
    相連區塊,把它換成乾淨的凸包面。
    回傳 (新 mesh, 平面清單[(z, 凸包頂點 XY, 是否朝上)])。"""
    from scipy.spatial import ConvexHull
    fn = mesh.face_normals
    horiz = np.abs(fn[:, 2]) > 0.999
    if not horiz.any():
        return mesh, []
    zf = np.round(mesh.vertices[mesh.faces][:, :, 2].mean(1), 6)
    keep = np.ones(len(mesh.faces), bool)
    new_v, new_f = [], []
    base = len(mesh.vertices)
    planes = []
    for z in np.unique(zf[horiz]):
        ids = np.where(horiz & (zf == z))[0]
        if len(ids) < 8:
            continue
        sub = mesh.submesh([ids], append=True)
        comps = trimesh.graph.connected_components(sub.face_adjacency, nodes=np.arange(len(sub.faces)))
        for comp in comps:
            fids = ids[comp]
            tri = mesh.vertices[mesh.faces[fids]]
            area = np.abs(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])[:, 2]).sum() / 2
            pts = tri[:, :, :2].reshape(-1, 2)
            try:
                hull = ConvexHull(pts)
            except Exception:
                continue
            if area <= hull.volume * ratio:      # 2D 凸包的 volume 就是面積
                continue
            up = (fn[fids, 2] > 0).sum() >= len(fids) / 2
            ring = pts[hull.vertices]             # 逆時針
            if not up:
                ring = ring[::-1]
            start = base + sum(len(v) for v in new_v)
            new_v.append(np.column_stack([ring, np.full(len(ring), z)]))
            new_f += [[start, start + i, start + i + 1] for i in range(1, len(ring) - 1)]
            keep[fids] = False
            planes.append((z, pts[hull.vertices], up))
    if not planes:
        return mesh, []
    verts = np.vstack([mesh.vertices] + new_v)
    faces = np.vstack([mesh.faces[keep], np.array(new_f)])
    return trimesh.Trimesh(verts, faces, process=True), planes


MARK_DEPTH = 0.1 * UNIT_MM / 1000   # 印字凹刻深度上限(約 0.025 mm)
MARK_LIFT = 0.02 / 1000             # 抬到平面外 0.02 mm


def lift_markings(mesh, planes):
    """印字(晶片型號等)在 OBJ 裡是比外殼表面略凹的水平面。外殼表面被
    clean_overlapping_planes 換成完整平面後會蓋住印字,這裡把它們抬到平面外側。"""
    from scipy.spatial import Delaunay
    if not planes:
        return 0
    fn = mesh.face_normals
    cz = mesh.triangles_center
    lifted = 0
    for z, ring, up in planes:
        tri = Delaunay(ring)
        depth = (z - cz[:, 2]) if up else (cz[:, 2] - z)
        sel = (np.abs(fn[:, 2]) > 0.999) & (depth > 0) & (depth < MARK_DEPTH)
        sel &= tri.find_simplex(cz[:, :2]) >= 0
        if not sel.any():
            continue
        vids = np.unique(mesh.faces[sel])
        verts = mesh.vertices.copy()
        verts[vids, 2] = z + (MARK_LIFT if up else -MARK_LIFT)
        mesh.vertices = verts
        lifted += int(sel.sum())
    return lifted


def fix_and_shade(mesh):
    """修正繞序使面朝外,再依角度做平滑著色(超過角度的邊保持銳利)。"""
    trimesh.repair.fix_normals(mesh, multibody=True)
    return trimesh.graph.smooth_shade(mesh, angle=SMOOTH_ANGLE)


def gltf_arrays(mesh, normals=None):
    """Z-up -> glTF Y-up:(x, y, z) -> (x, z, -y),這是旋轉,不會改變繞序。"""
    v = np.asarray(mesh.vertices)
    n = np.asarray(mesh.vertex_normals) if normals is None else normals
    rot = lambda a: np.stack([a[:, 0], a[:, 2], -a[:, 1]], 1)
    return rot(v), rot(n), np.asarray(mesh.faces)


def rounded_prism(x0, x1, y0, y1, z0, z1, r=0.0, seg=6):
    """垂直邊倒圓角的方塊(mm),回傳 trimesh(公尺)。底面、頂面是凸多邊形,用扇形三角化。"""
    r = min(r, (x1 - x0) / 2 - 1e-3, (y1 - y0) / 2 - 1e-3)
    if r <= 0:
        ring = np.array([(x0, y0), (x1, y0), (x1, y1), (x0, y1)])
    else:
        ring = []
        for cx, cy, a0 in ((x1 - r, y0 + r, -90), (x1 - r, y1 - r, 0), (x0 + r, y1 - r, 90), (x0 + r, y0 + r, 180)):
            a = np.radians(np.linspace(a0, a0 + 90, seg + 1))
            ring += list(zip(cx + r * np.cos(a), cy + r * np.sin(a)))
        ring = np.array(ring)
    n = len(ring)
    v = np.vstack([np.column_stack([ring, np.full(n, z0)]), np.column_stack([ring, np.full(n, z1)]),
                   [[(x0 + x1) / 2, (y0 + y1) / 2, z0], [(x0 + x1) / 2, (y0 + y1) / 2, z1]]])
    f = []
    for i in range(n):
        j = (i + 1) % n
        f += [[i, j, n + j], [i, n + j, n + i], [2 * n, j, i], [2 * n + 1, n + i, n + j]]
    m = trimesh.Trimesh(v / 1000.0, np.array(f), process=True)
    trimesh.repair.fix_normals(m)
    return m


def main():
    obj_zip, tex_dir, out = sys.argv[1:4]
    fix = json.loads(Path(sys.argv[4]).read_text()) if len(sys.argv) > 4 else {}
    tex_dir = Path(tex_dir)
    zf = zipfile.ZipFile(obj_zip)
    obj_name = next(n for n in zf.namelist() if n.lower().endswith('.obj'))
    mtl_name = next(n for n in zf.namelist() if n.lower().endswith('.mtl'))
    mats = parse_mtl(zf.read(mtl_name).decode('utf-8'))
    V, N, groups = parse_obj(zf.read(obj_name).decode('utf-8'))

    board_vi = np.unique(np.array(groups[BOARD_MTL])[:, :, 0])
    bmin, bmax = V[board_vi].min(0), V[board_vi].max(0)
    center = np.array([(bmin[0] + bmax[0]) / 2, (bmin[1] + bmax[1]) / 2, 0.0])
    size_mm = (bmax - bmin) * UNIT_MM
    print(f'板子中心(OBJ 座標) {center[:2]},外形 {size_mm[0]:.2f} × {size_mm[1]:.2f} mm,'
          f'厚 {size_mm[2]:.2f} mm')
    # 貼圖範圍:板子中心對到貼圖的 center_mm(make_textures.py --fit-outline 時不是原點)
    tex_half = TEX_HALF_MM
    if (tex_dir / 'texture.json').exists():
        tex_half = json.loads((tex_dir / 'texture.json').read_text())['half_mm']

    g = GLB()

    # ---- 板子:修正繞序後依面法線分成正面 / 背面 / 側邊 ----
    # 板子有開孔,trimesh 修正繞序不可靠,改用頂點高度分類:三點都在頂面/底面高度
    board = to_mesh(V, groups[BOARD_MTL], center)
    fzv = board.vertices[board.faces][:, :, 2]
    zt, zb = board.vertices[:, 2].max(), board.vertices[:, 2].min()
    is_top = np.isclose(fzv, zt).all(1)
    is_bot = np.isclose(fzv, zb).all(1)
    sub = {k: board.submesh([np.where(m)[0]], append=True)
           for k, m in dict(top=is_top, bottom=is_bot, edge=~(is_top | is_bot)).items()}
    # 正面三角形繞序朝 +Z、背面朝 -Z,側邊朝外(遠離板子中心軸)
    for k, sign in (('top', 1), ('bottom', -1)):
        m = sub[k]
        flip = m.face_normals[:, 2] * sign < 0
        m.faces[flip] = m.faces[flip][:, ::-1]
    e = sub['edge']
    radial = e.triangles_center.copy(); radial[:, 2] = 0
    # 外緣朝外;開孔的內壁中心離軸較近,但仍以「朝向孔外的板材」判斷較複雜,雙面材質可涵蓋
    flip = (e.face_normals * radial).sum(1) < 0
    # 外緣:離板子邊界 1 mm 內(圓板用半徑、方板用 bbox 都適用)
    half_xy = (bmax[:2] - bmin[:2]) * UNIT_MM / 2000
    r = np.linalg.norm(radial[:, :2], axis=1)
    outer = (r > half_xy.min() - 1e-3) if np.isclose(*half_xy) and r.max() < half_xy.min() * 1.05 \
        else (np.abs(radial[:, :2]) > half_xy - 1e-3).any(1)
    e.faces[flip & outer] = e.faces[flip & outer][:, ::-1]
    for side, nz in (('top', 1.0), ('bottom', -1.0)):
        mat = g.material(dict(
            name=f'PCB_{side}', doubleSided=True,
            pbrMetallicRoughness=dict(
                baseColorTexture=dict(index=g.image(tex_dir / f'{side}_basecolor.png')),
                metallicRoughnessTexture=dict(index=g.image(tex_dir / f'{side}_mr.png')),
                metallicFactor=1.0, roughnessFactor=1.0)))
        m = sub[side]
        nrm = np.tile([0.0, 0.0, nz], (len(m.vertices), 1))
        pos, nrm, idx = gltf_arrays(m, nrm)
        x_mm, y_mm = m.vertices[:, 0] * 1000, m.vertices[:, 1] * 1000
        uv = np.stack([(x_mm + tex_half) / (2 * tex_half),
                       (tex_half - y_mm) / (2 * tex_half)], 1)
        g.mesh(f'PCB_{side}', pos, nrm, idx, mat, uv)
    edge_mat = g.material(dict(name='PCB_edge', pbrMetallicRoughness=dict(
        baseColorFactor=[*srgb_to_linear([0.20, 0.22, 0.30]).tolist(), 1.0],
        metallicFactor=0.0, roughnessFactor=0.7), doubleSided=True))
    pos, nrm, idx = gltf_arrays(trimesh.graph.smooth_shade(sub['edge'], angle=SMOOTH_ANGLE))
    g.mesh('PCB_edge', pos, nrm, idx, edge_mat)

    # ---- 零件:每個 OBJ 材質一個網格 ----
    # 第一階段:清理重疊平面,收集清理過的平面
    meshes, all_planes = OrderedDict(), []
    for name, tris in groups.items():
        if name == BOARD_MTL or name in fix.get('drop', []):
            continue
        mesh, planes = clean_overlapping_planes(to_mesh(V, tris, center))
        if planes:
            print(f'  {name}: 清理 {len(planes)} 個重疊的平面')
        meshes[name] = mesh
        all_planes += planes
    # 第二階段:把被蓋住的印字抬出來,輸出
    for name, mesh in meshes.items():
        n = lift_markings(mesh, all_planes)
        if n:
            print(f'  {name}: 抬升 {n} 個印字三角形')
        m = mats.get(name, {})
        mat = g.material(pbr_for(name, m.get('Kd', [0.5] * 3), m.get('Ks')))
        pos, nrm, idx = gltf_arrays(fix_and_shade(mesh))
        g.mesh(f'parts_{name}', pos, nrm, idx, mat)

    # 修正檔:換掉壞掉的零件模型
    if fix.get('drop'):
        print('  刪除壞掉的零件材質群組:', ', '.join(fix['drop']))
    for a in fix.get('add', []):
        mat = g.material(dict(name=a['name'], pbrMetallicRoughness=dict(
            baseColorFactor=[*srgb_to_linear(a['color']).tolist(), 1.0],
            metallicFactor=float(a.get('metal', 0)), roughnessFactor=float(a.get('rough', 0.5)))))
        m = rounded_prism(*a['box'], r=a.get('round', 0.0))
        pos, nrm, idx = gltf_arrays(trimesh.graph.smooth_shade(m, angle=SMOOTH_ANGLE))
        g.mesh(f'parts_{a["name"]}', pos, nrm, idx, mat)
        print(f'  加入 {a["name"]}:{a["box"]}')

    g.save(out)
    print('輸出', out, f'{Path(out).stat().st_size / 1e6:.1f} MB')


if __name__ == '__main__':
    main()
