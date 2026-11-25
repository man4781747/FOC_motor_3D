# /// script
# requires-python = ">=3.10"
# ///
"""把 GLB 轉成單一檔案的 glTF JSON(二進位資料以 base64 data URI 內嵌)。

用法:
    uv run scripts/glb_to_gltf_json.py <輸入.glb> <輸出.gltf 或 .json>

用途:某些網頁託管只允許 .json 等文字格式時使用;three.js GLTFLoader 可直接載入。
"""
import base64
import json
import struct
import sys

src, dst = sys.argv[1:3]
data = open(src, 'rb').read()
magic, version, total = struct.unpack_from('<4sII', data, 0)
assert magic == b'glTF' and version == 2, '不是 glTF 2.0 的 GLB 檔'
off, gltf, binchunk = 12, None, b''
while off < total:
    length, kind = struct.unpack_from('<I4s', data, off)
    chunk = data[off + 8: off + 8 + length]
    if kind == b'JSON':
        gltf = json.loads(chunk)
    elif kind == b'BIN\0':
        binchunk = chunk
    off += 8 + length
gltf['buffers'][0]['uri'] = 'data:application/octet-stream;base64,' + base64.b64encode(binchunk).decode()
with open(dst, 'w') as f:
    json.dump(gltf, f, separators=(',', ':'))
print(dst, f'{len(open(dst, "rb").read()) / 1e6:.1f} MB')
