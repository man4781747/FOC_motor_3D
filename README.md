# EasyEDA PCB → 3D 模型(GLB)

把 EasyEDA 匯出的 **OBJ 3D 模型**與 **Gerber 製造檔**合成一個看起來像 EasyEDA 網頁 3D 預覽的模型,
輸出成標準 glTF(`.glb`)。產生的模型可以:

- 直接在瀏覽器裡轉動、縮放(`output/viewer/`,three.js)
- 匯入 Blender(File → Import → glTF 2.0)或任何支援 glTF 的工具

另外也用程式建了一顆 **2804 雲台無刷馬達**,並把它用銅柱裝到 PCB 背面,組成一個完整模型。

## 為什麼需要這個流程

EasyEDA 匯出的 OBJ 只有零件和板子的外形,而且有幾個問題:

| 問題 | 本專案的處理 |
|---|---|
| 板子只有一個單色材質,**沒有防焊、絲印、焊墊、走線、via** | 從 Gerber 產生正反兩面的板面貼圖,貼到板子上下表面 |
| 單位是 10 mil(0.254 mm),座標離原點很遠(約 X 4020, Y −3375) | 以板子中心為原點,換算成公尺(glTF 標準單位) |
| **法線不可靠**(零件擺放時頂點有轉、法線沒跟著轉),三角形繞序不一致 | 捨棄 OBJ 法線,用 trimesh 修正繞序並依角度重算平滑法線 |
| 含孔洞的平面(例如晶片頂面的印字)被三角化成**互相重疊的三角形**,渲染會閃爍 | 偵測並換成乾淨的凸包平面 |
| 晶片印字是比表面低 0.01 mm 的凹刻,清理平面後會被蓋住 | 把這類淺凹的印字抬到平面外 0.02 mm |
| 沒有物件分組,匯入後是一整塊 | 依材質拆成獨立網格(板子正面 / 背面 / 側邊 + 各零件材質) |

## 輸入檔案

放在專案根目錄,腳本直接讀 zip,不需要先解壓:

| 檔案 | 從 EasyEDA 匯出的方式 |
|---|---|
| `OBJ_PCB_*.zip` | 匯出 3D 模型 → OBJ(內含 `.obj` + `.mtl`) |
| `Gerber_*.zip` | 匯出 PCB 製造文件(Gerber),需含鑽孔檔 |

目前的板子:「小型獨立運作系統」,圓形 Ø40 mm、4 層板、板厚 2.54 mm(沿用 OBJ),主控 LQFP64。

## 使用方式

需要 [uv](https://docs.astral.sh/uv/)。每支腳本開頭都有 inline 相依宣告,`uv run` 會自動安裝。

```bash
# 1. Gerber -> 板面貼圖(output/textures/)
uv run scripts/make_textures.py Gerber_*.zip output/textures

# 2. OBJ + 貼圖 -> PCB 模型
uv run scripts/build_glb.py OBJ_PCB_*.zip output/textures output/pcb.glb

# 3. 2804 馬達本體
uv run scripts/build_motor.py output/motor_2804.glb

# 4. 馬達裝到 PCB 背面(需要第 1、2 步的結果)
uv run scripts/build_assembly.py output/pcb.glb output/textures output/pcb_motor_assembly.glb

# 5. 轉成網頁檢視器用的 .gltf.json(部分網頁託管不接受 .glb)
for m in pcb motor_2804 pcb_motor_assembly; do
  uv run scripts/glb_to_gltf_json.py output/$m.glb output/viewer/$m.gltf.json
done
```

換一版設計時,重新從 EasyEDA 匯出兩個 zip,再依序跑一次即可。

### 在瀏覽器檢視

GLTFLoader 需要透過 HTTP 讀檔,不能直接雙擊開 `index.html`:

```bash
uv run python -m http.server 8080 -d output/viewer
```

開啟 <http://localhost:8080>,左上角可切換三個模型,也可以直接用網址指定:`#pcb`、`#motor`、`#assembly`。

- 拖曳:自由旋轉(TrackballControls,不受 ±90° 限制)
- 滾輪 / 雙指:縮放;右鍵:平移
- 下方按鈕:斜視、正面、背面、側視、自動旋轉

## 檔案結構

```
scripts/
  make_textures.py      Gerber -> 正反面 basecolor 與 metallicRoughness 貼圖(gerbonara + cairosvg)
  build_glb.py          OBJ + 貼圖 -> pcb.glb;也提供共用的 GLB 寫入器 / 網格處理函式
  build_motor.py        用 manifold3d(CSG)建 2804 馬達,尺寸都是檔案開頭的常數
  build_assembly.py     讀入 pcb.glb,加上翻轉後的馬達、銅柱、螺絲、焊到 U/V/W 的出線
  build_gimbal.py       雙軸雲台組合 + 3 個列印件 STL
  build_print_step.py   3 個列印件 STEP(CadQuery),給 CAD 軟體修改
  glb_to_gltf_json.py   GLB -> 單檔 glTF JSON(base64 內嵌)
output/
  pcb.glb                    PCB(約 8.5 MB)
  motor_2804.glb             馬達本體(約 0.4 MB)
  pcb_motor_assembly.glb     組合模型(約 9.1 MB)
  textures/                  板面貼圖(4096 x 4096)
  gimbal_2axis.glb           雙軸雲台(約 10 MB)
  print/                     列印件 STL(已擺好列印方向)
  step/                      列印件 STEP(組裝位置座標)
  viewer/                    three.js 檢視頁 + 四個 .gltf.json 模型
```

## 模型規格

### PCB 貼圖顏色

防焊深藍、防焊下的銅較亮(看得出走線)、開窗處沉金、白色絲印、鑽孔深色。
顏色常數在 `make_textures.py` 開頭(`MASK`、`GOLD`、`SILK` 等),要換防焊顏色改那裡即可。

### 2804 雲台馬達(依賣家圖片)

| 部位 | 尺寸 |
|---|---|
| 外形 | Ø34.5 × 15 mm,上方轉子、下方 2.5 mm 定子座,中間一圈細縫 |
| 中空軸 | Ø6.5 孔,內有軸承環 |
| 頂部螺孔 | 4 × M2.5,孔位圓 Ø16;另有 2 個小孔 |
| 底部螺孔 | 4 × M2.5,孔位圓 Ø19 |
| 出線 | 黑 / 紅 / 白 3 條,約 30 mm |

定子座高度、小孔位置、出線方向是依照片目測。

### 組合方式

- PCB 的 4 個 NPTH 安裝孔 (0, ±9.5)、(±9.5, 0) 剛好是 M2.5、孔位圓 Ø19,與馬達底部螺孔相同
- 馬達翻轉朝下、轉 45° 對孔,裝在 **PCB 背面**
- 4 支 M2.5 六角銅柱(對邊 4 mm、長 6.5 mm),正面以 M2.5 盤頭螺絲固定
- 出線焊到背面的 U / V / W 焊墊;焊墊位置由背面貼圖自動偵測
- 線色對應(白 / 紅 / 黑)只是為了走線不交叉,**不代表實際相序**

## 已知限制

- 零件的金屬材質是依 MTL 顏色推估的(中性灰且反光強 → 金屬),少數零件可能不準
- `build_glb.py` 的重疊平面清理以凸包取代,若零件有凹形且重疊的平面,外形可能被填平
- 馬達只有外觀,沒有內部繞組與磁鐵

## 雙軸雲台(進行中)

`scripts/build_gimbal.py` 用兩組「PCB + 2804 馬達」模組組成水平(Pan)+ 俯仰(Tilt)雲台,
並輸出三個 3D 列印結構件的 STL(已擺好列印方向)到 `output/print/`。

```bash
uv run scripts/build_gimbal.py output/pcb_motor_assembly.glb output/gimbal_2axis.glb
uv run scripts/glb_to_gltf_json.py output/gimbal_2axis.glb output/viewer/gimbal_2axis.gltf.json
uv run scripts/build_print_step.py          # 列印件 STEP(output/step/),給 CAD 軟體修改
```

安裝方式(詳見腳本開頭說明):

1. **底座**:八角底板 + 4 根凸柱。水平模組 PCB 正面朝下放在凸柱上,M2.5×12 螺絲從底板下方沉頭孔鎖進模組銅柱。馬達朝上,定子與 PCB 固定,馬達線不會被扭到
2. **框架**(平台 + 立板一體):平台用 4 顆 M2.5 從上方鎖進水平馬達轉子面(Ø16)
3. **俯仰模組**:PCB 正面朝立板,以 M2.5 公母銅柱 8 mm 連接,螺絲從立板外側鎖入
4. **相機支架**:4 顆 M2.5 從外側鎖進俯仰馬達轉子面;相機用 4 顆 M2 鎖在前板

| 列印件 | 尺寸 (mm) | 列印方向 | 懸空 |
|---|---|---|---|
| `base.stl` | 62.8 × 62.8 × 12 | 底板朝下 | 72 mm²(沉頭孔頂的小橋接) |
| `frame.stl` | 89 × 64 × 66 | 平台朝下 | 0.2 mm² |
| `camera_bracket.stl` | 36 × 39 × 42 | 法蘭朝下 | 23 mm²(M2 小孔頂部) |

### 在 CAD 軟體修改列印件(STEP)

`output/print/*.stl` 是三角網格,CAD 軟體讀得到但幾乎不能改。要在 FreeCAD、SolidWorks、Fusion 裡修改,
請用 `output/step/*.step`:

- 由 `scripts/build_print_step.py` 用 CadQuery(OpenCASCADE)產生,是真正的實體:
  只有平面、圓柱面(孔)與圓錐面(凸柱 45° 段),沒有自由曲面
- 尺寸常數直接從 `build_gimbal.py` 匯入,改尺寸只要改那邊,兩支腳本重跑後 STL 與 STEP 一致
- 座標為**組裝位置**(mm、Z 朝上、水平軸 = Z 軸),方便和其他零件對位;列印方向請用 STL
- 腳本會把 STEP 與 manifold3d 版本互相比對(最大表面偏差 > 0.05 mm 就報錯),
  目前三件偏差 0.010–0.029 mm,只是多邊形近似圓的誤差
- 匯入後沒有特徵樹(dumb solid):FreeCAD 可直接加特徵、布林運算;SolidWorks 可用直接編輯
  (Move / Delete Face),Professional / Premium 版可用 FeatureWorks 辨識出孔、擠出等特徵
- 注意:在 CAD 裡做的修改不會回寫到 Python 腳本

## 交接:待辦事項

上一台電腦資源不足,以下工作尚未完成(2026-11-25):

- [x] **干涉檢查**(2026-11-25 完成):靜止、水平 60°、俯仰 ±60° / ±90°、水平 45° + 俯仰 90° 共 7 個姿勢,
  門檻 0.2 mm 下皆無干涉。做法:列印件先 `mesh.merge_vertices(merge_tex=True, merge_norm=True)` 恢復水密,
  再用 `trimesh.proximity.signed_distance` 檢查其他零件的頂點**加表面均勻取樣點**(只查頂點會漏掉穿過大平面
  中間的情況,例如凸柱頂進 PCB 板面)。最貼近處:水平模組一顆正面零件(`parts_mtl15`)陷入底座凸柱 0.16 mm,
  低於門檻;俯仰 -90° 時鏡頭朝下,離平台約 1.8 mm
- [x] **檢視頁雲台規格卡**(`output/viewer/index.html` 的 `MODELS.gimbal`):俯仰軸高 76 mm;
  水平軸說明「定子 + PCB 固定在底座,轉子帶動平台」;俯仰滑桿(`#tilt`)放寬到 ±90°
- [ ] 檢視頁改完後重新發布網頁(Artifact)
- [ ] 實際列印前,用切片軟體再確認一次 `output/print/` 的 STL

注意事項:

- 本機 headless Chrome 截圖時,10 MB 以上的模型常常來不及載入;改用 Python(trimesh + matplotlib)
  從 GLB 直接渲染來驗證比較可靠,渲染時記得套用節點變換(`scene.dump()` 或 `scene.graph`)
- 推送需使用 `man4781747` 帳號(`gh auth switch --user man4781747`)
