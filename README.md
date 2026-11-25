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

開啟 <http://localhost:8080>,左上角可切換模型,也可以直接用網址指定:`#pcb`、`#motor`、`#assembly`、`#gimbal`、`#car`。

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
  hardware.py           五金模型:有螺紋的盤頭 / 沉頭螺絲、公母銅柱、熱熔螺母、螺帽
  build_car.py          舵輪平台車組合 + 6 種列印件 STL / STEP(CadQuery)
  glb_to_gltf_json.py   GLB -> 單檔 glTF JSON(base64 內嵌)
output/
  pcb.glb                    PCB(約 8.5 MB)
  motor_2804.glb             馬達本體(約 0.4 MB)
  pcb_motor_assembly.glb     組合模型(約 9.1 MB)
  textures/                  板面貼圖(4096 x 4096)
  gimbal_2axis.glb           雙軸雲台(約 10 MB)
  print/                     列印件 STL(已擺好列印方向)
  step/                      列印件 STEP(組裝位置座標)
  car_swerve.glb             舵輪平台車(約 12 MB)
  car/print/, car/step/      平台車列印件 STL / STEP
  viewer/                    three.js 檢視頁 + 五個 .gltf.json 模型
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
- 4 支 M2.5 公母六角銅柱(對邊 4 mm、長 6.5 mm),公頭 3 mm 鎖進馬達底部螺孔;正面以 M2.5 × 6 盤頭螺絲
  穿過 PCB 鎖進銅柱母牙
- 出線焊到背面的 U / V / W 焊墊;焊墊位置由背面貼圖自動偵測
- 線色對應(白 / 紅 / 黑)只是為了走線不交叉,**不代表實際相序**

## 已知限制

- 零件的金屬材質是依 MTL 顏色推估的(中性灰且反光強 → 金屬),少數零件可能不準
- `build_glb.py` 的重疊平面清理以凸包取代,若零件有凹形且重疊的平面,外形可能被填平
- 馬達只有外觀,沒有內部繞組與磁鐵

## 水情盒主板(Main_PCB V2.4)

另一塊板子,輸入檔放在 `水情盒/`。100 × 100 mm 方形、2 層板(沒有背面絲印),
板子不在 Gerber 原點,所以貼圖要加 `--fit-outline` 依板框決定範圍:

```bash
uv run scripts/make_textures.py 水情盒/Gerber_*.zip output/water_box/textures --fit-outline
uv run scripts/build_glb.py 水情盒/OBJ_PCB_*.zip output/water_box/textures output/water_box/main_pcb_v2.4.glb
uv run scripts/glb_to_gltf_json.py output/water_box/main_pcb_v2.4.glb output/viewer/water_box_main_pcb.gltf.json
```

檢視頁網址加 `#waterbox`。

- `--fit-outline`:取板框中較長的線段算外框(板框層上有一個 1.2 × 8 mm 的小框不是板邊,會被略過),
  貼圖中心與半寬寫進 `textures/texture.json`,`build_glb.py` 讀它來算 UV;沒有這個檔時沿用 Ø40 板的 ±20.32 mm
- 鑽孔檔裡的 G85 槽孔(右下角 3 個橢圓孔)gerbonara 不支援,`make_textures.py` 先把它展開成密排圓孔
- 缺少的 Gerber 層(例如背面絲印)當成空白
- 已用貼圖疊零件頂點驗證對位:安裝孔、槽孔、焊墊都對齊

### RS485 擴展板(RS485_V1)

輸入檔放在 `rs485擴展板/`。100 × 100 mm、2 層板,同樣要加 `--fit-outline`:

```bash
uv run scripts/make_textures.py rs485擴展板/Gerber_*.zip output/rs485/textures --fit-outline
uv run scripts/build_glb.py rs485擴展板/OBJ_PCB_*.zip output/rs485/textures output/rs485/rs485_v1.glb \
    fixes/rs485_v1.json                    # 換掉 EasyEDA 匯出壞掉的電感模型
uv run scripts/glb_to_gltf_json.py output/rs485/rs485_v1.glb output/viewer/rs485_v1.gltf.json
```

檢視頁網址加 `#rs485`。

**電感模型修正**:EasyEDA 匯出的 OBJ 裡,功率電感(焊墊中心約 (18.3, 5.3))本身就是壞的
(頂面立起來、本體被拉成大三角形,端子偏離焊墊 1.8 mm),不是轉檔造成的。`build_glb.py` 第 4 個參數可以給
修正檔,`fixes/rs485_v1.json` 刪掉它的材質群組(mtl17 本體、mtl18 端子),依焊墊與絲印換成
10 × 9.1 × 4 mm 的遮蔽式電感(圓角黑色本體 + 兩端 L 形端子)。尺寸是依焊墊推估,
知道實際料號後改 JSON 的 `box` 即可。其他板子遇到壞掉的零件也可以用同樣方式處理。

### SEN0676 雷達液位感測器模組

`scripts/build_sen0676.py`:DFRobot SEN0676(80 GHz 雷達液位感測器)。輸入檔放在 `sen0676/`
(官方 STEP zip 與尺寸圖 PDF)。

```bash
uv run scripts/build_sen0676.py sen0676/SEN0676_*_stpfile_V1.0.zip output/sen0676/sen0676.glb \
    --step-out output/sen0676/sen0676_module.step
uv run scripts/glb_to_gltf_json.py output/sen0676/sen0676.glb output/viewer/sen0676.gltf.json
```

檢視頁網址加 `#sen0676`。

- **外殼**直接用官方 STEP(Pro/E 匯出,料號 LDRR04C-01):62 × 66.9 × 53.65 mm,長度與尺寸圖一致;
  顏色改成實品的黑色(STEP 原色是淺灰)。
  STEP 裡**只有外殼**,沒有電路板
- **電路板依尺寸圖第 1 頁重建**:35 × 35 mm、4 × Ø2.6 孔(孔距 31)、左側 5P 排針
  (TX / RX / DE / GND / VCC,VCC 離下緣 9 mm)、右側 4P / 3P 焊孔、雷達晶片、黑色晶片、被動元件、Ø30 絲印圓。
  元件位置依圖面比例量(約 ±0.3 mm),元件高度、板厚 1.6、板子顏色是估計值
- 外殼 35.28 方槽的內角約 R1.7,方角板子放不進去,所以板子四角做 **R2** 圓角
- **安裝**:雷達晶片朝鏡頭,元件面貼在 4 根凸柱頂(Ø2.8 × 0.78 沉孔、Ø2 × 5.5 導孔),
  **M2.5 × 6 盤頭自攻** × 4 從背面鎖入(鎖入 4.4 mm);5P 排針在背面,朝外殼開口
- 腳本試 0 / 90 / 180 / 270° 四個角度,用 OpenCASCADE 布林運算確認零件與外殼重疊體積為 0;
  四個角度都可以,採用 0°(尺寸圖沒有標板子在外殼裡的方向)
- `--step-out` 另存整個模組的 STEP(外殼 + 重建的電路板),可在 CAD 軟體使用

### KSS AG-12 電纜固定頭(M12 × 1.25)

`scripts/build_ag12.py`:KSS AG-12 尼龍電纜固定頭(緊迫頭),5 個零件:螺帽 AGL-12、墊片 AGW-12、
本體 AGB(與夾緊爪一體)、夾緊圈 AGS(TPE)、迫緊頭 AGN。官網沒有 CAD 檔,尺寸來源:

| 尺寸 | 數值 | 來源 |
|---|---|---|
| 螺紋 / 螺紋長 L1 / 開孔 | M12 × 1.25 / 7.9 / Ø12 +0.2 | 官網、型錄 0604 |
| 線徑 | 3 ~ 6.5 | 官網 |
| 螺帽 AGL-12 | 對邊 17.0、高 5.2 | 型錄 0611 |
| 墊片 AGW-12 | 外徑 17.0(厚度未公開,估 1.0) | 官網 AGW 頁 |
| 本體六角 / 迫緊頭外徑 / 總高 | 15.2 / Ø15.5 / 30(螺紋末端到頂) | 實測 |
| 凸緣、六角高、迫緊頭圓角與出線孔、內部的爪與夾緊圈 | — | 依產品照片估計 |

```bash
uv run scripts/build_ag12.py output/ag12/ag12.glb --panel 3.09    # --panel:面板厚度(螺帽位置)
uv run scripts/glb_to_gltf_json.py output/ag12/ag12.glb output/viewer/ag12.gltf.json
```

檢視頁網址加 `#ag12`(有爆炸圖與零件拆分)。注意:螺紋 L1 只有 7.9,裝在 3.09 mm 的水情盒牆上時,
扣掉墊片後面板內側只剩 3.8 mm 螺紋,AGL-12 螺帽(高 5.2)只咬合 3.8 mm。

### 水情盒組合(防水盒 + 底板)

`scripts/build_water_box.py` 把 `水情盒/底板_3.3mf` 放進 RS PRO 1862264 防水盒底部,用 4 顆螺絲鎖住:

```bash
uv run scripts/build_water_box.py 196136506-104-1862264.zip 水情盒/底板_3.3mf output/water_box/water_box.glb \
    --rs485 output/rs485/rs485_v1.glb \     # RS485 擴展板(選用,需先建好 rs485_v1.glb)
    --main output/water_box/main_pcb_v2.4.glb \ # 水情盒主板,裝在最上層(選用,需要 --rs485)
    --sen0676 output/sen0676/sen0676.glb       # SEN0676 雷達,裝在 -Y 牆的 Ø40 孔(選用)
uv run scripts/glb_to_gltf_json.py output/water_box/water_box.glb output/viewer/water_box.gltf.json
```

檢視頁網址加 `#waterboxasm`;拉爆炸圖或在零件拆分取消勾選盒蓋就能看到內部。

- 盒底內側有 4 根凸柱(Ø7.9、中心 Ø2.4 自攻導孔、比內底高 5.9 mm),位置 (±62.75, ±61.4);
  底板的 4 個 Ø4.0 孔剛好對準(誤差 0.02 mm),四角與兩側缺口避開盒蓋螺絲柱
- 螺絲:**ST3.0 × 6 盤頭自攻** × 4,穿過 2 mm 底板,鎖進凸柱 4 mm(導孔深 5.9 mm,不會頂到底)
- 腳本會自動找孔位與凸柱、檢查底板與盒身不干涉(截面重疊面積 0)、螺絲長度不超過導孔深度
- 底板頂面到盒緣約 49 mm
- 底板上有兩組 95.0 × 95.0 mm 的 Ø3.0 孔(中心離底板中心 40.0 / 50.7 mm),和 100 × 100 板子的四角孔
  (95.1 × 95.1)只差 0.1 mm
- **RS485 擴展板**(`--rs485`):裝在較置中的那組孔,**M3 × 12 尼龍母母六角柱** × 4(對邊 5.5)。
  下方 **M3 × 6 盤頭** × 4 從底板背面鎖入,要在底板放進盒子前先鎖好,螺絲頭離盒底 3.7 mm;
  上方用 **M3 × 30 + 6 金屬公母六角柱** × 4 取代螺絲:公頭穿過板子鎖進尼龍柱 3.5 mm,母牙端朝上
  (柱頂 z=56.9,離盒緣 4.3 mm)。
- **水情盒主板**(`--main`):裝在最上層,四角 Ø3.0 孔(同樣 95.1 × 95.1)放在金屬六角柱頂端,
  **M3 × 6 盤頭** × 4 從主板正面鎖進金屬柱母牙(鎖入 3.5 mm)。轉角與 RS485 板相同(`MAIN_YAW`)。
  主板正面 z=59.5,比盒緣(z=61.2)低 1.7 mm,較高的零件會伸進盒蓋裡;最高零件 z=70.5,離盒蓋內面 7.8 mm
- 由下往上的組裝順序:底板背面螺絲 + 尼龍柱 → 底板鎖進盒子 → RS485 板 + 金屬柱 → 主板 + 螺絲 → 盒蓋
- 兩塊板子都轉 270°(`PCB_YAW`,主板 `MAIN_YAW` 跟著),RS485 的 4 組端子朝 -X(盒子中央方向);
  這個方向金屬柱旁最近的零件離柱中心 3.8 mm(六角柱角半徑 3.18 mm),間隙約 0.6 mm
- **SEN0676 雷達**(`--sen0676`):盒身 **-X 短邊**的牆(厚 3.09)開 **Ø48.5 孔**(G1½ 螺紋外徑 47.8 要穿過),
  孔中心在 (-X, +Y) 角落盒蓋螺絲柱邊緣往 -Y **70 mm**、盒身上緣往下 **27 mm**(y=-3.77、z=35.85)。
  雷達從盒內裝、鏡頭朝外,六角法蘭一組對邊平行於盒底(平邊離盒內底 0.36 mm):
  1. 螺紋根部凹槽(Ø44.6 × 2.5)套 **O 型環 GB/T 3452.1 43.7 × 2.65**(拉伸 2.1%),法蘭把它壓在牆內面
     (平均壓縮 25%,法蘭面離牆 1.97 mm;O 型環外緣 Ø49.8,壓在孔邊外 0.67 mm 寬)
  2. 螺紋穿出牆外,**KSS AGJL-112 尼龍螺帽**(PF 1½",對邊 58.5、高 9.8,黑)從外面鎖緊,螺帽一組對邊也平行盒底;
     螺帽鎖在完整螺紋上,雷達凸出牆外 21.6 mm
  - 盒身牆有約 0.9° 拔模角:外牆在螺帽範圍內上下差 0.49 mm,螺帽貼在最凸處;
    內牆讓 O 型環壓縮率上緣 11%、下緣 39%(上緣偏低,要更確實可改用較粗的 O 型環或平墊片)
  - 孔頂(z=60.1)離密封條底只有 0.45 mm,開孔時要小心不要破到密封條溝
  - 離電路板六角柱 69.8 mm、雷達頂離盒蓋內面 11.4 mm;**唯一干涉:底板 -X 邊緣**(x=-102.5)壓在法蘭下方,
    需在底板開缺口。常數在 `SEN_HOLE_*`、`SEN_ORING*`、`SEN_NUT`
- 底板的 Ø3.0 孔剛好等於 M3 外徑,螺絲會直接咬進列印件;想要好穿可以把孔擴到 Ø3.2–3.4
- RS485 板最高的零件到盒緣約 26 mm;板子背面的插件腳離底板約 11 mm

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
4. **相機支架**:4 顆 M2.5 從外側鎖進俯仰馬達轉子面;前板正面有 4 個水滴形螺絲柱,背面埋 M3 熱熔螺母,
   相機用 4 顆 M3 × 6 從相機板背面鎖入(相機板孔需能過 M3)

檢視頁的雲台(`#gimbal`)有爆炸圖與零件拆分(14 種零件):底座往下、水平轉動部分往上,
俯仰模組與相機支架沿俯仰軸往外拆;爆炸圖下仍可拖水平 / 俯仰滑桿。

| 列印件 | 尺寸 (mm) | 列印方向 | 懸空 |
|---|---|---|---|
| `base.stl` | 62.8 × 62.8 × 12 | 底板朝下 | 72 mm²(沉頭孔頂的小橋接) |
| `frame.stl` | 89 × 64 × 66 | 平台朝下 | 0.2 mm² |
| `camera_bracket.stl` | 36 × 43 × 42 | 法蘭朝下 | 0.4 mm² |

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

## 舵輪平台車(8 × 2804 模組)

`scripts/build_car.py` 用 8 組「PCB + 2804 馬達」模組做 4 輪舵輪(swerve)平台車:每個角 1 顆驅動 + 1 顆轉向,
可前進、橫移、原地旋轉。輸出組合模型與 6 種列印件(STL 已擺好列印方向,STEP 可在 CAD 修改)。

```bash
uv run scripts/build_car.py output/pcb_motor_assembly.glb output/car_swerve.glb \
    --wheel 雲台輪胎.3MF --box 196136506-104-1862264.zip            # 防水盒(選用);列印件 -> output/car/{print,step}/
uv run scripts/glb_to_gltf_json.py output/car_swerve.glb output/viewer/car_swerve.gltf.json
```

檢視頁網址加 `#car`,可切換「前進 / 橫移」與「原地旋轉」模式、拖轉向與車輪滑桿,或按「播放」讓車輪轉動。

- **爆炸圖**:滑桿 0–100%,零件依組裝順序拆開(上蓋、橋架往上,舵輪往下,驅動模組與車輪沿輪軸往外)。
  每個模組拆成馬達、銅柱、PCB 三件;馬達本身不可拆,轉子與定子一起移動
- **零件拆分**:清單依列印件 / 馬達 / 電路板 / 電子 / 五金分組並列出數量;點名稱只看該零件(其他半透明),
  勾選框可隱藏零件;也可以直接在畫面上點零件選取,右上角顯示名稱、數量與 STL 檔名

零件資料與爆炸位移寫在 GLB 節點的 `extras`(各腳本的 `PARTS`、`EXPLODE`,由 `build_gimbal.make_extras` 產生),
檢視頁只讀這些資料;雙軸雲台(`#gimbal`)也有這兩個功能,其他模型加上同樣的 `extras` 即可使用。

每個角由下往上(詳見腳本開頭說明):

1. **車輪**:白色一體式杯狀輪框 `雲台輪胎.3MF`(`--wheel`,外部提供),Ø52 × 20,杯底 4 mm、杯內 Ø46 × 16,
   包住會轉的馬達轉子。外圈 5 道半圓凹槽(槽底 Ø50.86、圓弧半徑約 1.4)各套一條**黑色 O 型環**當胎皮:
   GB/T 3452.1 **48.7 × 2.65**(NBR),套上拉伸 4.4%、截面變細到約 2.6 mm,輪胎外徑約 **Ø56**;
   O 型環著地,輪框離地 2 mm。4 顆 M2.5 × 5 盤頭從外側 Ø6 沉頭孔鎖進轉子面 Ø16(鎖入約 3.4 mm,輪轂馬達直驅)。
   轉向軸通過輪寬中心面,原地轉向時輪胎不拖地。杯內壁離馬達出線 3.9 mm、輪頂離轉向叉法蘭 3 mm
2. **驅動模組**:PCB 正面朝轉向叉立板,M2.5 公母銅柱 8 mm 連接(同雲台俯仰模組)
3. **轉向叉**:L 形,法蘭用 4 顆 M2.5 從下方鎖進轉向馬達轉子面
4. **轉向模組**:轉子面朝下穿過底板 Ø40 孔,PCB 正面朝上,由**橋架**的 4 根凸柱頂住,M2.5 從橋架上方鎖進模組銅柱
5. **車身**:黑色底板 + 白色上蓋。橋架腳與上蓋螺絲柱埋 M3 熱熔螺母,M3 × 8 沉頭從底板下方鎖入;
   上蓋四角各 4 個孔對準橋架螺絲,不用拆蓋就能鎖。上蓋頂部中央有 4 個 M3 熱熔螺母,可直接鎖雙軸雲台底座;
   電池(3S 1300 mAh 佔位)放車尾

| 列印件 | 尺寸 (mm) | 數量 | 材料 / 列印方向 |
|---|---|---|---|
| `top_cover` | 250 × 200 × 32 | 1 | PLA,頂面朝下(頂緣 45° 倒角) |
| `bottom_plate` | 250 × 200 × 4 | 1 | PLA(黑),平放 |
| `bridge` | 64 × 64 × 28 | 4 | PLA,頂板朝下 |
| `fork` | 52 × 44 × 57 | 4 | PLA,法蘭頂面朝下 |
| 輪框 `雲台輪胎.3MF` | Ø52 × 20 | 4 | PLA 白,外部提供的列印檔 |
| O 型環 48.7 × 2.65 | — | 20 | 外購(每輪 5 條) |

整車 250 × 200 × 101 mm(車身 36 mm,下方是舵輪)。上蓋與底板需要 256 mm 以上的熱床。

**防水盒(選用)**:`--box` 指定 RS PRO 1862264 防水盒的 3MF(可直接給 TraceParts 下載的 zip),
盒子(225 × 175 × 81 mm)放在上蓋頂面、長邊沿車頭方向,**只是放著,沒有固定**;加上後整車高約 182 mm。
盒子會蓋住上蓋四角的橋架螺絲維修孔與中央的雲台孔位。盒子模型檔是外部下載的,不放進 repo。
上蓋倒過來印時,頂部 4 個熱熔螺母盲孔的孔頂是 Ø4 的小橋接(共約 50 mm²),不需要支撐。

### 五金規則(雲台與平台車共用)

| 鎖附位置 | 規格 |
|---|---|
| 馬達轉子面、馬達底部(本身有牙) | M2.5 |
| PCB 模組(安裝孔只能過 M2.5) | M2.5 公母銅柱:模組內 6.5 mm,鎖結構件的間隔銅柱 8 mm(公頭 6 mm 穿過 PCB) |
| 其他鎖進列印件的地方 | **M3 埋頭型熱熔螺母** M3 × 5、外徑 4.5,列印件開 Ø4.0 × 深 6 |

五金模型在 `scripts/hardware.py`,螺絲有螺紋、長度依實際鎖附計算(例如橋架螺絲 M2.5 × 14 鎖進模組銅柱約 3.5 mm)。

限制與待確認:

- **轉向的走線**:驅動模組的線要從會轉的叉接到車身,目前沒有滑環,轉向角度要限制在約 ±180° 內
  (或之後加滑環、改成無線供電通訊)
- 2804 是雲台馬達,扭力小,適合輕量車;整車重量建議控制在 1.5 kg 以內
- 電池只是尺寸佔位,實際尺寸不同時改 `BATTERY` 常數

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
