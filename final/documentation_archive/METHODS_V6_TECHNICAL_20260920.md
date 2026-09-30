# V6 方法詳細說明：病灶相關 AAL3 區域放射體學與融合預測

文件日期：2026-09-20。

本文件依據目前保存的 V6 程式、固定設定、特徵封存檔及已執行的模型紀錄撰寫。主體說明原始 V6 的影像處理、L／T／W 特徵、基準與融合模型、交叉驗證及不確定性估計；後段另述後續四任務比較及 SHAP 分析，以免把不同階段的操作混成同一個方法。本文提供詳細方法文字，不以既有高分倒推方法合理性，也不把未實作的程序寫成已完成。

依本次要求，所有圖與表均以 **pending** 保留位置；本文不嵌入既有圖檔或填入正式表格。已有可靠程式依據的設定仍完整寫在正文。影像取得、人工標註及倫理資訊若無法由本輪程式確認，另明確列為待研究團隊補充的資料。

## 1. 研究目的、方法定位與原稿的關係

V6 探討的問題是：在缺血性腦中風病人的擴散加權影像（diffusion-weighted imaging，DWI）中，使用病灶位置與 AAL3 解剖分區共同定義區域性影像特徵，是否能在既有臨床、病灶體積、解剖負荷與全病灶 radiomics 的基礎上，增加對結局的辨識能力。這個問題涉及「區域特徵是否提供額外資訊」，因此比較基準不能只有臨床變數，還必須納入已包含病灶大小、位置及全病灶特徵的模型。

V6 的出發點是，若直接在完整腦區中計算 radiomics，面積很小的病灶可能被大量非病灶組織稀釋。為檢驗這個假設，V6 同時建立三組互補特徵：L 著重病灶局部的一階強度統計及其與對側同源腦區的差異；T 在病灶與 AAL3 區域交集內計算 Original radiomics；W 則保留完整腦區 radiomics 的起點，再以該腦區的病灶占比加權。這些都是具體的特徵定義，不能僅用「融合特徵」一個名稱概括。

原稿的 F 是從完整 AAL3 腦區提取 radiomics，再將各腦區的特徵串接成病人向量；原稿也強調 F 的抽取不需先取得病灶遮罩。V6 的 L、T 與 W 均依賴病灶遮罩，因而属于病灶導向的區域放射體學延伸。V6 保留 DWI、AAL3、區域身分與可追溯特徵的概念，但改變了 ROI 或特徵加權方式。因此，V6 不能被描述為原稿 F 完全不變的實作，也不能承接原稿「不需病灶分割」的主張。

V6 中的「融合」還有兩種不同意思。第一種是影像與區域資訊共同決定特徵，例如病灶與 atlas 的交集；第二種是模型輸入或預測層級的融合，例如將臨床與多組影像特徵串接，或將兩個模型的機率加權。方法敘述需分別命名，避免把 ROI 的設計與分類模型的組合視為同一個步驟。

> 圖 1．V6 研究架構與原稿 F 的方法關係：**pending**。
>
> 表 1．原稿 F、V6 L／T／W 與後續 Fwhole 對照的定義：**pending**。

## 2. 分析單位、資料來源及研究樣本

### 2.1 病人層級的分析單位

本研究以病人為建模與評估單位。一位病人對應一列影像特徵及其臨床紀錄；模型不將同一病人的不同 voxel、腦區或影像切片當成獨立樣本。病人 ID 用於對齊資料、驗證唯一性及建立訓練／驗證分割，不作為模型輸入。

V6 延用 V5 凍結的 198 人影像 cohort。V5 又承接先前 V3 完成影像處理的個案。因此，198 人是這一條既有影像處理與品質檢查路線所保留下來的 cohort，不能直接視為原始招募總人數，也不能代表所有最初候選個案都通過了 V6 之前的檢查。完整納入與排除流程應由原始收案紀錄及 V3／V5 品質紀錄補足，不應由目前的 198 人特徵檔反推未保存的排除原因。

V6 不依結局或模型表現另行挑選病人；特徵抽取與矩陣組裝程式不讀取臨床結局。這表示自動特徵計算沒有使用 outcome 來挑 ROI 或挑病人，但不等於已證實原始人工分割、影像配準或人工 QC 的執行者對結局盲化。

### 2.2 結局標籤與有效樣本

原始 V6 的兩個任務為出院 mRS 與出院 NIHSS。mRS 使用資料表中的 `mRS_discharge`，大於 2 編碼為 1，否則為 0；NIHSS 使用 `NIHSS_discharge`，大於 4 編碼為 1，否則為 0。這是出院時點的二元分類，不能改寫為 90 天預後、長期功能恢復，或其他未在資料欄位中定義的時間點。

每個任務先保留同時有影像特徵與該任務非缺失標籤的病人，再進行分割。mRS 有效樣本為 191 人，其中 80 人為正類、111 人為負類；另外 7 位影像 cohort 病人缺失該標籤。NIHSS 有效樣本為 192 人，其中 78 人為正類、114 人為負類；另外 6 人缺失該標籤。標籤缺失不以訓練中位數或其他方式補值。不同任務保留的病人不一定完全相同，不能把四任務所有指標都寫成來自同一個無缺失的 198 人集合。

臨床資料沿用既有清理函式及原始資料表的歷史修正。V6 本次未重新定義 outcome，也未依新結果變更標籤。原始清理函式包含固定的資料欄位處理與少數歷史修正；這些操作應隨程式與來源檔保留，若正式稿要敘述其臨床理由，需由原始資料管理紀錄確認，而不能只從行列索引推定。

### 2.3 預測資訊的使用範圍

臨床 predictor 採保守的 core 集合，包括 `age`、`sex`、`AF`、`CHF`、`CAD`、`HTN`、`Lipid`、`DM`、`cancer`、`RecurSt`、`smoking` 與 `alcohol`，共 12 個原始欄位。程式將欄位轉為數值型態，保留資料來源的編碼，沒有新增 one-hot encoding。各二元欄位的 0／1 定義、缺失代碼及病史時點，仍應以臨床資料字典為準。

mRS／NIHSS 評分可作為結局，但不作為 predictor；排除範圍包括入院與出院總分、NIHSS 子項、分數變化及由其衍生的輸入。這項限制不能反向證明其餘欄位的時間一定正確：若要宣稱是在特定臨床時間點進行預測，仍須確認所有臨床欄位、DWI 與病灶遮罩在該時間點可取得。

### 2.4 尚須補充的臨床與影像來源資訊

收案機構與期間、原始納入排除條件、倫理審查及同意程序、掃描儀器與序列參數、DWI 的 b-value、原生 voxel spacing 與切片厚度、發病到影像檢查時間、病灶遮罩來源、標註者及一致性評估、配準演算法與人工 QC 規範，均需由原始研究紀錄提供。本文不以推測補上上述資訊。

> 圖 2．從原始收案至各任務有效 cohort 的流程圖：**pending**。
>
> 表 2．人口學、臨床特徵、缺失情形與結局分布：**pending**。
>
> 表 3．MRI 取得、病灶標註、配準與倫理資訊：**pending**；其中尚未確認的內容須由研究團隊補充。

## 3. 影像輸入、共同空間與品質檢查

### 3.1 V6 使用的影像檔

V6 使用既有的 MNI 空間 DWI `diff2MNI.nii.gz`、病灶遮罩 `lesion2MNI.nii.gz`、MNI152 2 mm brain support 及同格點的 AAL3 圖譜。這些檔案已在上游流程完成空間處理；V6 抽取程式本身沒有重新估計配準變換，也沒有以本次 outcome 指導影像配準。

抽取前檢查 DWI、病灶遮罩、brain support 與 AAL3 的影像尺寸、voxel spacing、原點與方向矩陣是否一致。尺寸必須完全相同，其他幾何量以絕對容許差 1e−5、相對容許差 0 檢查。這可以排除明顯的格點不一致，但不代表已證明個別病人的解剖配準正確；後者仍需影像品質與解剖對位證據。

本文所稱 voxel 數是此 2 mm MNI 網格上的計數，不是原生影像中獨立採樣點的數量。即使重採樣後 ROI 超過某個 voxel 閾值，也不能據此宣稱已克服原生厚切片、部分容積效應或小病灶纹理的不穩定。

### 3.2 遮罩與有效影像訊號

令第 i 位病人的 DWI 強度為 I_i(v)。病灶遮罩以輸入值大於 0.5 定義；AAL3 的背景 label 0 不列入腦區特徵。有效影像訊號定義為強度有限且大於 0 的 voxel，亦即同時排除 NaN、正負無窮及非正值。

V6 的局部病灶特徵使用「AAL3 腦區 ∩ 病灶 ∩ 有效影像訊號」作為實際計算支援。若某腦區有病灶，但交集內有效 voxel 數占原始病灶交集少於 90%，程式報錯，不將該區默默改當成無病灶。

V5 上游另要求 brain support 中至少 95% voxel 有有效 DWI 訊號；位於擴張 brain support 之外的病灶比例不可超過 50%。此處 brain support 以各軸半徑 2 voxel 的二值膨脹建立外圍容許帶。這些是實作中的固定 QC 門檻，不能被解讀為已完成完整的配準品質評分。

### 3.3 實際使用的 AAL3 區域

特徵欄位依實際 AAL3 影像中存在的正值 label 排序，並使用 label 對照檔保留解剖名稱。現有特徵封存檔包含 166 個區域。這個數字來自本次圖譜檔與特徵矩陣，不是直接沿用原稿互有差異的 160／166 敘述，也不表示所有區域在每位病人都可以算出完整 radiomics。

小區域或 coverage 不足的完整腦區仍保留欄位位置；特徵是否可用依固定規則處理，再由 training fold 的前處理決定哪些欄位在當次模型中能使用。程式沒有使用全 cohort 的 outcome 來選出「重要腦區」後才開始交叉驗證。

> 圖 3．共同影像格點、病灶與 AAL3 交集，以及品質檢查示意：**pending**。

## 4. 病人內 DWI 強度正規化

V6 固定採用 V5 的 `unclip_p99` 表徵。對每位病人，先在固定 brain support 與有效影像訊號的交集上計算 DWI 強度的第 1 與第 99 百分位數，記為 q1_i 與 q99_i，再以兩者間距作為尺度。這是每位病人自己影像內的正規化，不是用所有病人的平均或標準差建立一個跨病人的影像尺度。

對有效 voxel，其轉換可寫為：

$$
I_i^{*}(v)=100\max\left(\frac{I_i(v)-q1_i}{q99_i-q1_i},0\right).
$$

無效 voxel 在輸出陣列中設為 0。低於第 1 百分位數的值經下界截斷成 0；高於第 99 百分位數的值則可以大於 100。`unclip` 在此表示沒有把高值截到 100，並不代表完全不做任何截斷。

第 1／99 百分位數取自固定 brain support 中的有效訊號，`unclip_p99` 的數值尺度沒有用病灶位置排除參照組織。這與另一個歷史候選 `normal_tissue_unclipped` 不同，後者會排除膨脹病灶。V6 並未把兩種正規化混用，也未在各 outer test fold 看表現後改選正規化方式。

實作仍有一個需要保留的界線：共用的 V5 `normalize` 函式會同時計算多個候選影像與其 QC，包含排除膨脹病灶後的 normal-tissue 參照。函式要求 brain 參照及該 normal-tissue 參照各至少有 1,000 個有效 voxel，並檢查百分位距離不為 0。因此，雖然被 V6 採用的 `unclip_p99` 數值公式不依賴病灶定義，目前共用抽取程式仍會讀取病灶並執行相關檢查，不能宣稱整條實作已是無病灶遮罩的流程。

影像正規化與後續特徵表的 z-score 是不同步驟。前者只使用該病人的影像，不在病人間估計共同參數；後者則必須在每個 training fold 內估計。本方法同時使用這兩種處理，但作用層級與資料邊界不同。

> 圖 4．`unclip_p99` 的強度轉換與高強度訊號保留示意：**pending**。

## 5. 局部病灶、對側參照與共同可用性指標

### 5.1 局部病灶支援

令第 k 個 AAL3 區域為 A_k，第 i 位病人的病灶集合為 M_i，有效 voxel 集合為 V_i。原始病灶交集為 A_k∩M_i，實際供 L／T 計算的局部 ROI 為：

$$
\Omega_{ik}=A_k\cap M_i\cap V_i.
$$

「區域受到病灶影響」由 A_k∩M_i 是否非空決定；「能否抽取特徵」則還取決於有效 voxel 數、ROI 的空間分布，以及對側參照是否足夠。這兩個概念分開保留，不能把特徵算不出來的病人當成該區沒有病灶。

### 5.2 左右同源腦區的建立

V6 使用 AAL3 label 名稱末尾的 `_L` 與 `_R` 尋找左右同源區域。配對依解剖名稱進行，不以 label 編號的奇偶決定左右。沒有可對應名稱的中線區域或其他非成對區域，不強行指定對側。

對第 k 區的對側區域，先排除整位病人的病灶遮罩及其各軸半徑 2 MNI voxel 的二值膨脹範圍，再與有效 DWI 訊號取交集。左右兩側存在病灶時，排除操作仍使用完整病灶集合，因此對側病灶亦會被排除。2 voxel 是網格上的形態學設定，不應簡化成已驗證的固定生理距離或病灶周圍正常組織範圍。

此參照是「整個同源 AAL3 區域中的非病灶有效組織」，不是把原病灶形狀鏡射到對側後取一個同形 ROI，也不要求參照與原病灶具有相同體积或組織組成。稱為對側參照是準確的；稱為鏡像病灶控制則不符合目前程式。

### 5.3 四個共同 indicator

每個 AAL3 區域附加四個 0／1 指標：`affected` 表示原始病灶交集非空；`firstorder_available` 表示已算出局部八項統計；`contrast_available` 表示已算出八項對側差異；`texture_available` 表示已取得該局部 ROI 的完整 93 項 Original radiomics。

這四項指標描述病灶存在與計算支援，不是影像紋理本身。它們可能攜帶病灶位置、大小、是否達到最小 ROI 閾值及對側是否受累等資訊，因此本身就可能與結局相關。V6 將同一組 166×4＝664 個 indicator 加入 L、T 與 W，每個分支都不只是其名稱對應的連續 radiomics 特徵。

特別是 W 雖以完整腦區 radiomics 為起點，也附帶來自局部病灶支援的四個 indicator；不能把 W 描述成完全不含局部病灶可用性資訊的單純加權向量。

> 圖 5．病灶交集、同源對側參照與病灶膨脹排除：**pending**。
>
> 表 4．四個 indicator 的判定條件及其與缺失值的關係：**pending**。

## 6. L：局部一階統計與對側差異

### 6.1 八項局部統計

當局部 ROI 至少有 8 個有效 voxel 時，V6 對其正規化強度計算平均數、中央値、第 10 百分位數、第 90 百分位數、標準差、四分位距、均方根，以及離散強度熵，共八項數值。標準差使用 NumPy 的母體式定義，即分母為 voxel 數 n；四分位距為第 75 百分位數減去第 25 百分位數；均方根為強度平方平均的平方根。

L 的熵由固定寬度 5 的離散 bin 計算，先使用 floor(I*/5) 取得整數 bin，再由 ROI 內各 bin 的相對頻率 p_b 計算：

$$
H=-\sum_b p_b\log_2(p_b).
$$

這八項是 V6 自行計算的一階摘要，不能與 PyRadiomics 的 18 項 first-order features 視為同一組欄位。其百分位數與熵的實作應以 V6 函式及固定軟體版本為準。

### 6.2 對側差異

局部八項統計可用、左右配對存在，且對側參照至少含 64 個有效 voxel 時，對參照計算同樣八項統計。令局部某項統計為 a、對側參照為 b，其對比定義為：

$$
D(a,b)=\frac{a-b}{|a|+|b|+5}.
$$

分母中的常數 5 是固定設定，用以避免兩側數值很小時比例不穩定；不由 outcome 或外層效能估計。此公式保留相對大小的正負方向，同時使差異大小受兩側絕對值約束。它不是病灶值除以對側值，也不是個體內 z-score。

若原病灶不足 8 voxel、該區沒有左右配對，或對側有效支援少於 64 voxel，相關差異值不計算；由缺失值與 indicator 共同表達。對側參照是在原病灶統計可用的條件下才進一步判定，不能假設所有區域都各自建立一個對側比較。

### 6.3 L 的向量結構

每區有 8 個局部統計加 8 個對側差異，共 16 個連續特徵。166 區形成 2,656 欄，再附加 664 個 indicator，得到每位病人 3,320 欄的 L 向量。這是進入 training-fold 前處理前的原始欄位數，並非每個最終模型都使用 3,320 個 predictor。

> 圖 6．L 特徵計算與對側差異的示意：**pending**。

## 7. T：病灶與 AAL3 交集的 Original radiomics

### 7.1 可計算條件

T 使用與 L 相同的局部 ROI Ω_ik，但採較嚴格的支援要求：有效 voxel 數必須大於 64，且至少兩個影像座標軸上具有非零空間跨度。恰好 64 個 voxel 不符合此門檻；只有單一方向延伸的線狀 ROI，即使超過 64 voxel，也不通過預檢。

這個條件不是要求 ROI 一定有三維體積：例如單層但在另外兩軸都有範圍的 ROI 可能通過。PyRadiomics 設定為 `minimumROIDimensions=2`、`force2D=False`，因此不能把它描述成專門採用 native 2D in-plane texture 的方法，也不能將它與後續 V7 的 native 路線混稱。

### 7.2 抽取的特徵與參數

T 只啟用 Original 影像，不使用 wavelet、LoG 或其他影像濾波器。每個可用 ROI 抽取 93 項非 shape 特徵：18 項 first-order、24 項 GLCM、14 項 GLDM、16 項 GLRLM、16 項 GLSZM，以及 5 項 NGTDM。後五類合計 75 項纹理描述子。因此「T 有 93 個紋理特徵」並不精確，應寫為 93 項 Original radiomics，其中包含一階統計及紋理。

固定灰階離散化寬度為 5，距離設定為 `[1]`。PyRadiomics 內部 `normalize=False`，因為輸入已經做過前述病人內強度轉換；這個設定不表示原始 DWI 強度完全未經處理。`resampledPixelSpacing=None` 表示 extractor 不再額外做一次重採樣，但輸入本身已經是 MNI 2 mm 影像。其餘明列設定包括 `correctMask=True`、`minimumROISize=64`、`additionalInfo=False`。

每次抽取在保留完整 ROI 的前提下裁切運算區域，邊界盡可能保留一個 voxel 的背景，並拷貝對應的影像幾何資訊。這個裁切用於降低計算量，不是依預測標籤挑選影像 patch。對 Original-only 特徵而言，此處不加入卷積濾波影像所需要的額外背景。

計算前檢查 ROI 強度範圍除以 binWidth 是否超過 4,096，作為固定的資源與數值 QC 限制。抽取結果必須恰好有 93 項有效特徵；非預期錯誤會中止並留下紀錄，不會被靜默地當作有效結果。

### 7.3 T 的向量結構與解釋界線

166 區各 93 項連續特徵形成 15,438 欄，再加 664 個 indicator，T 共 16,102 欄。名字保留 AAL3 label、病灶角色與 radiomics descriptor，因此可追溯至區域與特徵種類。

當某區沒有病灶時，T 的連續值設為結構性 0；有病灶但 ROI 太小或不符合條件時，93 項連續值設為缺失。這使模型有機會利用局部計算是否可行的資訊，也意味著 T 的好表現不必然來自紋理。判讀時需要看實際入選內容，並與 flags-only 分支比較。

> 表 5．T 的特徵家族、數量與 PyRadiomics 參數：**pending**。
>
> 圖 7．局部 ROI 大小、空間支援與 T 可用性的示意：**pending**。

## 8. W：病灶占比加權的完整腦區 radiomics

W 沿用 V5 在完整 AAL3 區域內抽取的 93 項 Original radiomics。令完整第 k 區的第 j 個特徵為 F_ikj，該區病灶占比為：

$$
b_{ik}=\frac{|A_k\cap M_i|}{|A_k|}.
$$

這個占比的分子使用原始病灶遮罩與 atlas 的交集，分母為完整 atlas 區域的 voxel 數；它不是以有效病灶 ROI 的 voxel 數除以有效腦區 voxel 數，也不是病灶總量中落在該區的比例。

W 定義為 W_ikj＝b_ik F_ikj。若 b_ik＝0，程式明確將 W_ikj 設為 0，即使該完整腦區的原始 radiomics 不可用也一樣；若 b_ik＞0 而 F_ikj 缺失，乘积保留缺失狀態，再由 training-fold 補值處理。

完整腦區的原始抽取要求區域 voxel 數大於 64，且有效 DWI 覆蓋至少 90%。與 V6 的 L／T 不同，V5 的完整腦區 extractor 使用完整 atlas mask；通過 coverage 門檻後，ROI 內無效影像位置在正規化陣列中已設為 0，沒有再把整個 mask 改成有效 voxel 交集。因此，完整腦區與局部病灶兩条路線的有效訊號處理不完全相同，應照實說明。

W 的連續部分亦為 15,438 欄，再附加共同 664 個 indicator，總數為 16,102。它保留「先在完整腦區抽取 radiomics」的來源，但病灶占比使無病灶腦區的連續值歸零，已改變完整 F 的訊息內容。對非線性的紋理指標而言，乘上占比不等同於重新在病灶內計算該特徵；它是一個明確定義的加權表徵，不能被解釋為原特徵的數學等價形式。

> 圖 8．完整腦區 F 與病灶占比加權 W 的差異：**pending**。

## 9. 結構性零、缺失值與特徵矩陣組裝

V6 區分「沒有病灶」與「有病灶但無法計算」。前者的 L／T 連續特徵維持 0，四個 indicator 亦為 0；後者依可計算性將相應特徵設為 NaN，並以 indicator 呈現局部統計、對側差異及 T 的可用狀態。例如，一個 20 voxel 的有效病灶可有局部一階統計，但通常不符合 T 的大於 64 voxel 條件；其 L 與 T 不能一律視為完全缺失或完全可用。

L 的無病灶結構性 0 也不是在生物意義上宣稱病灶強度或左右差異恰好等於 0，而是模型編碼。模型可透過 indicator 區分部分情境，但不能因此保證模型只學到生物特性而非 ROI 大小或可計算性。

組裝時依同一組病人 ID、同一組區域順序與同一組 descriptor 名稱建立矩陣，並驗證 V5／V6 病人順序一致。該階段只進行固定規則的數值排列與 W 加權，不做 outcome-driven 特徵排名、不補值、不以全 cohort 估計 z-score，也不做 SMOTE。

資料量為 L：198×3,320；T：198×16,102；W：198×16,102。這些矩陣與來源 checkpoint、設定檔的 SHA-256 一起保存。重新讀取時，程式核對內容雜湊，避免在不知情下將不同抽取版本混入同一結果。

> 表 6．L／T／W 原始維度、結構性零與缺失處理：**pending**。

## 10. 強基準模型的組成

V6 的主要基準名稱為 `V1_base`，由三部分串接。第一部分為前述 12 個 core clinical predictor、病灶體積的自然對數轉換 log(1＋volume_mm³)，以及全病灶 radiomics availability。病灶體積由完整病灶 mask voxel 數乘上影像 spacing 的乘積得到；此體積與某一 AAL3 區域的 burden 不同。

第二部分為 166 區的 AAL3 病灶占比，在當次 training fold 內經前處理與 RFE 保留最多 20 項。第三部分為 V5 `unclip_p99` 的全病灶 Original radiomics，從 93 項候選中在 training fold 內選擇；mRS 最多保留 20 項，NIHSS 最多保留 35 項。

全病灶 radiomics 使用完整病灶 mask 作 ROI，與 T 的逐區病灶交集不同。上游 V5 同樣檢查病灶大小與有效覆蓋；無法計算時保留缺失。`global_available` 的實作是 V5 三種正規化候選的 global features 皆存在時為 1，不是依當次 outcome 判斷是否可用。這個指標與其餘小型基準欄位一樣，若在 training fold 為常數，會被移除。

core clinical、log-volume 與 global availability 只接受訓練內補值、常數剔除與縮放，不再做 supervised RFE；burden 與全病灶 radiomics 各自獨立做 RFE。串接後不再把全部基準特徵送進一次共同 RFE。若暫不考慮常數剔除及可用欄位不足，mRS 基準最多為 14＋20＋20＝54 欄，NIHSS 最多為 14＋20＋35＝69 欄；實際數量以各 fold 紀錄為準。

`V1_base` 與 V5 的 `B_R__unclip_p99` 使用對齊的規格。V6 在每個 fold 核對兩者的訓練／驗證 ID，並要求基準預測在絕對誤差 1e−12、相對誤差 0 下相同，以避免區域特徵比較時基準偷偷發生變動。

> 圖 9．強基準模型的三部分輸入與分支內選擇：**pending**。

## 11. 每個 training fold 內的前處理

### 11.1 有效欄位、補值與常數剔除

外層分割完成後，對每一組待處理特徵，先只用 training 病人判定哪些欄位至少有一個有限值。若 training fold 整欄皆缺失或非有限值，該欄在訓練與驗證兩邊共同刪除；即使驗證病人中剛好有有效值，也不能用其值救回此欄。

對保留欄位，以 training fold 的有限值計算中位數，再用同一個中位數補訓練與驗證矩陣的非有限值。結構性 0 是有限值，會參與中位數的計算，不被當成缺失而跳過。補值後再以 training 資料檢查 range 是否大於 0，移除訓練內常數欄位。驗證病人的變異不參與欄位保留判定。

### 11.2 標準化

對剩餘欄位，使用 training fold 補值後數據的平均 μ_j 與母體標準差 σ_j 做 z-score：

$$
z_{ij}=\frac{x_{ij}-\mu_j}{s_j},\qquad
s_j=\begin{cases}1,&\sigma_j<10^{-8}\\\sigma_j,&\text{其他情況。}\end{cases}
$$

相同的 μ_j 與 s_j 套用於 held-out 病人。統計量以 float64 計算，轉換後以 float32 保存。所有 indicator 亦隨該分支進行標準化，並非永遠以未縮放的 0／1 進入最終模型。

樹模型本身未必需要 z-score，但 V6 的 SMOTE 使用距離尋找鄰居，LR 也有 L2 正則化，因此此尺度處理是實際演算法的一部分。不能因主模型是 LightGBM，就在重現時任意省略而仍稱為完全相同的方法。

### 11.3 訓練與驗證的資訊邊界

補值中位數、欄位保留條件、標準化統計量與後續 supervised 選擇都只由 training fold 決定。驗證資料接受已決定的轉換，不回饋參數估計。這與抽取前對單一病人影像做固定正規化不同：前者是在多位病人間學習轉換，後者是每位病人獨立的 deterministic 計算。

> 圖 10．單一 fold 的資料邊界：訓練內估計與 held-out 套用：**pending**。

## 12. 分階段遞迴特徵消除

V6 的 `staged_rfe` 在 training fold 的原始病人資料上訓練 LightGBM，依特徵重要性排序，逐步移除較低排名的特徵。這一步發生在 SMOTE 之前，因此 RFE 的模型不使用合成病人。

刪除程序有四階段：先以目標 10,000、每次最多刪除 40,000 的規則縮減；再以目標 950、每次最多刪除 2,000 縮減；之後以目標 50、每次最多刪除 150 縮減；最後每次刪除一個，直到任務指定的保留數。若原本可用特徵已少於某階段目標，該階段直接跳過。每次保留數為當前數量減去該階段步長，但不得低於目標；不是把所有分支強行補到 10,000 欄。

burden 的目標固定為 20；全病灶 radiomics 與 L／T／W 的目標分別為 mRS 20、NIHSS 35。若 training-fold 前處理後的欄位數本來就少於目標，保留全部可用欄位，不生成虛構特徵來補足數量。

實際使用的 LightGBM importance 為 `importance_type='split'`，即分裂次數型重要性，不是 SHAP 排名，也不是明確指定的 gain importance。排名使用穩定排序，重要性相同時保留當前順序；同一 fold 的各次 RFE fit 使用相同 seed。這讓重現具有可追溯性，但相關特徵、相同或零重要性的欄位仍可能彼此替代，入選不代表已證實穩健的生物標記。

V6 只借用 V4 模組中的這些函式，沒有執行 V4 主程式的內層候選選模。原始 V6 使用固定候選與固定設定的 repeated outer CV，不是完整 nested CV，也沒有內層搜尋最好的 L／T／W 或 blend 權重。方法段不能因 import 的模組標題含有 nested 就寫成已完成 nested tuning。

> 表 7．RFE 階段、刪除步長與各分支保留數：**pending**。

## 13. LightGBM、SMOTE 與融合模型

### 13.1 主分類器

主模型為二元 LightGBM GBDT。已核對的設定為 learning rate 0.04、最大深度 3、最多 6 個 leaves、L1 正則化係數 0.1、L2 正則化係數 1.0、每棵樹使用 50% 特徵，以及 100 個 boosting estimators。100 來自本次固定 LightGBM 版本中的有效預設，並非重新搜尋得到。

其他有效參數包括 `min_child_samples=20`、`min_child_weight=0.001`、`min_split_gain=0`、`subsample=1.0` 與 `subsample_freq=0`。沒有另設 class weight，沒有提供 outer validation set 進行 early stopping，也沒有依 held-out 表現調整迭代數。執行時使用固定 random_state、`deterministic=True` 與 `force_col_wise=True`。

### 13.2 Training-only SMOTE

完成各分支的前處理與選擇，並形成該模型最後的訓練向量後，才在 training fold 內執行 SMOTE。實作是保存在專案中的 deterministic 二元線性插值函式，不是未記錄參數的外部重採樣呼叫。

首先辨識 training fold 少數類別，將其補到與多數類別相同樣本數。鄰居數最多為 3；若少數類原始病人較少，使用 min(3, n_minority−1)。需要補樣本但少數類不足兩人時直接報錯。每個合成點由少數類病人 x 與其一個少數類近鄰 x′，以均勻隨機係數 λ 線性插值得到：

$$
x_{\mathrm{syn}}=x+\lambda(x'-x),\qquad 0\leq\lambda<1.
$$

對融合模型，SMOTE 的鄰居是在已串接的完整模型向量中計算，不是先對每個分支分別產生合成病人後再拼接。held-out 病人不參與鄰居搜尋、合成點產生或類別平衡，也不被複製進訓練集合。

本實作沒有使用專門處理類別變數的 SMOTENC。臨床二元欄位及 availability indicator 經標準化後一併插值，合成資料可能出現介於兩個編碼之間的值。這是目前數值流程的實際性質，應列為解釋與模型設計限制，不能將所有合成向量都描述成一位具有完整臨床意義的真實病人。

### 13.3 單獨區域模型

`V1_L`、`V1_T` 與 `V1_W` 各自只使用對應分支在 training fold 選出的特徵，經該模型自己的 training-only SMOTE 後訓練 LightGBM。其目的在於觀察區域表徵單獨的辨識能力；這些模型仍包含分支中可被選入的共同 indicator。

### 13.4 特徵層級融合

`V1_base_L`、`V1_base_T` 與 `V1_base_W` 將強基準向量與對應區域分支選出的特徵直接串接，再對串接向量做 SMOTE 與 LightGBM 訓練。區域與基準分支先各自完成 RFE，串接後不再做共同的 supervised 篩選。

目前採作 V6 核心展示的 `V1_base_T` 可概括為：core clinical、log-volume、global availability、training-selected burden、training-selected global lesion radiomics，加上 training-selected T。若暫不計常數與可用性剔除，其輸入上限為 mRS 74 欄、NIHSS 104 欄。這不是只有區域 T，也不是原稿的 F-only。

### 13.5 機率層級融合

另固定建立 `V1_blend_L`、`V1_blend_T` 與 `V1_blend_W`，以同一 held-out 病人的兩個模型預測作加權：

$$
p_{\mathrm{blend}}=0.75p_{\mathrm{base}}+0.25p_{\mathrm{regional}}.
$$

兩個 component model 都只由該 fold 的 training 病人訓練。0.75／0.25 是設定中固定的權重，沒有根據 outer test AUC 調整，也沒有另訓練 stacking meta-classifier。機率平均與特徵串接是不同模型，不能互相替代其效能或 SHAP 解釋。

> 圖 11．單分支、特徵串接與固定機率融合：**pending**。
>
> 表 8．LightGBM 與 SMOTE 的固定參數：**pending**。

## 14. 次要分類器與比較矩陣

原始 V6 的次要模型為 ExtraTrees（ET）及 L2 logistic regression（LR）。ET 使用 400 棵樹、`min_samples_leaf=5`、`max_features=0.7`、`class_weight='balanced'` 與固定 random_state。LR 使用 `C=0.1`、L2 penalty、`solver='liblinear'`、`class_weight='balanced'`、`max_iter=5000` 與固定 random_state。

ET／LR 對輸入做相同的 training-only 缺失處理與標準化，但使用前處理後的全部可用欄位，不套用 LightGBM 的 RFE 結果，也不執行 SMOTE。其類別不平衡處理由 balanced class weight 完成。

ET／LR 的基準包含 core clinical、log-volume、global availability，以及前處理後的全部 burden 與全病灶 radiomics。因此，ET_base 與 LR_base 的輸入維度不同於經分支 RFE 的 V1_base。不同分類器之間的絕對分數不能被解釋成只改了 classifier 的嚴格單因素比較；若要判斷區域特徵對 ET 或 LR 的增量，應與該分類器自己的 base 比較。

每個結局的原始 V6 共 24 種方法：V1 基準一種、L／T／W 單分支三種、base 加各分支三種及固定 blend 三種，共十種；ET 與 LR 各有一種基準、三種單分支、三種基準融合，各七種。所有固定比較都報告，不只保留最高 AUC 的那一列。

> 表 9．原始 V6 的 24 種比較方法及其輸入／不平衡處理：**pending**。

## 15. 交叉驗證、OOF 預測與核心版本的選擇

### 15.1 重複分層五折交叉驗證

每個任務分別使用 repeated stratified five-fold cross-validation，重複三次，共 15 個分割。分割器 random_state 固定為 20260909。每次 repeat 中，一位病人只屬於一個 held-out fold；該次其餘四折作 training。不同 repeat 會重新分割，因此同一病人總共得到三次 OOF 預測。

同一任務的所有候選共用相同 train/test ID。分割以臨床表中符合条件的病人順序建立，影像矩陣再透過 ID 對齊，不假設影像檔與臨床表原本排序相同。程式檢查 ID 唯一性，並保存每個 split 的 train/test 名單。

令 split 編號 s 從 0 到 14，訓練 seed 為 20260909＋s＋1000×任務索引；原始 V6 的 mRS 索引為 0、NIHSS 為 1。此 seed 用於該 split 的模型及重採樣。多個方法用相同分割和對齊的 seed，旨在減少比較因隨機切分不同而產生的差異，不代表其模型參數或特徵矩陣完全相同。

### 15.2 OOF 的彙整

病人 i 在三次 repeat 的 OOF 機率為 p_i1、p_i2、p_i3，主要彙整使用：

$$
\bar p_i=\frac{p_{i1}+p_{i2}+p_{i3}}{3}.
$$

主要 AUROC 由各病人的單一平均機率 p̄_i 與真實標籤計算，不是把 15 個 fold 的 AUC 直接平均，也不是把同一病人的三次預測當成三位獨立病人。另保留每個 repeat 完整 OOF 的指標，再計算三個 repeat AUC 的平均與範圍。

「平均 OOF 機率後的 AUC」與「mean-repeat AUC」不一定相同，因為平均機率可能改變病人間排名。正文與表格必須寫出採用哪一個統計量，不能在不同模型間混用較好看的數字。

### 15.3 既有主要比較與後來採用的核心

V6 設定檔明列的主要比較是 mRS 的 `V1_base_L` 相對 `V1_base`。這是保存下來的該輪 protocol 內容，不等於整個研究從未看過這個 cohort，也不等於已有外部預註冊證明。

`V1_base_T` 是在看過既有比較後，因其分數而選作目前核心展示的分支。因此可稱為事後選出的 V6 核心或候選模型，不能改寫成原先預先指定的 primary endpoint／primary model，也不能用重新執行相同資料來消除這個選擇過程。

本次精確重跑重新訓練模型，但使用經雜湊核對的既有影像特徵檔，沒有再次由 MRI 抽取所有特徵。所有原始 V6 方法的 OOF 能重現，支持計算流程的可重現性；它不等同於外部驗證、時間外驗證或全新未接觸測試集。

> 圖 12．三次五折 CV 與每人三次 OOF 的彙整方式：**pending**。
>
> 表 10．主要比較、事後核心與次要比較的研究角色：**pending**。

## 16. 評估指標與不確定性

主要辨識指標為 AUROC。另報告 AUPRC、accuracy、sensitivity、specificity 與 Brier score。AUPRC 使用 scikit-learn `average_precision_score` 的 average precision 定義，應避免把它寫成已採其他插值方式計算的 PR 曲線面積。

分類標籤以機率大於或等於 0.5 判為正類，未使用 outer test label 選擇最佳閾值。accuracy 為正確分類比例；sensitivity 為正類中被辨識出的比例；specificity 為負類中被辨識出的比例。Brier score 為預測機率與二元標籤差的平方平均。Brier 兼具多種機率預測品質資訊，單獨報告它不代表已完成 calibration curve、校準斜率或臨床效益分析。

信賴區間採 2,000 次病人層級的分層 bootstrap。先將病人按真實正／負類分成兩組，再在每組內有放回抽取與原組相同數量的病人。每次抽樣同時套用到所有候選模型的固定 OOF 機率，使模型間差異保持病人配對。對 bootstrap AUC 分布取第 2.5 與第 97.5 百分位數，作為 conditional 95% CI。

兩模型增量定義為候選 AUC 減去 reference AUC；其區間由同一次抽樣中的 AUC 差值分布計算，不是把兩個模型各自 CI 端點相減。主要完整彙整以 V1_base 為 reference；另保存 ET／LR 各自對應 base 的 matched comparisons，判讀這兩類模型的區域增量時應使用後者。

bootstrap 的抽樣單位是病人平均 OOF 預測，不是 fold，也不是 SMOTE 合成點。bootstrap 中沒有重新做資料分割、補值、RFE、SMOTE 或模型訓練，因此得到的是固定預測條件下的病人抽樣不確定性。區間不涵蓋整條訓練流程的不穩定、跨版本研究選擇、最佳模型挑選或重複探索造成的樂觀偏差，也沒有進行多重比較校正。區間跨 0 的配對增量不能表述為已證实穩定改善。

> 表 11．評估指標、閾值與 bootstrap 規格：**pending**。
>
> 圖 13．OOF 評估與配對病人 bootstrap 流程：**pending**。

## 17. Availability／flags 的事後消融

V6 在分析入選特徵後，另增加一個針對 mRS 的事後機制檢查。該實驗取 T 中的 664 個 indicator，排除連續 radiomics，使用相同病人、相同 15 個分割、相同強基準，以及相同 training-only RFE／SMOTE／LightGBM 流程。indicator 分支在每個 training fold 重新選取最多 20 項，再與基準串接，形成 base＋flags 對照。

此對照用於檢驗：如果區域 T 的效能主要來自病灶是否存在或是否達到計算門檻，僅保留 flags 是否已能取得相近表現。比較仍以固定 OOF 的配對 bootstrap 進行。

這不是把某個已訓練 T 模型的連續欄位逐一刪掉、其餘權重保持不變的實驗。flags-only 重新執行 RFE，可能選到不同的 indicator，模型也重新訓練；因此差異不能被解釋為單一纹理描述子的精確因果貢獻。這個實驗是在看到 T 的入選構成後提出，也不是新的獨立驗證。

後續四任務補充比較另有「flags-only 單分支」模型，該模型不含上述 base。這兩種 flags 對照的輸入不同，不能只因名字都含 flags 就把結果互換。

> 圖 14．原始 T、base＋T、base＋flags 與單獨 flags 的消融關係：**pending**。
>
> 表 12．分支入選特徵的 radiomics／indicator 組成與消融設計：**pending**。

## 18. 後續四任務論文比較：與原始 V6 分開描述

### 18.1 任務與特徵集合

為補齊原稿的比較結構，後續實驗加入 T/E 與 LVO。T/E 使用原欄位 `T/E`，值為 E 時編碼為 1；LVO 使用原欄位 `LVO-index` 的二元值。T/E 有效樣本 196 人，正類 38、負類 158；LVO 有效樣本 196 人，正類 77、負類 119。兩者各有 2 位影像 cohort 病人缺失對應標籤。其臨床定義、判讀規範及時間點仍須查原始資料字典；不能一概將它們寫成與出院 mRS／NIHSS 相同意義的功能預後。

補充實驗使用五組特徵：C 為保守 core clinical；R 為全病灶 Original radiomics；S 為 AAL3 burden；T 為 V6 的局部病灶 radiomics 加 indicator；Fwhole 為 V5 完整 AAL3 的 Original radiomics、採 unclip_p99。T/E 的 C 另排除 AF 欄位，應在方法中明列，不與其他任務視為完全相同的臨床輸入。

Fwhole 不包含 V6 的共同局部 indicator，也不乘病灶占比；它與 W 不同。Fwhole 雖保留原稿完整區域的設計核心，但只有 Original 特徵，不能冒稱原稿全部 filters 的 F。S 只包含 AAL3 病灶占比，也不能冒稱原稿包含其他 atlas 的完整語意特徵集合。

### 18.2 分類器與選擇法

補充比較使用 Random Forest、XGBoost、CatBoost 與 LightGBM，而不是原始 V6 的 ET／LR 配置。Random Forest 設 100 棵樹及最大深度 4；XGBoost 設 100 個 estimators、learning rate 0.015、`eval_metric='logloss'`；CatBoost 設 Logloss、AUC evaluation metric、learning rate 0.01、30 iterations、depth 9、`od_type='Iter'`。沒有提供 outer test set 作 CatBoost 的 eval_set，因此不能將該參數寫成已由外層測試執行 early stopping。其餘未明列的參數依已記錄軟體版本的預設值。

四個 classifier 各自在 training fold 用自身的 feature importance 進行 staged RFE，再使用 training-only SMOTE 訓練最終模型；不是把 LightGBM 選好的特徵一律交給另外三個分類器。不同 classifier 的 importance 定義也不應當作相同量直接比較。

補充的 backward 選擇依原 notebook 真正操作重現：對候選特徵訓練一次分類器，再按 importance 取前 k 個。它不是逐步刪除每個候選並重新估計效能的 sequential backward selection。RFE 與這個 single-fit top-k 的差異必須在文字上說清楚。

### 18.3 特徵組合與模型數量

mRS、NIHSS、T/E、LVO 的分支保留數分別為最多 20、35、20、25。每個分支各自選擇，再串接；這包含 S，故補充 NIHSS 的 S 目標是 35，與原始 V6 基準固定 20 項 burden 不同。C 若可用欄位數不足目標，保留所有前處理後欄位。

每個任務有 36 種固定方法：四分類器乘五種單分支共 20 種；LightGBM 的五種 backward 單分支；T 與 Fwhole 各自再與 C、R、S、C＋R、C＋R＋S 形成五種組合，共十種；以及一種 LightGBM flags-only。每種方法各有 15 個 folds，四任務共 2,160 個 fold×方法組合。

補充的 T＋C＋R＋S 不含原始 V6 base 額外加入的 log-volume 與 global availability，S 的保留數在部分任務也不同，所以它不是 `V1_base_T` 的別名。不能把其 AUC 與其 SHAP 圖標示成原始 V6 核心的解釋。補充表格中的 ΔAUC 預設 reference 為 LightGBM＋RFE 的 C-only，亦不同於原始 V6 的強基準 reference。

> 表 13．四任務補充比較的模型、任務、特徵與參數：**pending**。
>
> 圖 15．原始 V6 核心與後續 T＋C＋R＋S 的輸入差異：**pending**。

## 19. SHAP 解釋與腦區映射

### 19.1 解釋模型與病人的範圍

本輪新增的 SHAP 由四任務補充實驗中的 LightGBM 模型產生。原始 `model_v6.py` 本身沒有保存相同形式的 held-out SHAP，因此目前圖檔不能直接當成 `V1_base_T`、ET、LR、RF、XGBoost 或 CatBoost 的解釋。此處說明的是已實作的補充 LightGBM 解釋流程。

每個 fold 只對該 fold 的 held-out 病人計算 TreeSHAP。使用 LightGBM `pred_contrib=True`，輸出各特徵貢獻及 expected value；其單位是 raw log odds，不是機率百分點。背景依 tree-path-dependent 方法，由訓練樹路徑資訊決定；訓練資料包含 training-only SMOTE 的合成樣本。因此它不是以未重採樣的原始病人另行建立的 interventional SHAP 背景。

程式檢查：每位 held-out 病人的特徵貢獻加上 expected value，可還原 raw score；對 raw score 套 sigmoid 後，也能還原保存的預測機率。這是數值一致性檢查，不能代替生物學有效性驗證。

### 19.2 跨 fold 的特徵重要性

不同 fold 可能選到不同特徵。彙整時建立該模型所有 fold 曾使用特徵的聯集；某病人在某個 held-out 模型沒有使用到的特徵，其貢獻設為 0。特徵重要性是所有 N×3 次 held-out 預測的平均絕對 SHAP：

$$
I_j=\frac{1}{3N}\sum_{i=1}^{N}\sum_{r=1}^{3}|\phi_{ijr}|.
$$

分母不只計入該特徵被選中的 folds，避免極少數入選的特徵被誤寫成在整個研究中普遍重要。另計算 15 個 folds 中該特徵被選入的比例。入選頻率與 SHAP 重要性是不同指標，不能把其中一項當成另一項。

summary 圖顯示前 20 個平均絕對 SHAP 較大的特徵，點代表 held-out 的特徵貢獻，顏色依 fold 內前處理後的特徵值作相對高低編碼。對未入選而沒有對應值的點，使用中性色；顏色不是原始 MRI 強度，也不能用來讀取臨床欄位的原始量綱。

### 19.3 AAL3 腦區圖

T 的特徵名稱含 `AAL3_區域編號`，Fwhole 的名稱保留區域編號尾碼，可直接對應回所用 AAL3 圖譜。對第 k 區，将屬於該區的各特徵平均絕對 SHAP 加總作為區域重要性，填入該 label 的所有 voxel，輸出 NIfTI 與切片圖。

radiomics 與 availability 指標分開聚合及繪圖。R 是全病灶描述子，沒有每項特徵對應到單一 AAL3 區域的定義，因此不強行把 R 的 SHAP 分配至某個腦區。現有製圖程式也沒有把 S 的 burden 圖與 T／Fwhole 的 radiomics 圖混成同一張輸出。

腦區圖以 canonical 方向顯示、標記左右與切片位置，色階表示平均絕對 SHAP 的加總。此聚合會受到各區入選特徵數、特徵相關性及模型背景的影響；它不是每個 voxel 的獨立 attribution，也不是因果損傷圖。使用絕對值後不保留正負方向，因此不能從色階直接讀出某區「改善」或「惡化」預後。

> 圖 16．held-out TreeSHAP、跨 fold 零補貢獻與入選頻率的彙整：**pending**。
>
> 圖 17．各特徵組的 SHAP summary：**pending**。
>
> 圖 18．radiomics 與 availability 分開呈現的 AAL3 重要性圖：**pending**。
>
> 表 14．重要特徵、腦區名稱、特徵種類與入選頻率：**pending**。

## 20. 資料洩漏控制與尚未消除的偏差

目前流程避免了幾項直接的驗證資料洩漏：未用全 cohort 的結局先選特徵；未用 held-out 病人決定中位數、常數欄位與標準化尺度；未在分割前產生跨訓練／驗證資料的 SMOTE 樣本；未讓同一病人的不同腦區當成獨立樣本分散到兩側；也未用 outer test label 選擇分類閾值、early-stopping 時點或固定 blend 權重。

但是，本研究已在相同 cohort 上進行多版本探索。研究者已看過候選的結果，再决定以 V6 及其中某條支線作為核心。這種方法與模型選擇的資訊回饋不會因單次 CV 的 RFE 改成 training-only 而自動消失，也不會因為重新訓練能得到一樣的分數就消失。現有結果應描述為反覆使用 cohort 的內部探索與可重現性評估。

特徵抽取不讀 outcome 並不足以證明全部上游資料都無洩漏。臨床欄位的取得時點、病灶標註是否使用後續資訊、配準與影像 QC 的人工決策、歷史標籤修正，以及原始排除個案的原因，仍需以來源紀錄檢查。這些不能從本次模型程式的 train/test ID 無交集直接推得。

目前也沒有全新外部測試集、獨立時間 cohort 或未碰觸保留集。病例數相對原始特徵數很小，RFE 的特徵集合可能不穩定，SMOTE 不能增加真正獨立的病人數；多重比較、固定候選取最高分與 outcome 缺失亦影響解釋。這些是界定推論範圍所需的限制，不應寫成已由本方法完全解決。

> 表 15．已執行的洩漏控制與仍需來源證據的風險：**pending**。

## 21. 品質驗證、重跑與軟體環境

抽取層面的既有測試涵蓋：按左右解剖名稱配對而非 label 奇偶；八項一階統計在已知輸入下的值；對側參照排除膨脹病灶；以及 64 voxel 邊界與退化 ROI 的處理。模型層面的測試則改變 held-out 數值，確認不會改變 training 前處理輸出與保留欄位，並確認含結構性零及缺失值的轉換後矩陣為有限值。這些是針對關鍵行為的檢查，不是對所有資料來源與所有可能洩漏的全面證明。

原始 V6 模型保存各 split 的病人 ID、seed、各分支入選特徵、RFE 剩餘數紀錄、前處理摘要及所有預測。獨立重跑核對每位病人三次 OOF、train/test 無交集、OOF 與分割預測一致、AUC 重新計算一致，並核對 V5 基準及新舊 V6 預測的精確再現。補充實驗另驗證全部四任務各 15 folds、36 方法，以及 SHAP 加總還原機率。

抽取 checkpoint 記錄 NumPy 1.26.4、SimpleITK 2.5.5 與 PyRadiomics v3.0.1。模型重跑記錄 NumPy 2.0.2、pandas 2.3.3、scikit-learn 1.6.1 與 LightGBM 4.6.0。補充分類器比較另記錄 XGBoost 2.1.4 與 CatBoost 1.2.8。影像抽取與建模位於不同環境，應分別記載，不能只給出其中一份版本清單。

抽取採病人層級平行處理，每個工作將 SimpleITK thread 數設為 1；模型依固定 jobs 執行，並記錄環境與來源程式的 hash。精確重跑的原始 V6 使用 jobs=8；補充實驗使用每個工作 jobs=2，並以多個 fold 工作平行處理。平行數不被當成模型超參數搜尋的一部分。

特徵封存檔、來源影像／遮罩、atlas、label、程式及 protocol 具有雜湊與 checkpoint 關聯。程式只在完整工作成功後寫入完成標記；報告產出與模型訓練各有各自的完成狀態，不能因訓練完成便假設圖表已成功生成。輸出保留逐人 OOF、所有方法指標、CI、特徵頻率、SHAP、腦區值、NIfTI、設定與稽核紀錄。

這套工作以 CV 評估與模型比較為目的，目前這些輸出本身不等於已封裝一個在全部病人上重訓、可直接部署到新病人的最終模型。若正式稿要描述臨床使用流程，尚需另外定義最終重訓、前處理物件保存、輸入 QC 與外部驗證程序。

> 表 16．軟體版本、測試、來源追溯與重現性檢查：**pending**。

## 22. 可供方法段使用的整合敘述

本研究使用既有 MNI 2 mm 空間 DWI、病灶遮罩與 AAL3 圖譜，建立病灶相關的區域放射體學表徵。影像在病人層級進行固定 brain-support 百分位強度正規化，以下界第 1 百分位數與第 99 百分位距離作為尺度，保留超過第 99 百分位數的高強度訊號。影像、病灶、腦區圖譜及 brain support 在特徵抽取前檢查共同格點與固定 coverage 條件。

區域表徵包含三個分支。L 在每個病灶與 AAL3 區域交集內計算八項強度摘要，並在支援足夠時與排除病灶及其膨脹範圍的對側完整同源腦區比較。T 在具有大於 64 個有效 voxel 及至少兩個非退化空間維度的病灶交集內，提取 93 項 Original radiomics。W 將完整 AAL3 腦區的 Original radiomics 乘上該區病灶占比。每個分支均附加區域病灶存在、一階統計可用、對側差異可用及局部 radiomics 可用四個指標。無病灶以結構性零表示，有病灶但無法計算則保留缺失並記錄可用性。

模型以臨床 core 變數、病灶體積、全病灶特徵可用性、AAL3 病灶占比及全病灶 radiomics 建立強基準，再分別加入 L、T 或 W。每個任務採三次重複分層五折交叉驗證。欄位有效性判定、補值、常數剔除、標準化與分階段 RFE 均限制於 training fold；選出的特徵完成串接後，僅在 training fold 內執行 SMOTE，再訓練固定參數的 LightGBM。另使用各自匹配基準的 ExtraTrees 與 L2 logistic regression 作次要比較，並保留固定權重的機率融合對照。

每位病人具有三次未參與該模型訓練的 OOF 機率，主要結果以其平均值計算 AUROC，另報告 mean-repeat AUROC 及固定 0.5 閾值下的分類指標。模型差異使用固定 OOF 的分層、配對病人 bootstrap 估計 conditional 95% CI。所有候選皆完整報告；後來採作核心的 V1_base_T 與原設定的主要比較 V1_base_L 明確區分。由於同一 cohort 已經過多版本探索，本研究結果定位為內部探索與可重現性評估，不作獨立外部驗證的宣稱。

後續為重現原稿比較結構而加入的四任務分類器、特徵組合與 held-out TreeSHAP 分析，需另列為補充方法。該分析中的 T＋C＋R＋S 與原始 V6 基準融合模型輸入不同，其結果與解釋圖不交互替代。重要腦區由有明確區域名稱的特徵映射回 AAL3，並將 radiomics 與 availability 重要性分開呈現，以免把計算支援或病灶位置資訊誤歸因於紋理。

## 23. 文件依據與後續編修位置

本文件直接依據以下專案檔案撰寫；它們是實作與執行紀錄來源，不取代正式論文中對演算法及軟體原始文獻的引用。

- [V6 固定 protocol](../../config/protocol.json)：分支、閾值、主要比較與 CV 設定。
- [V6 特徵抽取](../../extract_v6.py)：局部統計、左右配對、ROI 與對側參照。
- [V6 特徵組裝](../../assemble_v6.py)：L／T／W、共同 indicator、零與 NaN 的编码。
- [V6 建模](../../model_v6.py)：強基準、原始 24 方法、分割與匹配比較。
- [V5 抽取與正規化](<../../../ischemic stroke analyseV5/extract_v5.py>)：強度處理、完整腦區與全病灶 radiomics。
- [V5 特徵組裝](<../../../ischemic stroke analyseV5/assemble_v5.py>)：burden、體積與 global availability。
- [V4 共用函式](<../../../ischemic stroke analyseV4/run_v4.py>)：core clinical、staged RFE、縮放與 bootstrap。
- [基礎資料清理、分類器及 SMOTE](<../../../ischemic stroke analyseV2/leakage_free_experiment.py>)：標籤、訓練內前處理與分類指標。
- [V6 flags 事後消融](../../ablate_flags.py)：mRS base＋flags 的對照。
- [四任務補充實驗](../reproduce_tables.py)：36 方法、各 classifier 的 RFE 與 held-out SHAP。
- [四任務報告與腦區映射](../report_tables.py)：OOF 稽核、SHAP 彙整、腦區聚合與輸出。
- [原始 V6 重跑環境](core/config.json)及[補充比較環境](experiment_config.json)：實際版本與來源雜湊。
- [臨床 predictor 排除要求](../../../CLINICAL_PREDICTOR_EXCLUSIONS.md)：評分及衍生量的排除範圍。

圖 1 至圖 18、表 1 至表 16 均保留 **pending**。正式排版時可合併部分圖表，但不宜把方法定義不同的模型合併成同一名稱。臨床來源、倫理、掃描、分割與配準資訊需另依原始研究紀錄補齊；在補齊前，不應將本文件視為已涵蓋所有投稿所需的研究方法資訊。
