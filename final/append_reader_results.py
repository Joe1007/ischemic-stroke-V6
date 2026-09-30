"""Append numerical results to the standalone methods document, without images."""
from pathlib import Path
import json,hashlib,re
import pandas as pd
HERE=Path(__file__).resolve().parent;OUT=HERE/'output';DOC=OUT/'METHODS_V6_DETAILED.md'
tasks=['mRS','NIHSS','TE','LVO'];tasklabels={'mRS':'mRS','NIHSS':'NIHSS','TE':'T/E','LVO':'LVO'}
d=pd.read_csv(OUT/'all_metrics.csv');assert len(d)==144
core={t:pd.read_csv(OUT/'core'/f'{t}_summary.csv') for t in tasks[:2]}
assert all(len(x)==24 for x in core.values())
pairs={t:pd.read_csv(OUT/'core'/f'{t}_matched_comparisons.csv') for t in tasks[:2]}
used_sources=['all_metrics.csv','cohort_counts.csv']+[f'core/{t}_{s}.csv' for t in tasks[:2] for s in ['summary','matched_comparisons']]

def label(m):
    if m.startswith(('V1_','ET_','LR_')):
        model,block=m.split('_',1);model={'V1':'LightGBM','ET':'ET','LR':'LR'}[model]
        if block=='base':b='B' if model=='LightGBM' else model+' 基準'
        elif block.startswith('base_'):b=('B' if model=='LightGBM' else model+' 基準')+'＋'+block[5:]
        elif block.startswith('blend_'):b='機率融合 P_'+block[6:]
        else:b=block+' 單分支'
        return model+'：'+b
    model,selector,block=m.split('_',2)
    return {'RF':'隨機森林','LightGBM':'LightGBM','XGBoost':'XGBoost','CatBoost':'CatBoost'}[model]+'／'+{'RFE':'RFE','backward':'單次篩選'}[selector]+'／'+block.replace('_','＋').replace('flags','僅指標')

def num(x):return f'{float(x):.4f}'
def ci(r):return f'{r.auroc:.4f} [{r.auc_low:.4f}, {r.auc_high:.4f}]'
def delta(r):return f'{r.delta_auc:+.4f} [{r.delta_low:+.4f}, {r.delta_high:+.4f}]'
def table(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','|'+'|'.join(['---']*len(headers))+'|']+['| '+' | '.join(map(str,r))+' |' for r in rows])
lines=['## 24. 結果表：閱讀方式與分析樣本','',
'本節附上已完成的數值結果。表 R1 起為實際結果表，不是待填占位；方法章原有圖表占位仍保留 pending。所有模型使用相同任務內的病人分割，但不同任務可用病人數不同。此處不嵌入 SHAP、ROC 或腦區圖片。','',
'AUROC 以每位病人三次折外機率平均後計算；mean-repeat AUROC 是三個 repeat 分別計分後的平均，兩者分開列示。所有分類指標的閾值固定為 0.5。95% CI 是固定 OOF 的配對病人 bootstrap 條件區間，不包含重新訓練或研究選擇的不確定性。AUPRC 採 average precision 定義。','',
'B 是本研究的強基準，包含 C、log-volume、全病灶可用性、入選 S 與入選 R。B＋T 是主分析核心；補充的 T＋C＋R＋S 不含前兩項額外影像摘要，部分分支保留數也不同，不能混為同一個模型。ET／LR 使用各自的基準。','',
'**表 R1．各任務有效樣本與標籤分布。**','']
c=pd.read_csv(OUT/'cohort_counts.csv').set_index('task')
lines.append(table(['任務','影像分析集合','有效 n','正類','負類','標籤缺失'],[[tasklabels[t],198,int(c.loc[t,'n']),int(c.loc[t,'positive']),int(c.loc[t,'negative']),int(c.loc[t,'missing_target_in_image_cohort'])] for t in tasks]))
lines+=['','## 25. 主分析：強基準與區域特徵','',
'**表 R2．LightGBM 強基準及其區域特徵融合。** 下列 CI 為 AUROC 的條件區間。B＋L 是設定中列明的主要比较；B＋T 是檢視候選結果後指定的核心，沒有改寫原主要比較。','']
cs={t:x.set_index('method') for t,x in core.items()}
lines.append(table(['模型','mRS AUROC [95% CI]','NIHSS AUROC [95% CI]'],[[label(m),ci(cs['mRS'].loc[m]),ci(cs['NIHSS'].loc[m])] for m in ['V1_base','V1_base_L','V1_base_T','V1_base_W']]))
lines+=['','**表 R3．各模型相對自身基準的配對 ΔAUROC。** ET 與 LR 的增量使用各自基準，沒有混用 LightGBM 的 B。','']
pr={t:x.set_index('method') for t,x in pairs.items()}
methods=[m for m in pr['mRS'].index if m!=pr['mRS'].loc[m,'reference']]
lines.append(table(['方法','比較基準','mRS ΔAUROC [95% CI]','NIHSS ΔAUROC [95% CI]'],[[label(m),label(pr['mRS'].loc[m,'reference']),delta(pr['mRS'].loc[m]),delta(pr['NIHSS'].loc[m])] for m in methods]))
lines+=['','B＋T 的兩個配對增量區間皆跨 0，因此不能確認加入 T 帶來穩定改善。這個比較是完整融合模型的增量，不能直接解釋成紋理的獨立貢獻。','',
'**表 R4．主分析全部 24 種方法的 AUROC。**','']
lines.append(table(['方法','mRS','NIHSS'],[[label(m),num(cs['mRS'].loc[m,'auroc']),num(cs['NIHSS'].loc[m,'auroc'])] for m in cs['mRS'].index]))
lines+=['','## 26. 四任務分類器比較','',
'C＝臨床、R＝全病灶 radiomics、S＝AAL3 病灶負荷、T＝局部病灶 radiomics 加共同指標、Fwhole＝完整腦區 Original radiomics。下列比較各分類器使用自身的重要性做 RFE；並非所有分類器共用同一份入選特徵。','']
bytask={t:d[d.task.eq(t)].set_index('method') for t in tasks}
for n,t in enumerate(tasks,5):
    lines += [f'**表 R{n}．{tasklabels[t]}：分類器 × 單一特徵組的 AUROC。**','',table(['分類器','C','R','S','T','Fwhole'],[[{'RF':'隨機森林'}.get(m,m)]+[num(bytask[t].loc[f'{m}_RFE_{b}','auroc']) for b in ['C','R','S','T','Fwhole']] for m in ['RF','XGBoost','CatBoost','LightGBM']]),'']
lines+=['## 27. 特徵選擇與特徵組合','',
'**表 R9．LightGBM 的 RFE 與單次重要性篩選。** 原結果中的 backward 對應此處「單次篩選」，不表示逐步向後消除。','']
lines.append(table(['特徵','選擇法']+[tasklabels[t] for t in tasks],[[b,{'RFE':'RFE','backward':'單次篩選'}[sel]]+[num(bytask[t].loc[f'LightGBM_{sel}_{b}','auroc']) for t in tasks] for b in ['C','R','S','T','Fwhole'] for sel in ['RFE','backward']]))
for n,b in [(10,'T'),(11,'Fwhole')]:
    combinations=[b,b+'_C',b+'_R',b+'_S',b+'_C_R',b+'_C_R_S']
    lines += ['',f'**表 R{n}．{b} 與其他特徵組合的 AUROC，固定使用 LightGBM＋RFE。**','',table(['組合']+[tasklabels[t] for t in tasks],[[m.replace('_','＋')]+[num(bytask[t].loc['LightGBM_RFE_'+m,'auroc']) for t in tasks] for m in combinations])]
lines+=['','**表 R12．局部區域特徵與僅病灶／可用性指標的比較。** 此處 flags 是單分支，未加入 B，不能與 B＋flags 的消融混用。','',table(['特徵']+[tasklabels[t] for t in tasks],[[b]+[num(bytask[t].loc['LightGBM_RFE_'+m,'auroc']) for t in tasks] for b,m in [('T','T'),('僅病灶／可用性指標','flags')]]),'',
'T 與僅指標模型在多個任務上的 AUC 接近，但此處的點估計不能取代直接配對的不確定性分析。整體結果不支持單凭 T 的模型名稱就將效能歸因於紋理。','',
'## 28. 固定融合模型的完整指標','',
'**表 R13．四任務 LightGBM＋RFE、T＋C＋R＋S 的固定比較。**','',table(['任務','n','AUROC [95% CI]','AUPRC','Accuracy','Sensitivity','Specificity','Brier'],[[tasklabels[t],int(bytask[t].loc['LightGBM_RFE_T_C_R_S','n']),ci(bytask[t].loc['LightGBM_RFE_T_C_R_S'])]+[num(bytask[t].loc['LightGBM_RFE_T_C_R_S',k]) for k in ['auprc','accuracy','sensitivity','specificity','brier']] for t in tasks]),'',
'T/E 的 accuracy 與 sensitivity 差異明顯，應同時呈現；不因其 accuracy 較高就推論正類辨識良好。四任務數字均屬同一研究分析集合的內部評估。','',
'## 29. 完整指標附表：主分析全部方法','',
'以下保留每個方法的 AUROC 區間、mean-repeat AUROC、AUPRC、accuracy、sensitivity、specificity 與 Brier。方法名稱只描述分類器和輸入，不使用開發版本代號。','']
headers=['方法','AUROC [95% CI]','Mean-repeat AUROC','AUPRC','Accuracy','Sensitivity','Specificity','Brier']
keys=['repeat_auroc_mean','auprc','accuracy','sensitivity','specificity','brier']
for n,t in enumerate(tasks[:2],14):
    lines += [f'**表 R{n}．{tasklabels[t]} 主分析全部 24 種方法。**','',table(headers,[[label(r.method),ci(r)]+[num(getattr(r,k)) for k in keys] for r in core[t].itertuples()]),'']
lines+=['## 30. 完整指標附表：四任務全部補充方法','',
'每個任務完整列出 36 種固定方法。下一節另列各方法相對 C-only 的增量；不能把該 reference 誤認為主分析的 B。','']
for n,t in enumerate(tasks,16):
    lines += [f'**表 R{n}．{tasklabels[t]} 補充比較全部 36 種方法。**','',table(headers,[[label(r.method),ci(r)]+[num(getattr(r,k)) for k in keys] for r in d[d.task.eq(t)].itertuples()]),'']
lines+=['## 31. 完整配對增量附表：補充比較','',
'以下 reference 統一為該任務的 LightGBM＋RFE、C-only。其他分類器的 ΔAUROC 在此也是對此固定 reference，不是其自身 C-only 或 B 的增量。所有區間均為固定 OOF 的條件區間，未校正多重比較。','']
for n,t in enumerate(tasks,20):
    rows=d[d.task.eq(t)];assert rows.reference.eq('LightGBM_RFE_C').all()
    lines += [f'**表 R{n}．{tasklabels[t]}：全部方法相對 LightGBM 臨床模型的 ΔAUROC。**','',table(['方法','ΔAUROC [95% CI]'],[[label(r.method),delta(r)] for r in rows.itertuples()]),'']
lines+=['以上數值直接取自保存的模型評估資料，統一顯示至小數點後四位。圖像仍維持 pending；沒有在本附件嵌入 SHAP、ROC 或腦區圖。','']
s=DOC.read_text().split('## 24. 結果表：')[0].rstrip()
s=s.replace('# 研究方法：病灶相關腦區放射體學與臨床資訊融合預測','# 研究方法與結果：病灶相關腦區放射體學與臨床資訊融合預測',1)
s=s.replace('全文中的圖與表均以 **pending** 標示預定位置。方法參數寫在正文；尚待研究團隊提供的收案、倫理及影像取得資訊，集中列於資料來源與文末的待補項目。','方法章的示意圖及設計表保留 **pending** 位置。第 24 節起附上已完成的結果表 R1–R23，包含完整評估指標；不嵌入 SHAP、ROC 或腦區圖片。尚待研究團隊提供的收案、倫理及影像取得資訊，集中列於資料來源與第 23 節。')
s=s.replace('图','圖').replace('比較','比較')
s+='\n\n'+'\n'.join(lines)
s=s.replace('主要比较','主要比較').replace('單凭','單憑')
DOC.write_text(s)
auditpath=OUT/'METHODS_V6_DETAILED.audit.json';audit=json.loads(auditpath.read_text())
audit.update(sha256=hashlib.sha256(DOC.read_bytes()).hexdigest(),characters=len(s),lines=len(s.splitlines()),numbered_sections=31,filled_results_tables=23,filled_tables=23,embedded_figures=0,complete_main_metric_rows=48,complete_supplementary_metric_rows=144,quantitative_sources_sha256={n:hashlib.sha256((OUT/n).read_bytes()).hexdigest() for n in used_sources})
auditpath.write_text(json.dumps(audit,ensure_ascii=False,indent=2)+'\n')
assert not re.search(r'V[2-9]\b|V1_',s)
assert '![' not in s
assert set(re.findall(r'\*\*表 R(\d+)．',s))==set(map(str,range(1,24)))
print('Updated',DOC,'chars',len(s),'23 results tables; 48 main and 144 supplementary complete metric rows')
