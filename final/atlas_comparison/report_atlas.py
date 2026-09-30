#!/usr/bin/env python3
"""Audit the four-atlas runs, compute paired comparisons and write section 32 of METHODS_V6_DETAILED.md."""
import hashlib, json, re, sys
from pathlib import Path
import numpy as np
import pandas as pd
HERE=Path(__file__).resolve().parent;FINAL=HERE.parent;ROOT=FINAL.parents[1]
sys.path.insert(0,str(FINAL.parent))
import model_v6 as v6
mb=v6.mb
ATLASES=['AAL3','BA','JHU','JUELICH_25']
LABEL={'AAL3':'AAL3','BA':'Brodmann（BA）','JHU':'JHU 白質','JUELICH_25':'Jülich（25%）'}
TASKS=['mRS','NIHSS','TE','LVO'];TASKNAME={'mRS':'mRS','NIHSS':'NIHSS','TE':'T/E','LVO':'LVO'}
FLAGS={'affected','firstorder_available','contrast_available','texture_available'}
SEED=20260909;BOOT=2000
RES=HERE/'results';RES.mkdir(exist_ok=True)
DOC=FINAL/'output/METHODS_V6_DETAILED.md'

def tables_dir(a):return FINAL/'tables_20260917' if a=='AAL3' else HERE/f'tables_{a}'
def core_dir(a):return FINAL/'output/core' if a=='AAL3' else HERE/f'core_{a}'
def f4(x):return f'{x:.4f}'
def sd(x):return f'{x:+.4f}'

def collect(a,task):
    """Averaged OOF per method, with the same integrity checks as report_tables.py."""
    details=[];shap={}
    dirs=sorted(tables_dir(a).glob(task+'_[0-9][0-9][0-9]'));dirs=[d for d in dirs if d.is_dir()]
    assert len(dirs)==15,(a,task,len(dirs))
    for dest in dirs:
        s=json.loads((dest/'split.json').read_text());assert not set(s['train_ids'])&set(s['test_ids'])
        ms=[p for p in dest.glob('*.json') if p.name!='split.json'];assert len(ms)==36,(dest,len(ms))
        for fp in ms:
            m=json.loads(fp.read_text());p=np.array(m['probability'])
            assert len(p)==len(s['test_ids']) and np.isfinite(p).all()
            details.extend({'ID':i,'y':y,'repeat':s['split']//5,'method':m['method'],'probability':q} for i,y,q in zip(s['test_ids'],s['y_test'],p))
            z=dest/(m['method']+'_shap.npz')
            if m['method'] in ('LightGBM_RFE_T','LightGBM_RFE_Fwhole'):
                z=np.load(z,allow_pickle=True)
                np.testing.assert_allclose(1/(1+np.exp(-(z['shap'].sum(1)+z['expected']))),p,atol=1e-8)
                shap.setdefault(m['method'],[]).append(pd.DataFrame(z['shap'],columns=z['names'],index=z['ids']))
    d=pd.DataFrame(details);assert d.groupby(['ID','method']).size().eq(3).all()
    avg=d.groupby(['ID','method']).probability.mean().unstack();y=d.groupby('ID').y.first().loc[avg.index]
    return y,avg,shap

def region_of(name):
    block,rest=name.split('::',1)
    if block=='T':return int(re.match(r'[A-Z0-9_]+?_(\d{3})\|',rest).group(1)) if rest.rsplit('|',1)[-1] not in FLAGS else None
    if block=='Fwhole':return int(rest.rsplit('_',1)[-1])
    return None

def top_regions(parts,names,k=3):
    sh=pd.concat(parts).fillna(0);imp=sh.abs().mean()
    reg=pd.Series({n:region_of(n) for n in imp.index}).dropna().astype(int)
    agg=imp[reg.index].groupby(reg).sum().sort_values(ascending=False).head(k)
    return [(names.get(int(r),str(r)),v) for r,v in agg.items()]

def pretty_region(a,n):
    n=n.strip('_')
    if a=='BA':return 'BA '+str(int(n.rsplit('_',1)[-1]))
    return n.replace('GM_','').replace('WM_','')

def main():
    audit={'identical_atlas_independent_methods':{},'aal3_matches_published':None,'core_baseline_identical':{}}
    meta={a:json.loads((HERE/f'features_{a}.json').read_text()) for a in ATLASES}
    rnames={}
    for a in ATLASES:
        z=np.load(HERE/f'features_{a}.npz');rnames[a]={int(r):str(n) for r,n in zip(z['region_ids'],z['region_names'])}
    published=pd.read_csv(FINAL/'output/all_metrics.csv')
    rows=[];deltas=[];tops=[]
    for task in TASKS:
        ys={};preds={};shaps={}
        for a in ATLASES:ys[a],preds[a],shaps[a]=collect(a,task)
        for a in ATLASES[1:]:
            assert ys[a].index.equals(ys['AAL3'].index) and ys[a].equals(ys['AAL3'])
            same=[m for m in preds[a] if re.fullmatch(r'[A-Za-z]+_(RFE|backward)_(C|R)',m)]
            diff=max(float(np.abs(preds[a][m]-preds['AAL3'][m]).max()) for m in same)
            audit['identical_atlas_independent_methods'][f'{a}/{task}']={'methods':len(same),'max_abs_probability_difference':diff}
            assert diff<1e-9,(a,task,diff)
        y=ys['AAL3'].to_numpy()
        for a in ATLASES:
            summ=pd.DataFrame(mb.uncertainty(y,{m:preds[a][m].to_numpy() for m in preds[a]},'LightGBM_RFE_C',BOOT,SEED))
            summ['atlas']=a;summ['task']=task;summ['n']=len(y);rows.append(summ)
            if a=='AAL3':
                pub=published[published.task.eq(task)].set_index('method').auroc
                np.testing.assert_allclose(summ.set_index('method').auroc.loc[pub.index],pub,atol=1e-12)
                audit['aal3_matches_published']=True
        for m in preds['AAL3'].columns:
            if re.fullmatch(r'[A-Za-z]+_(RFE|backward)_(C|R)',m):continue
            r=pd.DataFrame(mb.uncertainty(y,{a:preds[a][m].to_numpy() for a in ATLASES},'AAL3',BOOT,SEED))
            r['method']=m;r['task']=task;r['atlas']=list(ATLASES);deltas.append(r)
        for a in ATLASES:
            for m in ('LightGBM_RFE_Fwhole','LightGBM_RFE_T'):
                for rank,(n,v) in enumerate(top_regions(shaps[a][m],rnames[a]),1):
                    tops.append({'task':task,'atlas':a,'method':m,'rank':rank,'region':pretty_region(a,n),'mean_abs_shap_sum':v})
    allm=pd.concat(rows);allm.to_csv(RES/'atlas_all_metrics.csv',index=False)
    dl=pd.concat(deltas)[['task','method','atlas','auroc','auc_low','auc_high','delta_auc','delta_low','delta_high']]
    dl.to_csv(RES/'atlas_paired_delta_vs_AAL3.csv',index=False)
    tp=pd.DataFrame(tops);tp.to_csv(RES/'atlas_top_shap_regions.csv',index=False)
    # Core (mRS / NIHSS): B+L/T/W with fixed B.
    core=[];coredelta=[]
    for task in ['mRS','NIHSS']:
        cp={a:pd.read_csv(core_dir(a)/f'{task}_oof_predictions.csv').set_index('ID') for a in ATLASES}
        for a in ATLASES[1:]:
            assert cp[a].index.equals(cp['AAL3'].index)
            d=float(np.abs(cp[a]['V1_base']-cp['AAL3']['V1_base']).max());audit['core_baseline_identical'][f'{a}/{task}']=d;assert d<1e-12
        for a in ATLASES:
            s=pd.read_csv(core_dir(a)/f'{task}_summary.csv');s['atlas']=a;core.append(s)
            m=pd.read_csv(core_dir(a)/f'{task}_matched_comparisons.csv');m['atlas']=a;m['task']=task;coredelta.append(m)
        y=cp['AAL3'].y.to_numpy()
        for m in ['V1_L','V1_T','V1_W','V1_base_L','V1_base_T','V1_base_W']:
            r=pd.DataFrame(mb.uncertainty(y,{a:cp[a][m].to_numpy() for a in ATLASES},'AAL3',BOOT,SEED))
            r['atlas']=list(ATLASES);r['method']=m;r['task']=task;deltas.append(r)
    core=pd.concat(core);core.to_csv(RES/'core_all_metrics.csv',index=False)
    coredelta=pd.concat(coredelta);coredelta.to_csv(RES/'core_matched_vs_B.csv',index=False)
    coreatlas=pd.concat(deltas[-12:])[['task','method','atlas','auroc','delta_auc','delta_low','delta_high']]
    coreatlas.to_csv(RES/'core_paired_delta_vs_AAL3.csv',index=False)
    (RES/'audit.json').write_text(json.dumps(audit,indent=2))
    write_doc(meta,allm,dl,tp,core,coredelta,coreatlas)

def write_doc(meta,allm,dl,tp,core,coredelta,coreatlas):
    A=ATLASES;L=[]
    get=lambda a,t,m:allm[(allm.atlas==a)&(allm.task==t)&(allm.method==m)].iloc[0]
    L+=['## 32. 四種腦區圖譜的比較','',
      '本節把前述流程中的 AAL3 分別替換為 Brodmann 分區（BA）、JHU 白質纖維束圖譜（JHU）與 Jülich 細胞構築機率圖譜 25% 閾值版本（Jülich），比較腦區定義對結果的影響。除了圖譜影像與腦區名稱外，影像、病灶遮罩、強度正規化、ROI 規則、特徵定義、病人分割、隨機種子、特徵選擇、分類器與評估方式全部與 AAL3 分析相同。',
      '',
      '### 32.1 替換範圍與前處理','',
      '會隨圖譜改變的特徵只有 S、L、T、W、Fwhole 與共同指標（flags）；C 與全病灶 R 不使用圖譜，四種圖譜完全相同。第 25 節的強基準 B 內含 AAL3 病灶負荷，本節刻意**保持 B 不變**，因此各圖譜的 B＋L／T／W 只差在加入的區域特徵；第 26–27 節的補充比較則以各圖譜自己的 S 取代 AAL3 的 S。',
      '',
      '通用抽取程式先以 AAL3 驗證：在 198 位病人、每人 166 個腦區上，重新計算的完整腦區特徵、病灶負荷與 L／T 原始數值，和主分析使用的矩陣逐位相同。BA 圖譜的 x 軸方向與 DWI 相反，以最近鄰重採樣對齊；此操作只是體素重新排列，各標籤的體素數完全不變。',
      '',
      '對側差異需要左右成對的腦區，依名稱中的 _L／_R 配對。BA 圖譜的每個分區都同時涵蓋左右半球，沒有可配對的區域，因此 BA 的 L 只有局部一階統計，沒有對側差異，contrast_available 指標也恆為 0。JHU 與 Jülich 的中線結構（例如胼胝體、穹窿）同樣沒有對側參照。',
      '',
      '**表 R24．四種圖譜的區域數與特徵可用性。** 可用區域數為 198 位病人的中位數；「無紋理病人」指沒有任何腦區達到 T 紋理門檻（大於 64 個有效體素）的病人數。',
      '',
      '| 圖譜 | 腦區數 | 左右配對腦區 | Fwhole 欄位 | T 欄位（含指標） | 受影響腦區 | 有一階統計 | 有對側差異 | 有紋理 | 無紋理病人 | 無一階統計病人 |',
      '|---|---|---|---|---|---|---|---|---|---|---|']
    for a in A:
        m=meta[a]
        L.append(f"| {LABEL[a]} | {m['regions']} | {m['paired_regions']} | {m['columns']['Fwhole']} | {m['columns']['T']} | {m['median_affected_regions']:g} | {m['median_firstorder_regions']:g} | {m['median_contrast_regions']:g} | {m['median_texture_regions']:g} | {m['patients_without_texture']} | {m['patients_without_firstorder']} |")
    L+=['','BA 的腦區最大（中位數約 3,100 體素），JHU 最小（約 290 體素）。腦區越大，完整腦區特徵混入的正常組織越多；腦區越小，病灶越容易被切成低於紋理門檻的碎片。四種圖譜中約 39%–43% 的病人沒有任何可計算的局部紋理，這個限制並非 AAL3 特有。','',
      '### 32.2 主分析：強基準加入各圖譜的區域特徵','',
      '**表 R25．LightGBM：B 加入各圖譜區域特徵的 AUROC 與相對 B 的配對 ΔAUROC。** B 在四種圖譜中完全相同（mRS 0.7581、NIHSS 0.7124）。',
      '','| 模型 | 圖譜 | mRS AUROC | mRS ΔAUROC vs B [95% CI] | NIHSS AUROC | NIHSS ΔAUROC vs B [95% CI] |','|---|---|---|---|---|---|']
    names={'V1_base_L':'B＋L','V1_base_T':'B＋T','V1_base_W':'B＋W','V1_L':'L 單分支','V1_T':'T 單分支','V1_W':'W 單分支'}
    for m in ['V1_base_L','V1_base_T','V1_base_W']:
        for a in A:
            cells=[]
            for t in ['mRS','NIHSS']:
                r=coredelta[(coredelta.atlas==a)&(coredelta.task==t)&(coredelta.method==m)].iloc[0]
                cells+= [f4(r.auroc),f'{sd(r.delta_auc)} [{sd(r.delta_low)}, {sd(r.delta_high)}]']
            L.append(f'| {names[m]} | {LABEL[a]} | '+' | '.join(cells)+' |')
    L+=['','**表 R26．LightGBM 區域特徵單分支的 AUROC。** 此表不含 B，只看區域特徵本身的辨識力。','',
        '| 模型 | 任務 | '+' | '.join(LABEL[a] for a in A)+' |','|---|---|'+'---|'*len(A)]
    for m in ['V1_L','V1_T','V1_W']:
        for t in ['mRS','NIHSS']:
            L.append(f'| {names[m]} | {t} | '+' | '.join(f4(core[(core.atlas==a)&(core.task==t)&(core.method==m)].auroc.iloc[0]) for a in A)+' |')
    L+=['','**表 R27．主分析各圖譜相對 AAL3 的配對 ΔAUROC。** 同一模型、同一批病人與同一組 OOF 分割，只比較圖譜。','',
        '| 模型 | 任務 | BA [95% CI] | JHU [95% CI] | Jülich [95% CI] |','|---|---|---|---|---|']
    for m in ['V1_base_L','V1_base_T','V1_base_W','V1_L','V1_T','V1_W']:
        for t in ['mRS','NIHSS']:
            s=coreatlas[(coreatlas.task==t)&(coreatlas.method==m)].set_index('atlas')
            L.append(f'| {names[m]} | {t} | '+' | '.join(f'{sd(s.loc[a,"delta_auc"])} [{sd(s.loc[a,"delta_low"])}, {sd(s.loc[a,"delta_high"])}]' for a in A[1:])+' |')
    L+=['','### 32.3 四任務補充比較','',
        '**表 R28．LightGBM（RFE）各圖譜的區域特徵 AUROC。** C 與 R 不受圖譜影響，列於最上方作參照；經核對，這些方法在四種圖譜下的逐人預測機率完全相同。','',
        '| 特徵 | 圖譜 | mRS | NIHSS | T/E | LVO |','|---|---|---|---|---|---|']
    for m,lab in [('LightGBM_RFE_C','C'),('LightGBM_RFE_R','R')]:
        L.append(f'| {lab} | （不適用） | '+' | '.join(f4(get('AAL3',t,m).auroc) for t in TASKS)+' |')
    feats=[('S','S'),('T','T'),('Fwhole','Fwhole'),('flags','僅共同指標'),('T_C','T＋C'),('T_R','T＋R'),('T_S','T＋S'),('T_C_R','T＋C＋R'),('T_C_R_S','T＋C＋R＋S'),
           ('Fwhole_C','Fwhole＋C'),('Fwhole_R','Fwhole＋R'),('Fwhole_S','Fwhole＋S'),('Fwhole_C_R','Fwhole＋C＋R'),('Fwhole_C_R_S','Fwhole＋C＋R＋S')]
    for key,lab in feats:
        for a in A:L.append(f'| {lab} | {LABEL[a]} | '+' | '.join(f4(get(a,t,'LightGBM_RFE_'+key).auroc) for t in TASKS)+' |')
    L+=['','**表 R29．其他分類器（RFE）的 S／T／Fwhole 單一特徵組 AUROC。**','',
        '| 分類器 | 特徵 | 圖譜 | mRS | NIHSS | T/E | LVO |','|---|---|---|---|---|---|---|']
    for kind,klab in [('RF','隨機森林'),('XGBoost','XGBoost'),('CatBoost','CatBoost')]:
        for key in ['S','T','Fwhole']:
            for a in A:L.append(f'| {klab} | {key} | {LABEL[a]} | '+' | '.join(f4(get(a,t,f'{kind}_RFE_{key}').auroc) for t in TASKS)+' |')
    L+=['','**表 R30．LightGBM 單次重要性篩選（backward）的 S／T／Fwhole AUROC。**','',
        '| 特徵 | 圖譜 | mRS | NIHSS | T/E | LVO |','|---|---|---|---|---|---|']
    for key in ['S','T','Fwhole']:
        for a in A:L.append(f'| {key} | {LABEL[a]} | '+' | '.join(f4(get(a,t,'LightGBM_backward_'+key).auroc) for t in TASKS)+' |')
    L+=['','**表 R31．主要 LightGBM 模型各圖譜相對 AAL3 的配對 ΔAUROC [95% CI]。**','',
        '| 特徵 | 任務 | BA | JHU | Jülich |','|---|---|---|---|---|']
    for key,lab in [('S','S'),('T','T'),('Fwhole','Fwhole'),('T_C_R_S','T＋C＋R＋S'),('Fwhole_C_R_S','Fwhole＋C＋R＋S')]:
        for t in TASKS:
            s=dl[(dl.task==t)&(dl.method=='LightGBM_RFE_'+key)].set_index('atlas')
            L.append(f'| {lab} | {TASKNAME[t]} | '+' | '.join(f'{sd(s.loc[a,"delta_auc"])} [{sd(s.loc[a,"delta_low"])}, {sd(s.loc[a,"delta_high"])}]' for a in A[1:])+' |')
    L+=['','**表 R32．各圖譜最佳 LightGBM 組合。** 每個任務、每種圖譜在 36 種補充方法中 AUROC 最高者。這是看過結果後的排序，只能作探索描述，不能當成預先指定的比較。','',
        '| 任務 | 圖譜 | 最高 AUROC 的方法 | AUROC [95% CI] |','|---|---|---|---|']
    for t in TASKS:
        for a in A:
            s=allm[(allm.atlas==a)&(allm.task==t)].sort_values('auroc',ascending=False).iloc[0]
            L.append(f'| {TASKNAME[t]} | {LABEL[a]} | {s.method} | {f4(s.auroc)} [{f4(s.auc_low)}, {f4(s.auc_high)}] |')
    L+=['','### 32.4 重要腦區','',
        '**表 R33．held-out SHAP 加總最高的前三個腦區。** 計算方式同第 19 節：同一腦區的 radiomics 特徵平均絕對 SHAP 加總，不含共同指標。四種圖譜的腦區定義不同，名稱不能逐一對應。','',
        '| 任務 | 圖譜 | Fwhole 前三名 | T 前三名 |','|---|---|---|---|']
    for t in TASKS:
        for a in A:
            cell=lambda m:'、'.join(tp[(tp.task==t)&(tp.atlas==a)&(tp.method==m)].sort_values('rank').region)
            L.append(f'| {TASKNAME[t]} | {LABEL[a]} | {cell("LightGBM_RFE_Fwhole")} | {cell("LightGBM_RFE_T")} |')
    L+=['','### 32.5 判讀與限制','']+interpretation(allm,dl,coredelta,coreatlas)+['',
        '本節四種圖譜使用相同病人、相同分割，並在看到結果前固定了全部方法；但同一 cohort 已反覆分析，且本節一次比較了四種圖譜與數十種方法，未作多重比較校正。任何圖譜的最高分都只能視為探索性結果；若要主張某一圖譜較佳，需要預先指定比較並以外部資料驗證。95% CI 與前文相同，是固定 OOF 的條件 bootstrap 區間。','',
        '數值檔：`atlas_comparison/results/` 內的 atlas_all_metrics.csv（全部方法完整指標）、atlas_paired_delta_vs_AAL3.csv、core_all_metrics.csv、core_matched_vs_B.csv、core_paired_delta_vs_AAL3.csv、atlas_top_shap_regions.csv 與 audit.json。']
    s=DOC.read_text();marker='\n## 32. 四種腦區圖譜的比較'
    if marker in s:s=s[:s.index(marker)]
    s=s.rstrip('\n')+'\n\n'+'\n'.join(L)+'\n'
    assert not re.search(r'V[2-9]\b|V1_',s) and '![' not in s
    DOC.write_text(s)
    ap=FINAL/'output/METHODS_V6_DETAILED.audit.json';au=json.loads(ap.read_text())
    au.update(sha256=hashlib.sha256(DOC.read_bytes()).hexdigest(),characters=len(s),lines=len(s.splitlines()),numbered_sections=32,
              filled_results_tables=33,filled_tables=33,atlas_comparison={'section':32,'tables':'R24-R33','atlases':ATLASES,
              'sources_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(RES.glob('*'))}})
    ap.write_text(json.dumps(au,ensure_ascii=False,indent=2)+'\n')
    mp=FINAL/'output/MANIFEST.json'
    if mp.exists():
        man=json.loads(mp.read_text())
        for n in ['METHODS_V6_DETAILED.md','METHODS_V6_DETAILED.audit.json']:
            if n in man:man[n]=hashlib.sha256((FINAL/'output'/n).read_bytes()).hexdigest()
        mp.write_text(json.dumps(man,ensure_ascii=False,indent=2))
    print('Wrote section 32 to',DOC)

def interpretation(allm,dl,coredelta,coreatlas):
    """Data-driven summary sentences; no claims beyond the paired intervals."""
    out=[]
    for t in ['mRS','NIHSS']:
        r=coredelta[(coredelta.task==t)&coredelta.method.isin(['V1_base_L','V1_base_T','V1_base_W'])&(coredelta.atlas!='AAL3')]
        pos=r[r.delta_low>0];neg=r[r.delta_high<0]
        out.append(f'- {t} 主分析：非 AAL3 圖譜的 9 個 B＋區域特徵組合中，ΔAUROC 區間完全大於 0 的有 {len(pos)} 個'+
                   (f'（{"、".join(LABEL[x.atlas]+" "+x.method.replace("V1_base_","B＋") for x in pos.itertuples())}）' if len(pos) else '')+
                   f'，完全小於 0 的有 {len(neg)} 個'+(f'（{"、".join(LABEL[x.atlas]+" "+x.method.replace("V1_base_","B＋") for x in neg.itertuples())}）' if len(neg) else '')+f'；最大 ΔAUROC 為 {sd(r.delta_auc.max())}。')
    for t in TASKS:
        s=dl[(dl.task==t)&dl.method.isin(['LightGBM_RFE_S','LightGBM_RFE_T','LightGBM_RFE_Fwhole'])&(dl.atlas!='AAL3')]
        sig=s[(s.delta_low>0)|(s.delta_high<0)]
        out.append(f'- {TASKNAME[t]}：S／T／Fwhole 單一特徵組相對 AAL3 的 9 個配對比較中，區間不含 0 的有 {len(sig)} 個'+
                   (('（'+'、'.join(LABEL[x.atlas]+' '+x.method.replace('LightGBM_RFE_','')+' '+sd(x.delta_auc) for x in sig.itertuples())+'）') if len(sig) else '')+'。')
    return out

if __name__=='__main__':main()
