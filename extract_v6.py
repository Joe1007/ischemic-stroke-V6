#!/usr/bin/env python3
"""Outcome-blind lesion/AAL3 re-extraction. No clinical data are read."""
import argparse
import json
import sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import SimpleITK as sitk
from radiomics import featureextractor

HERE=Path(__file__).resolve().parent; ROOT=HERE.parent; INPUTS=HERE/'inputs'
import extraction_base as eb

def paired_labels(labels):
    inverse={name:r for r,name in labels.items()}
    return {r:inverse[name[:-2]+('_R' if name.endswith('_L') else '_L')]
            for r,name in labels.items() if name.endswith(('_L','_R'))
            and name[:-2]+('_R' if name.endswith('_L') else '_L') in inverse}

def firstorder(x):
    x=np.asarray(x,dtype=np.float64)
    q10,q25,q50,q75,q90=np.percentile(x,[10,25,50,75,90])
    _,counts=np.unique(np.floor(x/5).astype(np.int64),return_counts=True)
    p=counts/counts.sum()
    return [float(v) for v in (x.mean(),q50,q10,q90,x.std(),q75-q25,
                               np.sqrt(np.mean(x*x)),-np.sum(p*np.log2(p)))]

def contrast(a,b):
    a=np.asarray(a);b=np.asarray(b)
    return ((a-b)/(np.abs(a)+np.abs(b)+5)).tolist()

def eligible_texture(mask):
    if mask.sum()<=64:return False
    pts=np.where(mask)
    return sum(p.max()>p.min() for p in pts)>=2

def subject(sid,output):
    sitk.ProcessObject.SetGlobalDefaultNumberOfThreads(1);eb.configure_radiomics_logging()
    p=json.loads((HERE/'config/protocol.json').read_text())
    vp=json.loads((INPUTS/'v5_protocol.json').read_text())
    source=INPUTS/'v5_checkpoints'/f'{sid}.json.gz'
    folder=ROOT/'data/origin data'/sid
    prov={'protocol':eb.canonical_hash(p),'source':eb.sha256(source),
          'code':eb.sha256(Path(__file__)),'helper':eb.sha256(Path(eb.__file__)),
          'image':eb.sha256(folder/vp['image']),'lesion':eb.sha256(folder/vp['lesion']),
          'atlas':eb.sha256(ROOT/vp['atlas']),'brain':eb.sha256(ROOT/vp['brain']),
          'labels':eb.sha256(ROOT/vp['labels']),
          'numpy':np.__version__,'simpleitk':sitk.Version_VersionString(),'pyradiomics':eb.radiomics.__version__}
    path=Path(output)/'checkpoints'/f'{sid}.json.gz'
    if path.exists():
        d=eb.read(path)
        if d['provenance']!=prov:raise ValueError('Resume mismatch '+sid)
        if d['status']!='complete':raise ValueError('Prior failed case '+sid)
        return {'ID':sid,'resumed':True}
    prior=eb.read(source)
    for key,oldkey in [('image','image_sha256'),('lesion','lesion_sha256'),('atlas','atlas_sha256'),('brain','brain_sha256')]:
        if prov[key]!=prior['provenance'][oldkey]:raise ValueError('V5 source changed '+sid)
    image=sitk.ReadImage(str(folder/vp['image']),sitk.sitkFloat32)
    lm=sitk.ReadImage(str(folder/vp['lesion']));atlas=sitk.ReadImage(str(ROOT/vp['atlas']))
    brain=sitk.ReadImage(str(ROOT/vp['brain']))
    assert all(eb.same_geometry(image,i) for i in (lm,atlas,brain))
    values=sitk.GetArrayFromImage(image);lesion=sitk.GetArrayFromImage(lm)>.5
    aa=sitk.GetArrayFromImage(atlas);ba=sitk.GetArrayFromImage(brain)>0
    valid=np.isfinite(values)&(values>0)
    arr=eb.normalize(values,ba,lesion,vp)[0]['unclip_p99']
    exclusion=sitk.GetArrayFromImage(sitk.BinaryDilate(sitk.GetImageFromArray(lesion.astype('uint8')),[2]*3))>0
    labels=eb.labels_from_file(ROOT/vp['labels']);pairs=paired_labels(labels)
    ext=featureextractor.RadiomicsFeatureExtractor(**vp['radiomics'])
    ext.disableAllImageTypes();ext.enableImageTypeByName('Original');ext.disableAllFeatures()
    for family in vp['feature_classes']:ext.enableFeatureClassByName(family)
    rows=[]
    for region in sorted(int(r) for r in np.unique(aa) if r>0):
        raw=(aa==region)&lesion;mask=raw&valid;n=int(mask.sum())
        if raw.sum() and n/raw.sum()<.9:raise ValueError(f'Insufficient lesion ROI coverage {sid}/{region}')
        row={'region':region,'name':labels[region],'lesion_voxels':int(raw.sum()),'valid_voxels':n,
             'firstorder':None,'contrast':None,'texture':None,'reference_voxels':0,
             'counterpart':pairs.get(region)}
        if n>=p['small_roi_min_voxels']:
            row['firstorder']=firstorder(arr[mask])
            if region in pairs:
                reference=(aa==pairs[region])&valid&~exclusion
                row['reference_voxels']=int(reference.sum())
                if reference.sum()>=p['reference_min_voxels']:
                    row['contrast']=contrast(row['firstorder'],firstorder(arr[reference]))
        if eligible_texture(mask):row['texture']=eb.radiomics_roi(arr,mask,image,ext,vp)
        rows.append(row)
    d={'ID':sid,'status':'complete','provenance':prov,'rows':rows}
    eb.write(path,d)
    return {'ID':sid,'affected':sum(r['lesion_voxels']>0 for r in rows),
            'firstorder':sum(r['firstorder'] is not None for r in rows),
            'texture':sum(r['texture'] is not None for r in rows)}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--jobs',type=int,default=6)
    parser.add_argument('--pilot',action='store_true');args=parser.parse_args()
    ids=np.load(INPUTS/'v5_features.npz')['ids'].tolist()
    if args.pilot:ids=[ids[i] for i in np.linspace(0,len(ids)-1,12).round().astype(int)]
    output=HERE/('pilot' if args.pilot else 'output');output.mkdir(parents=True,exist_ok=True)
    manifest={'ids':ids,'pilot':args.pilot,'protocol':eb.sha256(HERE/'config/protocol.json')}
    mp=output/'manifest.json'
    if mp.exists() and json.loads(mp.read_text())!=manifest:raise ValueError('Manifest mismatch')
    mp.write_text(json.dumps(manifest,indent=2))
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures=[pool.submit(subject,s,output) for s in ids]
        for i,f in enumerate(as_completed(futures)):print(i+1,len(ids),f.result(),flush=True)
    (output/'COMPLETE.json').write_text(json.dumps(manifest))

if __name__=='__main__':main()
