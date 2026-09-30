#!/usr/bin/env python3
"""Outcome-blind V5 whole-region + V6 lesion-region extraction for an arbitrary MNI atlas.

Per region this reproduces exactly what V5 (unclip_p99 whole-region Original radiomics,
lesion fraction) and V6 (lesion first-order, homologous contrast, lesion-intersection
texture) compute for AAL3; only the label image and label names change.
No clinical data are read.
"""
import argparse
import json
import re
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import numpy as np
import SimpleITK as sitk
from radiomics import featureextractor

HERE=Path(__file__).resolve().parent;V6=HERE.parents[1];ROOT=V6.parent;INPUTS=V6/'inputs'
sys.path.insert(0,str(V6))
import extract_v6 as v6
eb=v6.eb
ATLASES={'AAL3':('data/MNI/AAL3_MNI.nii.gz','data/MNI/AAL3.txt'),
         'BA':('data/MNI/BA_MNI.nii.gz','data/MNI/BA.txt'),
         'JHU':('data/MNI/JHU_MNI.nii.gz','data/MNI/JHU.txt'),
         'JUELICH_25':('data/MNI/JUELICH_25_MNI.nii.gz','data/MNI/JUELICH_25.txt')}


def label_names(atlas,present):
    labels=eb.labels_from_file(ROOT/ATLASES[atlas][1])
    return {r:labels.get(r,f'{atlas}_{r:03d}') for r in present}


def paired_labels(labels):
    """Pair by anatomical name; tolerate JUELICH trailing '_' and JHU space/underscore typos."""
    key={r:re.sub(r'[\s_]+','_',n).strip('_') for r,n in labels.items()}
    inverse={k:r for r,k in key.items()}
    out={}
    for r,k in key.items():
        if k.endswith(('_L','_R')):
            other=k[:-2]+('_R' if k.endswith('_L') else '_L')
            if other in inverse:out[r]=inverse[other]
    return out


def load_atlas(atlas,reference):
    """Resample to the DWI grid by nearest neighbour; must be an exact voxel permutation."""
    img=sitk.ReadImage(str(ROOT/ATLASES[atlas][0]))
    if eb.same_geometry(img,reference):return sitk.GetArrayFromImage(img).astype(np.int32)
    out=sitk.Resample(img,reference,sitk.Transform(),sitk.sitkNearestNeighbor,0,img.GetPixelID())
    a=sitk.GetArrayFromImage(img);b=sitk.GetArrayFromImage(out)
    ua,ca=np.unique(a,return_counts=True);ub,cb=np.unique(b,return_counts=True)
    if not (np.array_equal(ua,ub) and np.array_equal(ca,cb)):raise ValueError(atlas+' resampling is not lossless')
    return b.astype(np.int32)


def subject(sid,atlas,output):
    sitk.ProcessObject.SetGlobalDefaultNumberOfThreads(1);eb.configure_radiomics_logging()
    p=json.loads((V6/'config/protocol.json').read_text())
    vp=json.loads((INPUTS/'v5_protocol.json').read_text())
    folder=ROOT/'data/origin data'/sid
    prov={'atlas':atlas,'atlas_sha256':eb.sha256(ROOT/ATLASES[atlas][0]),
          'labels_sha256':eb.sha256(ROOT/ATLASES[atlas][1]),
          'v6_protocol':eb.canonical_hash(p),'v5_protocol':eb.canonical_hash(vp),
          'code':eb.sha256(Path(__file__)),'v6_extractor':eb.sha256(Path(v6.__file__)),
          'helper':eb.sha256(Path(eb.__file__)),
          'image':eb.sha256(folder/vp['image']),'lesion':eb.sha256(folder/vp['lesion']),
          'brain':eb.sha256(ROOT/vp['brain']),
          'numpy':np.__version__,'simpleitk':sitk.Version_VersionString(),'pyradiomics':eb.radiomics.__version__}
    path=Path(output)/'checkpoints'/f'{sid}.json.gz'
    if path.exists():
        d=eb.read(path)
        if d['provenance']!=prov:raise ValueError('Resume mismatch '+sid)
        return {'ID':sid,'resumed':True}
    image=sitk.ReadImage(str(folder/vp['image']),sitk.sitkFloat32)
    lm=sitk.ReadImage(str(folder/vp['lesion']));brain=sitk.ReadImage(str(ROOT/vp['brain']))
    assert all(eb.same_geometry(image,i) for i in (lm,brain))
    aa=load_atlas(atlas,image)
    values=sitk.GetArrayFromImage(image);lesion=sitk.GetArrayFromImage(lm)>.5
    ba=sitk.GetArrayFromImage(brain)>0
    valid=np.isfinite(values)&(values>0)
    arr=eb.normalize(values,ba,lesion,vp)[0]['unclip_p99']
    exclusion=sitk.GetArrayFromImage(sitk.BinaryDilate(sitk.GetImageFromArray(lesion.astype('uint8')),[2]*3))>0
    regions=sorted(int(r) for r in np.unique(aa) if r>0)
    labels=label_names(atlas,regions);pairs=paired_labels(labels)
    ext=featureextractor.RadiomicsFeatureExtractor(**vp['radiomics'])
    ext.disableAllImageTypes();ext.enableImageTypeByName('Original');ext.disableAllFeatures()
    for family in vp['feature_classes']:ext.enableFeatureClassByName(family)
    rows=[]
    for region in regions:
        whole=aa==region;n=int(whole.sum());hit=int((whole&lesion).sum())
        row={'region':region,'name':labels[region],'roi_voxels':n,'lesion_fraction':hit/max(n,1),
             'whole':None,'whole_skip':None}
        # V5 whole-region rule (unclip_p99 only).
        if n>vp['minimum_roi_voxels_exclusive'] and float((whole&valid).sum()/max(n,1))>=vp['minimum_roi_coverage']:
            row['whole']=eb.radiomics_roi(arr,whole,image,ext,vp)
        else:row['whole_skip']='prespecified_size_or_coverage_skip'
        # V6 lesion-region rule.
        raw=whole&lesion;mask=raw&valid;nv=int(mask.sum())
        if raw.sum() and nv/raw.sum()<.9:raise ValueError(f'Insufficient lesion ROI coverage {sid}/{region}')
        row.update({'lesion_voxels':int(raw.sum()),'valid_voxels':nv,'firstorder':None,'contrast':None,
                    'texture':None,'reference_voxels':0,'counterpart':pairs.get(region)})
        if nv>=p['small_roi_min_voxels']:
            row['firstorder']=v6.firstorder(arr[mask])
            if region in pairs:
                reference=(aa==pairs[region])&valid&~exclusion
                row['reference_voxels']=int(reference.sum())
                if reference.sum()>=p['reference_min_voxels']:
                    row['contrast']=v6.contrast(row['firstorder'],v6.firstorder(arr[reference]))
        if v6.eligible_texture(mask):row['texture']=eb.radiomics_roi(arr,mask,image,ext,vp)
        rows.append(row)
    d={'ID':sid,'status':'complete','provenance':prov,'pairs':{str(k):v for k,v in pairs.items()},'rows':rows}
    eb.write(path,d)
    return {'ID':sid,'regions':len(rows),'affected':sum(r['lesion_voxels']>0 for r in rows),
            'texture':sum(r['texture'] is not None for r in rows)}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--atlas',required=True,choices=ATLASES)
    ap.add_argument('--jobs',type=int,default=6);ap.add_argument('--limit',type=int)
    ap.add_argument('--output',type=Path);args=ap.parse_args()
    ids=np.load(INPUTS/'v5_features.npz')['ids'].tolist()
    if args.limit:ids=[ids[i] for i in np.linspace(0,len(ids)-1,args.limit).round().astype(int)]
    out=args.output or HERE/f'extract_{args.atlas}';out.mkdir(parents=True,exist_ok=True)
    manifest={'atlas':args.atlas,'ids':ids,'limit':args.limit}
    mp=out/'manifest.json'
    if mp.exists() and json.loads(mp.read_text())!=manifest:raise ValueError('Manifest mismatch')
    mp.write_text(json.dumps(manifest,indent=2))
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures=[pool.submit(subject,s,args.atlas,out) for s in ids]
        for i,f in enumerate(as_completed(futures)):print(i+1,len(ids),f.result(),flush=True)
    (out/'COMPLETE.json').write_text(json.dumps(manifest))

if __name__=='__main__':main()
