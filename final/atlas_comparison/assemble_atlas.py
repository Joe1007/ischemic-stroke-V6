#!/usr/bin/env python3
"""Assemble per-atlas S/Fwhole/L/T/W matrices in the exact V5/V6 layout. No outcomes read."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;V6=HERE.parents[1];ROOT=V6.parent;INPUTS=V6/'inputs'
FLAGS=['affected','firstorder_available','contrast_available','texture_available']

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def assemble(atlas):
    directory=HERE/f'extract_{atlas}';manifest=json.loads((directory/'COMPLETE.json').read_text())
    assert manifest['atlas']==atlas and manifest['limit'] is None
    old=np.load(INPUTS/'v5_features.npz')
    assert manifest['ids']==old['ids'].tolist()
    features=old['global_names'].tolist();stats=json.loads((V6/'config/protocol.json').read_text())['firstorder']
    patients=[]
    for sid in manifest['ids']:
        with gzip.open(directory/'checkpoints'/f'{sid}.json.gz','rt') as h:d=json.load(h)
        assert d['ID']==sid and d['status']=='complete' and d['provenance']['atlas']==atlas
        patients.append(d)
    regions=[r['region'] for r in patients[0]['rows']];names={r['region']:r['name'] for r in patients[0]['rows']}
    n=len(patients);R=len(regions)
    F=np.full((n,R*93),np.nan,dtype='float32');burden=np.zeros((n,R),dtype='float32')
    l=np.zeros((n,R*16),dtype='float32');t=np.zeros((n,R*93),dtype='float32');avail=np.zeros((n,R*4),dtype='float32')
    for i,d in enumerate(patients):
        assert [r['region'] for r in d['rows']]==regions
        for j,row in enumerate(d['rows']):
            burden[i,j]=row['lesion_fraction']
            if row['whole'] is not None:
                assert sorted(row['whole'])==features;F[i,j*93:(j+1)*93]=[row['whole'][f] for f in features]
            a=row['lesion_voxels']>0;fo=row['firstorder'];co=row['contrast'];tex=row['texture']
            avail[i,j*4:(j+1)*4]=[a,fo is not None,co is not None,tex is not None]
            if a:
                l[i,j*16:j*16+8]=fo if fo is not None else np.nan
                l[i,j*16+8:(j+1)*16]=co if co is not None else np.nan
                t[i,j*93:(j+1)*93]=[tex[f] for f in features] if tex is not None else np.nan
    fnames=[f'{f}_{r:03d}' for r in regions for f in features]
    ln=[f'{atlas}_{r:03d}|{kind}|{f}' for r in regions for kind in ['lesion','contrast'] for f in stats]
    tn=[f'{atlas}_{r:03d}|lesion|{f}' for r in regions for f in features]
    fn=[f'{atlas}_{r:03d}|{f}' for r in regions for f in FLAGS]
    w=np.where(np.repeat(burden>0,93,axis=1),F*np.repeat(burden,93,axis=1),0)
    out=HERE/f'features_{atlas}.npz'
    np.savez_compressed(out,ids=old['ids'],F_unclip_p99=F,names=np.array(fnames),burden=burden,
        burden_names=np.array([f'{atlas}_fraction_{r:03d}' for r in regions]),
        L=np.hstack([l,avail]),T=np.hstack([t,avail]),W=np.hstack([w,avail]),
        L_names=np.array(ln+fn),T_names=np.array(tn+fn),W_names=np.r_[np.char.add('occupancy_weighted|',np.array(fnames)),fn],
        region_ids=np.array(regions),region_names=np.array([names[r] for r in regions]))
    pairs=patients[0]['pairs']
    meta={'atlas':atlas,'regions':R,'paired_regions':len(pairs),
          'whole_region_computable':int(np.isfinite(F[0,::93]).sum()),
          'median_affected_regions':float(np.median(avail[:,0::4].sum(1))),
          'median_firstorder_regions':float(np.median(avail[:,1::4].sum(1))),
          'median_contrast_regions':float(np.median(avail[:,2::4].sum(1))),
          'median_texture_regions':float(np.median(avail[:,3::4].sum(1))),
          'patients_without_texture':int((avail[:,3::4].sum(1)==0).sum()),
          'patients_without_firstorder':int((avail[:,1::4].sum(1)==0).sum()),
          'columns':{'Fwhole':F.shape[1],'S':R,'L':l.shape[1]+avail.shape[1],'T':t.shape[1]+avail.shape[1]},
          'archive_sha256':sha(out),'assembler_sha256':sha(Path(__file__)),'outcomes_read':False}
    (HERE/f'features_{atlas}.json').write_text(json.dumps(meta,indent=2))
    print(json.dumps(meta))

def verify_aal3():
    """Assembled AAL3 must equal the historical V5/V6 arrays used by every AAL3 result."""
    a=np.load(HERE/'features_AAL3.npz');o5=np.load(INPUTS/'v5_features.npz')
    o6=np.load(V6/'output/features.npz')
    for k in ['ids','names','burden_names','F_unclip_p99','burden']:np.testing.assert_array_equal(a[k],o5[k])
    for k in ['L','T','W','L_names','T_names','W_names']:np.testing.assert_array_equal(a[k],o6[k])
    print('AAL3 assembled arrays identical to V5/V6 archives')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('atlases',nargs='+');args=ap.parse_args()
    for atlas in args.atlases:
        assemble(atlas)
        if atlas=='AAL3':verify_aal3()
