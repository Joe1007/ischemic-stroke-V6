#!/usr/bin/env python3
import gzip
import json
import hashlib
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;INPUTS=HERE/'inputs'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def assemble():
    directory=HERE/'output';manifest=json.loads((directory/'COMPLETE.json').read_text())
    assert not manifest['pilot']
    oldpath=INPUTS/'v5_features.npz'
    oldmanifest=json.loads((oldpath.parent/'features_manifest.json').read_text())
    assert sha(oldpath)==oldmanifest['archive_sha256']
    old=np.load(oldpath,allow_pickle=False)
    assert manifest['ids']==old['ids'].tolist()
    regions=[int(n.rsplit('_',1)[1]) for n in old['burden_names']]
    feature_names=old['global_names'].tolist()
    assert old['names'].tolist()==[f'{f}_{r:03d}' for r in regions for f in feature_names]
    protocol=json.loads((HERE/'config/protocol.json').read_text())
    ph=hashlib.sha256(json.dumps(protocol,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    stats=protocol['firstorder']
    ln=[f'AAL3_{r:03d}|{kind}|{f}' for r in regions for kind in ['lesion','contrast'] for f in stats]
    tn=[f'AAL3_{r:03d}|lesion|{f}' for r in regions for f in feature_names]
    flags=['affected','firstorder_available','contrast_available','texture_available']
    fn=[f'AAL3_{r:03d}|{f}' for r in regions for f in flags]
    l=np.zeros((len(old['ids']),len(ln)),dtype='float32');t=np.zeros((len(l),len(tn)),dtype='float32')
    avail=np.zeros((len(l),len(fn)),dtype='float32');qc=[];hashes={}
    for i,sid in enumerate(old['ids']):
        path=directory/'checkpoints'/f'{sid}.json.gz';hashes[str(sid)]=sha(path)
        with gzip.open(path,'rt') as h:d=json.load(h)
        assert d['ID']==sid and d['status']=='complete'
        assert d['provenance']['protocol']==ph
        assert d['provenance']['code']==sha(HERE/'extract_v6.py')
        assert sha(INPUTS/'v5_checkpoints'/f'{sid}.json.gz')==d['provenance']['source']
        assert [r['region'] for r in d['rows']]==regions
        for j,row in enumerate(d['rows']):
            a=row['lesion_voxels']>0;fo=row['firstorder'];co=row['contrast'];tex=row['texture']
            avail[i,j*4:(j+1)*4]=[a,fo is not None,co is not None,tex is not None]
            if a:
                l[i,j*16:j*16+8]=fo if fo is not None else np.nan
                l[i,j*16+8:(j+1)*16]=co if co is not None else np.nan
                t[i,j*93:(j+1)*93]=[tex[f] for f in feature_names] if tex is not None else np.nan
        qc.append({'ID':str(sid),**{f:int(avail[i,k::4].sum()) for k,f in enumerate(flags)}})
    w=np.where(np.repeat(old['burden']>0,93,axis=1),old['F_unclip_p99']*np.repeat(old['burden'],93,axis=1),0)
    np.savez_compressed(directory/'features.npz',ids=old['ids'],L=np.hstack([l,avail]),T=np.hstack([t,avail]),
        W=np.hstack([w,avail]),L_names=np.array(ln+fn),T_names=np.array(tn+fn),
        W_names=np.r_[np.char.add('occupancy_weighted|',old['names']),fn])
    metadata={'archive_sha256':sha(directory/'features.npz'),'assembler_sha256':sha(Path(__file__)),
        'v5_archive_sha256':sha(oldpath),'protocol_sha256':sha(HERE/'config/protocol.json'),
        'checkpoint_sha256':hashes,'qc':qc,'outcomes_read':False}
    (directory/'features_manifest.json').write_text(json.dumps(metadata,indent=2))
    print('Assembled',len(l),'subjects; median availability',np.median([[q[f] for f in flags] for q in qc],axis=0))
if __name__=='__main__':assemble()
