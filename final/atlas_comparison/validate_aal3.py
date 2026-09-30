#!/usr/bin/env python3
"""The generic extractor run with AAL3 must reproduce the V5 and V6 AAL3 checkpoints exactly."""
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;V6=HERE.parents[1];ROOT=V6.parent;INPUTS=V6/'inputs'
sys.path.insert(0,str(HERE))
import extract_atlas as ea
eb=ea.eb

def main(directory):
    labels=eb.labels_from_file(ROOT/'data/MNI/AAL3.txt')
    assert ea.paired_labels(labels)==ea.v6.paired_labels(labels),'Generic pairing changed AAL3 pairs'
    ids=json.loads((directory/'manifest.json').read_text())['ids'];checked=0
    for sid in ids:
        new=eb.read(directory/'checkpoints'/f'{sid}.json.gz')
        old5=eb.read(INPUTS/'v5_checkpoints'/f'{sid}.json.gz')
        old6=eb.read(V6/'output/checkpoints'/f'{sid}.json.gz')
        r5={r['region']:r for r in old5['rows']};r6={r['region']:r for r in old6['rows']}
        assert [r['region'] for r in new['rows']]==sorted(r5)==sorted(r6)
        for row in new['rows']:
            a=r5[row['region']];b=r6[row['region']]
            assert row['roi_voxels']==a['roi_voxels'] and row['lesion_fraction']==a['lesion_fraction']
            assert (row['whole'] is None)==('unclip_p99' not in a['features'])
            if row['whole'] is not None:assert row['whole']==a['features']['unclip_p99']
            for k in ['lesion_voxels','valid_voxels','firstorder','contrast','texture','reference_voxels','counterpart']:
                assert row[k]==b[k],(sid,row['region'],k)
            checked+=1
    print(json.dumps({'subjects':len(ids),'region_rows_checked':checked,'exact_match_v5_unclip_p99_and_v6':True}))

if __name__=='__main__':main(Path(sys.argv[1]))
