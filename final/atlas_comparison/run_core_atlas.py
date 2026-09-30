#!/usr/bin/env python3
"""V6 core comparison (model_v6.task_run, 24 methods, mRS/NIHSS) with L/T/W from another atlas.

The strong baseline B is deliberately left as in V6 (AAL3 burden + global lesion R), so that
B is identical for every atlas and each atlas only changes the added regional block.
task_run asserts per fold that B reproduces the V5 predictions exactly.
"""
import argparse, json, sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;FINAL=HERE.parent;ROOT=FINAL.parents[1];INPUTS=FINAL.parent/'inputs'
sys.path.insert(0,str(FINAL.parent))
import model_v6 as v6
mb=v6.mb

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--atlas',required=True);ap.add_argument('--task',required=True)
    ap.add_argument('--jobs',type=int,default=2);args=ap.parse_args()
    args.output=HERE/f'core_{args.atlas}';args.output.mkdir(exist_ok=True)
    p=json.loads((FINAL.parent/'config/protocol.json').read_text())
    data=np.load(HERE/f'features_{args.atlas}.npz');old=np.load(INPUTS/'v5_features.npz')
    np.testing.assert_array_equal(data['ids'],old['ids'])
    clinical=mb.load_clinical();mb.valid_ids(clinical.ID.tolist())
    mb.write_json(args.output/f'config_{args.task}.json',{'atlas':args.atlas,'task':args.task,'jobs':args.jobs,
        'baseline':'unchanged V6/V5 AAL3 baseline','regional_blocks':p['blocks'],
        'sources':{str(f.relative_to(ROOT)):mb.digest(f) for f in [Path(__file__),Path(v6.__file__),HERE/f'features_{args.atlas}.npz']}})
    v6.task_run(args.task,data,old,clinical,p,args)
    mb.write_json(args.output/f'COMPLETE_{args.task}.json',{'task':args.task})

if __name__=='__main__':main()
