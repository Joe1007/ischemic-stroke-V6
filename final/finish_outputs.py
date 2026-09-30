#!/usr/bin/env python3
"""Wait for successful training, then validate and package the requested outputs."""
import fcntl,json,os,subprocess,sys,time,shutil
from pathlib import Path
from datetime import datetime
HERE=Path(__file__).resolve().parent
OUT=HERE/'output';OUT.mkdir(exist_ok=True)
lock=(OUT/'.report.lock').open('w')
fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
training=HERE/'tables_20260917'
def status(state,count):
    content=f'# 輸出狀態\n\n更新：{datetime.now().isoformat(timespec="seconds")}\n\n{state}\n\n已完成模型組合：{count}/2160。\n\n[已完成的 V6 核心報告](core/RESULTS.md)\n\n完整圖表生成與稽核成功後，才會產生 RESULTS.md 與 COMPLETE.json。\n'
    p=OUT/'STATUS.md.tmp';p.write_text(content);p.replace(OUT/'STATUS.md')
subprocess.run([sys.executable,str(HERE/'audit_core.py')],check=True)  # writes into output/core
last=None
while not (training/'COMPLETE.json').exists():
    count=sum(p.name!='split.json' for p in training.glob('*_[0-9][0-9][0-9]/*.json'))
    status('四任務訓練進行中；完整表圖尚未完成。',count)
    if count!=last:print(datetime.now().isoformat(),count,'/2160',flush=True);last=count
    time.sleep(30)
status('訓練完成，正在稽核、彙整指標及繪製圖表。',2160)
try:
    subprocess.run([sys.executable,str(HERE/'report_tables.py'),'--input',str(training),'--output',str(OUT)],check=True)
except BaseException:
    status('報告生成失敗或中斷，請檢查 finish_outputs.log；尚未宣稱完成。',2160)
    raise
status('已完成；請閱讀 [RESULTS.md](RESULTS.md)。所有產物與 hash 清單保存在本目錄。',2160)
print('Completed output:',OUT,flush=True)
