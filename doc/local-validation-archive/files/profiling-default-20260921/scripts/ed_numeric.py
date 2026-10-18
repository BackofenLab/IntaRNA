#!/usr/bin/env python3
"""Separate untimed follow-up locating the observed ED-round-trip energy delta."""
import csv,hashlib,io,json,subprocess
from benchmark import ROOT,BIN,ENV,dump
DEST=ROOT/'diagnostics'
commands=[json.loads(x) for x in (DEST/'commands.jsonl').read_text().splitlines()]
columns='id1,id2,start1,end1,start2,end2,E,bpList,ED1,ED2,E_init,E_loops,E_dangleL,E_dangleR,E_endL,E_endR,E_hybrid'
result={};records=[]
for name in ['ed_matched_computed','ed_matched_reloaded']:
    original=next(r['command'] for r in commands if r['name']==name)
    command=[str(BIN)]+[('--outCsvCols='+columns) if a.startswith('--outCsvCols=') else a for a in original[1:]]
    p=subprocess.run(command,env=ENV,capture_output=True,timeout=60);p.check_returncode()
    (DEST/(name+'_components.stdout')).write_bytes(p.stdout);(DEST/(name+'_components.stderr')).write_bytes(p.stderr)
    result[name]=list(csv.DictReader(io.StringIO(p.stdout.decode()),delimiter=';'))
    assert len(result[name])==1
    records.append({'name':name+'_components','command':command,'exit_code':p.returncode,'stdout_sha256':hashlib.sha256(p.stdout).hexdigest()})
dump(DEST/'energy-components.json',result);dump(DEST/'numeric-commands.json',records)
print(json.dumps(result,indent=2))
