#!/usr/bin/env python3
"""Untimed default-equivalence checks and an ED writer/reader reproduction."""
import json,subprocess,re,hashlib,csv,io
from pathlib import Path
from benchmark import ROOT,BIN,ENV,DEFAULTS,canonical,validate,dump,sha
DEST=ROOT/'diagnostics'
def compare_outputs(first,second):
    def rows(name):
        text='\n'.join(line for line in (DEST/(name+'.stdout')).read_text().splitlines() if line and not line.startswith('#'))
        return sorted(csv.DictReader(io.StringIO(text),delimiter=';'),key=lambda row:(row['id1'],row['id2'],row['start1'],row['start2']))
    a,b=rows(first),rows(second)
    differences=[{'row':i,'fields':{key:{'computed':x[key],'loaded':y[key]} for key in x if x[key]!=y[key]}} for i,(x,y) in enumerate(zip(a,b)) if x!=y]
    structural=['id1','id2','start1','end1','start2','end2','bpList']
    return {'full_output_equal':a==b,'same_row_count':len(a)==len(b),'coordinates_and_bpList_equal':len(a)==len(b) and all(all(x[k]==y[k] for k in structural) for x,y in zip(a,b)),'differences':differences}
def execute(name,args):
    command=[str(BIN)]+args
    p=subprocess.run(command,env=ENV,capture_output=True,timeout=240)
    (DEST/(name+'.stdout')).write_bytes(p.stdout);(DEST/(name+'.stderr')).write_bytes(p.stderr)
    valid,rows=validate(p.stdout)
    row={'name':name,'command':command,'exit_code':p.returncode,'csv_valid':valid,'reported_rows':rows,'stdout_sha256':hashlib.sha256(p.stdout).hexdigest(),'canonical_sha256':hashlib.sha256(canonical(p.stdout)).hexdigest()}
    with (DEST/'commands.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
    p.check_returncode();assert valid,name
    return row

def main():
    DEST.mkdir(exist_ok=True)
    cases={x['name']:x for x in json.loads((ROOT/'cases.json').read_text())}
    checks=[]
    for name in ['bio_fhla','bio_phob','bio_ilve']:
        row=execute(name+'_implicit',[a for a in cases[name]['args'] if a not in DEFAULTS])
        reference=ROOT/'results'/(name+'.measure1.stdout')
        equal=row['canonical_sha256']==hashlib.sha256(canonical(reference.read_bytes())).hexdigest()
        checks.append({'case':name,'explicit_defaults_equal_implicit_defaults_bpList':equal})
        assert equal,name
    cache=DEST/'ed';cache.mkdir(exist_ok=True)
    base=cases['length_t1000_q100']['args']
    generated=execute('ed_write',base+['--out=tAcc:'+str(cache/'target.ed'),'--out=qAcc:'+str(cache/'query.ed'),'--verbose'])
    load_base=[a for a in base if a!='--acc=C']+['--tAcc=E','--tAccFile='+str(cache/'target.ed'),'--qAcc=E','--qAccFile='+str(cache/'query.ed')]
    loaded=execute('ed_reload',load_base+['--verbose'])
    rows={}
    for name in ['target','query']:
        p=cache/(name+'.ed');lines=p.read_text().splitlines();header=lines[1];maxlen=int(header.split()[-1]);data=[line.split() for line in lines[2:] if line.strip()]
        rows[name]={'sha256':sha(p),'header':header,'last_header_length':maxlen,'data_rows':len(data),'fields_per_row':sorted({len(r) for r in data}),'numeric_final_column_rows':sum(r[-1] not in ['NA','nan','NaN'] for r in data),'first_valid_last_column':next({'row':int(r[0]),'value':r[-1]} for r in data if r[-1] not in ['NA','nan','NaN']),'last_row_final_column':data[-1][-1]}
    matchedflags=['--tIntLenMax=149','--qIntLenMax=99']
    computed=execute('ed_matched_computed',[a for a in base if a!='--intLenMax=0']+matchedflags)
    reloaded=execute('ed_matched_reloaded',[a for a in load_base if a!='--intLenMax=0']+matchedflags)
    messages=[line for line in (DEST/'ed_reload.stdout').read_text().splitlines() if 'reducing maximal interaction length' in line]
    assert len(messages)==2 and any('to 149' in line for line in messages) and any('to 99' in line for line in messages)
    dump(DEST/'checks.json',{'implicit_default_checks':checks,'ed_files':rows,'reader_messages':messages,'unmatched_comparison':compare_outputs('ed_write','ed_reload'),'matched_149_99_comparison':compare_outputs('ed_matched_computed','ed_matched_reloaded'),'note':'ED loading is a separate diagnostic requested in the review, not part of default-mode performance timings; coordinates/base pairs and energy equality are checked separately.'})
    print(json.dumps(json.loads((DEST/'checks.json').read_text()),indent=2))
if __name__=='__main__':main()
