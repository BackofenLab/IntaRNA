#!/usr/bin/env python3
"""Independently audit saved evidence; never runs IntaRNA or rewrites measurements."""
import csv,hashlib,io,json,math,statistics
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
RESULT=ROOT/'results'
errors=[];warnings=[];findings=[]
def expect(condition,message):
    if not condition:errors.append(message)
def digest(data):return hashlib.sha256(data).hexdigest()
def load(path):return json.loads(path.read_text())
def records(path):return [json.loads(line) for line in path.read_text().splitlines()]
def normalized(data):return b'\n'.join(sorted(x for x in data.splitlines() if x and not x.startswith(b'#') and not x.startswith(b'Creating experiment directory ')))+b'\n'
def equal(a,b):return math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-10) if isinstance(a,(int,float)) and isinstance(b,(int,float)) else a==b

def main():
    plan=load(ROOT/'study-plan.json');cases=load(ROOT/'cases.json');env=load(ROOT/'environment.json');runs=records(RESULT/'runs.jsonl');summary=load(RESULT/'summary.json')
    definitions={x['name']:x for x in cases};lookup={x['case']:x for x in summary}
    expect(len(definitions)==plan['cases']==18,'Case coverage')
    expect(env['revision']==plan['revision'],'Revision provenance')
    binary=ROOT/'source/src/bin/IntaRNA'
    if binary.exists():expect(digest(binary.read_bytes())==env['binary_sha256'],'Executable SHA-256')
    else:warnings.append('Compiled executable omitted; recorded per-run binary hashes are checked.')
    expect(env['environment_controlled']==plan['controlled_environment'],'Controlled environment')
    expect(set(plan['scientific_defaults'])=={'--model=X','--mode=H','--acc=C','--accW=150','--accL=100','--windowWidth=0','--seedBP=7','--seedMaxUP=0','--intLenMax=0'},'Requested default-mode scope')
    allowed={'target','query','model','mode','acc','accW','accL','windowWidth','seedBP','seedMaxUP','intLenMax','threads','outMode','outCsvCols'}
    for case in cases:
        flags=[a.split('=',1)[0][2:] for a in case['args']]
        expect(len(flags)==len(set(flags)) and set(flags)==allowed,case['name']+' parameter set')
        expect(set(plan['scientific_defaults']).issubset(case['args']),case['name']+' scientific defaults')
    checked_inputs=0
    for name,entry in load(ROOT/'inputs/manifest.json').items():
        expect(digest((ROOT/'inputs'/name).read_bytes())==entry['sha256'],'Input hash '+name);checked_inputs+=1
    expected={(c['name'],phase,rep) for c in cases for phase,rep in [('warmup',0)]+[('measure',n) for n in range(1,6)]}
    observed=Counter((r['case'],r['phase'],r['repetition']) for r in runs)
    expect(set(observed)==expected and all(v==1 for v in observed.values()),'Exactly one warmup plus five measured records per case')
    end=None;out_hash={};structure_hashes={}
    for r in runs:
        name=f"{r['case']}.{r['phase']}{r['repetition']}";data=(RESULT/(name+'.stdout')).read_bytes();key=(r['case'],r['phase'],r['repetition'])
        expect(r['exit_code']==0 and not r['timeout'],name+' exit status')
        expect(r['binary_sha256']==env['binary_sha256'],name+' binary provenance')
        expect(r['command'][6:]==definitions[r['case']]['args'],name+' command definition')
        expect(digest(data)==r['output_sha256'],name+' stdout hash')
        norm=digest(normalized(data));out_hash[key]=norm
        expect(norm==r['canonical_sha256'],name+' canonical hash')
        cells=list(csv.reader(io.StringIO('\n'.join(x for x in data.decode().splitlines() if x and not x.startswith('#'))),delimiter=';'))
        expect(bool(cells) and cells[0]==['id1','id2','start1','end1','start2','end2','E','bpList'] and all(len(c)==8 for c in cells),name+' full-structure CSV')
        expect(r['csv_valid'] and len(cells)-1==r['output_rows'],name+' derived metadata')
        structure_hashes[key]=digest(json.dumps(sorted(tuple(row[:6]+row[7:]) for row in cells[1:])).encode())
        expect(not (RESULT/(name+'.stderr')).read_bytes(),name+' empty stderr')
        vals=list(map(float,(RESULT/(name+'.time')).read_text().strip().split('\t')))
        for field,value in zip(['wall_s','user_s','system_s','max_rss_kib','major_faults','minor_faults','voluntary_context_switches','involuntary_context_switches','fs_inputs','fs_outputs','time_exit_code'],vals):expect(equal(value,r[field]),name+' GNU time '+field)
        expect(r['elapsed_s']>0 and -0.02<=r['elapsed_s']-r['wall_s']<0.1,name+' monotonic/GNU wall agreement')
        start=datetime.fromisoformat(r['started_utc']);stop=datetime.fromisoformat(r['ended_utc'])
        expect(stop>=start and (end is None or start>=end),name+' sequential timing');end=stop
    for case in cases:
        group=[r for r in runs if r['case']==case['name']];measured=[r for r in group if r['phase']=='measure'];times=[r['elapsed_s'] for r in measured];saved=lookup[case['name']]
        stable_all=len({out_hash[(r['case'],r['phase'],r['repetition'])] for r in group})==1
        if not stable_all:findings.append({'kind':'within_case_output_instability','case':case['name'],'records':[f"{r['case']}.{r['phase']}{r['repetition']}.stdout" for r in group]})
        expected_values={'n':5,'median_s':statistics.median(times),'min_s':min(times),'max_s':max(times),'range_over_median_pct':100*(max(times)-min(times))/statistics.median(times),'median_user_s':statistics.median(r['user_s'] for r in measured),'median_system_s':statistics.median(r['system_s'] for r in measured),'median_rss_mib':statistics.median(r['max_rss_kib']/1024 for r in measured),'peak_rss_mib':max(r['max_rss_kib']/1024 for r in measured),'stable_output':len({r['canonical_sha256'] for r in measured})==1,'stable_bpList':len({structure_hashes[(r['case'],r['phase'],r['repetition'])] for r in measured})==1,'reported_rows':measured[0]['output_rows'],'case':case['name'],'category':case['category'],'threads':case['threads'],'predictions':case['predictions']}
        expect(set(saved)==set(expected_values),case['name']+' summary fields')
        for k,v in expected_values.items():expect(equal(saved[k],v),case['name']+' summary '+k)
    csv_rows=list(csv.DictReader((RESULT/'summary.csv').open()))
    expect(len(csv_rows)==len(summary),'CSV summary coverage')
    for row in csv_rows:
        for k,v in lookup[row['case']].items():expect(row[k]==str(v),row['case']+' CSV '+k)
    thread_equality={};thread_structure_equality={}
    for rna in ['gcvb','chix']:
        keys=[key for key in out_hash if key[0] in {f'batch_{rna}_{t}t' for t in [1,6,12]}]
        thread_equality[rna]=len({out_hash[key] for key in keys})==1
        thread_structure_equality[rna]=len({structure_hashes[key] for key in keys})==1
        if not thread_equality[rna]:findings.append({'kind':'thread_count_output_mismatch','query':rna,'coordinates_and_bpList_equal':thread_structure_equality[rna]})
    diag=load(ROOT/'diagnostics/checks.json')
    expect(all(r['explicit_defaults_equal_implicit_defaults_bpList'] for r in diag['implicit_default_checks']),'Explicit vs implicit defaults')
    for r in records(ROOT/'diagnostics/commands.jsonl'):
        data=(ROOT/'diagnostics'/(r['name']+'.stdout')).read_bytes()
        expect(digest(data)==r['stdout_sha256'] and digest(normalized(data))==r['canonical_sha256'],r['name']+' diagnostic hashes')
    for name,entry in diag['ed_files'].items():expect(digest((ROOT/'diagnostics/ed'/(name+'.ed')).read_bytes())==entry['sha256'],name+' ED hash')
    profile_checks=load(ROOT/'profiles/checks.json')
    for name,entry in load(ROOT/'profiles/input-manifest.json').items():
        expect(digest((ROOT/'inputs'/name).read_bytes())==entry['sha256'],'Profile input hash '+name)
    expect(Counter(r['kind'] for r in profile_checks)=={'cpu':4,'heap':2,'phases':3},'Profile coverage')
    expect(all(r['output_equal_bpList'] for r in profile_checks),'All profiler/phase outputs retain bpList')
    for r in records(ROOT/'profiles/commands.jsonl'):
        expect(r['exit_code']==0 and datetime.fromisoformat(r['started_utc'])>=end,r['name']+' profile after clean timing')
    # Independently re-check all nine instrumentation outputs against preserved raw references.
    for r in profile_checks:
        name=r['profile'];kind=r['kind'];suffix={'cpu':'cpu','heap':'heap','phases':'phases'}[kind]
        if name=='chix384':reference=(ROOT/'profiles/chix384_reference.stdout').read_bytes()
        else:
            case={'gcvb':'batch_gcvb_1t','gcvb_repeat':'batch_gcvb_1t'}.get(name,name)
            reference=(RESULT/(case+'.measure1.stdout')).read_bytes()
        actual=(ROOT/'profiles'/(name+'_'+suffix+'.stdout')).read_bytes()
        expect(normalized(actual)==normalized(reference),name+' independent '+kind+' output check')
    parallel=load(ROOT/'diagnostics/parallel/summary.json')
    repeat_rows=records(ROOT/'diagnostics/parallel/checks.jsonl')
    expect(Counter(r['threads'] for r in repeat_rows)=={1:3,6:5,12:20},'Follow-up repeat coverage')
    repeat_reference=normalized((RESULT/'batch_chix_1t.measure1.stdout').read_bytes())
    for r in repeat_rows:
        data=(ROOT/'diagnostics/parallel'/(r['name']+'.stdout')).read_bytes()
        expect(digest(data)==r['stdout_sha256'] and r['binary_sha256']==env['binary_sha256'],r['name']+' repeat provenance')
        expect(r['output_equal']==(normalized(data)==repeat_reference and r['exit_code']==0),r['name']+' repeat comparison')
    for threads in [1,6,12]:
        expected_counts={'runs':sum(r['threads']==threads for r in repeat_rows),'different':sum(r['threads']==threads and not r['output_equal'] for r in repeat_rows)}
        expect(parallel['counts'][str(threads)]==expected_counts,'Repeat counts '+str(threads))
    for r in load(ROOT/'diagnostics/numeric-commands.json'):
        expect(r['exit_code']==0 and digest((ROOT/'diagnostics'/(r['name']+'.stdout')).read_bytes())==r['stdout_sha256'],r['name']+' energy-component evidence')
    if parallel['observed_differences']:findings.append({'kind':'parallel_output_instability_reproduced','counts':parallel['counts'],'observed_differences':[{k:r[k] for k in ['name','threads','differences']} for r in parallel['observed_differences']]})
    if not diag['matched_149_99_comparison']['full_output_equal']:findings.append({'kind':'matched_ED_output_mismatch','comparison':diag['matched_149_99_comparison']})
    report={'audited_utc':datetime.now(timezone.utc).isoformat(),'status':('pass_with_findings' if findings else 'pass') if not errors else 'fail','integrity_status':'pass' if not errors else 'fail','scientific_output_status':'differences_observed' if findings else 'all_checked_equal','scientific_findings':findings,'errors':errors,'warnings':warnings,'cases':len(cases),'measured_runs':sum(r['phase']=='measure' for r in runs),'warmups':sum(r['phase']=='warmup' for r in runs),'input_hashes_checked':checked_inputs,'profile_checks':len(profile_checks),'thread_equality_full_output':thread_equality,'thread_equality_coordinates_and_bpList':thread_structure_equality,'default_equivalence_cases':len(diag['implicit_default_checks']),'additional_repeat_checks':len(repeat_rows)}
    (RESULT/'audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
    raise SystemExit(bool(errors))
if __name__=='__main__':main()
