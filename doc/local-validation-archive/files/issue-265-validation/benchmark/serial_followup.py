#!/usr/bin/env python3
import hashlib,json,os,pathlib,statistics,subprocess
root=pathlib.Path('/home/alex/CodexProjects/IntaRNA/issue-265-validation/benchmark/run')
original=json.loads((root/'results.json').read_text())
env=dict(os.environ,**original['environment'])
reference=original['warmups'][0]['hashes']
results=[]
for rep in range(4):
 for label in (['base','candidate'] if rep%2==0 else ['candidate','base']):
  template=next(r for r in original['results'] if r['label']==label and r['threads']==1)
  directory=root/('%s-t1-idle%d' %(label,rep+1));directory.mkdir()
  source_directory=root/template['name']
  command=[arg.replace(str(source_directory),str(directory)) for arg in template['command']]
  command[command.index('-c')+1]='8'
  with (directory/'stdout.txt').open('w') as stdout,(directory/'stderr.txt').open('w') as stderr:
   subprocess.run(command,env=env,stdout=stdout,stderr=stderr,check=True)
  wall,user,system,rss=(directory/'time.txt').read_text().split()
  lines=(directory/'predictions.csv').read_text().splitlines();norm=lines[0]+'\n'+'\n'.join(sorted(lines[1:]))+'\n'
  hashes={'predictions.csv':hashlib.sha256(norm.encode()).hexdigest()}
  for f in sorted(directory.glob('query*.ed')):hashes[f.name]=hashlib.sha256(f.read_bytes()).hexdigest()
  result={'label':label,'repetition':rep+1,'wall_s':float(wall),'user_s':float(user),'system_s':float(system),'peak_rss_kib':int(rss),'matches_base_serial':hashes==reference,'command':command}
  results.append(result);print(label,rep+1,wall,user,system,rss,result['matches_base_serial'],flush=True)
summary=[{'label':label,'threads':1,'median_wall_s':statistics.median(r['wall_s'] for r in results if r['label']==label),'min_wall_s':min(r['wall_s'] for r in results if r['label']==label),'max_wall_s':max(r['wall_s'] for r in results if r['label']==label),'median_peak_rss_kib':statistics.median(r['peak_rss_kib'] for r in results if r['label']==label)} for label in ['base','candidate']]
(root/'serial-followup.json').write_text(json.dumps({'cpu_affinity':'8','note':'All four repetitions occurred after build and test processes completed. Same immutable binary snapshots, inputs/options and output checks as original benchmark.','results':results,'summary':summary},indent=2)+'\n')
print(json.dumps(summary,indent=2),flush=True)
if not all(r['matches_base_serial'] for r in results):raise SystemExit('Output mismatch')
