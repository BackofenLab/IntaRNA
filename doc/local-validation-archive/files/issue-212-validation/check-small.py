import csv, io, subprocess, json
b='/home/alex/CodexProjects/IntaRNA/pr-worktrees/issue-212/src/bin/IntaRNA'
def run(args):
 p=subprocess.run([b,'--threads=1','--outMode=C','--outCsvCols=start1,end1,start2,end2,E','-n','100','--outDeltaE=100']+args,capture_output=True,text=True)
 return p.returncode,p.stdout,p.stderr
for t,q in [('CCCCAAUU','CCGGGGUU'),('CCCAAU','CGGGUU'),('CCCAAUU','CGGGUU'),('CCCAAU','GGGUU'),('CCAAU','GGUU')]:
 for mode in ('H','M'):
  for ov in ('N','T','Q'):
   code,out,err=run(['-t',t,'-q',q,'--seedBP=2','--energy=B','--acc=N','--mode='+mode,'--outOverlap='+ov])
   print(t,q,mode,ov,out.replace('\n',' / '),err)
