from pathlib import Path
import hashlib,json,subprocess
r=Path('/home/alex/CodexProjects/IntaRNA');src=r/'pr-worktrees/issue-274';v=r/'issue-274-validation'
base='21b9d8fc33aadd3c44d788c4da10a91cf6694ce1'
names=subprocess.check_output(['git','-C',str(src),'ls-files','src','tests'],text=True).splitlines()
names=[name for name in names if Path(name).suffix in ('.h','.cpp','.cc','.icc')]
result={'baseline_revision':base,'source_verification':{},'modified_source_sha256':{},'build_flags':{},'binaries':{}}
for name,tree in [('base',r/'issue-270-anchored-validation/release'),('release',v/'release'),('debug',v/'debug')]:
 mismatch=[]
 for f in names:
  expected=subprocess.check_output(['git','-C',str(src),'show',f'{base}:{f}']) if name=='base' else (src/f).read_bytes()
  if (tree/f).read_bytes()!=expected:mismatch.append(f)
 result['source_verification'][name]={'checked_files':len(names),'mismatches':mismatch}
 if mismatch:raise RuntimeError((name,mismatch))
 binary=tree/'src/bin/IntaRNA';result['binaries'][name]={'sha256':hashlib.sha256(binary.read_bytes()).hexdigest()}
 flags={}
 for file in ('Makefile','src/IntaRNA/Makefile','src/bin/Makefile'):
  flags[file]=[line for line in (tree/file).read_text().splitlines() if line.startswith(('CXX =','CXXFLAGS =','CPPFLAGS =','AM_CXXFLAGS =','LDFLAGS ='))]
 result['build_flags'][name]=flags
for name in subprocess.check_output(['git','-C',str(src),'diff','--name-only'],text=True).splitlines():
 if name.startswith(('src/','tests/')):result['modified_source_sha256'][name]=hashlib.sha256((src/name).read_bytes()).hexdigest()
result['compiler']=subprocess.check_output([str(r/'issue-246-validation/toolchain-gcc14/bin/x86_64-conda-linux-gnu-g++'),'--version'],text=True).splitlines()[0]
(v/'source-verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(result['source_verification'])
