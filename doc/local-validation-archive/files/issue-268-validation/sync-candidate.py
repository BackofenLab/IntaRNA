from pathlib import Path
import subprocess, shutil
root=Path('/home/alex/CodexProjects/IntaRNA')
source=root/'pr-worktrees/issue-268'
paths=subprocess.check_output(['git','-C',str(source),'ls-files','-c','-o','--exclude-standard','-z']).decode().split('\0')
for mode in ('release','debug'):
    changed=[]
    for path in filter(None,paths):
        src=source/path
        if not src.is_file(): continue
        dst=root/'issue-268-validation'/mode/path
        if not dst.exists() or src.read_bytes()!=dst.read_bytes():
            dst.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(src,dst)
            changed.append(path)
    print(mode,'synced',changed)
