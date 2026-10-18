from pathlib import Path
import subprocess, shutil
root=Path(__file__).resolve().parent
source=root.parent/'pr-worktrees/issue-270'
paths=subprocess.check_output(['git','-C',str(source),'ls-files','-c','-o','--exclude-standard','-z']).decode().split('\0')
for mode in ('release','debug'):
    changed=[]
    for path in filter(None,paths):
        src=source/path
        if not src.is_file(): continue
        dst=root/mode/path
        if not dst.exists() or src.read_bytes()!=dst.read_bytes():
            dst.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(src,dst)
            changed.append(path)
    print(mode,'synced',changed)
