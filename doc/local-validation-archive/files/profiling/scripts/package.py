#!/usr/bin/env python3
"""Bundle report and portable evidence, leaving large build/profile binaries local."""
import hashlib,json,shutil,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
for name in ['runApiTests.log','runIntaRNA.log','test-suite.log']:
    shutil.copyfile(ROOT/'source/tests'/name,ROOT/'logs'/('tests-'+name))
files=[]
for name in ['README.md','REPORT.md','environment.json','cases.json']:
    files.append(ROOT/name)
for dirname in ['scripts','inputs','results','council','logs']:
    for p in (ROOT/dirname).rglob('*'):
        if p.is_file() and '__pycache__' not in p.parts and p.name!='stream_cleanup_probe':files.append(p)
files += [p for p in (ROOT/'profiles').iterdir() if p.is_file() and not p.name.startswith('probe')]
files += [p for p in (ROOT/'source/src').rglob('*') if p.is_file() and p.suffix in {'.cpp','.cc','.h','.icc'}]
manifest={str(p.relative_to(ROOT)):dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(set(files))}
dest=ROOT/'intarna-profile-evidence.zip'
with zipfile.ZipFile(dest,'w',compression=zipfile.ZIP_DEFLATED) as z:
    for p in sorted(set(files)):z.write(p,'profiling/'+str(p.relative_to(ROOT)))
    z.writestr('profiling/bundle-manifest.json',json.dumps(manifest,indent=2)+'\n')
(ROOT/'bundle-info.json').write_text(json.dumps(dict(file=dest.name,bytes=dest.stat().st_size,
    sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),files=len(manifest),
    note='Large source build artifacts, compiled probe, gprofng experiment directories and old-binary capability probes remain local and are excluded.'),indent=2)+'\n')
print(dest, dest.stat().st_size, 'bytes')
