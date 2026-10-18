#!/usr/bin/env python3
"""Create a checksummed evidence archive; exclude binaries and raw experiments."""
import hashlib,json,shutil,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
for name in ['runApiTests.log','runIntaRNA.log','test-suite.log']:
    shutil.copyfile(ROOT/'source/tests'/name,ROOT/'logs'/('tests-'+name))
selected=[]
for name in ['README.md','REPORT.md','ED-READER-EVIDENCE.md','OUTPUT-STABILITY.md','environment.json','cases.json','study-plan.json','review-requirements.json']:
    selected.append(ROOT/name)
for folder in ['scripts','inputs','results','diagnostics','logs']:
    selected.extend(p for p in (ROOT/folder).rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='ed_conversion_probe')
selected.extend(p for p in (ROOT/'profiles').iterdir() if p.is_file())
selected=sorted(set(selected));manifest={str(p.relative_to(ROOT)):{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in selected}
archive=ROOT/'default-mode-evidence.zip'
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED) as z:
    for p in selected:z.write(p,'default-mode/'+str(p.relative_to(ROOT)))
    z.writestr('default-mode/manifest.json',json.dumps(manifest,indent=2)+'\n')
meta={'file':archive.name,'bytes':archive.stat().st_size,'sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),'files':len(selected),'source_revision':'0a4568ad6e3da52f935221ff0da98a4a171e9833','exclusions':'Compiled source/build artifacts and raw gprofng experiment directories; sources remain pinned in upstream Git.'}
(ROOT/'bundle-info.json').write_text(json.dumps(meta,indent=2)+'\n');(ROOT/'SHA256SUMS').write_text(meta['sha256']+'  '+archive.name+'\n')
print(json.dumps(meta,indent=2))
