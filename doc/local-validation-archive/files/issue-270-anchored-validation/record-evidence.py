from pathlib import Path
import hashlib,json,shutil,subprocess
root=Path(__file__).resolve().parent
repo=root.parent/'pr-worktrees/issue-270-anchored'
dest=repo/'doc/benchmarks/seed-extension'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
for label,folder in [('native','native-runtime'),('pairs','native-pairs'),('exact','kernel-exact'),('heuristic','kernel-heuristic')]:
    for name in ('results.json','samples.tsv'):
        shutil.copy2(root/folder/name,dest/f'anchored-{label}-{name}')
source_names=subprocess.check_output(['git','diff','--name-only','--','src','tests'],cwd=repo,text=True).splitlines()
flags={}
for variant,folder in [('base',root.parent/'issue-268-validation/release'),('candidate',root/'release')]:
    flags[variant]={}
    for line in (folder/'src/IntaRNA/Makefile').read_text().splitlines():
        if any(line.startswith(key+' =') for key in ('AM_CXXFLAGS','CPPFLAGS','CXX','CXXFLAGS','LDFLAGS')):
            key,value=line.split(' =',1);flags[variant][key]=value.strip()
    flags[variant]['configure']=subprocess.check_output([str(folder/'config.status'),'--config'],cwd=folder,text=True).strip()
metadata={
    'baseline':'3866c62d10846fab1adc38b937d0c9b96cc06e51',
    'date':'2026-10-16',
    'candidate_source_sha256':{name:digest(repo/name) for name in source_names},
    'native_build_flags':flags,
    'source_verification':read(root/'source-verification.json'),
    'kernel_build':read(root/'kernel/manifest.json'),
    'kernel_binary_sha256':{v:digest(root/'kernel'/v/'benchmark') for v in ('base','candidate')},
    'kernel_dump_sha256':{str(p.relative_to(root)):digest(p) for mode in ('exact','heuristic') for p in sorted((root/f'kernel-{mode}-verify').glob('*.txt')) if p.name!='last-rss.txt'},
    'native_manifests':{label:read(root/folder/'manifest.json') for label,folder in [('ordinary','native-runtime'),('pair_output','native-pairs')]},
    'native_verification':{'commands':read(root/'native-comparison/commands.json'),'checks':read(root/'native-comparison/results.json')},
    'validation':{'release':'make tests -j2: 7/7 drivers; 115 API cases; 457827 assertions','debug':'make tests -j2: 7/7 drivers; 115 API cases; 457828 assertions','kernel_exact':read(root/'kernel-exact-verify/verification.json'),'kernel_heuristic':read(root/'kernel-heuristic-verify/verification.json')},
    'measurement_conditions':{'repetitions':5,'threads':1,'ordering':'alternating baseline/candidate fresh processes','concurrent_project_build_test_install':False,'other_workstation_activity_cpu_affinity_frequency':'uncontrolled','platform_validation':'Linux GCC14.4; macOS/Apple Clang not run'},
}
(dest/'anchored-source-hashes.json').write_text(json.dumps(metadata,indent=2)+'\n')
print('Recorded benchmark summaries, samples and provenance')
