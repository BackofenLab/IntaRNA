from pathlib import Path
import difflib
import hashlib
import json
import argparse
parser=argparse.ArgumentParser(description="Instrument a disposable configured IntaRNA build to attribute seeded partition work.")
parser.add_argument('--source',type=Path,required=True,help='unmodified baseline source checkout')
parser.add_argument('--build',type=Path,required=True,help='disposable configured build copy; must contain source files')
parser.add_argument('--output',type=Path,required=True,help='destination for patch, source snapshots, and hashes')
args=parser.parse_args()
base=args.source.resolve()
build=args.build.resolve()
output=args.output.resolve()
if base==build:
    parser.error('--build must differ from --source; never instrument production source')
if not (build/'config.status').is_file():
    parser.error('--build must be a configured disposable build copy')
output.mkdir(parents=True,exist_ok=True)
def replace_once(text, old, new, count=1):
    if text.count(old)!=count:
        raise RuntimeError(f'expected {count} replacement target(s), found {text.count(old)}: {old!r}')
    return text.replace(old,new,count)
files=['src/IntaRNA/PredictorMfeEns2dSeedExtension.cpp','src/IntaRNA/SeededPartitionFunction.cpp']
original={name:(base/name).read_text() for name in files}
text=original[files[0]]
text='#include <chrono>\n#include <iostream>\n'+text
text=replace_once(text,'\tprepareStackSeedExtensions(n,m,seeds);','\tstd::cerr << "#profile_seeds";\n\tsize_t profileSeedCount=0;\n\tfor(size_t i=0;i<n;++i) for(size_t j=0;j<m;++j) {\n\t\tconst size_t length=seeds.seedLength(i,j);\n\t\tif(length) { ++profileSeedCount; std::cerr << " [" << i << "," << j << "," << length << "]"; }\n\t}\n\tstd::cerr << " admitted_count=" << profileSeedCount << \'\\n\';\n\tconst auto profileBegin=std::chrono::steady_clock::now();\n\tprepareStackSeedExtensions(n,m,seeds);\n\tconst auto profilePrepared=std::chrono::steady_clock::now();',1)
text=replace_once(text,'\tconst auto result=K::compute(domain,seeds,weights,pairProbabilities!=nullptr);','\tconst auto result=K::compute(domain,seeds,weights,pairProbabilities!=nullptr);\n\tconst auto profileComputed=std::chrono::steady_clock::now();',1)
needle='\toutput.setExactPartition(true);\n\treportOptima();'
replacement=needle+'''\n\tconst auto profileReported=std::chrono::steady_clock::now();
\tstd::cerr << "#profile_predictor prepare_s=" << std::chrono::duration<double>(profilePrepared-profileBegin).count()
\t\t<< " compute_s=" << std::chrono::duration<double>(profileComputed-profilePrepared).count()
\t\t<< " report_s=" << std::chrono::duration<double>(profileReported-profileComputed).count() << '\\n';'''
assert text.count(needle)==1
text=replace_once(text,needle,replacement,1)
modified={files[0]:text}
text=original[files[1]]
text='#include <chrono>\n#include <iostream>\n#include <cstdint>\n'+text
needle='\tseeds.forEachStart([&](size_t p1,size_t p2) {'
assert text.count(needle)==1
text=replace_once(text,needle,'''\tuint64_t profileSources=0,profileLeftCells=0,profileBoundedLeftCells=0,profileGroups=0,profileAllGroups=0,profileRightCells=0,profileCroppedCells=0,profileAnchors=0;
\tdouble profileLeftSeconds=0,profileRightSeconds=0;
'''+needle,1)
needle='\t\tconst size_t n=std::min(d.span1,d.n-p1),m=std::min(d.span2,d.m-p2);'
assert text.count(needle)==1
text=replace_once(text,needle,needle+'''\n\t\t++profileSources; profileLeftCells+=n*m;
\t\tconst auto profileLeftBegin=std::chrono::steady_clock::now();''',1)
needle='\t\tfor(size_t i=0;i<n;++i) for(size_t j=0;j<m;++j) {\n\t\t\tvalid(i,j)=w.valid({p1+i,p2+j});'
text=replace_once(text,needle,'\t\tsize_t profilePrefixN=0,profilePrefixM=0;\n'+needle,1)
needle='\t\t\tending(i,j)=seeds.seedEndingLength(p1+i,p2+j);'
text=replace_once(text,needle,needle+'''\n\t\t\tconst size_t profileLength=ending(i,j);
\t\t\tif(profileLength && profileLength<=i+1 && profileLength<=j+1) {
\t\t\t\tprofilePrefixN=std::max(profilePrefixN,i+1); profilePrefixM=std::max(profilePrefixM,j+1);
\t\t\t}''',1)
needle='\t\tauto accepts=[&](size_t i,size_t j,size_t len)'
text=replace_once(text,needle,'\t\tprofileBoundedLeftCells+=profilePrefixN*profilePrefixM;\n'+needle,1)
needle='\t\t// Retained right extensions share first-hit injections by linearity:'
text=replace_once(text,needle,'''\t\tprofileLeftSeconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-profileLeftBegin).count();
'''+needle,1)
needle='\t\tauto fillRight=[&](const Group & group) {'
text=replace_once(text,needle,needle+'''\n\t\t\tconst auto profileRightBegin=std::chrono::steady_clock::now();
\t\t\t++profileGroups; if(group.all) ++profileAllGroups;
\t\t\tsize_t profileMin1=group.end.first,profileMin2=group.end.second;
\t\t\tforAnchors(group,[&](Pair q) { ++profileAnchors; profileMin1=std::min(profileMin1,q.first); profileMin2=std::min(profileMin2,q.second); });
\t\t\tprofileRightCells+=(group.end.first+1)*(group.end.second+1);
\t\t\tprofileCroppedCells+=(group.end.first-profileMin1+1)*(group.end.second-profileMin2+1);''',1)
needle='\t\t};\n\t\tfor(const auto & group:groups) {'
assert text.count(needle)==1
text=replace_once(text,needle,'''\t\t\tprofileRightSeconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-profileRightBegin).count();
'''+needle,1)
needle='\t});\n\treturn out;\n}\n}'
assert text.count(needle)==1
text=replace_once(text,needle,'''\t});
\tstd::cerr << "#profile_kernel sources=" << profileSources << " left_cells=" << profileLeftCells
\t\t<< " bounded_prefix_cells=" << profileBoundedLeftCells << " fillRight_calls=" << profileGroups << " unrestricted_groups=" << profileAllGroups
\t\t<< " right_rectangle_cells=" << profileRightCells << " minima_cropped_cells=" << profileCroppedCells
\t\t<< " firsthit_anchors=" << profileAnchors << " prefix_s=" << profileLeftSeconds
\t\t<< " fillRight_s=" << profileRightSeconds << '\\n';
\treturn out;
}
}''',1)
modified[files[1]]=text
patch=[]
hashes={}
for name in files:
    (output/(Path(name).name+'.base')).write_text(original[name])
    (build/name).write_text(modified[name])
    (output/(Path(name).name+'.instrumented')).write_text(modified[name])
    patch.extend(difflib.unified_diff(original[name].splitlines(True), modified[name].splitlines(True), fromfile='a/'+name,tofile='b/'+name))
    hashes[name]={'base_sha256':hashlib.sha256(original[name].encode()).hexdigest(),'instrumented_sha256':hashlib.sha256(modified[name].encode()).hexdigest()}
(output/'instrumentation.patch').write_text(''.join(patch))
(output/'source-hashes.json').write_text(json.dumps(hashes,indent=2)+'\n')
