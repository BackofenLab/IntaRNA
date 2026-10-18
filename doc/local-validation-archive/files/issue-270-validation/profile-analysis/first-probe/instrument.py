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
files=['src/IntaRNA/PredictorMfeEns2dSeedExtension.cpp','src/IntaRNA/SeededPartitionFunction.cpp']
original={name:(base/name).read_text() for name in files}
text=original[files[0]]
text='#include <chrono>\n#include <iostream>\n'+text
text=text.replace('\tprepareStackSeedExtensions(n,m,seeds);','\tconst auto profileBegin=std::chrono::steady_clock::now();\n\tprepareStackSeedExtensions(n,m,seeds);\n\tconst auto profilePrepared=std::chrono::steady_clock::now();',1)
text=text.replace('\tconst auto result=K::compute(domain,seeds,weights,pairProbabilities!=nullptr);','\tconst auto result=K::compute(domain,seeds,weights,pairProbabilities!=nullptr);\n\tconst auto profileComputed=std::chrono::steady_clock::now();',1)
needle='\toutput.setExactPartition(true);\n\treportOptima();'
replacement=needle+'''\n\tconst auto profileReported=std::chrono::steady_clock::now();
\tstd::cerr << "#profile_predictor prepare_s=" << std::chrono::duration<double>(profilePrepared-profileBegin).count()
\t\t<< " compute_s=" << std::chrono::duration<double>(profileComputed-profilePrepared).count()
\t\t<< " report_s=" << std::chrono::duration<double>(profileReported-profileComputed).count() << '\\n';'''
assert text.count(needle)==1
text=text.replace(needle,replacement,1)
modified={files[0]:text}
text=original[files[1]]
text='#include <chrono>\n#include <iostream>\n#include <cstdint>\n'+text
needle='\tseeds.forEachStart([&](size_t p1,size_t p2) {'
assert text.count(needle)==1
text=text.replace(needle,'''\tuint64_t profileSources=0,profileLeftCells=0,profileGroups=0,profileAllGroups=0,profileRightCells=0,profileCroppedCells=0,profileAnchors=0;
\tdouble profileLeftSeconds=0,profileRightSeconds=0;
'''+needle,1)
needle='\t\tconst size_t n=std::min(d.span1,d.n-p1),m=std::min(d.span2,d.m-p2);'
assert text.count(needle)==1
text=text.replace(needle,needle+'''\n\t\t++profileSources; profileLeftCells+=n*m;
\t\tconst auto profileLeftBegin=std::chrono::steady_clock::now();''',1)
needle='\t\t// Retained right extensions share first-hit injections by linearity:'
text=text.replace(needle,'''\t\tprofileLeftSeconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-profileLeftBegin).count();
'''+needle,1)
needle='\t\tauto fillRight=[&](const Group & group) {'
text=text.replace(needle,needle+'''\n\t\t\tconst auto profileRightBegin=std::chrono::steady_clock::now();
\t\t\t++profileGroups; if(group.all) ++profileAllGroups;
\t\t\tsize_t profileMin1=group.end.first,profileMin2=group.end.second;
\t\t\tforAnchors(group,[&](Pair q) { ++profileAnchors; profileMin1=std::min(profileMin1,q.first); profileMin2=std::min(profileMin2,q.second); });
\t\t\tprofileRightCells+=(group.end.first+1)*(group.end.second+1);
\t\t\tprofileCroppedCells+=(group.end.first-profileMin1+1)*(group.end.second-profileMin2+1);''',1)
needle='\t\t};\n\t\tfor(const auto & group:groups) {'
assert text.count(needle)==1
text=text.replace(needle,'''\t\t\tprofileRightSeconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-profileRightBegin).count();
'''+needle,1)
needle='\t});\n\treturn out;\n}\n}'
assert text.count(needle)==1
text=text.replace(needle,'''\t});
\tstd::cerr << "#profile_kernel sources=" << profileSources << " left_cells=" << profileLeftCells
\t\t<< " fillRight_calls=" << profileGroups << " unrestricted_groups=" << profileAllGroups
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
