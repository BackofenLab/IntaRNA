from pathlib import Path
import difflib
p=Path('src/IntaRNA/SeededPartitionFunction.cpp')
s=p.read_text()
s=s.replace('#include <array>','#include <array>\n#include <unordered_map>\n#include <vector>')
s=s.replace('\t\tMatrix<Z_type> h(matrix_detail::product(n,m),states),dh;', '''\t\t// Store all suffix states together for each admitted local base pair. The
\t\t// state count depends on the runtime seed length, so std::array cannot fit.
\t\t// Pair output appends reverse states to the same allocation; partition-only
\t\t// calls allocate only the forward states. Sequence 2 remains reversed.
\t\tstd::unordered_map<Pair,std::vector<Z_type>,Pair::Hash> state;
\t\tconst size_t valuesPerPair=matrix_detail::product(states,pairs?2:1);''')
s=s.replace('\t\tif (pairs) { dh.resize(n*m,states); dt.resize(n,m); }', '\t\tif (pairs) dt.resize(n,m);')
s=s.replace('// Cache vertex admission and seed endings in this start\'s cropped domain.\n', '// Cache vertex admission and seed endings in this start\'s cropped domain.\n\t\tsize_t admitted=0;\n')
s=s.replace('\t\t\tvalid(i,j)=w.valid({p1+i,p2+j});','\t\t\tvalid(i,j)=w.valid({p1+i,p2+j});\n\t\t\tadmitted+=valid(i,j)!=0;')
s=s.replace('\t\tauto accepts=', '''\t\t// Reserve exactly the admitted count and finish insertion before using any
\t\t// state pointers. Traverse coordinates, never hash iteration order.
\t\tstate.reserve(admitted);
\t\tfor (size_t i=0;i<n;++i) for(size_t j=0;j<m;++j)
\t\t\tif(valid(i,j)) state.try_emplace(Pair{i,j},valuesPerPair,Z_type(0));
\t\tauto accepts=''',1)
s=s.replace('h(0,initial)=init;', 'state.at({0,0})[initial]=init;')
s=s.replace('\t\t\tconst size_t q=i*m+j;','\t\t\tZ_type * const h=state.at({i,j}).data();')
s=s.replace('''\t\t\t\tfor(size_t c=0;c<states;++c) {
\t\t\t\t\tauto & v=h(q,dest(i,j,c,true));v=A::add(v,A::multiply(h(q-m-1,c),t));''','''\t\t\t\tconst Z_type * const previous=state.at({i-1,j-1}).data();
\t\t\t\tfor(size_t c=0;c<states;++c) {
\t\t\t\t\tauto & v=h[dest(i,j,c,true)];v=A::add(v,A::multiply(previous[c],t));''')
s=s.replace('h(q,dest(i,j,c? s2:0,false))','h[dest(i,j,c? s2:0,false)]')
s=s.replace('h(q,c)', 'h[c]').replace('h(q,s2)', 'h[s2]').replace('h(q,s1)', 'h[s1]')
s=s.replace('if(pairs) { dh[s2]=b;if(!d.noLP) dh[s1]=b; }', 'if(pairs) { h[states+s2]=b;if(!d.noLP) h[states+s1]=b; }')
s=s.replace('const size_t i=ii-1,j=jj-1,q=i*m+j; if(!valid(i,j)) continue;', 'const size_t i=ii-1,j=jj-1; if(!valid(i,j)) continue;\n\t\t\tZ_type * const dh=state.at({i,j}).data()+states;')
s=s.replace('dh(q,dest(i,j,c,true))','dh[dest(i,j,c,true)]')
s=s.replace('''\t\t\t\tfor(size_t c=0;c<states;++c) {
\t\t\t\t\tconst Z_type a=A::multiply(dh[dest(i,j,c,true)],t);
\t\t\t\t\tdh(q-m-1,c)=A::add(dh(q-m-1,c),a);
\t\t\t\t\tmass=A::add(mass,A::multiply(a,h(q-m-1,c)));''','''\t\t\t\tZ_type * const previous=state.at({i-1,j-1}).data();
\t\t\t\tfor(size_t c=0;c<states;++c) {
\t\t\t\t\tconst Z_type a=A::multiply(dh[dest(i,j,c,true)],t);
\t\t\t\t\tprevious[states+c]=A::add(previous[states+c],a);
\t\t\t\t\tmass=A::add(mass,A::multiply(a,previous[c]));''')
s=s.replace('dh(q,dest(i,j,c?s2:0,false))','dh[dest(i,j,c?s2:0,false)]')
s=s.replace('dh(0,initial)','state.at({0,0})[states+initial]')
assert 'h(q' not in s and 'q=' not in s
out=Path('/home/alex/CodexProjects/IntaRNA/issue-262-validation/hash-prototype/SeededPartitionFunction.cpp')
out.write_text(s)
Path('doc/seed-partition-storage/hash-storage.patch').write_text(''.join(difflib.unified_diff(p.read_text().splitlines(keepends=True),s.splitlines(keepends=True),fromfile='a/src/IntaRNA/SeededPartitionFunction.cpp',tofile='b/src/IntaRNA/SeededPartitionFunction.cpp')))
