#include "IntaRNA/SeededPartitionFunction.h"
#include "IntaRNA/PartitionArithmetic.h"
#include <array>
#include <unordered_map>
#include <vector>

// Seed-free suffix lengths plus two absorbing seeded states (singleton/run).
// Nonstack edges operate on aggregate totals, independently of seed length.
namespace IntaRNA {
SeededPartitionFunction::Result 
SeededPartitionFunction::
compute(
		const IntaRNA::SeededPartitionFunction::Domain & d,
		const IntaRNA::StackSeedDomain & seeds,
		const IntaRNA::SeededPartitionFunction::Weights & w, bool pairs)
{
	using namespace IntaRNA;
	using A=PartitionArithmetic;
	using Pair=SeededPartitionFunction::Pair;
	SeededPartitionFunction::Result out;
	// Empty seed families or zero spans admit no interaction, even with pair output.
	if (pairs) out.mass=Matrix<Z_type>(d.n,d.m,0);
	if (!seeds.maxSeedLength() || !d.span1 || !d.span2) return out;
	// At least two length states are necessary for noLP even with singleton seeds.
	const size_t cap=std::max(size_t(2),seeds.maxSeedLength()), s1=cap, s2=cap+1, states=cap+2;
	seeds.forEachStart([&](size_t p1,size_t p2) {
		// Each eligible left boundary owns a separate forward and reverse problem.
		const Pair p{p1,p2}; if (!w.valid(p)) return;
		const size_t n=std::min(d.span1,d.n-p1),m=std::min(d.span2,d.m-p2);
		// Store all suffix states together for each admitted local base pair. The
		// state count depends on the runtime seed length, so std::array cannot fit.
		// Pair output appends reverse states to the same allocation; partition-only
		// calls allocate only the forward states. Sequence 2 remains reversed.
		std::unordered_map<Pair,std::vector<Z_type>,Pair::Hash> state;
		const size_t valuesPerPair=matrix_detail::product(states,pairs?2:1);
		Matrix<std::array<Z_type,2>> total(n,m),dt;
		Matrix<unsigned char> valid(n,m);
		Matrix<size_t> ending(n,m);
		if (pairs) dt.resize(n,m);
		// Cache vertex admission and seed endings in this start's cropped domain.
		size_t admitted=0;
		for (size_t i=0;i<n;++i) for(size_t j=0;j<m;++j) {
			valid(i,j)=w.valid({p1+i,p2+j});
			admitted+=valid(i,j)!=0;
			ending(i,j)=seeds.seedEndingLength(p1+i,p2+j);
		}
		// Reserve exactly the admitted count and finish insertion before using any
		// state pointers. Traverse coordinates, never hash iteration order.
		state.reserve(admitted);
		for (size_t i=0;i<n;++i) for(size_t j=0;j<m;++j)
			if(valid(i,j)) state.try_emplace(Pair{i,j},valuesPerPair,Z_type(0));
		auto accepts=[&](size_t i,size_t j,size_t len) { auto e=ending(i,j);return e && e<=len; };
		// A stack extends the current suffix; a loop resets it to one pair. Once
		// seeded, retain only whether the final pair still needs a noLP partner.
		auto dest=[&](size_t i,size_t j,size_t c,bool stacked) {
			if (c>=cap) return stacked || !d.noLP?s2:s1;
			const size_t len=stacked?std::min(cap,c+2):1;
			return accepts(i,j,len)?(len>=2 || !d.noLP?s2:s1):len-1;
		};
		// Enumerate loop predecessors separately from the direct stacking edge.
		auto nonstack=[&](size_t i,size_t j,const auto & visit) {
			for(size_t u=1;u<=i && u-1<=d.loop1;++u) for(size_t v=1;v<=j && v-1<=d.loop2;++v)
				if((u!=1 || v!=1) && valid(i-u,j-v)) visit(i-u,j-v);
		};
		// Initiation belongs to the leftmost pair exactly once, including singletons.
		const size_t initial=accepts(0,0,1)?(d.noLP?s1:s2):0;
		const Z_type init=A::check(w.initiation(p)); state.at({0,0})[initial]=init;
		// Forward pass: extend prefixes and finalize every admitted right boundary.
		for(size_t i=0;i<n;++i) for(size_t j=0;j<m;++j) {
			if(!valid(i,j)) continue;
			Z_type * const h=state.at({i,j}).data();
			// Stacking transitions preserve the suffix length needed for seed admission.
			if(i && j && valid(i-1,j-1)) {
				const Z_type t=A::check(w.transition({p1+i-1,p2+j-1},{p1+i,p2+j}));
				const Z_type * const previous=state.at({i-1,j-1}).data();
				for(size_t c=0;c<states;++c) {
					auto & v=h[dest(i,j,c,true)];v=A::add(v,A::multiply(previous[c],t));
				}
			}
			// Loop transitions need only the unseeded/seeded predecessor totals.
			nonstack(i,j,[&](size_t u,size_t v) {
				if(total(u,v)[0]==0 && total(u,v)[1]==0) return;
				const Z_type t=A::check(w.transition({p1+u,p2+v},{p1+i,p2+j}));
				for(size_t c=0;c<2;++c) {
					auto & value=h[dest(i,j,c? s2:0,false)];value=A::add(value,A::multiply(total(u,v)[c],t));
				}
			});
			// Exclude lonely terminal pairs from loop predecessors and completed paths.
			for(size_t c=d.noLP?1:0;c<cap;++c) total(i,j)[0]=A::add(total(i,j)[0],h[c]);
			total(i,j)[1]=A::add(h[s2],d.noLP?Z_type(0):h[s1]);
			if(total(i,j)[1]>0) {
				// Share the same boundary coefficient with reporting and the outside pass.
				const Z_type b=A::check(w.boundary(p,{p1+i,p2+j}));
				out.z=A::add(out.z,A::multiply(b,total(i,j)[1]));
				if(pairs) { h[states+s2]=b;if(!d.noLP) h[states+s1]=b; }
				if(w.complete) w.complete(p,{p1+i,p2+j},total(i,j)[1],b);
			}
		}
		if(!pairs) return;
		// Reverse pass: propagate boundary derivatives back through the same edges.
		for(size_t ii=n;ii>0;--ii) for(size_t jj=m;jj>0;--jj) {
			const size_t i=ii-1,j=jj-1; if(!valid(i,j)) continue;
			Z_type * const dh=state.at({i,j}).data()+states;
			// Distribute derivatives of the two aggregate totals to their suffix states.
			for(size_t c=d.noLP?1:0;c<cap;++c) dh[c]=A::add(dh[c],dt(i,j)[0]);
			dh[s2]=A::add(dh[s2],dt(i,j)[1]);
			if(!d.noLP) dh[s1]=A::add(dh[s1],dt(i,j)[1]);
			auto & mass=out.mass(p1+i,p2+j);
			// Every incoming edge owns its destination pair. Its forward weight times
			// the propagated derivative contributes that pair's unnormalized mass.
			if(i && j && valid(i-1,j-1)) {
				const Z_type t=A::check(w.transition({p1+i-1,p2+j-1},{p1+i,p2+j}));
				Z_type * const previous=state.at({i-1,j-1}).data();
				for(size_t c=0;c<states;++c) {
					const Z_type a=A::multiply(dh[dest(i,j,c,true)],t);
					previous[states+c]=A::add(previous[states+c],a);
					mass=A::add(mass,A::multiply(a,previous[c]));
				}
			}
			nonstack(i,j,[&](size_t u,size_t v) {
				const Z_type t=A::check(w.transition({p1+u,p2+v},{p1+i,p2+j}));
				for(size_t c=0;c<2;++c) {
					const Z_type a=A::multiply(dh[dest(i,j,c?s2:0,false)],t);
					dt(u,v)[c]=A::add(dt(u,v)[c],a);
					mass=A::add(mass,A::multiply(a,total(u,v)[c]));
				}
			});
		}
		// The initial pair has no incoming edge, so collect its initiation term here.
		out.mass(p1,p2)=A::add(out.mass(p1,p2),A::multiply(state.at({0,0})[states+initial],init));
	});
	return out;
}
}
