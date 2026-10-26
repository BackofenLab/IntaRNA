#include "IntaRNA/SeededPartitionFunction.h"
#include "IntaRNA/PartitionArithmetic.h"
#include <algorithm>
#include <array>
#include <unordered_map>
#include <vector>

namespace IntaRNA {
namespace {
using Kernel=SeededPartitionFunction;
using Pair=Kernel::Pair;
using A=PartitionArithmetic;
using Support=std::array<Z_type,2>;
using RightStates=std::array<Z_type,4>;

/** Shared left-most-anchor factorization, with bounded boundary publication. */
Kernel::Result
computeAnchored(const Kernel::Domain & d,const StackSeedDomain & seeds,
		const Kernel::Weights & w,bool pairs)
{
	Kernel::Result out;
	if(pairs) out.mass=Matrix<Z_type>(d.n,d.m,0);
	if(!d.n || !d.m || !seeds.maxSeedLength() || !d.span1 || !d.span2) return out;

	// A containing seed cannot be the first completed owner: its inner seed
	// finishes earlier, or finishes together and is shorter. Among the remaining
	// occurrences on a path, start and end order agree; reversing their order
	// would require overlapping diagonal stacks and hence containment. Therefore
	// the left-most minimal anchor is also the heuristic's first completed seed,
	// including its shortest co-ending tie break. Both modes exclude earlier
	// minimal seeds while extending left and leave the right side unrestricted.
	const auto anchors=seeds.minimalOccurrences();
	std::unordered_map<Pair,size_t,Pair::Hash> anchorIndex;
	anchorIndex.reserve(anchors.size());
	std::vector<Pair> retainedRight;
	if(w.selectRight) retainedRight.reserve(anchors.size());
	size_t cap=2;
	for(size_t a=0;a<anchors.size();++a) {
		anchorIndex.emplace(Pair(anchors[a].i,anchors[a].j),a);
		cap=std::max(cap,anchors[a].length);
		// Freeze the endpoint once per canonical owner, before batches or the
		// outside pass can recompute its reusable extension tables.
		if(w.selectRight) retainedRight.push_back(w.selectRight(
				{anchors[a].i,anchors[a].j},
				{anchors[a].i+anchors[a].length-1,anchors[a].j+anchors[a].length-1}));
	}
	const size_t maxN=std::min(d.n,d.span1),maxM=std::min(d.m,d.span2);
	const size_t cells=matrix_detail::product(maxN,maxM);
	// The extension tables belong to the workspace, not individual seeds. Left
	// states record a capped initial stack length and anchor support; right
	// states record anchor support and support of their current frontier pair.
	Matrix<Support> left(cells,cap),leftTotal(maxN,maxM),rightTotal(maxN,maxM);
	Matrix<RightStates> right(maxN,maxM);
	Matrix<Support> dLeft,dLeftTotal;
	Matrix<RightStates> dRight;
	if(pairs) {
		dLeft.resize(cells,cap,false); dLeftTotal.resize(maxN,maxM,false);
		dRight.resize(maxN,maxM,false);
	}
	Matrix<unsigned char> valid(d.n,d.m);
	for(size_t i=0;i<d.n;++i) for(size_t j=0;j<d.m;++j) valid(i,j)=w.valid({i,j});

	// Publish complete boundaries once, in row-major order, without a global
	// four-coordinate table. At most 256 left boundaries share a batch. Exact
	// and outside computations target 2 MiB per table; the ordinary heuristic
	// may use 16 MiB because it stores no boundary/adjoint table and benefits
	// directly from recomputing each anchor across fewer batches.
	constexpr size_t batchLimit=256;
	const size_t batchBytes=w.selectRight && !pairs ? 16*1024*1024 : 2*1024*1024;
	const size_t batchSize=std::max(size_t(1),std::min(batchLimit,batchBytes/sizeof(Z_type)/cells));
	Matrix<Z_type> hybrid(batchSize,cells),boundary;
	if(pairs) boundary.resize(batchSize,cells,false);
	std::vector<Pair> starts;
	std::vector<Z_type> initiation;
	std::vector<size_t> candidates;
	std::vector<unsigned char> included(anchors.size(),0);
	starts.reserve(batchSize); initiation.reserve(batchSize);
	candidates.reserve(anchors.size());
	auto reset=[]<class T>(Matrix<T> & matrix,size_t rows,size_t columns) {
		if(columns==matrix.size2()) std::fill_n(&matrix(0,0),rows*columns,T{});
		else for(size_t i=0;i<rows;++i) std::fill_n(&matrix(i,0),columns,T{});
	};
	auto resetLeft=[&](Matrix<Support> & matrix,size_t rows,size_t columns) {
		// Left states use a fixed maxM row stride and a contiguous stack axis.
		for(size_t i=0;i<rows;++i) std::fill_n(&matrix(i*maxM,0),columns*cap,Support{});
	};

	// These dimensions and endpoints describe the anchor currently in the shared
	// tables. Offset zero is a stored pair, so singleton seeds need a full span.
	Pair s{0,0},t{0,0};
	Pair rightOffset{0,0};
	size_t length=0,nL=0,mL=0,nR=0,mR=0,supportBegin=0;
	bool singleton=false,selectedRight=false;
	Z_type seedWeight=0;
	auto leftPair=[&](size_t i,size_t j) { return Pair(s.first-i,s.second-j); };
	auto rightPair=[&](size_t i,size_t j) { return Pair(t.first+i,t.second+j); };
	auto earlierLength=[&](Pair p) {
		const auto found=anchorIndex.find(p);
		return found==anchorIndex.end()?0:anchors[found->second].length;
	};
	auto nonstack=[&](size_t i,size_t j,const auto & visit) {
		for(size_t u=1;u<=i && u-1<=d.loop1;++u)
		for(size_t v=1;v<=j && v-1<=d.loop2;++v)
			if(u!=1 || v!=1) visit(i-u,j-v);
	};
	auto fill=[&](size_t a) {
		const auto anchor=anchors[a];
		s={anchor.i,anchor.j}; length=anchor.length;
		t={s.first+length-1,s.second+length-1};
		singleton=d.noLP && length==1;
		supportBegin=singleton?0:1;
		// Crop both extension rectangles to the batch's possible left boundaries.
		// Paths to these sources cannot visit positions beyond their leftmost
		// coordinates; their mandatory left offsets also bound the right side.
		nL=0; mL=0;
		size_t minLeft1=maxN,minLeft2=maxM;
		bool unrestrictedRight=!w.selectRight;
		for(const Pair p:starts) {
			if(p.first>s.first || p.second>s.second) continue;
			const size_t u=s.first-p.first,v=s.second-p.second;
			if(u>maxN-length || v>maxM-length) continue;
			nL=std::max(nL,u+1); mL=std::max(mL,v+1);
			minLeft1=std::min(minLeft1,u); minLeft2=std::min(minLeft2,v);
			if(p==s) unrestrictedRight=true;
		}
		if(!nL || !mL) { seedWeight=0; return; }
		nR=std::min(d.n-t.first,maxN-length-minLeft1+1);
		mR=std::min(d.m-t.second,maxM-length-minLeft2+1);
		selectedRight=false;
		if(w.selectRight) {
			const Pair end=retainedRight[a];
			// Only a strictly later endpoint denotes a two-sided extension. Check
			// coordinates before subtracting, including user-supplied selectors.
			if(end.first>t.first && end.second>t.second && end.first<d.n && end.second<d.m) {
				rightOffset={end.first-t.first,end.second-t.second};
				selectedRight=rightOffset.first<nR && rightOffset.second<mR;
			}
			// A no-left source retains all right endpoints. Other batches need
			// only paths to the selected endpoint and the seed-only cell (0,0).
			if(!unrestrictedRight) {
				nR=selectedRight?rightOffset.first+1:1;
				mR=selectedRight?rightOffset.second+1:1;
			}
		}
		resetLeft(left,nL,mL); reset(leftTotal,nL,mL);
		reset(right,nR,mR); reset(rightTotal,nR,mR);
		// Seed energy consists of its actual stack transitions. Initiation is
		// attached later to the complete interaction's left boundary exactly once.
		seedWeight=1;
		for(size_t k=0;k<length;++k) {
			if(!valid(s.first+k,s.second+k)) { seedWeight=0; break; }
			if(k) seedWeight=A::multiply(seedWeight,A::check(w.transition(
					{s.first+k-1,s.second+k-1},{s.first+k,s.second+k})));
		}
		if(seedWeight==0) return;
		left(0,length-1)[singleton?0:1]=1;
		for(size_t i=0;i<nL;++i) for(size_t j=0;j<mL;++j) {
			const Pair p=leftPair(i,j); if(!valid(p.first,p.second)) continue;
			const size_t forbidden=earlierLength(p);
			if(i && j) {
				const Pair q=leftPair(i-1,j-1);
				if(valid(q.first,q.second)) {
					const Z_type edge=A::check(w.transition(p,q));
					for(size_t c=0;c<cap;++c) {
						const size_t run=std::min(cap,c+2); if(forbidden && forbidden<=run) continue;
						for(size_t h=supportBegin;h<2;++h) {
							const Z_type previous=left((i-1)*maxM+j-1,c)[h];
							if(previous==0) continue;
							const size_t support=singleton && i==1 && j==1?1:h;
							auto & value=left(i*maxM+j,run-1)[support];
							value=A::add(value,A::multiply(previous,edge));
						}
					}
				}
			}
			// An initial singleton may obtain support from the other extension.
			// Other unsupported frontier pairs cannot be left by a nonstack edge.
			if(forbidden!=1) nonstack(i,j,[&](size_t u,size_t v) {
				const Pair q=leftPair(u,v); if(!valid(q.first,q.second)) return;
				if(leftTotal(u,v)[0]==0 && leftTotal(u,v)[1]==0) return;
				const Z_type edge=A::check(w.transition(p,q));
				for(size_t h=supportBegin;h<2;++h) left(i*maxM+j,0)[h]=A::add(
						left(i*maxM+j,0)[h],A::multiply(leftTotal(u,v)[h],edge));
			});
			for(size_t c=d.noLP && !(singleton && i==0 && j==0)?1:0;c<cap;++c)
				for(size_t h=supportBegin;h<2;++h) leftTotal(i,j)[h]=A::add(leftTotal(i,j)[h],left(i*maxM+j,c)[h]);
		}
		// An anchor whose seed-free left extension reaches no source in this
		// batch cannot contribute, so its right recurrence need not be evaluated.
		bool reachable=false;
		for(const Pair p:starts) {
			if(p.first>s.first || p.second>s.second) continue;
			const size_t u=s.first-p.first,v=s.second-p.second;
			if(u<nL && v<mL && (leftTotal(u,v)[0]>0 || leftTotal(u,v)[1]>0)) { reachable=true; break; }
		}
		if(!reachable) { seedWeight=0; return; }
		right(0,0)[singleton?0:3]=1;
		for(size_t i=0;i<nR;++i) for(size_t j=0;j<mR;++j) {
			const Pair q=rightPair(i,j); if(!valid(q.first,q.second)) continue;
			if(i && j) {
				const Pair p=rightPair(i-1,j-1);
				if(valid(p.first,p.second)) {
					const Z_type edge=A::check(w.transition(p,q));
					for(size_t h=supportBegin;h<2;++h) {
						const size_t support=singleton && i==1 && j==1?1:h;
						const auto & prev=right(i-1,j-1);
						auto & value=right(i,j)[support*2+1];
						value=A::add(value,A::multiply(A::add(prev[h*2],prev[h*2+1]),edge));
					}
				}
			}
			nonstack(i,j,[&](size_t u,size_t v) {
				const Pair p=rightPair(u,v); if(!valid(p.first,p.second)) return;
				if(rightTotal(u,v)[0]==0 && rightTotal(u,v)[1]==0) return;
				const Z_type edge=A::check(w.transition(p,q));
				for(size_t h=supportBegin;h<2;++h) {
					auto & value=right(i,j)[h*2+(d.noLP?0:1)];
					value=A::add(value,A::multiply(rightTotal(u,v)[h],edge));
				}
			});
			for(size_t h=supportBegin;h<2;++h) rightTotal(i,j)[h]=singleton && i==0 && j==0
					?A::add(right(i,j)[h*2],right(i,j)[h*2+1]):right(i,j)[h*2+1];
		}
	};

	auto contract=[&](bool outside) {
		if(seedWeight==0) return;
		if(outside) {
			resetLeft(dLeft,nL,mL); reset(dLeftTotal,nL,mL); reset(dRight,nR,mR);
		}
		Z_type seedMass=0;
		for(size_t b=0;b<starts.size();++b) {
			const Pair p=starts[b];
			if(p.first>s.first || p.second>s.second) continue;
			const size_t u=s.first-p.first,v=s.second-p.second;
			if(u>=nL || v>=mL || (leftTotal(u,v)[0]==0 && leftTotal(u,v)[1]==0)) continue;
			const Z_type source=A::multiply(initiation[b],seedWeight);
			Support sourceLeft{};
			for(size_t h=supportBegin;h<2;++h) sourceLeft[h]=A::multiply(source,leftTotal(u,v)[h]);
			auto joinRight=[&](size_t i,size_t j) {
				const size_t endpoint=(u+length-1+i)*maxM+v+length-1+j;
				for(size_t h=supportBegin;h<2;++h) for(size_t k=supportBegin;k<2;++k) {
					if(singleton && h==0 && k==0) continue;
					const Z_type l=leftTotal(u,v)[h],r=rightTotal(i,j)[k];
					if(l==0 || r==0) continue;
					const Z_type value=A::multiply(sourceLeft[h],r);
					if(!outside) {
						hybrid(b,endpoint)=A::add(hybrid(b,endpoint),value);
						continue;
					}
					const Z_type coefficient=boundary(b,endpoint);
					if(coefficient==0) continue;
					seedMass=A::add(seedMass,A::multiply(value,coefficient));
					const Z_type adj=A::multiply(source,coefficient);
					dLeftTotal(u,v)[h]=A::add(dLeftTotal(u,v)[h],A::multiply(adj,r));
					const size_t frontier=singleton && i==0 && j==0?0:1;
					dRight(i,j)[k*2+frontier]=A::add(dRight(i,j)[k*2+frontier],A::multiply(adj,l));
				}
			};
			if(!w.selectRight || (u==0 && v==0)) {
				// Exact prediction and a heuristic seed with no left extension
				// retain every right endpoint in the source's inclusive span.
				for(size_t i=0;i<nR && u+length+i<=maxN;++i)
				for(size_t j=0;j<mR && v+length+j<=maxM;++j) joinRight(i,j);
			} else {
				// True left extensions retain the seed ending plus at most one
				// chosen right endpoint. Visit the same two roots in both passes.
				joinRight(0,0);
				if(selectedRight && rightOffset.first<=maxN-length-u && rightOffset.second<=maxM-length-v)
					joinRight(rightOffset.first,rightOffset.second);
			}
		}
		if(!outside || seedMass==0) return;
		// Every forced seed pair belongs to the anchor exactly once. Extension
		// marking below excludes the two initial frontier cells shared with it.
		for(size_t k=0;k<length;++k) out.mass(s.first+k,s.second+k)=A::add(out.mass(s.first+k,s.second+k),seedMass);
		for(size_t ii=nL;ii>0;--ii) for(size_t jj=mL;jj>0;--jj) {
			const size_t i=ii-1,j=jj-1; const Pair p=leftPair(i,j);
			if(!valid(p.first,p.second)) continue;
			const size_t forbidden=earlierLength(p);
			for(size_t c=d.noLP && !(singleton && i==0 && j==0)?1:0;c<cap;++c)
				for(size_t h=supportBegin;h<2;++h) dLeft(i*maxM+j,c)[h]=A::add(dLeft(i*maxM+j,c)[h],dLeftTotal(i,j)[h]);
			// Nonnegative forward zeros have no contributing ancestors. Discard
			// their adjoints before evaluating transitions in the reverse pass.
			for(size_t c=0;c<cap;++c) for(size_t h=supportBegin;h<2;++h)
				if(left(i*maxM+j,c)[h]==0) dLeft(i*maxM+j,c)[h]=0;
			if(i || j) for(size_t c=0;c<cap;++c) for(size_t h=supportBegin;h<2;++h)
				out.mass(p.first,p.second)=A::add(out.mass(p.first,p.second),A::multiply(left(i*maxM+j,c)[h],dLeft(i*maxM+j,c)[h]));
			if(i && j) {
				const Pair q=leftPair(i-1,j-1);
				if(valid(q.first,q.second)) {
					const Z_type edge=A::check(w.transition(p,q));
					for(size_t c=0;c<cap;++c) {
						const size_t run=std::min(cap,c+2); if(forbidden && forbidden<=run) continue;
						for(size_t h=supportBegin;h<2;++h) {
							const size_t support=singleton && i==1 && j==1?1:h;
							if(left((i-1)*maxM+j-1,c)[h]==0 || dLeft(i*maxM+j,run-1)[support]==0) continue;
							auto & adj=dLeft((i-1)*maxM+j-1,c)[h];
							adj=A::add(adj,A::multiply(edge,dLeft(i*maxM+j,run-1)[support]));
						}
					}
				}
			}
			if(forbidden!=1) nonstack(i,j,[&](size_t u,size_t v) {
				const Pair q=leftPair(u,v); if(!valid(q.first,q.second)) return;
				if((dLeft(i*maxM+j,0)[0]==0 && dLeft(i*maxM+j,0)[1]==0)
						|| (leftTotal(u,v)[0]==0 && leftTotal(u,v)[1]==0)) return;
				const Z_type edge=A::check(w.transition(p,q));
				for(size_t h=supportBegin;h<2;++h) dLeftTotal(u,v)[h]=A::add(dLeftTotal(u,v)[h],A::multiply(edge,dLeft(i*maxM+j,0)[h]));
			});
		}
		for(size_t ii=nR;ii>0;--ii) for(size_t jj=mR;jj>0;--jj) {
			const size_t i=ii-1,j=jj-1; const Pair q=rightPair(i,j);
			if(!valid(q.first,q.second)) continue;
			for(size_t c=supportBegin*2;c<4;++c) if(right(i,j)[c]==0) dRight(i,j)[c]=0;
			if(i || j) for(size_t c=supportBegin*2;c<4;++c) out.mass(q.first,q.second)=A::add(
					out.mass(q.first,q.second),A::multiply(right(i,j)[c],dRight(i,j)[c]));
			if(i && j) {
				const Pair p=rightPair(i-1,j-1);
				if(valid(p.first,p.second)) {
					const Z_type edge=A::check(w.transition(p,q));
					for(size_t h=supportBegin;h<2;++h) {
						const size_t support=singleton && i==1 && j==1?1:h;
						if(dRight(i,j)[support*2+1]==0) continue;
						const Z_type adj=A::multiply(edge,dRight(i,j)[support*2+1]);
						for(size_t c=0;c<2;++c) if(right(i-1,j-1)[h*2+c]>0)
							dRight(i-1,j-1)[h*2+c]=A::add(dRight(i-1,j-1)[h*2+c],adj);
					}
				}
			}
			nonstack(i,j,[&](size_t u,size_t v) {
				const Pair p=rightPair(u,v); if(!valid(p.first,p.second)) return;
				if((dRight(i,j)[d.noLP?0:1]==0 && dRight(i,j)[2+(d.noLP?0:1)]==0)
						|| (rightTotal(u,v)[0]==0 && rightTotal(u,v)[1]==0)) return;
				const Z_type edge=A::check(w.transition(p,q));
				for(size_t h=supportBegin;h<2;++h) {
					const size_t frontier=singleton && u==0 && v==0?0:1;
					auto & adj=dRight(u,v)[h*2+frontier];
					adj=A::add(adj,A::multiply(edge,dRight(i,j)[h*2+(d.noLP?0:1)]));
				}
			});
		}
	};

	auto processBatch=[&]() {
		reset(hybrid,starts.size(),cells); if(pairs) reset(boundary,starts.size(),cells);
		candidates.clear(); initiation.clear();
		// Query only anchor starts inside the batch's span boxes. Sparse seed
		// domains are cheaper to scan directly than probing every possible pair
		// in every source span; dense domains retain the bounded grid probe.
		for(const Pair p:starts) {
			initiation.push_back(A::check(w.initiation(p)));
		}
		if(anchors.size()<cells) {
			for(size_t a=0;a<anchors.size();++a) {
				const auto & anchor=anchors[a];
				for(const Pair p:starts) {
					if(p.first>anchor.i || p.second>anchor.j) continue;
					const size_t i=anchor.i-p.first,j=anchor.j-p.second;
					if(i>=maxN || j>=maxM || anchor.length>maxN-i || anchor.length>maxM-j) continue;
					included[a]=1; candidates.push_back(a); break;
				}
			}
		} else {
			for(const Pair p:starts)
			for(size_t i=p.first;i<d.n && i-p.first<maxN;++i)
			for(size_t j=p.second;j<d.m && j-p.second<maxM;++j) {
				const auto found=anchorIndex.find({i,j});
				if(found==anchorIndex.end() || included[found->second]) continue;
				const auto & anchor=anchors[found->second];
				if(anchor.length>maxN-(i-p.first) || anchor.length>maxM-(j-p.second)) continue;
				included[found->second]=1; candidates.push_back(found->second);
			}
			std::sort(candidates.begin(),candidates.end());
		}
		for(const size_t a:candidates) { fill(a); contract(false); }
		for(size_t b=0;b<starts.size();++b) {
			const Pair p=starts[b];
			for(size_t i=0;i<maxN && i<d.n-p.first;++i)
			for(size_t j=0;j<maxM && j<d.m-p.second;++j) {
				const Z_type value=hybrid(b,i*maxM+j); if(value==0) continue;
				const Pair q{p.first+i,p.second+j};
				const Z_type coefficient=A::check(w.boundary(p,q));
				if(pairs) boundary(b,i*maxM+j)=coefficient;
				out.z=A::add(out.z,A::multiply(value,coefficient));
				if(w.complete) w.complete(p,q,value,coefficient);
			}
		}
		// Recompute the same anchor factors and differentiate using the exact
		// published coefficients, including zero-valued excluded boundaries.
		if(pairs) for(const size_t a:candidates) { fill(a); contract(true); }
		for(const size_t a:candidates) included[a]=0;
		starts.clear();
	};
	seeds.forEachStart([&](size_t i,size_t j) {
		if(!valid(i,j)) return;
		starts.emplace_back(i,j);
		if(starts.size()==batchSize) processBatch();
	});
	if(!starts.empty()) processBatch();
	return out;
}
} // namespace

SeededPartitionFunction::Result
SeededPartitionFunction::
compute(const Domain & d, const StackSeedDomain & seeds, const Weights & w, bool pairs)
{
	return computeAnchored(d,seeds,w,pairs);
}
}
