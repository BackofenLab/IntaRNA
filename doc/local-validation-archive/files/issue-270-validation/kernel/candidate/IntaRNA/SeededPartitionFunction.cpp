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

/** Exact left-most-anchor factorization, with bounded boundary publication. */
Kernel::Result
computeAnchored(const Kernel::Domain & d,const StackSeedDomain & seeds,
		const Kernel::Weights & w,bool pairs)
{
	Kernel::Result out;
	if(pairs) out.mass=Matrix<Z_type>(d.n,d.m,0);
	if(!d.n || !d.m || !seeds.maxSeedLength() || !d.span1 || !d.span2) return out;

	// Containing seeds add no paths to the exact union. In the remaining family,
	// start order equals end order along each diagonal, including overlapping
	// seeds. Thus excluding an earlier seed while extending left is sufficient;
	// no earlier seed can first finish beyond the anchor's right end.
	const auto anchors=seeds.minimalOccurrences();
	std::unordered_map<Pair,size_t,Pair::Hash> anchorIndex;
	anchorIndex.reserve(anchors.size());
	size_t cap=2;
	for(size_t a=0;a<anchors.size();++a) {
		anchorIndex.emplace(Pair(anchors[a].i,anchors[a].j),a);
		cap=std::max(cap,anchors[a].length);
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
	// four-coordinate table. At most 256 left boundaries share a batch, with a
	// 2 MiB budget per boundary table unless a single span already exceeds it.
	// Anchors crossing batches are recomputed in the same reusable workspace.
	constexpr size_t batchLimit=256,batchBytes=2*1024*1024;
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
	size_t length=0,nL=0,mL=0,nR=0,mR=0,supportBegin=0;
	bool singleton=false;
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
		for(const Pair p:starts) {
			if(p.first>s.first || p.second>s.second) continue;
			const size_t u=s.first-p.first,v=s.second-p.second;
			if(u>maxN-length || v>maxM-length) continue;
			nL=std::max(nL,u+1); mL=std::max(mL,v+1);
			minLeft1=std::min(minLeft1,u); minLeft2=std::min(minLeft2,v);
		}
		if(!nL || !mL) { seedWeight=0; return; }
		nR=std::min(d.n-t.first,maxN-length-minLeft1+1);
		mR=std::min(d.m-t.second,maxM-length-minLeft2+1);
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
			for(size_t i=0;i<nR && u+length+i<=maxN;++i)
			for(size_t j=0;j<mR && v+length+j<=maxM;++j) {
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
		// Query only anchor starts inside the batch's span boxes, then restore
		// deterministic seed order independently of hash-table iteration order.
		for(const Pair p:starts) {
			initiation.push_back(A::check(w.initiation(p)));
			for(size_t i=p.first;i<d.n && i-p.first<maxN;++i)
			for(size_t j=p.second;j<d.m && j-p.second<maxM;++j) {
				const auto found=anchorIndex.find({i,j});
				if(found==anchorIndex.end() || included[found->second]) continue;
				const auto & anchor=anchors[found->second];
				if(anchor.length>maxN-(i-p.first) || anchor.length>maxM-(j-p.second)) continue;
				included[found->second]=1; candidates.push_back(found->second);
			}
		}
		std::sort(candidates.begin(),candidates.end());
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
	// Exact prediction extends each left-most irredundant seed independently.
	// Frozen heuristic choices retain their established first-completed owner.
	if(!w.selectRight) return computeAnchored(d,seeds,w,pairs);
	using A=PartitionArithmetic;
	using Ends=std::array<Z_type,2>; // pending lonely pair / stack-supported end
	Result out;
	if (pairs) out.mass=Matrix<Z_type>(d.n,d.m,0);
	if (!seeds.maxSeedLength() || !d.span1 || !d.span2) return out;

	// Freeze heuristic choices before either pass. A first hit is owned by the
	// shortest admitted seed ending here, including coincident variable lengths.
	std::unordered_map<Pair,Pair,Pair::Hash> selected;
	if (w.selectRight) for(size_t i=0;i<d.n;++i) for(size_t j=0;j<d.m;++j) {
		const size_t length=seeds.seedEndingLength(i,j);
		if(length) selected.emplace(Pair(i,j),w.selectRight({i+1-length,j+1-length},{i,j}));
	}
	// Seed-free heuristic prefixes retain stack length; right extensions need
	// only noLP, together with the frozen right-endpoint selection.
	const size_t cap=std::max(size_t(2),seeds.maxSeedLength());
	// Allocate the maximal inclusive span once. Offset zero is a stored pair,
	// so an extension of at most span-1 positions needs span cells, even for
	// singleton seeds. Keep these shapes fixed and clear only active cells;
	// this also prevents stale values when active rectangles shrink/regrow.
	const size_t maxN=std::min(d.span1,d.n),maxM=std::min(d.span2,d.m);
	const size_t cells=matrix_detail::product(maxN,maxM);
	Matrix<Z_type> left(cells,cap),leftTotal(maxN,maxM),hybrid(maxN,maxM);
	Matrix<Ends> hit(maxN,maxM),right(maxN,maxM);
	Matrix<unsigned char> valid(maxN,maxM);
	Matrix<size_t> ending(maxN,maxM),byEnd;
	Matrix<Z_type> boundary,dLeft,dTotal;
	Matrix<Ends> dHit,adjoint;
	if(w.selectRight) byEnd.resize(maxN,maxM,false);
	if(pairs) {
		boundary.resize(maxN,maxM,false);
		dLeft.resize(cells,cap,false); dTotal.resize(maxN,maxM,false);
		dHit.resize(maxN,maxM,false); adjoint.resize(maxN,maxM,false);
	}
	auto reset=[]<class T>(Matrix<T> & matrix,size_t rows,size_t columns) {
		// Matrix storage is row-major. Clear contiguous full-width regions in
		// one pass; narrower rectangles leave the inactive row tails untouched.
		if(columns==matrix.size2()) std::fill_n(&matrix(0,0),rows*columns,T{});
		else for(size_t i=0;i<rows;++i) std::fill_n(&matrix(i,0),columns,T{});
	};

	// Flat linked anchor lists and a dense endpoint index avoid allocating a
	// vector and hash nodes for every heuristic group at every left boundary.
	// Zero denotes an empty list or an endpoint without a restricted group.
	struct Group { Pair end; bool all; size_t first=0; Pair begin{0,0}; };
	struct Anchor { Pair pair; size_t next; };
	std::vector<Group> groups;
	std::vector<Anchor> anchors;
	groups.reserve(w.selectRight?cells+1:1); anchors.reserve(cells);
	auto addAnchor=[&](Group & group,Pair q) {
		if(!group.first) group.begin=q;
		anchors.push_back({q,group.first}); group.first=anchors.size();
		group.begin.first=std::min(group.begin.first,q.first);
		group.begin.second=std::min(group.begin.second,q.second);
	};
	auto forAnchors=[&](const Group & group,const auto & visit) {
		for(size_t a=group.first;a;a=anchors[a-1].next) visit(anchors[a-1].pair);
	};
	seeds.forEachStart([&](size_t p1,size_t p2) {
		const Pair p{p1,p2}; if (!w.valid(p)) return;
		const size_t n=std::min(d.span1,d.n-p1),m=std::min(d.span2,d.m-p2);
		reset(left,n*m,cap); reset(leftTotal,n,m); reset(hit,n,m);
		reset(hybrid,n,m);
		if(w.selectRight) reset(byEnd,n,m);
		if(pairs) { reset(boundary,n,m); reset(dHit,n,m); }
		for(size_t i=0;i<n;++i) for(size_t j=0;j<m;++j) {
			valid(i,j)=w.valid({p1+i,p2+j});
			ending(i,j)=seeds.seedEndingLength(p1+i,p2+j);
		}
		auto accepts=[&](size_t i,size_t j,size_t len) { return ending(i,j) && ending(i,j)<=len; };
		auto nonstack=[&](size_t i,size_t j,const auto & visit) {
			for(size_t u=1;u<=i && u-1<=d.loop1;++u) for(size_t v=1;v<=j && v-1<=d.loop2;++v)
				if((u!=1 || v!=1) && valid(i-u,j-v)) visit(i-u,j-v);
		};
		auto transition=[&](size_t i,size_t j,size_t u,size_t v) {
			return A::check(w.transition({p1+u,p2+v},{p1+i,p2+j}));
		};
		auto extendLeft=[&](size_t i,size_t j,size_t len,Z_type weight) {
			// An accepting edge leaves the seed-free recurrence permanently. The
			// prefix already owns the seed's preceding pairs and their energies.
			auto & value=accepts(i,j,len)?hit(i,j)[d.noLP && len==1?0:1]:left(i*m+j,len-1);
			value=A::add(value,weight);
		};
		// N: nonaccepting prefixes. F: first-hit injections, including singleton
		// source hits. Initiation is owned by the source exactly once.
		extendLeft(0,0,1,A::check(w.initiation(p)));
		for(size_t i=0;i<n;++i) for(size_t j=0;j<m;++j) {
			if(!valid(i,j)) continue;
			if(i && j && valid(i-1,j-1)) {
				const Z_type t=transition(i,j,i-1,j-1);
				for(size_t c=0;c<cap;++c) {
					const Z_type prefix=left((i-1)*m+j-1,c);
					if(prefix>0) extendLeft(i,j,std::min(cap,c+2),A::multiply(prefix,t));
				}
			}
			nonstack(i,j,[&](size_t u,size_t v) {
				if(leftTotal(u,v)>0) extendLeft(i,j,1,A::multiply(leftTotal(u,v),transition(i,j,u,v)));
			});
			// A nonstack edge may leave only a supported end under noLP.
			for(size_t c=d.noLP?1:0;c<cap;++c) leftTotal(i,j)=A::add(leftTotal(i,j),left(i*m+j,c));
		}

		// Retained right extensions share first-hit injections by linearity:
		// T[p,r] = sum_q F[p,q] B[q,r]. Injections share work only when they
		// have the same retained endpoint (or retain all endpoints).
		groups.clear(); anchors.clear();
		groups.push_back({{n-1,m-1},true});
		for(size_t i=0;i<n;++i) for(size_t j=0;j<m;++j) {
			if(hit(i,j)[0]==0 && hit(i,j)[1]==0) continue;
			const size_t length=ending(i,j);
			if(!w.selectRight || (i+1==length && j+1==length)) {
				addAnchor(groups[0],{i,j});
				continue;
			}
			// No right extension is always retained, even if the chosen endpoint
			// is beyond this left boundary's inclusive span limits.
			hybrid(i,j)=A::add(hybrid(i,j),hit(i,j)[1]);
			const Pair end=selected.at({p1+i,p2+j});
			if(end.first<=p1+i || end.second<=p2+j || end.first-p1>=n || end.second-p2>=m) continue;
			const Pair local{end.first-p1,end.second-p2};
			auto & groupIndex=byEnd(local.first,local.second);
			if(!groupIndex) {
				groupIndex=groups.size(); groups.push_back({local,false});
			}
			addAnchor(groups[groupIndex],{i,j});
		}
		// A group has no paths above or left of all its injected anchors. Keep
		// both recurrences inside that rectangle, including predecessor lookups,
		// so reused storage outside it cannot leak into a later group.
		auto resetGroup=[&](Matrix<Ends> & matrix,const Group & group) {
			for(size_t i=group.begin.first;i<=group.end.first;++i)
				std::fill_n(&matrix(i,group.begin.second),group.end.second-group.begin.second+1,Ends{});
		};
		auto nonstackRight=[&](const Group & group,size_t i,size_t j,const auto & visit) {
			for(size_t u=1;u<=i-group.begin.first && u-1<=d.loop1;++u)
			for(size_t v=1;v<=j-group.begin.second && v-1<=d.loop2;++v)
				if((u!=1 || v!=1) && valid(i-u,j-v)) visit(i-u,j-v);
		};
		auto fillRight=[&](const Group & group) {
			resetGroup(right,group);
			forAnchors(group,[&](Pair q) { right(q.first,q.second)=hit(q.first,q.second); });
			// B is seed-unconstrained: no seed lookup or subtraction occurs here.
			for(size_t i=group.begin.first;i<=group.end.first;++i)
			for(size_t j=group.begin.second;j<=group.end.second;++j) {
				if(!valid(i,j)) continue;
				if(i>group.begin.first && j>group.begin.second && valid(i-1,j-1)) {
					const auto & prev=right(i-1,j-1);
					if(prev[0]>0 || prev[1]>0)
						right(i,j)[1]=A::add(right(i,j)[1],A::multiply(A::add(prev[0],prev[1]),transition(i,j,i-1,j-1)));
				}
				nonstackRight(group,i,j,[&](size_t u,size_t v) {
					if(right(u,v)[1]>0) {
						auto & value=right(i,j)[d.noLP?0:1];
						value=A::add(value,A::multiply(right(u,v)[1],transition(i,j,u,v)));
					}
				});
			}
		};
		for(const auto & group:groups) {
			if(!group.first) continue;
			fillRight(group);
			if(group.all) {
				for(size_t i=group.begin.first;i<n;++i) for(size_t j=group.begin.second;j<m;++j)
					hybrid(i,j)=A::add(hybrid(i,j),right(i,j)[1]);
			} else {
				auto & value=hybrid(group.end.first,group.end.second);
				value=A::add(value,right(group.end.first,group.end.second)[1]);
			}
		}
		// Join full boundaries before accessibility, reporting and differentiation.
		// Each callback receives the complete sum over canonical seed owners once.
		for(size_t i=0;i<n;++i) for(size_t j=0;j<m;++j) if(hybrid(i,j)>0) {
			const Pair q{p1+i,p2+j};
			const Z_type b=A::check(w.boundary(p,q));
			if(pairs) boundary(i,j)=b;
			out.z=A::add(out.z,A::multiply(hybrid(i,j),b));
			if(w.complete) w.complete(p,q,hybrid(i,j),b);
		}
		if(!pairs) return;

		// Reverse the retained right contractions. Recompute one group at a time
		// to bound memory independently of the number of heuristic choices.
		for(size_t g=0;g<groups.size();++g) {
			const auto & group=groups[g]; if(!group.first) continue;
			if(groups.size()!=1) fillRight(group);
			resetGroup(adjoint,group);
			for(size_t ii=group.end.first+1;ii>group.begin.first;--ii)
			for(size_t jj=group.end.second+1;jj>group.begin.second;--jj) {
				const size_t i=ii-1,j=jj-1; if(!valid(i,j)) continue;
				if(right(i,j)[0]==0 && right(i,j)[1]==0) continue;
				auto & a=adjoint(i,j);
				if(group.all || (i==group.end.first && j==group.end.second)) a[1]=A::add(a[1],boundary(i,j));
				// Nonnegative forward zeros have no contributing ancestors. Only
				// differentiate states reached by an injection in this group.
				for(size_t c=0;c<2;++c) if(right(i,j)[c]==0) a[c]=0;
				// A vertex is owned by its most recent pair: T*dT counts every
				// right-side pair, including the accepting pair, exactly once.
				auto & mass=out.mass(p1+i,p2+j);
				for(size_t c=0;c<2;++c) mass=A::add(mass,A::multiply(right(i,j)[c],a[c]));
				if(a[1]>0 && i>group.begin.first && j>group.begin.second && valid(i-1,j-1)
						&& (right(i-1,j-1)[0]>0 || right(i-1,j-1)[1]>0)) {
					const Z_type value=A::multiply(a[1],transition(i,j,i-1,j-1));
					for(size_t c=0;c<2;++c) if(right(i-1,j-1)[c]>0)
						adjoint(i-1,j-1)[c]=A::add(adjoint(i-1,j-1)[c],value);
				}
				if(a[d.noLP?0:1]>0) nonstackRight(group,i,j,[&](size_t u,size_t v) {
					if(right(u,v)[1]>0)
						adjoint(u,v)[1]=A::add(adjoint(u,v)[1],A::multiply(a[d.noLP?0:1],transition(i,j,u,v)));
				});
			}
			forAnchors(group,[&](Pair q) {
				for(size_t c=0;c<2;++c)
					dHit(q.first,q.second)[c]=A::add(dHit(q.first,q.second)[c],adjoint(q.first,q.second)[c]);
			});
		}
		if(w.selectRight) for(size_t i=0;i<n;++i) for(size_t j=0;j<m;++j) {
			if(hit(i,j)[0]==0 && hit(i,j)[1]==0) continue;
			const size_t length=ending(i,j);
			if(i+1==length && j+1==length) continue;
			// The no-right branch was not injected into an unrestricted group.
			dHit(i,j)[1]=A::add(dHit(i,j)[1],boundary(i,j));
			out.mass(p1+i,p2+j)=A::add(out.mass(p1+i,p2+j),A::multiply(hit(i,j)[1],boundary(i,j)));
		}

		// Reverse F accepting edges and N seed-free edges. F itself is not marked
		// again: its accepting pair was already counted in the right contraction.
		reset(dLeft,n*m,cap); reset(dTotal,n,m);
		auto destinationAdjoint=[&](size_t i,size_t j,size_t len) {
			return accepts(i,j,len)?dHit(i,j)[d.noLP && len==1?0:1]:dLeft(i*m+j,len-1);
		};
		for(size_t ii=n;ii>0;--ii) for(size_t jj=m;jj>0;--jj) {
			const size_t i=ii-1,j=jj-1,q=i*m+j; if(!valid(i,j)) continue;
			for(size_t c=d.noLP?1:0;c<cap;++c) dLeft(q,c)=A::add(dLeft(q,c),dTotal(i,j));
			auto & mass=out.mass(p1+i,p2+j);
			for(size_t c=0;c<cap;++c) mass=A::add(mass,A::multiply(left(q,c),dLeft(q,c)));
			if(i && j && valid(i-1,j-1)) {
				const Z_type t=transition(i,j,i-1,j-1);
				for(size_t c=0;c<cap;++c) {
					auto & value=dLeft(q-m-1,c);
					value=A::add(value,A::multiply(t,destinationAdjoint(i,j,std::min(cap,c+2))));
				}
			}
			nonstack(i,j,[&](size_t u,size_t v) {
				dTotal(u,v)=A::add(dTotal(u,v),A::multiply(transition(i,j,u,v),destinationAdjoint(i,j,1)));
			});
		}
	});
	return out;
}
}
