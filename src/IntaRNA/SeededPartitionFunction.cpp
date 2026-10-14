#include "IntaRNA/SeededPartitionFunction.h"
#include "IntaRNA/PartitionArithmetic.h"
#include <algorithm>
#include <array>
#include <unordered_map>
#include <vector>

namespace IntaRNA {
SeededPartitionFunction::Result
SeededPartitionFunction::
compute(const Domain & d, const StackSeedDomain & seeds, const Weights & w, bool pairs)
{
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
	// Seed-free prefixes retain stack length; right extensions need only noLP.
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
	struct Group { Pair end; bool all; size_t first=0; };
	struct Anchor { Pair pair; size_t next; };
	std::vector<Group> groups;
	std::vector<Anchor> anchors;
	groups.reserve(w.selectRight?cells+1:1); anchors.reserve(cells);
	auto addAnchor=[&](Group & group,Pair q) {
		anchors.push_back({q,group.first}); group.first=anchors.size();
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

		// Exact right extensions share all first-hit injections by linearity:
		// T[p,r] = sum_q F[p,q] B[q,r]. Heuristic injections share work only
		// when they have the same retained endpoint (or all endpoints).
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
		auto fillRight=[&](const Group & group) {
			reset(right,group.end.first+1,group.end.second+1);
			forAnchors(group,[&](Pair q) { right(q.first,q.second)=hit(q.first,q.second); });
			// B is seed-unconstrained: no seed lookup or subtraction occurs here.
			for(size_t i=0;i<=group.end.first;++i) for(size_t j=0;j<=group.end.second;++j) {
				if(!valid(i,j)) continue;
				if(i && j && valid(i-1,j-1)) {
					const auto & prev=right(i-1,j-1);
					if(prev[0]>0 || prev[1]>0)
						right(i,j)[1]=A::add(right(i,j)[1],A::multiply(A::add(prev[0],prev[1]),transition(i,j,i-1,j-1)));
				}
				nonstack(i,j,[&](size_t u,size_t v) {
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
				for(size_t i=0;i<n;++i) for(size_t j=0;j<m;++j)
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
			reset(adjoint,group.end.first+1,group.end.second+1);
			for(size_t ii=group.end.first+1;ii>0;--ii) for(size_t jj=group.end.second+1;jj>0;--jj) {
				const size_t i=ii-1,j=jj-1; if(!valid(i,j)) continue;
				auto & a=adjoint(i,j);
				if(group.all || (i==group.end.first && j==group.end.second)) a[1]=A::add(a[1],boundary(i,j));
				// A vertex is owned by its most recent pair: T*dT counts every
				// right-side pair, including the accepting pair, exactly once.
				auto & mass=out.mass(p1+i,p2+j);
				for(size_t c=0;c<2;++c) mass=A::add(mass,A::multiply(right(i,j)[c],a[c]));
				if(i && j && valid(i-1,j-1)) {
					const Z_type value=A::multiply(a[1],transition(i,j,i-1,j-1));
					for(size_t c=0;c<2;++c) adjoint(i-1,j-1)[c]=A::add(adjoint(i-1,j-1)[c],value);
				}
				nonstack(i,j,[&](size_t u,size_t v) {
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
