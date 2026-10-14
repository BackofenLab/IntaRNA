#include "IntaRNA/SeededPartitionFunction.h"
#include "IntaRNA/PartitionArithmetic.h"
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
	seeds.forEachStart([&](size_t p1,size_t p2) {
		const Pair p{p1,p2}; if (!w.valid(p)) return;
		const size_t n=std::min(d.span1,d.n-p1),m=std::min(d.span2,d.m-p2);
		Matrix<Z_type> left(matrix_detail::product(n,m),cap),leftTotal(n,m);
		Matrix<Ends> hit(n,m);
		Matrix<unsigned char> valid(n,m);
		Matrix<size_t> ending(n,m);
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
		struct Group { Pair end; bool all; std::vector<Pair> anchors; };
		std::vector<Group> groups;
		groups.push_back({{n-1,m-1},true,{}});
		std::unordered_map<Pair,size_t,Pair::Hash> byEnd;
		Matrix<Z_type> hybrid(n,m),boundary(n,m);
		Matrix<Ends> dHit;
		if(pairs) dHit.resize(n,m);
		for(size_t i=0;i<n;++i) for(size_t j=0;j<m;++j) {
			if(hit(i,j)[0]==0 && hit(i,j)[1]==0) continue;
			const size_t length=ending(i,j);
			if(!w.selectRight || (i+1==length && j+1==length)) {
				groups[0].anchors.emplace_back(i,j);
				continue;
			}
			// No right extension is always retained, even if the chosen endpoint
			// is beyond this left boundary's inclusive span limits.
			hybrid(i,j)=A::add(hybrid(i,j),hit(i,j)[1]);
			const Pair end=selected.at({p1+i,p2+j});
			if(end.first<=p1+i || end.second<=p2+j || end.first-p1>=n || end.second-p2>=m) continue;
			const Pair local{end.first-p1,end.second-p2};
			auto [at,fresh]=byEnd.emplace(local,groups.size());
			if(fresh) groups.push_back({local,false,{}});
			groups[at->second].anchors.emplace_back(i,j);
		}
		auto fillRight=[&](const Group & group) {
			Matrix<Ends> right(group.end.first+1,group.end.second+1);
			for(auto q:group.anchors) right(q.first,q.second)=hit(q.first,q.second);
			// B is seed-unconstrained: no seed lookup or subtraction occurs here.
			for(size_t i=0;i<right.size1();++i) for(size_t j=0;j<right.size2();++j) {
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
			return right;
		};
		Matrix<Ends> right;
		for(const auto & group:groups) {
			if(group.anchors.empty()) continue;
			right=fillRight(group);
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
			const Z_type b=A::check(w.boundary(p,q)); boundary(i,j)=b;
			out.z=A::add(out.z,A::multiply(hybrid(i,j),b));
			if(w.complete) w.complete(p,q,hybrid(i,j),b);
		}
		if(!pairs) return;

		// Reverse the retained right contractions. Recompute one group at a time
		// to bound memory independently of the number of heuristic choices.
		for(size_t g=0;g<groups.size();++g) {
			const auto & group=groups[g]; if(group.anchors.empty()) continue;
			if(groups.size()!=1) right=fillRight(group);
			Matrix<Ends> adjoint(right.size1(),right.size2());
			for(size_t ii=right.size1();ii>0;--ii) for(size_t jj=right.size2();jj>0;--jj) {
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
			for(auto q:group.anchors) for(size_t c=0;c<2;++c)
				dHit(q.first,q.second)[c]=A::add(dHit(q.first,q.second)[c],adjoint(q.first,q.second)[c]);
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
		Matrix<Z_type> dLeft(n*m,cap),dTotal(n,m);
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
