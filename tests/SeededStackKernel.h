#ifndef INTARNA_TEST_SEEDEDSTACKKERNEL_H_
#define INTARNA_TEST_SEEDEDSTACKKERNEL_H_
// Maximal-stack candidate retained only for oracle and timing comparisons.
#include "IntaRNA/SeededPartitionFunction.h"
#include "IntaRNA/PartitionArithmetic.h"
#include <algorithm>
#include <vector>

inline IntaRNA::SeededPartitionFunction::Result seededStack(const IntaRNA::SeededPartitionFunction::Domain & d, const IntaRNA::StackSeedDomain & seeds,
		const IntaRNA::SeededPartitionFunction::Weights & w, bool pairs)
{
	using namespace IntaRNA;
	using Pair=SeededPartitionFunction::Pair;
	using Result=SeededPartitionFunction::Result;
	using A=PartitionArithmetic;
	using State=std::array<Z_type,2>;
	Result result;
	if (pairs) result.mass=Matrix<Z_type>(d.n,d.m,0);
	if (!seeds.maxSeedLength() || !d.span1 || !d.span2) return result;
	seeds.forEachStart([&](size_t p1,size_t p2) {
		const Pair p{p1,p2};
		if (!w.valid(p)) return;
		const size_t n=std::min(d.span1,d.n-p1), m=std::min(d.span2,d.m-p2);
		Matrix<State> h(n,m), entry(n,m), dh, de;
		Matrix<unsigned char> valid(n,m);
		if (pairs) { dh.resize(n,m); de.resize(n,m); }
		for (size_t i=0;i<n;++i) for (size_t j=0;j<m;++j) valid(i,j)=w.valid({p1+i,p2+j});
		// Iterate incoming NONSTACK edges only. Completed maximal runs uniquely
		// delimit these entries; cutting a stack here would duplicate paths.
		auto incoming=[&](size_t i,size_t j,const auto & visit) {
			for (size_t u=1;u<=i && u-1<=d.loop1;++u)
			for (size_t v=1;v<=j && v-1<=d.loop2;++v) {
				if ((u==1 && v==1) || !valid(i-u,j-v)) continue;
				visit(i-u,j-v);
			}
		};
		for (size_t i=0;i<n;++i) for (size_t j=0;j<m;++j) {
			if (!valid(i,j)) continue;
			const Pair q{p1+i,p2+j};
			if (i==0 && j==0) entry(i,j)[0]=A::check(w.initiation(p));
			incoming(i,j,[&](size_t u,size_t v) {
				if (h(u,v)[0]==0 && h(u,v)[1]==0) return;
				const Z_type t=A::check(w.transition({p1+u,p2+v},q));
				for (size_t c=0;c<2;++c) entry(i,j)[c]=A::add(entry(i,j)[c],A::multiply(h(u,v)[c],t));
			});
			bool hit=false;
			Z_type product=1;
			for (size_t k=0;k<=std::min(i,j);++k) {
				const size_t a=i-k,b=j-k;
				if (!valid(a,b)) break;
				if (k) {
					const Z_type t=A::check(w.transition({p1+a,p2+b},{p1+a+1,p2+b+1}));
					if (t==0) break;
					product=A::multiply(product,t);
				}
				const size_t length=seeds.seedLength(p1+a,p2+b);
				hit |= length>0 && length<=k+1;
				if (d.noLP && k==0) continue;
				for (size_t c=0;c<2;++c) {
					const size_t accepted=c || hit;
					h(i,j)[accepted]=A::add(h(i,j)[accepted],A::multiply(entry(a,b)[c],product));
				}
			}
			if (h(i,j)[1]>0) {
				const Z_type b=A::check(w.boundary(p,q));
				result.z=A::add(result.z,A::multiply(b,h(i,j)[1]));
				if (pairs) dh(i,j)[1]=b;
				if (w.complete) w.complete(p,q,h(i,j)[1],b);
			}
		}
		if (!pairs) return;
		// Recompute each diagonal product chain with O(run length) scratch, then
		// reverse it once. Each entry marks the first pair, product edges the rest.
		std::vector<Z_type> product, edge, g;
		for (size_t ii=n;ii>0;--ii) for (size_t jj=m;jj>0;--jj) {
			const size_t i=ii-1,j=jj-1;
			if (!valid(i,j)) continue;
			product.clear(); edge.clear(); g.clear();
			bool hit=false;
			for (size_t k=0;k<=std::min(i,j);++k) {
				const size_t a=i-k,b=j-k;
				if (!valid(a,b)) break;
				Z_type t=1;
				if (k) { t=A::check(w.transition({p1+a,p2+b},{p1+a+1,p2+b+1})); if (t==0) break; }
				edge.push_back(t);
				product.push_back(k?A::multiply(product.back(),t):Z_type(1));
				g.push_back(0);
				const size_t length=seeds.seedLength(p1+a,p2+b);
				hit |= length>0 && length<=k+1;
				if (d.noLP && k==0) continue; // retain the empty product node
				for (size_t c=0;c<2;++c) {
					const Z_type adj=dh(i,j)[c || hit];
					de(a,b)[c]=A::add(de(a,b)[c],A::multiply(adj,product.back()));
					g.back()=A::add(g.back(),A::multiply(adj,entry(a,b)[c]));
				}
			}
			for (size_t k=product.size();k>1;) {
				--k;
				auto & mass=result.mass(p1+i-k+1,p2+j-k+1);
				mass=A::add(mass,A::multiply(g[k],product[k]));
				g[k-1]=A::add(g[k-1],A::multiply(g[k],edge[k]));
			}
			for (size_t c=0;c<2;++c) {
				auto & mass=result.mass(p1+i,p2+j);
				mass=A::add(mass,A::multiply(de(i,j)[c],entry(i,j)[c]));
			}
			incoming(i,j,[&](size_t u,size_t v) {
				if (de(i,j)[0]==0 && de(i,j)[1]==0) return;
				const Z_type t=A::check(w.transition({p1+u,p2+v},{p1+i,p2+j}));
				for (size_t c=0;c<2;++c) dh(u,v)[c]=A::add(dh(u,v)[c],A::multiply(de(i,j)[c],t));
			});
		}
	});
	return result;
}
#endif
