#include "IntaRNA/SeededPartitionFunction.h"
#include "IntaRNA/PartitionArithmetic.h"

// Seed-free suffix lengths plus two absorbing seeded states (singleton/run).
// Nonstack edges operate on aggregate totals, independently of seed length.
namespace IntaRNA {
SeededPartitionFunction::Result SeededPartitionFunction::compute(
		const IntaRNA::SeededPartitionFunction::Domain & d,
		const IntaRNA::StackSeedDomain & seeds,
		const IntaRNA::SeededPartitionFunction::Weights & w, bool pairs)
{
	using namespace IntaRNA;
	using A=PartitionArithmetic;
	using Pair=SeededPartitionFunction::Pair;
	SeededPartitionFunction::Result out;
	if (pairs) out.mass=Matrix<Z_type>(d.n,d.m,0);
	if (!seeds.maxSeedLength() || !d.span1 || !d.span2) return out;
	// At least two length states are necessary for noLP even with singleton seeds.
	const size_t cap=std::max(size_t(2),seeds.maxSeedLength()), s1=cap, s2=cap+1, states=cap+2;
	seeds.forEachStart([&](size_t p1,size_t p2) {
		const Pair p{p1,p2}; if (!w.valid(p)) return;
		const size_t n=std::min(d.span1,d.n-p1),m=std::min(d.span2,d.m-p2);
		Matrix<Z_type> h(n*m,states),dh;
		Matrix<std::array<Z_type,2>> total(n,m),dt;
		Matrix<unsigned char> valid(n,m);
		Matrix<size_t> ending(n,m);
		if (pairs) { dh.resize(n*m,states); dt.resize(n,m); }
		for (size_t i=0;i<n;++i) for(size_t j=0;j<m;++j) {
			valid(i,j)=w.valid({p1+i,p2+j});
			ending(i,j)=seeds.seedEndingLength(p1+i,p2+j);
		}
		auto accepts=[&](size_t i,size_t j,size_t len) { auto e=ending(i,j);return e && e<=len; };
		auto dest=[&](size_t i,size_t j,size_t c,bool stacked) {
			if (c>=cap) return stacked || !d.noLP?s2:s1;
			const size_t len=stacked?std::min(cap,c+2):1;
			return accepts(i,j,len)?(len>=2 || !d.noLP?s2:s1):len-1;
		};
		auto nonstack=[&](size_t i,size_t j,const auto & visit) {
			for(size_t u=1;u<=i && u-1<=d.loop1;++u) for(size_t v=1;v<=j && v-1<=d.loop2;++v)
				if((u!=1 || v!=1) && valid(i-u,j-v)) visit(i-u,j-v);
		};
		const size_t initial=accepts(0,0,1)?(d.noLP?s1:s2):0;
		const Z_type init=A::check(w.initiation(p)); h(0,initial)=init;
		for(size_t i=0;i<n;++i) for(size_t j=0;j<m;++j) {
			if(!valid(i,j)) continue;
			const size_t q=i*m+j;
			if(i && j && valid(i-1,j-1)) {
				const Z_type t=A::check(w.transition({p1+i-1,p2+j-1},{p1+i,p2+j}));
				for(size_t c=0;c<states;++c) {
					auto & v=h(q,dest(i,j,c,true));v=A::add(v,A::multiply(h(q-m-1,c),t));
				}
			}
			nonstack(i,j,[&](size_t u,size_t v) {
				if(total(u,v)[0]==0 && total(u,v)[1]==0) return;
				const Z_type t=A::check(w.transition({p1+u,p2+v},{p1+i,p2+j}));
				for(size_t c=0;c<2;++c) {
					auto & value=h(q,dest(i,j,c? s2:0,false));value=A::add(value,A::multiply(total(u,v)[c],t));
				}
			});
			for(size_t c=d.noLP?1:0;c<cap;++c) total(i,j)[0]=A::add(total(i,j)[0],h(q,c));
			total(i,j)[1]=A::add(h(q,s2),d.noLP?Z_type(0):h(q,s1));
			if(total(i,j)[1]>0) {
				const Z_type b=A::check(w.boundary(p,{p1+i,p2+j}));
				out.z=A::add(out.z,A::multiply(b,total(i,j)[1]));
				if(pairs) { dh(q,s2)=b;if(!d.noLP) dh(q,s1)=b; }
				if(w.complete) w.complete(p,{p1+i,p2+j},total(i,j)[1],b);
			}
		}
		if(!pairs) return;
		for(size_t ii=n;ii>0;--ii) for(size_t jj=m;jj>0;--jj) {
			const size_t i=ii-1,j=jj-1,q=i*m+j; if(!valid(i,j)) continue;
			for(size_t c=d.noLP?1:0;c<cap;++c) dh(q,c)=A::add(dh(q,c),dt(i,j)[0]);
			dh(q,s2)=A::add(dh(q,s2),dt(i,j)[1]);
			if(!d.noLP) dh(q,s1)=A::add(dh(q,s1),dt(i,j)[1]);
			auto & mass=out.mass(p1+i,p2+j);
			if(i && j && valid(i-1,j-1)) {
				const Z_type t=A::check(w.transition({p1+i-1,p2+j-1},{p1+i,p2+j}));
				for(size_t c=0;c<states;++c) {
					const Z_type a=A::multiply(dh(q,dest(i,j,c,true)),t);
					dh(q-m-1,c)=A::add(dh(q-m-1,c),a);
					mass=A::add(mass,A::multiply(a,h(q-m-1,c)));
				}
			}
			nonstack(i,j,[&](size_t u,size_t v) {
				const Z_type t=A::check(w.transition({p1+u,p2+v},{p1+i,p2+j}));
				for(size_t c=0;c<2;++c) {
					const Z_type a=A::multiply(dh(q,dest(i,j,c?s2:0,false)),t);
					dt(u,v)[c]=A::add(dt(u,v)[c],a);
					mass=A::add(mass,A::multiply(a,total(u,v)[c]));
				}
			});
		}
		out.mass(p1,p2)=A::add(out.mass(p1,p2),A::multiply(dh(0,initial),init));
	});
	return out;
}
}
