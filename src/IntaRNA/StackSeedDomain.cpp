#include "IntaRNA/StackSeedDomain.h"
#include <algorithm>
#include <stdexcept>

namespace IntaRNA {
StackSeedDomain::StackSeedDomain(const SeedHandler & handler, size_t n, size_t m, size_t w1, size_t w2) {
	if (!handler.guaranteesStackOnlySeeds()) throw std::invalid_argument("stack seed domain requires a stack-only handler");
	if (n==0 || m==0) return;
	// Use the range end as an invalid initial cursor, avoiding lastPos+offset wrap.
	size_t i=n, j=m;
	while (handler.updateToNextSeed(i,j,0,n-1,0,m-1)) {
		const size_t len=handler.getSeedLength1(i,j);
		if (len != handler.getSeedLength2(i,j)) throw std::logic_error("stack seed has unequal strand lengths");
		add({i,j,len},n,m,w1,w2);
	}
}
StackSeedDomain::StackSeedDomain(const std::vector<Occurrence> & seeds, size_t n, size_t m, size_t w1, size_t w2) {
	for (auto s:seeds) add(s,n,m,w1,w2);
}
void StackSeedDomain::add(Occurrence s, size_t n, size_t m, size_t w1, size_t w2) {
	w1=std::min(w1,n); w2=std::min(w2,m);
	if (s.length==0 || s.i>=n || s.j>=m || s.length>n-s.i || s.length>m-s.j
			|| s.length>w1 || s.length>w2) return;
	auto [p,inserted]=lengths.emplace(std::make_pair(s.i,s.j),s.length);
	if (!inserted && p->second!=s.length) throw std::invalid_argument("multiple retained seed lengths at one start");
	maximum=std::max(maximum,s.length);
	const size_t end1=s.i+s.length, end2=s.j+s.length;
	boxes.push_back({end1>w1?end1-w1:0,s.i,end2>w2?end2-w2:0,s.j});
}
void StackSeedDomain::forEachStart(const std::function<void(size_t,size_t)> & visit) const {
	if (boxes.empty()) return;
	size_t first=boxes.front().i0, last=0;
	for (auto b:boxes) { first=std::min(first,b.i0); last=std::max(last,b.i1); }
	std::vector<std::pair<size_t,size_t>> intervals;
	for (size_t i=first;i<=last;++i) {
		intervals.clear();
		for (auto b:boxes) if (b.i0<=i && i<=b.i1) intervals.emplace_back(b.j0,b.j1);
		std::sort(intervals.begin(),intervals.end());
		size_t next=0;
		for (auto [lo,hi]:intervals) {
			for (size_t j=std::max(lo,next);j<=hi;++j) visit(i,j);
			next=std::max(next,hi+1);
		}
	}
}
}
