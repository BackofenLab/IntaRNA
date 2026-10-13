#include "IntaRNA/StackSeedDomain.h"
#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <utility>

namespace IntaRNA {

StackSeedDomain::
StackSeedDomain(const SeedHandler & handler, size_t n, size_t m, size_t w1, size_t w2) {
	// Reject unsupported input before interpreting seed lengths as diagonal runs.
#if INTARNA_IN_DEBUG_MODE
	if (!handler.guaranteesStackOnlySeeds()) throw std::invalid_argument("stack seed domain requires a stack-only handler");
#endif
	if (n==0 || m==0) return;

	// Use the range end as an invalid initial cursor, avoiding lastPos+offset wrap.
	size_t i=n, j=m;
	while (handler.updateToNextSeed(i,j,0,n-1,0,m-1)) {
		// Stack-only admission guarantees equal strand lengths internally.
		const size_t len=handler.getSeedLength1(i,j);
		assert(len == handler.getSeedLength2(i,j));
		// Retain only occurrences fully contained in the active interaction domain.
		add({i,j,len},n,m,w1,w2);
	}
}


StackSeedDomain::
StackSeedDomain(const std::vector<Occurrence> & seeds, size_t n, size_t m, size_t w1, size_t w2) {
	// Apply the same containment and span filters to externally supplied seeds.
	for (auto s:seeds) add(s,n,m,w1,w2);
}


void 
StackSeedDomain::
add(Occurrence s, size_t n, size_t m, size_t w1, size_t w2) {
	// Check containment using remaining lengths so endpoint arithmetic cannot wrap.
	w1=std::min(w1,n); w2=std::min(w2,m);
	if (s.length==0 || s.i>=n || s.j>=m || s.length>n-s.i || s.length>m-s.j
			|| s.length>w1 || s.length>w2) return;

	// A start must identify one retained length, including for explicit input.
	auto [p,inserted]=lengths.emplace(Interaction::BasePair(s.i,s.j),s.length);
#if INTARNA_IN_DEBUG_MODE
	if (!inserted && p->second!=s.length) throw std::invalid_argument("multiple retained seed lengths at one start");
#endif

	// The shortest seed ending here suffices to accept a stacked path suffix.
	auto [end, fresh]=endingLengths.emplace(Interaction::BasePair(s.i+s.length-1,s.j+s.length-1),s.length);
	if (!fresh) end->second=std::min(end->second,s.length);
	maximum=std::max(maximum,s.length);

	// A left boundary can precede this seed only while its end fits both spans.
	const size_t end1=s.i+s.length, end2=s.j+s.length;
	boxes.push_back({end1>w1?end1-w1:0,s.i,end2>w2?end2-w2:0,s.j});
}

void 
StackSeedDomain::
forEachStart(const std::function<void(size_t,size_t)> & visit) const {
	if (boxes.empty()) return;
	// Restrict the row sweep to rows covered by at least one eligible-start box.
	size_t first=boxes.front().i0, last=0;
	for (auto b:boxes) { first=std::min(first,b.i0); last=std::max(last,b.i1); }
	// These pairs are inclusive intervals on strand two, not base-pair coordinates.
	std::vector<std::pair<size_t,size_t>> intervals;
	for (size_t i=first;i<=last;++i) {
		// Sort this row's column intervals independently of seed insertion order.
		intervals.clear();
		for (auto b:boxes) if (b.i0<=i && i<=b.i1) intervals.emplace_back(b.j0,b.j1);
		std::sort(intervals.begin(),intervals.end());
		// Skip already visited columns where boxes overlap; emit row-major order.
		size_t next=0;
		for (auto [lo,hi]:intervals) {
			for (size_t j=std::max(lo,next);j<=hi;++j) visit(i,j);
			next=std::max(next,hi+1);
		}
	}
}
}
