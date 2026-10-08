#ifndef INTARNA_STACKSEEDDOMAIN_H_
#define INTARNA_STACKSEEDDOMAIN_H_

#include "IntaRNA/SeedHandler.h"
#include <functional>
#include <map>
#include <utility>
#include <vector>

namespace IntaRNA {
/** Admitted stack occurrences and the deduplicated union of eligible-start boxes.
 * Coordinates are local, zero-based and increasing on both strands. Seed energies
 * have already been filtered by the handler and are never reused as path weights.
 */
class StackSeedDomain {
public:
	/** A consecutive diagonal run, including singleton explicit seeds. */
	struct Occurrence { size_t i, j, length; };
	/** Compile a filled handler, requiring both endpoints inside the active range.
	 * @param handler filled stack-only handler, including any index offset
	 * @param n first range length
	 * @param m second range length
	 * @param span1 inclusive first-strand interaction span limit
	 * @param span2 inclusive second-strand interaction span limit
	 */
	StackSeedDomain(const SeedHandler & handler, size_t n, size_t m, size_t span1, size_t span2);
	/** Compile already admitted occurrences (also useful for independent weight models).
	 * @param occurrences retained occurrences; out-of-range/overlong ones are excluded
	 * @param n first range length
	 * @param m second range length
	 * @param span1 inclusive first span limit
	 * @param span2 inclusive second span limit
	 */
	StackSeedDomain(const std::vector<Occurrence> & occurrences, size_t n, size_t m, size_t span1, size_t span2);
	/** @return retained length at this start, or zero if no occurrence starts here */
	size_t seedLength(size_t i, size_t j) const;
	/** @return maximum retained occurrence length, zero for an empty family */
	size_t maxSeedLength() const;
	/** Enumerate each eligible left boundary once, in row-major order.
	 * The caller intersects these boxes with accessible complementary vertices.
	 * @param visit callback receiving local first/second coordinates
	 */
	void forEachStart(const std::function<void(size_t,size_t)> & visit) const;
private:
	//! admitted occurrences keyed by local start
	std::map<std::pair<size_t,size_t>, size_t> lengths;
	//! eligible inclusive rectangles [i0,i1] x [j0,j1]
	struct Box { size_t i0, i1, j0, j1; };
	std::vector<Box> boxes;
	size_t maximum = 0;
	void add(Occurrence seed, size_t n, size_t m, size_t span1, size_t span2);
};
inline size_t StackSeedDomain::seedLength(size_t i, size_t j) const {
	auto p=lengths.find({i,j}); return p==lengths.end()?0:p->second;
}
inline size_t StackSeedDomain::maxSeedLength() const { return maximum; }
}
#endif
