
#ifndef INTARNA_INTERACTIONRANGE_H_
#define INTARNA_INTERACTIONRANGE_H_

#include "IntaRNA/general.h"
#include "IntaRNA/IndexRange.h"
#include "IntaRNA/RnaSequence.h"

namespace IntaRNA {

// dummy declaration to avoid inclusion loop
class Interaction;

/**
 * Represents two ranges that form an interaction
 */
class InteractionRange {
public:

	//! the first interaction partner
	const RnaSequence * s1;

	//! the second interaction partner
	const RnaSequence * s2;

	//! the index range of s1 interacting with r2
	IndexRange r1;

	//! the index range of s2 interacting with r1
	IndexRange r2;

	//! the energy value associated with this interaction range
	E_type energy;

public:

	/**
	 * Construction of an interaction range given two index ranges
	 * @param rna1 the first RNA interacting
	 * @param rna2 the second RNA interacting
	 * @param r1 the range of rna1 interacting
	 * @param r2 the range of rna2 interacting
	 * @param energy the energy value assigned to the range
	 */
	InteractionRange(
				const RnaSequence & rna1
				, const RnaSequence & rna2
				, const IndexRange & r1 = IndexRange(0,RnaSequence::lastPos)
				, const IndexRange & r2 = IndexRange(RnaSequence::lastPos,0)
				, const E_type energy = E_INF )
		: s1(&rna1)
		, s2(&rna2)
		, r1(r1.from, (r1.to > rna1.size()?rna1.size()-1:r1.to))
		, r2((r2.from > rna2.size()?rna2.size()-1:r2.from), r2.to)
		, energy(energy)
	{
#if INTARNA_IN_DEBUG_MODE
		if (!isSane())
			throw std::runtime_error("InteractionRange("+toString(*this)+") not sane!");
#endif
	}

	/**
	 * Construction of an interaction range given two index ranges
	 */
	InteractionRange( const Interaction& interaction )
		: s1(NULL), s2(NULL), r1(), r2(), energy(E_INF)
	{
		// copy data
		this->operator =(interaction);
	}

	/**
	 * Destruction
	 */
	virtual ~InteractionRange()
	{}


	/**
	 * Checks whether or not the range encoding is reasonable.
	 * @return from <= to
	 */
	bool isSane() const
	{
		// range of first sequence ascending
		// range of second sequence descending
		return r1.isAscending() && r2.isDescending();
	}

	/**
	 * Prints the range's boundaries to stream
	 * @param out the ostream to write to
	 * @param r the InteractionRange object to add
	 * @return the altered stream out
	 */
	friend std::ostream& operator<<(std::ostream& out, const InteractionRange& r)
	{
		return (out <<"("<<r.r1<<"#"<<r.r2<<")");
	}

	/**
	 * Defines a less-than order on interaction ranges.
	 * @param r the range to compare to
	 * @return (r1 < r.r1) || (!(r1>r.r1) && r2<r.r2)
	 */
	const bool operator < ( const InteractionRange &r ) const
	{
		return ( r1 < r.r1 || (!(r.r1<r1) && r2<r.r2) );
	}

	/**
	 * Copies the range of the given interaction into this.
	 * @param interaction the interaction to get the range from
	 * @return the altered range object (*this)
	 */
	InteractionRange & operator= ( const Interaction & interaction );

};

} // namespace

#endif /* INTERACTIONRANGE_H_ */
