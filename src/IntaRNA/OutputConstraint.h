
#ifndef OUTPUTCONSTRAINT_H_
#define OUTPUTCONSTRAINT_H_

#include "IntaRNA/general.h"

namespace IntaRNA {

/**
 *
 * Data structure that contains all constraints to be applied to (suboptimal)
 * output generation.
 *
 * @author Martin Mann
 *
 */
class OutputConstraint
{

public:

	//! different possibilities to en-/disable overlapping of interaction sites
	//! if suboptimal solutions are enumerated
	enum ReportOverlap {
		OVERLAP_NONE = 0,
		OVERLAP_SEQ1 = 1,
		OVERLAP_SEQ2 = 2,
		OVERLAP_BOTH = 3
	};


public:

	//! the maximal number of (sub)optimal interactions to be reported to the output handler
	const size_t reportMax;

	//! defines whether and where overlapping interaction sites are allowed for reporting
	const ReportOverlap reportOverlap;

	//! the maximal energy of a reported interaction (<= 0.0)
	const E_type maxE;

	//! the maximal energy difference to the mfe of a reported interaction
	const E_type deltaE;

public:

	/**
	 * Construction of an output constraint
	 *
	 * @param reportMax the maximal number of (sub)optimal interactions to be
	 *            reported to the output handler
	 * @param reportOverlap defines whether and where overlapping interaction
	 *            sites are allowed for reporting
	 * @param maxE maximal energy of a reported interaction (<= 0.0)
	 * @param deltaE maximal energy difference of a reported interaction to mfe
	 */
	OutputConstraint(	  const size_t reportMax = 1
						, const ReportOverlap reportOverlap = OVERLAP_BOTH
						, const E_type maxE = 0.0
						, const E_type deltaE = E_INF );

	//! destruction
	virtual ~OutputConstraint();
};

} // namespace

#endif /* OUTPUTCONSTRAINT_H_ */
