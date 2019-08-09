
#ifndef INTARNA_HELIXCONSTRAINT_H
#define INTARNA_HELIXCONSTRAINT_H

#include "IntaRNA/general.h"
#include "IntaRNA/IndexRangeList.h"

#include <cstddef>
#include <iostream>

namespace IntaRNA {

/**
 * Encodes helix constraints to be used for interaction prediction.
 *
 * @author Rick Gelhausen
 *
 */
class HelixConstraint {

public:

	/**
	 *  Constructor
	 *
	 * @param minBP minimal number of base pairs a helix is allowed to have (>= 2)
	 * @param maxBP maximal number of base pairs a helix is allowed to have (>= bpMin)
	 * @param maxIL maximal size for each internal loop during helix computation
	 * @param maxED maximal allowed energy for ED-values
	 * @param maxE maximal allowed energy during helix computation
	 * @parem noED whether or not ED values are added in the computation of helices
	 */
	HelixConstraint( const size_t minBP
				, const size_t maxBP
				, const size_t maxIL
				, const E_type maxED
				, const E_type maxE
			    , const bool noED
				);

	virtual ~HelixConstraint();

	/**
	 *  Provides the minimum number of base pairs allowed within an helix (>=2)
	 *
	 * @return the minimum number of base pairs an helix is allowed to have (>=2)
	 */
	size_t
	getMinBasePairs() const;

	/**
	 *  Provides the maximum number of base pairs allowed within an helix (>=minBP)
	 *
	 * @return the maximum number of base pairs an helix is allowed to have (>=maxBP)
	 */
	size_t
	getMaxBasePairs() const;

	/**
	 * Provides the maximally allowed size of each internal loop considered in an helix
	 * @return  the maximally allowed size of each internal loop considered in an helix
	 */
	size_t
	getMaxIL() const;

	/**
	 * Provides the maximally allowed ED value (per sequence) for a helix to be considered
	 * @return the maximally allowed ED value (per sequence) for a helix
	 */
	E_type
	getMaxED() const;

	/**
	 * Provides the maximally allowed energy for a helix to be considered
	 * @return the maximally allowed energy for a helix
	 */
	E_type
	getMaxE() const;

	/**
	 * Whether or not only the full interaction energy (including E_init, ED,
	 * dangling end contributions, etc) are to be used for
	 * energy check (getMaxE()) or the loop-terms only.
	 * @return true if the full helix energy is to be considered
	 */
	bool
	evalFullE() const;

	/**
	 * Provides the minimum interior loop length depending on the maximum allowed number of unpaired bases
	 * @return minimum interior loop length
	 */
	size_t
	getMinInternalLoopSize() const;

	/**
     * Provides the maximal length of the helix in seq1
     * @return the maximal length of the helix in seq1
     */
	size_t
	getMaxLength1() const;

	/**
	 * Provides the maximal length of the helix in seq2
	 * @return the maximal length of the helix in seq2
	 */
	size_t
	getMaxLength2() const;

	/**
	 * Prints the helix constraint details to stream
	 * @param out the ostream to write to
	 * @param c the object to add
	 * @return the altered stream out
	 */
	friend std::ostream& operator<<(std::ostream& out, const HelixConstraint& c);

protected:

	//! the minimal number of base pairs allowed in the helix (>=2)
	size_t minBP;

	//! the maximal number of base pairs allowed in the helix (>=maxBP)
	size_t maxBP;

	//! the maximally allowed size of internal loops in the helix
	size_t maxIL;

	//! the maximal ED value (per sequence) allowed for a helix
	E_type maxED;

	//! the maximal energie allowed for a helix
	E_type maxE;

	//! decision variable, whether the full helix energy is to be used for
	//! energy checks (true) or the loop-terms only (false)
	bool checkFullE;
};


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

inline
HelixConstraint::HelixConstraint(
		const size_t minBP_
		, const size_t maxBP_
		, const size_t maxIL_
		, const E_type maxED_
		, const E_type maxE_
		, const bool checkFullE_)
	:
		minBP(minBP_)
	  , maxBP(maxBP_)
	  , maxIL(maxIL_)
	  , maxED(maxED_)
	  , maxE(maxE_)
	  , checkFullE(checkFullE_)
{
	if (minBP < 2) throw std::runtime_error("HelixConstraint() : minimal base pair number ("+toString(minBP)+") < 2");
}

/////////////////////////////////////////////////////////////////////////////

inline
HelixConstraint::~HelixConstraint() {
}

/////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixConstraint::
getMinBasePairs() const {
	return minBP;
}

/////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixConstraint::
getMaxBasePairs() const {
	return maxBP;
}

/////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixConstraint::
getMaxIL() const {
	return maxIL;
}

/////////////////////////////////////////////////////////////////////////////

inline
E_type
HelixConstraint::
getMaxED() const {
	return maxED;
}

/////////////////////////////////////////////////////////////////////////////

inline
E_type
HelixConstraint::
getMaxE() const {
	return maxE;
}

/////////////////////////////////////////////////////////////////////////////

inline
bool
HelixConstraint::
evalFullE() const {
	return checkFullE;
}

/////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixConstraint::
getMaxLength1() const {
	return getMaxBasePairs() + (getMaxBasePairs()-1) * getMaxIL();
}

/////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixConstraint::
getMaxLength2() const {
	return getMaxBasePairs() + (getMaxBasePairs()-1) * getMaxIL();
}

/////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixConstraint::
getMinInternalLoopSize() const {
	return getMaxIL();
}

/////////////////////////////////////////////////////////////////////////////

inline
std::ostream&
operator<<(std::ostream& out, const HelixConstraint& c)
{
	out <<"HelixConstraint( minBP="<<c.getMinBasePairs()
			<<", maxBP="<<c.getMaxBasePairs()
			<<", maxIL="<<c.getMaxIL()
			<<", maxED="<<E_2_Ekcal(c.getMaxED())
			<<", maxE="<<E_2_Ekcal(c.getMaxE())
			<<", evalFullE="<<(c.evalFullE()?"true":"false")
		    <<")";
	return out;

}

/////////////////////////////////////////////////////////////////////////////

} // namespace

#endif //INTARNA_HELIXCONSTRAINT_H
