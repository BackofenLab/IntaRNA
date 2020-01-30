
#ifndef INTARNA_SEEDCONSTRAINT_H_
#define INTARNA_SEEDCONSTRAINT_H_


#include "IntaRNA/general.h"
#include "IntaRNA/IndexRangeList.h"

#include <cstddef>
#include <iostream>

namespace IntaRNA {

/**
 * Encodes seed constraints to be used for interaction prediction.
 *
 * @author Martin Mann
 *
 */
class SeedConstraint {

public:

	/**
	 * Construction
	 *
	 * @param bp of base pairs a seed has to have (>= 2)
	 * @param maxUnpairedOverall the maximal (summed) number of unpaired bases
	 *        within both seq1 and seq2 allowed within a seed
	 * @param maxUnpaired1 the maximal number of unpaired bases within seq1
	 *        allowed within a seed
	 * @param maxUnpaired2 the maximal number of unpaired bases within seq2
	 *        allowed within a seed
	 * @param maxE maximal energy a seed is allowed to have
	 * @param maxED maximal ED value a seed is allowed to have for each subsequence
	 * @param maxE maximal hybridization energy (including E_init) a seed is allowed to have
	 * @param ranges1 the index ranges of seq1 to be searched for seeds
	 * @param ranges2reversed the index reversed ranges of seq2 to be searched for seeds
	 * @param explicitSeeds the encodings of explicit seed interactions to be used
	 * @param noGUallowed whether or not GU base pairs are allowed within seeds
	 * @param noGUendAllowed whether or not GU base pairs are allowed at the ends of seeds
	 * @param noLP whether or not lonely base pairs are allowed
	 */
	SeedConstraint(  const size_t bp
				, const size_t maxUnpairedOverall
				, const size_t maxUnpaired1
				, const size_t maxUnpaired2
				, const E_type maxE
				, const E_type maxED
				, const E_type maxEhybrid
				, const IndexRangeList & ranges1
				, const IndexRangeList & ranges2reversed
				, const std::string & explicitSeeds
				, const bool noGUallowed
				, const bool noGUendAllowed
				, const bool noLP
				);

	virtual ~SeedConstraint();


	/**
	 * Provides the number of base pairs to be present within a seed
	 *
	 * @return the number of base pairs a seed has to have
	 */
	size_t
	getBasePairs() const;

	/**
	 * Provides the overall maximal number of unpaired bases within a seed
	 *
	 * @return the overall maximal number of unpaired bases within
	 *         a seed is allowed to have
	 */
	size_t
	getMaxUnpairedOverall() const;

	/**
	 * Provides the maximal number of unpaired bases within the first sequence
	 * within a seed
	 *
	 * @return the maximal number of unpaired bases within the first sequence
	 *         a seed is allowed to have
	 */
	size_t
	getMaxUnpaired1() const;

	/**
	 * Provides the maximal number of unpaired bases within the second sequence
	 * within a seed
	 *
	 * @return the maximal number of unpaired bases within the second sequence
	 *         a seed is allowed to have
	 */
	size_t
	getMaxUnpaired2() const;

	/**
	 * Provides the maximally allowed energy for seeds to be considered
	 *
	 * @return the maximally allowed energy for a seed
	 */
	E_type
	getMaxE() const;

	/**
	 * Provides the maximally allowed ED value for each seed subsequence to be
	 * considered
	 *
	 * @return the maximally allowed ED value (per sequence) for a seed
	 */
	E_type
	getMaxED() const;

	/**
	 * Provides the maximally allowed hybridization energy (including E_init)
	 * for seeds to be considered
	 *
	 * @return the maximally allowed hybridization energy for a seed
	 */
	E_type
	getMaxEhybrid() const;

	/**
	 * Provides the maximal length of the seed in seq1
	 * @return the maximal length of the seed in seq1
	 */
	size_t
	getMaxLength1() const;

	/**
	 * Provides the maximal length of the seed in seq2
	 * @return the maximal length of the seed in seq2
	 */
	size_t
	getMaxLength2() const;

	/**
	 * Whether or not GU base pairs are allowed within seeds
	 * @return true if GU base pairs are allowed; false otherwise
	 */
	bool
	isGUallowed() const;

	/**
	 * Whether or not GU base pairs are allowed at seed ends
	 * @return true if GU base pairs are allowed at seed ends; false otherwise
	 */
	bool
	isGUendAllowed() const;

	/**
	 * Whether or not lonely base pairs are allowed
	 * @return true if lonely base pairs are allowed; false otherwise
	 */
	bool
	isLpAllowed() const;

	/**
	 * Index ranges in seq1 to be searched for seeds or empty if all indices
	 * are to be considered.
	 * @return index ranges to be searched or empty list if all indices relevant
	 */
	const IndexRangeList &
	getRanges1() const;

	/**
	 * Index ranges in seq2 to be searched for seeds or empty if all indices
	 * are to be considered.
	 * @return index ranges to be searched or empty list if all indices relevant
	 */
	const IndexRangeList &
	getRanges2() const;

	/**
	 * Access to the explicit seed encodings to be used for prediction.
	 * @return the provided explicit seed encodings
	 */
	const std::string &
	getExplicitSeeds() const;

	/**
	 * Index ranges in seq1 to be searched for seeds or empty if all indices
	 * are to be considered.
	 * @return index ranges to be searched or empty list if all indices relevant
	 */
	IndexRangeList &
	getRanges1();

	/**
	 * Index ranges in seq2 to be searched for seeds or empty if all indices
	 * are to be considered.
	 * @return index ranges to be searched or empty list if all indices relevant
	 */
	IndexRangeList &
	getRanges2();

	/**
	 * Access to the explicit seed encodings to be used for prediction.
	 * @return the provided explicit seed encodings
	 */
	std::string &
	getExplicitSeeds();

	/**
	 * Prints the seed constraint details to stream
	 * @param out the ostream to write to
	 * @param c the object to add
	 * @return the altered stream out
	 */
	friend std::ostream& operator<<(std::ostream& out, const SeedConstraint& c);

protected:

	//! the number of base pairs to be present in a seed
	size_t bp;

	//! the overall summed maximally allowed number of unpaired bases in a seed
	size_t maxUnpairedOverall;

	//! the maximally allowed number of unpaired bases in seed seq1
	size_t maxUnpaired1;

	//! the maximally allowed number of unpaired bases in seed seq2
	size_t maxUnpaired2;

	//! the maximal energy allowed for a seed
	E_type maxE;

	//! the maximal ED value (per sequence) allowed for a seed
	E_type maxED;

	//! the maximal hybridization energy (incl E_init) allowed for a seed
	E_type maxEhybrid;

	//! the index ranges of seq1 to be searched for seeds
	IndexRangeList ranges1;

	//! the index ranges of seq2 to be searched for seeds
	IndexRangeList ranges2;

	//! the string encoding of the explicit seed interactions to be used
	std::string explicitSeeds;

	//! whether or not GU base pairs are allowed within seeds
	bool bpGUallowed;

	//! whether or not GU base pairs are allowed at seed ends
	bool bpGUendAllowed;

	//! whether or not lonely base pairs are allowed
	bool lpAllowed;

};


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

inline
SeedConstraint::SeedConstraint(
		 const size_t bp_
		, const size_t maxUnpairedOverall_
		, const size_t maxUnpaired1_
		, const size_t maxUnpaired2_
		, const E_type maxE_
		, const E_type maxED_
		, const E_type maxEhybrid
		, const IndexRangeList & ranges1
		, const IndexRangeList & ranges2reversed
		, const std::string & explicitSeeds
		, const bool noGUallowed
		, const bool noGUendAllowed
		, const bool noLP
		)
 :
	  bp(bp_)
	, maxUnpairedOverall(maxUnpairedOverall_)
	, maxUnpaired1(std::min(maxUnpaired1_,maxUnpairedOverall_)) // exclude too large boundaries
	, maxUnpaired2(std::min(maxUnpaired2_,maxUnpairedOverall_)) // exclude too large boundaries
	, maxE(maxE_)
	, maxED(maxED_)
	, maxEhybrid(maxEhybrid)
	, ranges1(ranges1)
	, ranges2(ranges2reversed)
	, explicitSeeds(explicitSeeds)
	, bpGUallowed(!noGUallowed)
	, bpGUendAllowed(!noGUendAllowed)
	, lpAllowed(!noLP)
{
	if (bp < 2) throw std::runtime_error("SeedConstraint() : base pair number ("+toString(bp)+") < 2");
}

/////////////////////////////////////////////////////////////////////////////

inline
SeedConstraint::~SeedConstraint() {
}

/////////////////////////////////////////////////////////////////////////////

inline
size_t
SeedConstraint::
getBasePairs() const {
	return bp;
}

/////////////////////////////////////////////////////////////////////////////

inline
E_type
SeedConstraint::
getMaxE() const {
	return maxE;
}

/////////////////////////////////////////////////////////////////////////////

inline
E_type
SeedConstraint::
getMaxED() const {
	return maxED;
}

/////////////////////////////////////////////////////////////////////////////

inline
E_type
SeedConstraint::
getMaxEhybrid() const {
	return maxEhybrid;
}

/////////////////////////////////////////////////////////////////////////////

inline
size_t
SeedConstraint::
getMaxUnpaired1() const {
	return maxUnpaired1;
}

/////////////////////////////////////////////////////////////////////////////

inline
size_t
SeedConstraint::
getMaxUnpaired2() const {
	return maxUnpaired2;
}

/////////////////////////////////////////////////////////////////////////////

inline
size_t
SeedConstraint::
getMaxUnpairedOverall() const {
	return maxUnpairedOverall;
}

/////////////////////////////////////////////////////////////////////////////

inline
size_t
SeedConstraint::
getMaxLength1() const {
	return getBasePairs() + getMaxUnpaired1();
}

/////////////////////////////////////////////////////////////////////////////

inline
size_t
SeedConstraint::
getMaxLength2() const {
	return getBasePairs() + getMaxUnpaired2();
}

/////////////////////////////////////////////////////////////////////////////

inline
const IndexRangeList &
SeedConstraint::
getRanges1() const {
	return ranges1;
}

/////////////////////////////////////////////////////////////////////////////

inline
const IndexRangeList &
SeedConstraint::
getRanges2() const {
	return ranges2;
}

/////////////////////////////////////////////////////////////////////////////

inline
IndexRangeList &
SeedConstraint::
getRanges1() {
	return ranges1;
}

/////////////////////////////////////////////////////////////////////////////

inline
IndexRangeList &
SeedConstraint::
getRanges2() {
	return ranges2;
}

/////////////////////////////////////////////////////////////////////////////

inline
const std::string &
SeedConstraint::
getExplicitSeeds() const {
	return explicitSeeds;
}

/////////////////////////////////////////////////////////////////////////////

inline
std::string &
SeedConstraint::
getExplicitSeeds() {
	return explicitSeeds;
}

/////////////////////////////////////////////////////////////////////////////

inline
bool
SeedConstraint::
isGUallowed() const {
	return bpGUallowed;
}

/////////////////////////////////////////////////////////////////////////////

inline
bool
SeedConstraint::
isGUendAllowed() const {
	return bpGUendAllowed;
}

/////////////////////////////////////////////////////////////////////////////

inline
bool
SeedConstraint::
isLpAllowed() const {
	return lpAllowed;
}

/////////////////////////////////////////////////////////////////////////////

inline
std::ostream&
operator<<(std::ostream& out, const SeedConstraint& c)
{
	out <<"SeedConstraint( bp="<<c.getBasePairs()
			<<", up="<<c.getMaxUnpairedOverall()
			<<", up1="<<c.getMaxUnpaired1()
			<<", up2="<<c.getMaxUnpaired2()
			<<", E="<<E_2_Ekcal(c.getMaxE())
			<<", ED="<<E_2_Ekcal(c.getMaxED())
			<<", Ehybrid="<<E_2_Ekcal(c.getMaxEhybrid())
			<<", noGU="<<(c.isGUallowed() ? "false":"true")
			<<", noGUend="<<(c.isGUendAllowed() ? "false":"true")
			<<")";
	return out;
}

/////////////////////////////////////////////////////////////////////////////

} // namespace

#endif /* SEEDCONSTRAINT_H_ */
