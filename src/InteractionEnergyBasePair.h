
#ifndef INTERACTIONENERGYBASEPAIR_H_
#define INTERACTIONENERGYBASEPAIR_H_

#include "InteractionEnergy.h"

/**
 * Implements a simple energy interface that is based on base pair counts only.
 *
 * @author Martin Mann 2014
 */
class InteractionEnergyBasePair: public InteractionEnergy {


public:

	/**
	 * Construct energy utility object given the accessibility ED values for
	 * two sequences.
	 *
	 * @param accS1 accessibility of the first sequence
	 * @param accS2 accessibility of the second sequence
	 * @param maxInternalLoopSize1 maximal number of enclosed unpaired positions
	 *          between two intermolecular base pairs in sequence 1, ie it holds
	 *          for an intermolecular loop closed by base pairs (i1,i2) and
	 *          (j1,j2) : (j1-i1+1) <= maxInternalLoopSize
	 * @param maxInternalLoopSize2 maximal number of enclosed unpaired positions
	 *          between two intermolecular base pairs in sequence 2, ie it holds
	 *          for an intermolecular loop closed by base pairs (i1,i2) and
	 *          (j1,j2) : (j2-i2+1) <= maxInternalLoopSize
	 * @param initES whether or not to compute and initialize ES values
	 *
	 */
	InteractionEnergyBasePair( const Accessibility & accS1
					, const ReverseAccessibility & accS2
					, const size_t maxInternalLoopSize1 = 16
					, const size_t maxInternalLoopSize2 = 16
					, const bool initES = false
				);

	virtual ~InteractionEnergyBasePair();


	/**
	 * Provides the ensemble energy (ES) of all intramolecular substructures
	 * that can be formed within a given region of sequence 1 under the
	 * assumption that the region is part of an (intermolecular) multiloop,
	 * i.e. at least one base pair is formed by each substructure.
	 *
	 * If no structure can be formed within the region, E_INF is returned.
	 *
	 * @param i1 the start of the structured region of seq1
	 * @param j1 the end of the structured region of seq1
	 * @return the ES value for [i1,j1] or E_INF if no intramolecular
	 *         structure can be formed
	 */
	virtual
	E_type
	getES1( const size_t i1, const size_t j1 ) const;

	/**
	 * Provides the ensemble energy (ES) of all intramolecular substructures
	 * that can be formed within a given region of sequence 2 under the
	 * assumption that the region is part of an (intermolecular) multiloop,
	 * i.e. at least one base pair is formed by each substructure.
	 *
	 * If no structure can be formed within the region, E_INF is returned.
	 *
	 * @param i2 the start of the structured region of seq2
	 * @param j2 the end of the structured region of seq2
	 * @return the ES value for [i2,j2] or E_INF if no intramolecular
	 *         structure can be formed
	 */
	virtual
	E_type
	getES2( const size_t i2, const size_t j2 ) const;

	/**
	 * Provides the energy contribution for a given number of unpaired
	 * nucleotides under the
	 * assumption that the region is part of an (intermolecular) multiloop.
	 *
	 * @param numUnpaired the number of unpaired bases
	 * @return the energy contribution of the given number of unpaired bases
	 *         within an intramolecular multiloop, which is always 0
	 */
	virtual
	E_type
	getE_multiUnpaired( const size_t numUnpaired ) const;


	/**
	 * Provides the energy contribution/penalty of the helix repesented by the
	 * interaction right of a multi-site gap starting with base pair (j1,j2)
	 *
	 * @param j1 the end of the gap in seq1, ie the first base paired in the
	 *           interaction site to the right of the gap
	 * @param j2 the end of the gap in seq2, ie the first base paired in the
	 *           interaction site to the right of the gap
	 *
	 * @return 0.0
	 */
	virtual
	E_type
	getE_multiHelix( const size_t j1, const size_t j2 ) const;

	/**
	 * Provides the energy contribution/penalty for closing an intermolecular
	 * multiloop on the left of a multi-site gap.
	 *
	 * @return 0.0
	 */
	virtual
	E_type
	getE_multiClosing() const;


	/**
	 * Provides the duplex initiation energy.
	 *
	 * @return the energy for duplex initiation
	 */
	virtual
	E_type
	getE_init() const;


	/**
	 * Computes the energy estimate for the interaction loop region closed by
	 * the intermolecular base pairs (i1,i2) and (j1,j2) where the regions
	 * [i1,j1] and [i2,j2] are considered unpaired.
	 * The energy estimate is the negated number of gained base pairs by
	 * closing this loop, i.e. -1 or E_INF is the internal loop size exceeds
	 * the allowed maximum (see constructor).
	 *
	 * @param i1 the index of the first sequence (<j1) interacting with i2
	 * @param j1 the index of the first sequence (>i1) interacting with j2
	 * @param i2 the index of the second sequence (<j2) interacting with i1
	 * @param j2 the index of the second sequence (>i2) interacting with j1
	 *
	 * @return -1 or E_INF if the allowed loop size is exceeded or no valid internal loop boundaries
	 */
	virtual
	E_type
	getE_interLeft( const size_t i1, const size_t j1, const size_t i2, const size_t j2 ) const;

	/**
	 * Computes the dangling end energy penalty estimate for the left side of
	 * an interaction loop region closed on the left by the intermolecular
	 * base pair (i1,i2).
	 *
	 * This penalty is always zero for this base pair based energy function.
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param i2 the index of the second sequence interacting with i1
	 *
	 * @return 0
	 */
	virtual
	E_type
	getE_danglingLeft( const size_t i1, const size_t i2 ) const;


	/**
	 * Computes the dangling end energy penalty estimate for the right side of
	 * an interaction loop region closed on the right by the intermolecular
	 * base pair (j1,j2).
	 *
	 * @param j1 the index of the first sequence interacting with j2
	 * @param j2 the index of the second sequence interacting with j1
	 *
	 * @return the dangling end penalty for the right side of the interaction
	 */
	virtual
	E_type
	getE_danglingRight( const size_t j1, const size_t j2 ) const;

	/**
	 * Provides the penalty for closing an interaction with the given
	 * base pair on the "left side" (i1 = 5' end of seq1 of the interaction)
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param i2 the index of the second sequence interacting with i1
	 *
	 * @return 0.0
	 */
	virtual
	E_type
	getE_endLeft( const size_t i1, const size_t i2 ) const {
		return 0.0;
	}

	/**
	 * Provides the penalty for closing an interaction with the given
	 * base pair on the "right side" (j1 = 3' end of seq1 of the interaction)
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param i2 the index of the second sequence interacting with i1
	 *
	 * @return 0.0
	 */
	virtual
	E_type
	getE_endRight( const size_t j1, const size_t j2 ) const {
		return 0.0;
	}

	/**
	 * Returns always RT=1 due to the lack of reasonable values for this energy
	 * function.
	 * @return 1.0
	 */
	virtual
	E_type
	getRT() const {
		return 1.0;
	}

	/**
	 * Provides the best energy gain via stacking possible for this energy
	 * model
	 * @return -1
	 */
	virtual
	E_type
	getBestE_interLoop() const;

	/**
	 * Provides the best energy gain possible for left/right dangle
	 * for this energy model
	 * @return 0
	 */
	virtual
	E_type
	getBestE_dangling() const;

	/**
	 * Provides the best energy gain possible for left/right interaction ends
	 * for this energy model
	 * @return 0
	 */
	virtual
	E_type
	getBestE_end() const {
		return 0;
	}

};



////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

inline
InteractionEnergyBasePair::InteractionEnergyBasePair(
		const Accessibility & accS1
		, const ReverseAccessibility & accS2
		, const size_t maxInternalLoopSize1
		, const size_t maxInternalLoopSize2
		, const bool initES
	)
 :
	InteractionEnergy(accS1, accS2, maxInternalLoopSize1, maxInternalLoopSize2)
{
	if (initES) {
		NOTIMPLEMENTED("InteractionEnergyVrna() : ES computation missing");
	}
}

////////////////////////////////////////////////////////////////////////////

inline
InteractionEnergyBasePair::~InteractionEnergyBasePair()
{
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyBasePair::
getES1( const size_t i1, const size_t j1 ) const
{
#if IN_DEBUG_MODE
	// sanity check
	if (i1>j1) throw std::runtime_error("InteractionEnergy::getES1(i1="+toString(i1)+" > j1="+toString(j1));
	if (j1>=size1()) throw std::runtime_error("InteractionEnergy::getES1() : j1="+toString(j1)+" >= size1()="+toString(size1()));
#endif

	NOTIMPLEMENTED("InteractionEnergyVrna::getES2() : ES computation missing");
	// return computed value
	return E_INF;
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyBasePair::
getES2( const size_t i2, const size_t j2 ) const
{
#if IN_DEBUG_MODE
	// sanity check
	if (i2>j2) throw std::runtime_error("InteractionEnergy::getES2(i2="+toString(i2)+" > j2="+toString(j2));
	if (j2>=size2()) throw std::runtime_error("InteractionEnergy::getES2() : j2="+toString(j2)+" >= size2()="+toString(size2()));
#endif

	NOTIMPLEMENTED("InteractionEnergyVrna::getES2() : ES computation missing");
	// return computed value
	return E_INF;
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyBasePair::
getE_multiUnpaired( const size_t numUnpaired ) const
{
	return 0.0;
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyBasePair::
getE_multiHelix( const size_t j1, const size_t j2 ) const
{
	return 0.0;
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyBasePair::
getE_multiClosing() const
{
	return 0.0;
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyBasePair::
getE_init() const
{
	return -1.0;
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyBasePair::
getE_interLeft( const size_t i1, const size_t j1, const size_t i2, const size_t j2 ) const
{
	// if valid internal loop
	if ( isValidInternalLoop(i1,j1,i2,j2) ) {
		// return negated number of gained base pairs by closing this loop = -1
		return getBestE_interLoop();
	} else {
		return E_INF;
	}
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyBasePair::
getE_danglingLeft( const size_t i1, const size_t i2 ) const
{
	// no dangling end contribution
	return (E_type)0;
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyBasePair::
getE_danglingRight( const size_t j1, const size_t j2 ) const
{
	// no dangling end contribution
	return (E_type)0;
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyBasePair::
getBestE_interLoop() const
{
	return getE_init();
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyBasePair::
getBestE_dangling() const
{
	return (E_type)0.0;
}

////////////////////////////////////////////////////////////////////////////


#endif /* INTERACTIONENERGYBASEPAIR_H_ */
