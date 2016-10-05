
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
	 *
	 */
	InteractionEnergyBasePair( const Accessibility & accS1
					, const ReverseAccessibility & accS2
					, const size_t maxInternalLoopSize1 = 16
					, const size_t maxInternalLoopSize2 = 16
				);

	virtual ~InteractionEnergyBasePair();


	/**
	 * Computes the energy estimate for the interaction loop region closed by
	 * the intermolecular base pairs (i1,i2) and (j1,j2) where the regions
	 * [i1,j1] and [i2,j2] are considered unpaired.
	 * The energy estimate is the negated number of gained base pairs by
	 * closing this loop, i.e. -1 or E_INF is the internal loop size exceeds
	 * the allowed maximum (see constructor).
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param j1 the index of the first sequence interacting with j2
	 * @param i2 the index of the second sequence interacting with i1
	 * @param j2 the index of the second sequence interacting with j1
	 *
	 * @return -1 or E_INF if the allowed loop size is
	 *         exceeded
	 */
	virtual
	E_type
	getInterLoopE( const size_t i1, const size_t j1, const size_t i2, const size_t j2 ) const;

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
	getDanglingLeft( const size_t i1, const size_t i2 ) const;


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
	getDanglingRight( const size_t j1, const size_t j2 ) const;

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
	getBestStackingEnergy() const;

	/**
	 * Provides the best energy gain possible for interaction initiation
	 * for this energy model
	 * @return -1
	 */
	virtual
	E_type
	getBestInitEnergy() const;

	/**
	 * Provides the best energy gain possible for left/right dangle
	 * for this energy model
	 * @return 0
	 */
	virtual
	E_type
	getBestDangleEnergy() const;
};

#endif /* INTERACTIONENERGYBASEPAIR_H_ */
