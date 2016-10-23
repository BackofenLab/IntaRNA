
#ifndef INTERACTIONENERGYIDXOFFSET_H_
#define INTERACTIONENERGYIDXOFFSET_H_

#include "InteractionEnergy.h"

/**
 * Wrapper for a given InteractionEnergy object where indices are shifted by
 * a given offset. This is useful for local interaction computations.
 *
 * @author Martin Mann
 *
 */
class InteractionEnergyIdxOffset : public InteractionEnergy
{
public:

	/**
	 * construction
	 *
	 * @param energyOriginal wrapped energy object used for computations
	 * @param offset1 the index offset for sequence 1
	 * @param offset2 the index offset for sequence 2
	 */
	InteractionEnergyIdxOffset( const InteractionEnergy & energyOriginal
								, const size_t offset1 = 0
								, const size_t offset2 = 0 );


	virtual ~InteractionEnergyIdxOffset();

	/**
	 * Access to the currently used index offset for sequence 1
	 * @return the index offset for sequence 1 used
	 */
	size_t getOffset1() const;

	/**
	 * Sets the index offset to be used for sequence 1
	 * @param offset1 the index offset for sequence 1 to be used
	 */
	void setOffset1(size_t offset1);

	/**
	 * Access to the currently used index offset for sequence 2
	 * @return the index offset for sequence 2 used
	 */
	size_t getOffset2() const;

	/**
	 * Sets the index offset to be used for sequence 2
	 * @param offset1 the index offset for sequence 2 to be used
	 */
	void setOffset2(size_t offset2);


	///////////////  OVERWRITTEN MEMBERS USING OFFSET  /////////////////

	/**
	 * Provides the ED penalty for making a region with sequence 1 accessible
	 *
	 * Note, the indices are shifted by an offset for computation.
	 *
	 * @param i1 the start of the accessible region
	 * @param j1 the end of the accessible region
	 * @return the ED value for [i1,j1]
	 */
	virtual
	E_type
	getED1( const size_t i1, const size_t j1 ) const;

	/**
	 * Provides the ED penalty for making a region with (the reversed)
	 * sequence 2 accessible
	 *
	 * Note, the indices are shifted by an offset for computation.
	 *
	 * @param i2 the start of the accessible region
	 * @param j2 the end of the accessible region
	 * @return the ED value for [i2,j2]
	 */
	virtual
	E_type
	getED2( const size_t i2, const size_t j2 ) const;

	/**
	 * Whether or not position i is accessible for interaction in sequence 1
	 *
	 * Note, the index is shifted by an offset for computation.
	 *
	 * @param i the position of interest in sequence 1
	 * @return true if the position can partake in an interaction; false otherwise
	 */
	virtual
	bool
	isAccessible1( const size_t i ) const;

	/**
	 * Whether or not position i is accessible for interaction in sequence 2
	 *
	 * Note, the index is shifted by an offset for computation.
	 *
	 * @param i the position of interest in sequence 2
	 * @return true if the position can partake in an interaction; false otherwise
	 */
	virtual
	bool
	isAccessible2( const size_t i ) const;

	/**
	 * Checks whether or not two positions (shifted by offset) can form a base pair
	 * @param i1 index in first sequence
	 * @param i2 index in second sequence
	 * @return true if seq1(i1) can form a base pair with seq2(i2)
	 */
	virtual
	bool
	areComplementary( const size_t i1, const size_t i2 ) const;


	/**
	 * Length of sequence 1 excluding the index offset
	 * @return length of sequence 1 excluding index offset
	 */
	virtual
	size_t
	size1() const;

	/**
	 * Length of sequence 2 excluding index offset
	 * @return length of sequence 2 excluding index offset
	 */
	virtual
	size_t
	size2() const;

	/**
	 * Provides the duplex initiation energy.
	 *
	 * @return the energy for duplex initiation
	 */
	virtual
	E_type
	getE_init( ) const;


	/**
	 * Computes the energy estimate for the 'left side' interaction loop region
	 * closed by the intermolecular base pairs (i1,i2) and enclosing (j1,j2)
	 * where the regions [i1,j1] and [i2,j2] are considered unpaired or E_INF
	 * is the internal loop size exceeds the allowed maximum (see constructor).
	 *
	 * Note, the indices are shifted by an offset for computation.
	 *
	 * Note, the right interaction base pair (j1,j2) is not included in the
	 * returned energy value.
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param j1 the index of the first sequence interacting with j2 with i1<=j1
	 * @param i2 the index of the second sequence interacting with i1
	 * @param j2 the index of the second sequence interacting with j1 with i2<=j2
	 *
	 * @return the energy for the loop
	 *         or E_INF if the allowed loop size is exceeded or no valid internal loop boundaries
	 */
	virtual
	E_type
	getE_interLeft( const size_t i1, const size_t j1, const size_t i2, const size_t j2 ) const;


	/**
	 * Computes the dangling end energy penalties for the left side
	 * (i1-1 and i2-1) of the interaction closed by the intermolecular
	 * base pair (i1,i2).
	 *
	 * Note, the indices are shifted by an offset for computation.
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param i2 the index of the second sequence interacting with i1
	 *
	 * @return the dangling end penalty for the left side of the interaction
	 */
	virtual
	E_type
	getE_danglingLeft( const size_t i1, const size_t i2 ) const;


	/**
	 * Computes the dangling end energy penalties for the right side
	 * (j1+1 and j2+1) of the interaction closed by the intermolecular
	 * base pair (j1,j2).
	 *
	 * Note, the indices are shifted by an offset for computation.
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
	 * Note, the indices are shifted by an offset for computation.
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param i2 the index of the second sequence interacting with i1
	 *
	 * @return the loop closure penalty for the left side of the interaction
	 */
	virtual
	E_type
	getE_endLeft( const size_t i1, const size_t i2 ) const;

	/**
	 * Provides the penalty for closing an interaction with the given
	 * base pair on the "right side" (j1 = 3' end of seq1 of the interaction)
	 *
	 * Note, the indices are shifted by an offset for computation.
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param i2 the index of the second sequence interacting with i1
	 *
	 * @return the loop closure penalty for the right side of the interaction
	 */
	virtual
	E_type
	getE_endRight( const size_t j1, const size_t j2 ) const;

	/**
	 * Computes the probability of the dangling ends for the left side
	 * (i1-1 and i2-1) of the interaction closed by the intermolecular
	 * base pair (i1,i2) for an interaction of [i1,j1] with [i2,j2].
	 *
	 * Note, the indices are shifted by an offset for computation.
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param j1 the index of the first sequence interacting with j2 with i1<=j1
	 * @param i2 the index of the second sequence interacting with i1
	 * @param j2 the index of the second sequence interacting with j1 with i2<=j2
	 *
	 * @return the dangling end probability for the left side of the interaction
	 */
	virtual
	E_type
	getPr_danglingLeft( const size_t i1, const size_t j1, const size_t i2, const size_t j2 ) const;

	/**
	 * Computes the probability of the dangling ends for the right side
	 * (j1+1 and j2+1) of the interaction closed by the intermolecular
	 * base pair (j1,j2) for an interaction of [i1,j1] with [i2,j2].
	 *
	 * Note, the indices are shifted by an offset for computation.
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param j1 the index of the first sequence interacting with j2 with i1<=j1
	 * @param i2 the index of the second sequence interacting with i1
	 * @param j2 the index of the second sequence interacting with j1 with i2<=j2
	 *
	 * @return the dangling end probability for the right side of the interaction
	 */
	virtual
	E_type
	getPr_danglingRight( const size_t i1, const size_t j1, const size_t i2, const size_t j2 ) const;


	/**
	 * Provides the base pair encoding for the given indices after shifting by
	 * the used offset
	 * @param i1 the index in the first sequence
	 * @param i2 the index in the (reversed) second sequence
	 * @return the according base pair (i1+offset1,reverseIdx(i2+offset2))
	 */
	virtual
	Interaction::BasePair
	getBasePair( const size_t i1, const size_t i2 ) const;


	/**
	 * Provides the index within the first sequence of the given base pair
	 * shifted by the offset.
	 * @return the shifted index of the first sequence within the base pair encoding
	 */
	virtual
	size_t
	getIndex1( const Interaction::BasePair & bp ) const;


	/**
	 * Provides the index within the second sequence of the given base pair
	 * shifted by the offset.
	 * @return the shifted index of the second sequence within the base pair encoding
	 */
	virtual
	size_t
	getIndex2( const Interaction::BasePair & bp ) const;


	/**
	 * Access to the normalized temperature for Boltzmann weight computation
	 */
	virtual
	E_type
	getRT() const;


	/**
	 * Provides the best energy gain via stacking possible for this energy
	 * model
	 * @return the best stacking energy gain produced by getE_interLoop()
	 */
	virtual
	E_type
	getBestE_interLoop() const;

	/**
	 * Provides the best energy gain possible for left/right dangle
	 * for this energy model
	 * @return the best initiation energy gain produced by getE_danglingLeft() or
	 *          getE_danglingRight()
	 */
	virtual
	E_type
	getBestE_dangling() const;

	/**
	 * Provides the best energy gain possible for left/right interaction ends
	 * for this energy model
	 * @return the best end energy gain produced by getE_endLeft() or
	 *          getE_endRight()
	 */
	virtual
	E_type
	getBestE_end() const;

protected:

	/** the wrapped energy computation handler */
	const InteractionEnergy & energyOriginal;

	/** the index offset in sequence 1 */
	size_t offset1;

	/** the index offset in sequence 2 */
	size_t offset2;

};

#endif /* INTERACTIONENERGYIDXOFFSET_H_ */
