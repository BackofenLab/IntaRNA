
#ifndef INTARNA_INTERACTIONENERGYIDXOFFSET_H_
#define INTARNA_INTERACTIONENERGYIDXOFFSET_H_

#include "IntaRNA/InteractionEnergy.h"

namespace IntaRNA {

/**
 * Wrapper for a given InteractionEnergy object where indices are shifted by
 * a given positive offset (shifted towards infinity).
 * This is useful for local interaction computations.
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
	 * @param offset2 the index offset for sequence 2 to be used
	 */
	void setOffset2(size_t offset2);


	/**
	 * Provides the overall energy for an interaction from [i1,j1] in the first
	 * sequence and [i2,j2] in the second sequence given the hybridization
	 * energy contribution.
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param j1 the index of the first sequence interacting with j2 with i1<=j1
	 * @param i2 the index of the second sequence interacting with i1
	 * @param j2 the index of the second sequence interacting with j1 with i2<=j2
	 * @param hybridE the hybridization energy for the interaction
	 *
	 * @return E = hybridE
	 * 				+ ED1(i1,j1) + ED2(i2,j2)
	 * 				+ Edangle(i1,i2) + Edangle(j1,j2)
	 * 				+ Eend(i1,i2) + Eend(j1,j2)
	 */
	virtual
	E_type
	getE( const size_t i1, const size_t j1
			, const size_t i2, const size_t j2
			, const E_type hybridE ) const;

	/**
	 * Provides the ensemble energy for a given partition function Z.
	 *
	 * @param Z the ensemble's partition function to convert
	 *
	 * @return E = -RT * log( Z )
	 */
	virtual
	E_type
	getE( const Z_type Z ) const;

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
	 * Checks whether or not two positions can form a GU base pair
	 * @param i1 index in first sequence
	 * @param i2 index in second sequence
	 * @return true if seq1(i1) can form a GU base pair with seq2(i2)
	 */
	virtual
	bool
	isGU( const size_t i1, const size_t i2 ) const;

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
	 * Provides the energy contribution of an interaction site gap, i.e. the
	 * provided regions are without intermolecular base pairs but are considered
	 * to be involved in intramolecular base pairs only. The multi-site gap is
	 * scored according to a multiloop in a single structure prediction model.
	 * The ends of the two regions are supposed to form an intermolecular base
	 * pair each, i.e. (i1,i2) and (j1,j2) have to be complementary.
	 *
	 * @param i1 the start of the structured region of seq1
	 * @param j1 the end of the structured region of seq1
	 * @param i2 the start of the structured region of seq2
	 * @param j2 the end of the structured region of seq2
	 * @param ES_mode defines for which sequence intramolecular structure
	 *          contributions are to be considered
	 * @return the energy contribution of a multi-site interaction gap
	 */
	virtual
	E_type
	getE_multi(  const size_t i1, const size_t j1
				, const size_t i2, const size_t j2
				, const ES_multi_mode ES_mode ) const;

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
	 *         within an intramolecular multiloop
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
	 *
	 * @return the energy contribution/penalty of the intermolecular helix
	 *         within an intramolecular multiloop
	 */
	virtual
	E_type
	getE_multiHelix( const size_t j1, const size_t j2 ) const;

	/**
	 * Provides the energy contribution/penalty for closing an intermolecular
	 * interaction left of a multi-site gap.
	 *
	 * @return the energy contribution/penalty of the intermolecular helix
	 *         within an intramolecular multiloop
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
	 * @param j1 the index of the first sequence interacting with j2
	 * @param j2 the index of the second sequence interacting with j1
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
	Z_type
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
	Z_type
	getPr_danglingRight( const size_t i1, const size_t j1, const size_t i2, const size_t j2 ) const;



	/**
	 * Access to the normalized temperature for Boltzmann weight computation
	 */
	virtual
	Z_type
	getRT() const;

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
	 * Provides the (reversed) index within the second sequence of the given base pair
	 * shifted by the offset.
	 * @return the shifted index of the second sequence within the base pair encoding
	 */
	virtual
	size_t
	getIndex2( const Interaction::BasePair & bp ) const;


	/**
	 * Checks whether or not the given indices mark valid internal loop
	 * boundaries, i.e.
	 *  - (i1,i2) and (j1,j2) are complementary
	 *  - i1..j1 and i2..j2 are allowed loop regions
	 *  - no boundary overlap ( (j1-i1==0 && j2-i2==0) || (j1-i1>0 && j2-i2>0) )
	 *  - if !internalLoopGU : both ends are no GU base pairs
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param j1 the index of the first sequence interacting with j2 with i1<=j1
	 * @param i2 the index of the second sequence interacting with i1
	 * @param j2 the index of the second sequence interacting with j1 with i2<=j2
	 *
	 * @return true if the boundaries are sound for internal loop calculation;
	 *         false otherwise
	 */
	virtual
	bool
	isValidInternalLoop( const size_t i1, const size_t j1, const size_t i2, const size_t j2 ) const;

	/**
	 * Provides the overall ensemble energy for sequence 1
	 * given its accessibility constraints
	 * @return Eall(constraint-conform intra-molecular structures for seq1)
	 */
	virtual
	E_type
	getEall1() const;

	/**
	 * Provides the overall ensemble energy for sequence 2
	 * given its accessibility constraints
	 * @return Eall(constraint-conform intra-molecular structures for seq2)
	 */
	virtual
	E_type
	getEall2() const;


protected:

	/** the wrapped energy computation handler */
	const InteractionEnergy & energyOriginal;

	/** the index offset in sequence 1 */
	size_t offset1;

	/** the index offset in sequence 2 */
	size_t offset2;

};



//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


inline
InteractionEnergyIdxOffset::
InteractionEnergyIdxOffset( const InteractionEnergy & energyOriginal
			, const size_t offset1
			, const size_t offset2 )

 :
	InteractionEnergy( energyOriginal )
	, energyOriginal(energyOriginal)
	, offset1(offset1)
	, offset2(offset2)
{
#if INTARNA_IN_DEBUG_MODE
	// input sanity checks
	setOffset1(offset1);
	setOffset2(offset2);
#endif
}

//////////////////////////////////////////////////////////////////////////

inline
InteractionEnergyIdxOffset::~InteractionEnergyIdxOffset()
{
}

//////////////////////////////////////////////////////////////////////////

inline
size_t
InteractionEnergyIdxOffset::
getOffset1() const
{
	return offset1;
}

//////////////////////////////////////////////////////////////////////////

inline
void
InteractionEnergyIdxOffset::
setOffset1(size_t offset1)
{
#if INTARNA_IN_DEBUG_MODE
	if (offset1 >= energyOriginal.size1()) {
		throw std::runtime_error("InteractionEnergyIdxOffset : offset1 "+toString(offset1)
				+" > seq1.length "+toString(energyOriginal.size1()));
	}
#endif
	this->offset1 = offset1;
}

//////////////////////////////////////////////////////////////////////////

inline
size_t
InteractionEnergyIdxOffset::
getOffset2() const
{
	return offset2;
}

//////////////////////////////////////////////////////////////////////////

inline
void
InteractionEnergyIdxOffset::
setOffset2(size_t offset2)
{
#if INTARNA_IN_DEBUG_MODE
	if (offset2 >= energyOriginal.size2()) {
		throw std::runtime_error("InteractionEnergyIdxOffset : offset2 "+toString(offset2)
				+" > seq2.length "+toString(energyOriginal.size2()));
	}
#endif
	this->offset2 = offset2;
}

//////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyIdxOffset::
getED1( const size_t i1, const size_t j1 ) const
{
	return energyOriginal.getED1( i1+offset1, j1+offset1 );
}

//////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyIdxOffset::
getED2( const size_t i2, const size_t j2 ) const
{
	return energyOriginal.getED2( i2+offset2, j2+offset2 );
}

////////////////////////////////////////////////////////////////////////////

inline
bool
InteractionEnergyIdxOffset::
isAccessible1( const size_t i ) const
{
	return energyOriginal.isAccessible1(i+offset1);
}

////////////////////////////////////////////////////////////////////////////

inline
bool
InteractionEnergyIdxOffset::
isAccessible2( const size_t i ) const
{
	return energyOriginal.isAccessible2(i+offset2);
}

////////////////////////////////////////////////////////////////////////////

inline
bool
InteractionEnergyIdxOffset::
isGU( const size_t i1, const size_t i2 ) const
{
	return energyOriginal.isGU( i1+offset1, i2+offset2);
}

//////////////////////////////////////////////////////////////////////////

inline
bool
InteractionEnergyIdxOffset::
areComplementary( const size_t i1, const size_t i2 ) const
{
	return energyOriginal.areComplementary( i1+offset1, i2+offset2 );
}

//////////////////////////////////////////////////////////////////////////

inline
size_t
InteractionEnergyIdxOffset::
size1() const
{
	return energyOriginal.size1()-offset1;
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
InteractionEnergyIdxOffset::
size2() const
{
	return energyOriginal.size2()-offset2;
}

//////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyIdxOffset::
getES1( const size_t i1, const size_t j1 ) const
{
	return energyOriginal.getES1( i1+offset1, j1+offset1 );
}

//////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyIdxOffset::
getES2( const size_t i2, const size_t j2 ) const
{
	return energyOriginal.getES2( i2+offset2, j2+offset2 );
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyIdxOffset::
getE_multiUnpaired( const size_t numUnpaired ) const
{
	return energyOriginal.getE_multiUnpaired( numUnpaired );
}


//////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyIdxOffset::
getE_multiHelix( const size_t j1, const size_t j2 ) const
{
	return energyOriginal.getE_multiHelix(j1+offset1, j2+offset2);
}

//////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyIdxOffset::
getE_multiClosing() const
{
	return energyOriginal.getE_multiClosing();
}

//////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyIdxOffset::
getE_init( ) const
{
	return energyOriginal.getE_init();
}

//////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyIdxOffset::
getE_interLeft( const size_t i1, const size_t j1, const size_t i2, const size_t j2 ) const
{
	return energyOriginal.getE_interLeft(i1+offset1, j1+offset1, i2+offset2, j2+offset2);
}


//////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyIdxOffset::
getE_danglingLeft( const size_t i1, const size_t i2 ) const
{
	return energyOriginal.getE_danglingLeft( i1+offset1, i2+offset2 );
}

//////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyIdxOffset::
getE_danglingRight( const size_t j1, const size_t j2 ) const
{
	return energyOriginal.getE_danglingRight( j1+offset1, j2+offset2 );
}

//////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyIdxOffset::
getE_endLeft( const size_t i1, const size_t i2 ) const
{
	return energyOriginal.getE_endLeft( i1+offset1, i2+offset2 );
}

//////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyIdxOffset::
getE_endRight( const size_t j1, const size_t j2 ) const
{
	return energyOriginal.getE_endRight( j1+offset1, j2+offset2 );
}

//////////////////////////////////////////////////////////////////////////

inline
Z_type
InteractionEnergyIdxOffset::
getPr_danglingLeft( const size_t i1, const size_t j1, const size_t i2, const size_t j2 ) const
{
	return energyOriginal.getPr_danglingLeft(i1+offset1, j1+offset1, i2+offset2, j2+offset2);
}

//////////////////////////////////////////////////////////////////////////

inline
Z_type
InteractionEnergyIdxOffset::
getPr_danglingRight( const size_t i1, const size_t j1, const size_t i2, const size_t j2 ) const
{
	return energyOriginal.getPr_danglingRight(i1+offset1, j1+offset1, i2+offset2, j2+offset2);
}

//////////////////////////////////////////////////////////////////////////

inline
Interaction::BasePair
InteractionEnergyIdxOffset::
getBasePair( const size_t i1, const size_t i2 ) const
{
	return energyOriginal.getBasePair( i1+offset1, i2+offset2 );
}

//////////////////////////////////////////////////////////////////////////

inline
size_t
InteractionEnergyIdxOffset::
getIndex1( const Interaction::BasePair & bp ) const
{
#if INTARNA_IN_DEBUG_MODE
	if (energyOriginal.getIndex1(bp)<offset1) throw std::runtime_error("InteractionEnergyIdxOffset::getIndex1("+toString(energyOriginal.getIndex1(bp))+") < offset1 = "+toString(offset1));
#endif
	return energyOriginal.getIndex1(bp)-offset1;
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
InteractionEnergyIdxOffset::
getIndex2( const Interaction::BasePair & bp ) const
{
#if INTARNA_IN_DEBUG_MODE
	if (energyOriginal.getIndex2(bp)<offset2) throw std::runtime_error("InteractionEnergyIdxOffset::getIndex2("+toString(energyOriginal.getIndex2(bp))+") < offset2 = "+toString(offset2));
#endif
	return energyOriginal.getIndex2(bp)-offset2;
}

////////////////////////////////////////////////////////////////////////////

inline
Z_type
InteractionEnergyIdxOffset::
getRT() const
{
	return energyOriginal.getRT();
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyIdxOffset::
getE( const Z_type Z ) const
{
	return energyOriginal.getE( Z );
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyIdxOffset::
getE( const size_t i1, const size_t j1
		, const size_t i2, const size_t j2
		, const E_type hybridE ) const
{
	return energyOriginal.getE( i1+offset1, j1+offset1, i2+offset2, j2+offset2, hybridE );
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyIdxOffset::
getE_multi(  const size_t i1, const size_t j1
			, const size_t i2, const size_t j2
			, const ES_multi_mode ES_mode ) const
{
	return energyOriginal.getE_multi( i1+offset1, j1+offset1, i2+offset2, j2+offset2, ES_mode );
}

////////////////////////////////////////////////////////////////////////////

inline
bool
InteractionEnergyIdxOffset::
isValidInternalLoop( const size_t i1, const size_t j1, const size_t i2, const size_t j2 ) const
{
	return energyOriginal.isValidInternalLoop( i1+offset1, j1+offset1, i2+offset2, j2+offset2 );
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyIdxOffset::
getEall1() const
{
	return energyOriginal.getEall1();
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
InteractionEnergyIdxOffset::
getEall2() const
{
	return energyOriginal.getEall2();
}

////////////////////////////////////////////////////////////////////////////


} // namespace


#endif /* INTERACTIONENERGYIDXOFFSET_H_ */
