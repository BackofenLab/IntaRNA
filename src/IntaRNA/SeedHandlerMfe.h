
#ifndef INTARNA_SEEDHANDLERMFE_H_
#define INTARNA_SEEDHANDLERMFE_H_

#include "IntaRNA/InteractionEnergy.h"
#include "IntaRNA/SeedConstraint.h"
#include "IntaRNA/SeedHandler.h"

#include <boost/multi_array.hpp>

#include <boost/numeric/ublas/matrix.hpp>

namespace IntaRNA {

/**
 * Handler to provide mfe seed interaction information for each intermolecular
 * index combination (= left end of seed).
 *
 */
class SeedHandlerMfe : public SeedHandler
{
public:

	//! 5D matrix type to hold the mfe energies for seed interactions
	//! of the ranges i1..(i1+bp+u1-1) with i2..(i2+bp+u2-1), with
	//! i1,i2 = the start index of the seed in seq1/2
	//! bp = the number of base pairs within the seed
	//! bpInbetween = the number of base pairs enclosed by left and right base pair, ie. == (bp-2)
	//! u1/u2 = the number of unpaired positions within the seed,
	//! using the index [i1][i2][bpInbetween][u1][u2] or a SeedIndex object
	typedef boost::multi_array<E_type,5> SeedRecMatrix;

	//! defines the seed data {{ i1, i2, bpInbetween, u1, u2 }} to access elements of
	//! the SeedRecMatrix
	typedef boost::array<SeedRecMatrix::index, 5> SeedIndex;

	//! matrix to store the seed information for each seed left side (i1,i2);
	//! it holds both the energy (first) as well as the length of the seed using
	//! the length combination using encodeSeedLength()
	typedef boost::numeric::ublas::matrix< std::pair<E_type, size_t> > SeedMatrix;


public:

	/**
	 * Construction
	 * @param energy the energy function to be used for seed prediction
	 * @param seedConstraint the seed constraint to be applied
	 */
	SeedHandlerMfe(
			const InteractionEnergy & energy
			, const SeedConstraint & seedConstraint
			);

	/**
	 * destruction
	 */
	virtual ~SeedHandlerMfe();

	/**
	 * Computes the seed matrix for the given interval boundaries
	 * @param i1 the first index of seq1 that might interact
	 * @param j1 the last index of seq1 that might interact
	 * @param i2 the first index of seq2 that might interact
	 * @param j2 the last index of seq2 that might interact
	 * @return the number of potential seed interactions
	 */
	virtual
	size_t
	fillSeed(const size_t i1, const size_t j1, const size_t i2, const size_t j2);

	/**
	 * Identifies the base pairs of the mfe seed interaction starting at i1,i2
	 * and writes them to the provided container
	 *
	 * NOTE: the left- and right-most base pairs are excluded!
	 *
	 * @param interaction the container to add the base pairs too
	 * @param i1 the start of the seed in seq1
	 * @param i2 the start of the seed in seq2
	 */
	virtual
	void
	traceBackSeed( Interaction & interaction, const size_t i1, const size_t i2) const;


	/**
	 * Access to the mfe of any seed with left-most base pair (i1,i2)
	 * @param i1 the left most interacting base of seq1
	 * @param i2 the left most interacting base of seq2
	 * @return the mfe of any seed starting at (i1,i2) or E_INF if none possible
	 */
	virtual
	E_type
	getSeedE( const size_t i1, const size_t i2 ) const;

	/**
	 * Checks whether or not a given base pair is the left-most base pair of
	 * any seed
	 * @param i1 the interacting base of seq1
	 * @param i2 the interacting base of seq2
	 * @return true if (i1,i2) is the left most base pair of some seed; false
	 *         otherwise
	 */
	virtual
	bool
	isSeedBound( const size_t i1, const size_t i2 ) const;

	/**
	 * Access to the length in seq1 of the mfe seed with left-most base pair (i1,i2)
	 * @param i1 the left most interacting base of seq1
	 * @param i2 the left most interacting base of seq2
	 * @return the length in seq1 of the mfe seed starting at (i1,i2) or 0 if none possible
	 */
	virtual
	size_t
	getSeedLength1( const size_t i1, const size_t i2 ) const;

	/**
	 * Access to the length in seq2 of the mfe seed with left-most base pair (i1,i2)
	 * @param i1 the left most interacting base of seq1
	 * @param i2 the left most interacting base of seq2
	 * @return the length in seq2 of the mfe seed starting at (i1,i2) or 0 if none possible
	 */
	virtual
	size_t
	getSeedLength2( const size_t i1, const size_t i2 ) const;



protected:

	//! the recursion data for the computation of a seed interaction
	//! i1..(i1+bpInbetween+u1-1) with i2..(i2+bpInbetween+u2-1)
	//! using the indexing [i1][i2][bpInbetween][u1][u2]
	SeedRecMatrix seedE_rec;

	//! the seed mfe information for seeds starting at (i1,i2)
	//! TODO replace with sparse data structure
	SeedMatrix seed;

	//! offset for seq1 indices for the current (restricted) matrices
	size_t offset1;

	//! offset for seq2 indices for the current (restricted) matrices
	size_t offset2;

	/**
	 * Provides the seed energy during recursion.
	 *
	 * @param i1 the seed left end in seq 1
	 * @param i2 the seed left end in seq 2
	 * @param bpInbetween the number of seed base pairs enclosed by the left-
	 *        and right-most base pair, ie. bpSeed-2
	 * @param u1 the number of unpaired bases within seq 1
	 * @param u2 the number of unpaired bases within seq 2
	 *
	 * @return the energy of the according (sub)seed
	 */
	E_type
	getSeedE( const size_t i1, const size_t i2, const size_t bpInbetween, const size_t u1, const size_t u2 ) const;

	/**
	 * Fills the seed energy during recursion.
	 *
	 * NOTE: internally a ring-list data structure is used which reuses memory
	 * instead of allocating mem for all possible parameter combinations. Thus,
	 * you have to call the method in appropriate order depending on your seed
	 * recursion.
	 *
	 * @param i1 the seed left end in seq 1
	 * @param i2 the seed left end in seq 2
	 * @param bpInbetween the number of seed base pairs enclosed by the left-
	 *        and right-most base pair, ie. bpSeed-2
	 * @param u1 the number of unpaired bases within seq 1
	 * @param u2 the number of unpaired bases within seq 2
	 * @param E the energy value to be set
	 */
	void
	setSeedE( const size_t i1, const size_t i2, const size_t bpInbetween, const size_t u1, const size_t u2, const E_type E );

	/**
	 * Encodes the seed lengths into one number
	 * @param l1 the length of the seed in seq1
	 * @param l2 the length of the seed in seq2
	 * @return the combined encoding = (l1 + l2*(max_l1+1))
	 */
	size_t
	encodeSeedLength( const size_t l1, const size_t l2 ) const;

	/**
	 * Decodes the length of the seed within sequence 1 from an encoding
	 * generated with encodeSeedLength()
	 * @param code the lengths encoding
	 * @return the length of the seed in seq1
	 */
	size_t
	decodeSeedLength1( const size_t code ) const;

	/**
	 * Decodes the length of the seed within sequence 2 from an encoding
	 * generated with encodeSeedLength()
	 * @param code the lengths encoding
	 * @return the length of the seed in seq2
	 */
	size_t
	decodeSeedLength2( const size_t code ) const;

	/**
	 * Fills a given interaction with the according
	 * hybridizing base pairs of the provided seed interaction (excluding
	 * the right-most seed base pair)
	 *
	 * @param interaction IN/OUT the interaction to fill
	 * @param i1 the seed left end in seq 1 (index including offset)
	 * @param i2 the seed left end in seq 2 (index including offset)
	 * @param bp the number of base pairs (bp+2) within the seed so far
	 * @param u1 the number of unpaired bases within seq 1
	 * @param u2 the number of unpaired bases within seq 2
	 */
	void
	traceBackSeed( Interaction & interaction
			, const size_t i1, const size_t i2, const size_t bp
			, const size_t u1, const size_t u2 ) const;

};


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


inline
SeedHandlerMfe::SeedHandlerMfe(
		const InteractionEnergy & energy
		, const SeedConstraint & seedConstraint
		)
	:
		SeedHandler(energy,seedConstraint)
		, seedE_rec( SeedIndex({{ 0,0,0,0,0 }}))
		, seed()
		, offset1(0)
		, offset2(0)
{
#if INTARNA_IN_DEBUG_MODE
	if ( ! seedConstraint.getExplicitSeeds().empty()) {
		LOG(WARNING) <<"explicit seeds definitions not supported by mfe-seed handler (and thus ignored)";
	}
#endif
}

////////////////////////////////////////////////////////////////////////////

inline
SeedHandlerMfe::~SeedHandlerMfe()
{
}

//////////////////////////////////////////////////////////////////////////

inline
void
SeedHandlerMfe::
traceBackSeed( Interaction & interaction
		, const size_t i1
		, const size_t i2
		) const
{
#if INTARNA_IN_DEBUG_MODE
	if ( i1 < offset1 ) throw std::runtime_error("SeedHandlerMfe::traceBackSeed(i1="+toString(i1)+") is out of range (>"+toString(offset1)+")");
	if ( i1-offset1 >= seed.size1() ) throw std::runtime_error("SeedHandlerMfe::traceBackSeed(i1="+toString(i1)+") is out of range (<"+toString(seed.size1()+offset1)+")");
	if ( i2 < offset2 ) throw std::runtime_error("SeedHandlerMfe::traceBackSeed(i2="+toString(i2)+") is out of range (>"+toString(offset2)+")");
	if ( i2-offset2 >= seed.size2() ) throw std::runtime_error("SeedHandlerMfe::traceBackSeed(i2="+toString(i2)+") is out of range (<"+toString(seed.size2()+offset2)+")");
	if ( !( isSeedBound(i1,i2) ) ) throw std::runtime_error("SeedHandlerMfe::traceBackSeed(i1="+toString(i1)+",i2="+toString(i2)+") no seed known (E_INF)");
	if ( i1+getSeedLength1(i1,i2)-1-offset1 >= seed.size1() ) throw std::runtime_error("SeedHandlerMfe::traceBackSeed(i1="+toString(i1)+") seed length ("+toString(getSeedLength1(i1,i2))+") exceeds of range (<"+toString(seed.size1()+offset1)+")");
	if ( i2+getSeedLength2(i1,i2)-1-offset2 >= seed.size2() ) throw std::runtime_error("SeedHandlerMfe::traceBackSeed(i2="+toString(i2)+") seed length ("+toString(getSeedLength2(i1,i2))+") exceeds of range (<"+toString(seed.size2()+offset2)+")");
#endif

	// get number of base pairs within the seed
	const size_t seedBps = getConstraint().getBasePairs();

	// trace back the according seed
	traceBackSeed( interaction, i1, i2
			, seedBps-2
			, getSeedLength1(i1,i2)-seedBps
			, getSeedLength2(i1,i2)-seedBps );
}

//////////////////////////////////////////////////////////////////////////

inline
E_type
SeedHandlerMfe::
getSeedE( const size_t i1, const size_t i2 ) const
{
	return seed(i1-offset1,i2-offset2).first;
}

//////////////////////////////////////////////////////////////////////////

inline
bool
SeedHandlerMfe::
isSeedBound( const size_t i1, const size_t i2 ) const
{
	return E_isNotINF( getSeedE(i1,i2) );
}

//////////////////////////////////////////////////////////////////////////

inline
size_t
SeedHandlerMfe::
getSeedLength1( const size_t i1, const size_t i2 ) const
{
	return decodeSeedLength1(seed(i1-offset1,i2-offset2).second);
}

//////////////////////////////////////////////////////////////////////////

inline
size_t
SeedHandlerMfe::
getSeedLength2( const size_t i1, const size_t i2 ) const
{
	return decodeSeedLength2(seed(i1-offset1,i2-offset2).second);
}

//////////////////////////////////////////////////////////////////////////

inline
E_type
SeedHandlerMfe::
getSeedE( const size_t i1, const size_t i2, const size_t bpInbetween, const size_t u1, const size_t u2 ) const
{
#if INTARNA_IN_DEBUG_MODE
	if ( i1 < offset1 ) throw std::runtime_error("SeedHandlerMfe::getSeedE(i1="+toString(i1)+") is out of range (>"+toString(offset1)+")");
	if ( i1-offset1 >= seed.size1() ) throw std::runtime_error("SeedHandlerMfe::getSeedE(i1="+toString(i1)+") is out of range (<"+toString(seed.size1()+offset1)+")");
	if ( i2 < offset2 ) throw std::runtime_error("SeedHandlerMfe::getSeedE(i2="+toString(i2)+") is out of range (>"+toString(offset2)+")");
	if ( i2-offset2 >= seed.size2() ) throw std::runtime_error("SeedHandlerMfe::getSeedE(i2="+toString(i2)+") is out of range (<"+toString(seed.size2()+offset2)+")");
#endif

	return seedE_rec( SeedIndex({{
		  (SeedRecMatrix::index) (i1-offset1)
		, (SeedRecMatrix::index) (i2-offset2)
		, (SeedRecMatrix::index) bpInbetween
		, (SeedRecMatrix::index) u1
		, (SeedRecMatrix::index) u2 }}) );
}

//////////////////////////////////////////////////////////////////////////

inline
void
SeedHandlerMfe::
setSeedE( const size_t i1, const size_t i2, const size_t bpInbetween, const size_t u1, const size_t u2, const E_type E )
{
#if INTARNA_IN_DEBUG_MODE
	if ( i1 < offset1 ) throw std::runtime_error("SeedHandlerMfe::setSeedE(i1="+toString(i1)+") is out of range (>"+toString(offset1)+")");
	if ( i1-offset1 >= seed.size1() ) throw std::runtime_error("SeedHandlerMfe::setSeedE(i1="+toString(i1)+") is out of range (<"+toString(seed.size1()+offset1)+")");
	if ( i2 < offset2 ) throw std::runtime_error("SeedHandlerMfe::setSeedE(i2="+toString(i2)+") is out of range (>"+toString(offset2)+")");
	if ( i2-offset2 >= seed.size2() ) throw std::runtime_error("SeedHandlerMfe::setSeedE(i2="+toString(i2)+") is out of range (<"+toString(seed.size2()+offset2)+")");
#endif

	seedE_rec( SeedIndex({{
		  (SeedRecMatrix::index) (i1-offset1)
		, (SeedRecMatrix::index) (i2-offset2)
		, (SeedRecMatrix::index) bpInbetween
		, (SeedRecMatrix::index) u1
		, (SeedRecMatrix::index) u2 }}) ) = E;
}

//////////////////////////////////////////////////////////////////////////

inline
size_t
SeedHandlerMfe::
encodeSeedLength( const size_t l1, const size_t l2 ) const
{
	return l1 + l2*(seedConstraint.getMaxLength1()+1);
}

//////////////////////////////////////////////////////////////////////////

inline
size_t
SeedHandlerMfe::
decodeSeedLength1( const size_t code ) const
{
	return code % (seedConstraint.getMaxLength1()+1);
}

//////////////////////////////////////////////////////////////////////////

inline
size_t
SeedHandlerMfe::
decodeSeedLength2( const size_t code ) const
{
	return code / (seedConstraint.getMaxLength1()+1);
}

//////////////////////////////////////////////////////////////////////////

} // namespace

#endif /* SEEDHANDLERMFE_H_ */
