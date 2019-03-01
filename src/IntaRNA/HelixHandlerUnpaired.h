
#ifndef INTARNA_HELIXHANDLERUNPAIRED_H_
#define INTARNA_HELIXHANDLERUNPAIRED_H_

#include "IntaRNA/InteractionEnergy.h"
#include "IntaRNA/HelixConstraint.h"
#include "IntaRNA/HelixHandler.h"

#include <boost/multi_array.hpp>

#include <boost/numeric/ublas/matrix.hpp>

namespace IntaRNA {

/**
 * Handler to provide helix interaction information that allows for minor
 * interior loops with a constraint maximal number of unpaired bases (per loop)
 *
 * @author Rick Gelhausen
 */
class HelixHandlerUnpaired : public HelixHandler {

public:

	//! 3D matrix type to hold the mfe energies for helix interactions
	//! of the ranges i1..(i1+bpMax-1) and i2..(i2+bpMax-1), with
	//! i1,i2 = the start index of the helix in seq1/2
	//! bp = the maximal number of base pairs within the helix (>=bpMin)
	//! using the index [i1][i2][bp] or a HelixIndex object
	typedef boost::multi_array<std::pair<E_type, size_t>,3> HelixRecMatrix;

	//! defines the helix data {{ i1, i2, bp }} to acces elements of
	//! the HelixRecMatrix
	typedef boost::array<HelixRecMatrix::index, 3> HelixIndex;

	//! matrix to store the helix information for each helix left side (i1, i2)
	//! it holds both the energy (first) as well as the length of the helix using
	//! the length combination of encodeHelixLength()
	//! The third entry is the bestBP, i.e. the optimal number of bases for this left boundary
	typedef boost::numeric::ublas::matrix< std::tuple<E_type, size_t, size_t> > HelixMatrix;
	typedef boost::numeric::ublas::matrix< std::pair<E_type, size_t> > HelixSeedMatrix;


public:

	/**
	 * Constructor
	 * @param energy the energy function to be used
	 */
	HelixHandlerUnpaired(
			const InteractionEnergy & energy
			, const HelixConstraint & helixConstraint
			, SeedHandler * const seedHandler = NULL
	);

	/**
	 * destructor
	 */
	virtual ~HelixHandlerUnpaired();

	/**
	 * Access to the underlying helix constraint
	 * @return the used helix constraint
	 */
	virtual
	const HelixConstraint&
	getConstraint() const;

	/**
	 * Access to the underlying interaction energy function
	 * @return the used energy function
	 */
	virtual
	const InteractionEnergy&
	getInteractionEnergy() const;

	/**
	 * Compute the helix matrix for the given interval boundaries
	 * @param i1 the first index of seq1 that might interact
	 * @param j1 the last index of seq1 that might interact
	 * @param i2 the first index of seq2 that might interact
	 * @param j2 the last index of seq2 that might interact
	 * @return
	 */
	virtual
	size_t
	fillHelix( const size_t i1, const size_t j1, const size_t i2, const size_t j2 );

	/**
	 * Computes the helix information for the given interval boundaries for helices containing a seed
	 * @param i1 the first index of seq1 that might interact
	 * @param j1 the last index of seq1 that might interact
	 * @param i2 the first index of seq2 that might interact
	 * @param j2 the last index of seq2 that might interact
	 * @return the number of potential helix interactions
	 */
	virtual
	size_t
	fillHelixSeed(const size_t i1, const size_t j1, const size_t i2, const size_t j2);

	/**
	 * Identifies the base pairs of the mfe helix interaction starting at i1,i2
	 * and writes them to the provided container
	 *
	 * Note: the right most base pair is excluded!
	 *
	 * @param interaction the container to add the base pairs too
	 * @param i1 the start of the helix in seq1
	 * @param i2 the start of the helix in seq2
	 */
	virtual
	void
	traceBackHelix( Interaction & interaction, const size_t i1, const size_t i2);

	/**
	 * Identifies the base pairs of the mfe helix interaction starting at i1,i2
	 * and writes them to the provided container
	 *
	 * Note: the right most base pair is excluded!
	 *
	 * @param interaction the container to add the base pairs too
	 * @param i1 the start of the helix in seq1
	 * @param i2 the start of the helix in seq2
	 */
	virtual
	void
	traceBackHelixSeed( Interaction & interaction, const size_t i1, const size_t i2);

	/**
	 * Access to the mfe of any helix with left-most base pair (i1, i2)
	 * @param i1 the left most interaction base of seq1
	 * @param i2 the left most interaction base of seq2
	 * @return the mfe of any helix starting at (i1,i2) or E_INF if none possible
	 */
	virtual
	E_type
	getHelixE( const size_t i1, const size_t i2 ) const;

	/**
     * Access to the mfe of any helix with left-most base pair (i1,i2)
	 * @param i1 the left most interacting base of seq1
	 * @param i2 the left most interaction base of seq2
	 * @return the mfe of any heli starting at (i1,i2) or E_INF if none possible
	 */
	virtual
	E_type
	getHelixSeedE( const size_t i1, const size_t i2 ) const;

	/**
	 * Access to the length in seq1 of the mfe helix with left-most base pair (i1,i2)
	 * @param i1 the left most interacting base of seq1
	 * @param i2 the left most interacting base of seq2
	 * @return the length in seq1 of the mfe helix starting at (i1,i2) or 0 if none possible
	 */
	virtual
	size_t
	getHelixLength1( const size_t i1, const size_t i2 ) const;

	/**
	 * Access to the length in seq2 of the mfe helix with left-most base pair (i1,i2)
	 * @param i1 the left most interacting base of seq1
	 * @param i2 the left most interacting base of seq2
	 * @return the length in seq2 of the mfe helix starting at (i1,i2) or 0 if none possible
	 */
	virtual
	size_t
	getHelixLength2( const size_t i1, const size_t i2 ) const;

	/**
	 * Access to the length in seq1 of the mfe helix with left-most base pair (i1,i2)
	 * @param i1 the left most interacting base of seq1
	 * @param i2 the left most interacting base of seq2
	 * @param bp the number of base pairs allowed
	 * @return the length in seq1 of the mfe helix starting at (i1,i2) or 0 if none possible
	 */
	virtual
	size_t
	getHelixLength1( const size_t i1, const size_t i2, const size_t bp ) const;

	/**
	 * Access to the length in seq2 of the mfe helix with left-most base pair (i1,i2)
	 * @param i1 the left most interacting base of seq1
	 * @param i2 the left most interacting base of seq2
	 * @param bp the number of base pairs allowed
	 * @return the length in seq2 of the mfe helix starting at (i1,i2) or 0 if none possible
	 */
	virtual
	size_t
	getHelixLength2( const size_t i1, const size_t i2, const size_t bp ) const;

	/**
	 * Access to the length in seq1 of the mfe helix with left-most base pair (i1,i2) containing a seed
	 * @param i1 the left most interacting base of seq1
	 * @param i2 the left most interacting base of seq2
	 * @return the length in seq1 of the mfe helix starting at (i1,i2) or 0 if none possible
	 */
	virtual
	size_t
	getHelixSeedLength1( const size_t i1, const size_t i2 ) const;

	/**
	 * Access to the length in seq2 of the mfe helix with left-most base pair (i1,i2) containing a seed
	 * @param i1 the left most interacting base of seq1
	 * @param i2 the left most interacting base of seq2
	 * @return the length in seq2 of the mfe helix starting at (i1,i2) or 0 if none possible
	 */
	virtual
	size_t
	getHelixSeedLength2( const size_t i1, const size_t i2 ) const;

	/**
	 * Set the seedHandler in order to compute helixSeed
	 * @param seedHandler seedHandler to be used in the helix computation
	 */
	void setSeedHandler(SeedHandler & seedHandler);

protected:

	/**
	 * @param i1 the helix left end in seq 1 (index including offset)
	 * @param i2 the helix left end in seq 2 (index including offset)
 	 * @return get optimal number of bases for this left boundary
	 */
	size_t
	getBestBP( const size_t i1, const size_t i2) const;

	/**
	 * Provides the helix energy during recursion
	 *
	 * @param i1 the helix left end in seq 1 (index including offset)
	 * @param i2 the helix left end in seq 2 (index including offset)
	 * @param bp the number of base pairs
	 * @param u1 the number of unpaired bases within seq 1
	 * @param u2 the number of unpaired bases within seq 2
	 *
	 * @return the energy of the according helix
	 */
	E_type
	getHelixE( const size_t i1, const size_t i2, const size_t bp);


	/**
	 * Fills the helix energy during recursion
	 *
	 * @param i1 the helix left end in seq 1 (index including offset)
	 * @param i2 the helix left end in seq 2 (index including offset)
	 * @param bp the number of base pairs
	 * @param u1 the number of unpaired bases within seq 1
	 * @param u2 the number of unpaired bases within seq 2
	 * @param E the energy value to be set
	 */
	void
	setHelixPair( const size_t i1, const size_t i2, const size_t bp
				, const E_type E, const size_t length );


	/**
	 * Encodes the helix lengths into one number
	 * @param l1 the length of the helix in seq1
	 * @param l2 the length of the helix in seq2
	 * @return the combined encoding = (l1 + l2*(max_l1+1))
	 */
	size_t
	encodeHelixLength( const size_t l1, const size_t l2 ) const;

	/**
	 * Decodes the length of the helix within sequence 1 from an encoding
	 * generated with encodeHelixLength()
	 * @param code the lengths encoding
	 * @return the length of the helix in seq1
	 */
	size_t
	decodeHelixLength1( const size_t code ) const;


	/**
	 * Decodes the length of the helix within sequence 2 from an encoding
	 * generated with encodeHelixLength()
	 * @param code the lengths encoding
	 * @return the length of the helix in seq2
	 */
	size_t
	decodeHelixLength2( const size_t code ) const;

	/**
	 * Encodes the seed lengths into one number
	 * @param l1 the length of the seed in seq1
	 * @param l2 the length of the seed in seq2
	 * @return the combined encoding = (l1 + l2*(max_l1+1))
	 */
	size_t
	encodeHelixSeedLength( const size_t l1, const size_t l2 ) const;

	/**
	 * Decodes the length of the seed within sequence 1 from an encoding
	 * generated with encodeHelixSeedLength()
	 * @param code the lengths encoding
	 * @return the length of the seed in seq1
	 */
	size_t
	decodeHelixSeedLength1( const size_t code ) const;

	/**
	 * Decodes the length of the seed within sequence 2 from an encoding
	 * generated with encodeHelixSeedLength()
	 * @param code the lengths encoding
	 * @return the length of the seed in seq2
	 */
	size_t
	decodeHelixSeedLength2( const size_t code ) const;

	/**
	 * Fills a given interaction with the according
	 * hybridizing base pairs of the provided helix interaction
	 *
	 * @param interaction IN/OUT the interaction to fill
	 * @param i1 the helix left end in seq 1 (index including offset)
	 * @param i2 the helix left end in seq 2 (index including offset)
	 * @param bp the number of base pairs (bestBP)
	 * @param u1 the number of unpaired bases within seq 1
	 * @param u2 the number of unpaired bases within seq 2
	 */
	void
	traceBackHelix( Interaction & interaction
			, const size_t i1, const size_t i2, const size_t bp);

protected:

	//! the used energy function
	const InteractionEnergy & energy;

	//! the helix constraint to be applied
	const HelixConstraint & helixConstraint;

	//! the recursion data for the computation of a helix interaction
	//! bp: the number of bases
	//! i1..(i1+bp-1) and i2..(i2+bp-1)
	//! using the indexing [i1][i2][bp]
	HelixRecMatrix helixE_rec;

	//! the helix mfe information for helix starting at (i1, i2)
	HelixMatrix helix;

	//! the helix mfe information for helix with seed starting at (i1, i2)
	HelixSeedMatrix helixSeed;

	//! offset for seq1 indices for the current matrices
	size_t offset1;

	//! offset for seq2 indices for the current matrices
	size_t offset2;

	//! used seedHandler
	SeedHandler * seedHandler;
};


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

inline
HelixHandlerUnpaired::HelixHandlerUnpaired(
		const InteractionEnergy & energy
		, const HelixConstraint & helixConstraint
		, SeedHandler * const seedHandler
)
		:
		energy(energy)
		, helixConstraint(helixConstraint)
		, seedHandler(seedHandler)
		, helixE_rec( HelixIndex({{ 0, 0, 0 }}))
		, helix()
		, helixSeed()
		, offset1(0)
		, offset2(0)
{
	if (seedHandler != NULL) {
		setSeedHandler(*seedHandler);
	}
}

////////////////////////////////////////////////////////////////////////////

inline
HelixHandlerUnpaired::~HelixHandlerUnpaired()
{
}

////////////////////////////////////////////////////////////////////////////

inline
const InteractionEnergy&
HelixHandlerUnpaired::
getInteractionEnergy() const
{
	return energy;
}

////////////////////////////////////////////////////////////////////////////

inline
const HelixConstraint&
HelixHandlerUnpaired::
getConstraint() const
{
	return helixConstraint;
}
////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerUnpaired::
getBestBP(const size_t i1, const size_t i2) const
{
	return std::get<2>(helix(i1-offset1, i2-offset2));
}

////////////////////////////////////////////////////////////////////////////

inline
void
HelixHandlerUnpaired::traceBackHelix(Interaction &interaction
		, const size_t i1
		, const size_t i2
)
{
#if INTARNA_IN_DEBUG_MODE
	if ( i1 < offset1 ) throw std::runtime_error("HelixHandlerUnpaired::traceBackHelix(i1="+toString(i1)+") is out of range (>"+toString(offset1)+")");
	if ( i1-offset1 >= helix.size1() ) throw std::runtime_error("HelixHandlerUnpaired::traceBackHelix(i1="+toString(i1)+") is out of range (<"+toString(helix.size1()+offset1)+")");
	if ( i2 < offset2 ) throw std::runtime_error("HelixHandlerUnpaired::traceBackHelix(i2="+toString(i2)+") is out of range (>"+toString(offset2)+")");
	if ( i2-offset2 >= helix.size2() ) throw std::runtime_error("HelixHandlerUnpaired::traceBackHelix(i2="+toString(i2)+") is out of range (<"+toString(helix.size2()+offset2)+")");
	if ( E_isINF( getHelixE(i1,i2) ) ) throw std::runtime_error("HelixHandlerUnpaired::traceBackHelix(i1="+toString(i1)+",i2="+toString(i2)+") no helix known (E_INF)");
	if ( i1+getHelixLength1(i1,i2)-1-offset1 >= helix.size1() ) throw std::runtime_error("HelixHandlerUnpaired::traceBackHelix(i1="+toString(i1)+") helix length ("+toString(getHelixLength1(i1,i2))+") exceeds of range (<"+toString(helix.size1()+offset1)+")");
	if ( i2+getHelixLength2(i1,i2)-1-offset2 >= helix.size2() ) throw std::runtime_error("HelixHandlerUnpaired::traceBackHelix(i2="+toString(i2)+") helix length ("+toString(getHelixLength2(i1,i2))+") exceeds of range (<"+toString(helix.size2()+offset2)+")");
#endif
	// get number of base pairs allowed within the helix
	const size_t bestBP = getBestBP(i1,i2);

	// trace back the according helix
	traceBackHelix( interaction, i1-offset1, i2-offset2, bestBP);
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
HelixHandlerUnpaired::
getHelixE(const size_t i1, const size_t i2) const
{
	return std::get<0>(helix(i1-offset1, i2-offset2));
}


////////////////////////////////////////////////////////////////////////////

inline
E_type
HelixHandlerUnpaired::
getHelixSeedE(const size_t i1, const size_t i2) const
{
	return helixSeed(i1-offset1, i2-offset2).first;
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerUnpaired::
getHelixLength1(const size_t i1, const size_t i2) const
{
	return decodeHelixLength1(std::get<1>(helix(i1-offset1, i2-offset2)));
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerUnpaired::
getHelixLength2(const size_t i1, const size_t i2) const
{
	return decodeHelixLength2(std::get<1>(helix(i1-offset1, i2-offset2)));
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerUnpaired::
getHelixLength1(const size_t i1, const size_t i2, const size_t bp) const {
// if no base pair is given return 0 in order to simplify helixHandlerSeed computation.
	if (bp <= 1) {
		return 0;
	} else {
		return decodeHelixLength1(helixE_rec(HelixIndex(
				{{(HelixRecMatrix::index) i1, (HelixRecMatrix::index) i2, (HelixRecMatrix::index) bp}})).second);
	}
}
////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerUnpaired::
getHelixLength2(const size_t i1, const size_t i2, const size_t bp) const
{
	// if no base pair is given return 0 in order to simplify helixHandlerSeed computation.
	if (bp <= 1) {
		return 0;
	} else {
		return decodeHelixLength2(helixE_rec(HelixIndex(
				{{(HelixRecMatrix::index) i1, (HelixRecMatrix::index) i2, (HelixRecMatrix::index) bp}})).second);
	}
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerUnpaired::
getHelixSeedLength1(const size_t i1, const size_t i2) const
{
	return decodeHelixSeedLength1(helixSeed(i1-offset1, i2-offset2).second);
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerUnpaired::
getHelixSeedLength2(const size_t i1, const size_t i2) const
{
	return decodeHelixSeedLength2(helixSeed(i1-offset1, i2-offset2).second);
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
HelixHandlerUnpaired::
getHelixE(const size_t i1, const size_t i2, const size_t bp)
{
	// if no base pair is given return 0 in order to simplify helixHandlerSeed computation.
	if (bp <= 1) {
		return 0;
	} else {
		return helixE_rec(HelixIndex({{(HelixRecMatrix::index) i1, (HelixRecMatrix::index) i2, (HelixRecMatrix::index) bp}})).first;
	}
}

////////////////////////////////////////////////////////////////////////////

inline
void
HelixHandlerUnpaired::
setHelixPair(const size_t i1, const size_t i2, const size_t bp, const E_type E, const size_t length)
{
	helixE_rec( HelixIndex({{
			(HelixRecMatrix::index) i1
			, (HelixRecMatrix::index) i2
			, (HelixRecMatrix::index) bp}}) ) = std::make_pair(E, E_isINF(E) ? 0 : length);
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerUnpaired::
encodeHelixLength(const size_t l1, const size_t l2) const
{
	return l1 + l2 * (helixConstraint.getMaxLength1()+1);
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerUnpaired::
decodeHelixLength1(const size_t code) const
{
	return code % (helixConstraint.getMaxLength1()+1);
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerUnpaired::
decodeHelixLength2(const size_t code) const
{
	return code / (helixConstraint.getMaxLength1()+1);
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerUnpaired::
encodeHelixSeedLength(const size_t l1, const size_t l2) const
{
	return l1 + l2 * (helixConstraint.getMaxLength1()+ seedHandler->getConstraint().getMaxLength1()+1);
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerUnpaired::
decodeHelixSeedLength1(const size_t code) const
{
	return code % (helixConstraint.getMaxLength1()+ seedHandler->getConstraint().getMaxLength1()+1);
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerUnpaired::
decodeHelixSeedLength2(const size_t code) const
{
	return code / (helixConstraint.getMaxLength1()+ seedHandler->getConstraint().getMaxLength1()+1);
}

////////////////////////////////////////////////////////////////////////////

inline
void
HelixHandlerUnpaired::
setSeedHandler(SeedHandler & seedHandler) {
	this->seedHandler = &seedHandler;
}



///////////////////////////////////////////////////////////////////////////

} // namespace

#endif //INTARNA_HELIXHANDLERUNPAIRED_H_
