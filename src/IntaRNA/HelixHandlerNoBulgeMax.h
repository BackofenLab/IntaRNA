
#ifndef INTARNA_HELIXHANDLERNOBULGEMAX_H_
#define INTARNA_HELIXHANDLERNOBULGEMAX_H_

#include "IntaRNA/InteractionEnergy.h"
#include "IntaRNA/HelixConstraint.h"
#include "IntaRNA/HelixHandler.h"

#include <boost/unordered_map.hpp>

namespace IntaRNA {

/**
 * Handler to provide maximal canonical helix information that is based on
 * stacked base pairs only.
 *
 * @author Martin Raden
 */
class HelixHandlerNoBulgeMax : public HelixHandler {


public:

	//! type used to store helix information (energy and length)
	typedef std::pair<E_type, size_t> HelixData;

	//! container type for sparse helix information for a given left-most
	//! base pair (i1,i2)
	//! it holds both the energy (first) as well as the length of the helix
	typedef boost::unordered_map< Interaction::BasePair, HelixData > HelixHash;

protected:

	//! shortening typedef
	typedef HelixHash::key_type BP;

public:

	/**
	 * Constructor
	 * @param energy the energy function to be used
	 */
	HelixHandlerNoBulgeMax(
			const InteractionEnergy & energy
			, const HelixConstraint & helixConstraint
			, SeedHandler * const seedHandler = NULL
	);

	/**
	 * destructor
	 */
	virtual ~HelixHandlerNoBulgeMax();

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
	 * @return number of valid helices
	 */
	virtual
	size_t
	fillHelix( const size_t i1, const size_t j1, const size_t i2, const size_t j2 );

	/**
	 * Compute the helix matrix, containing a seed, for the given interval boundaries
	 * @param i1 the first index of seq1 that might interact
	 * @param j1 the last index of seq1 that might interact
	 * @param i2 the first index of seq2 that might interact
	 * @param j2 the last index of seq2 that might interact
	 * @return
	 */
	virtual
	size_t
	fillHelixSeed( const size_t i1, const size_t j1, const size_t i2, const size_t j2 );

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
	 * Identifies the base pairs of the mfe helix interaction, containing a seed, starting at i1,i2
	 * and writes them to the provided container
	 *
	 * @param interaction the container to add the base pairs too
	 * @param i1 the start of the helix in seq1
	 * @param i2 the start of the helix in seq2
	 */
	virtual
	void
	traceBackHelixSeed( Interaction & interaction, const size_t i1, size_t i2);

	/**
	 * Access to the mfe of any helix with left-most base pair (i1, i2)
	 * @param i1 the left most interaction base of seq1
	 * @param i2 the left most interaction base of seq2
	 * @return the mfe of any helix starting at (i1,i2) or E_INF if none possible
	 * 			or return 0 if bp is 0 or 1 (used in helixHandlerSeed computation to avoid too many conditions)
	 */
	virtual
	E_type getHelixE( const size_t i1, const size_t i2 ) const;

	/**
	 * Access to the mfe of any helix, containing a seed, with left-most base pair (i1, i2)
	 * @param i1 the left most interaction base of seq1
	 * @param i2 the left most interaction base of seq2
	 * @return the mfe of any helix starting at (i1,i2) or E_INF if none possible
	 */
	virtual
	E_type getHelixSeedE( const size_t i1, const size_t i2 ) const;

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
	 * Access to the length in seq1 of the mfe helix, containing a seed, with left-most base pair (i1,i2)
	 * @param i1 the left most interacting base of seq1
	 * @param i2 the left most interacting base of seq2
	 * @return the length in seq1 of the mfe helix starting at (i1,i2) or 0 if none possible
	 */
	virtual
	size_t
	getHelixSeedLength1( const size_t i1, const size_t i2 ) const;

	/**
	 * Access to the length in seq2 of the mfe helix, containing a seed, with left-most base pair (i1,i2)
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
	 * Encodes the seed lengths into one number
	 * @param l1 the length of the seed in seq1
	 * @param l2 the length of the seed in seq2
	 * @return the combined encoding = (l1 + l2*(max_l1+1))
	 */
	size_t
	encodeHelixSeedLength( const size_t l1, const size_t l2 ) const;

	/**
	 * Decodes the length of the seed within sequence 1 from an encoding
	 * generated with encodeSeedLength()
	 * @param code the lengths encoding
	 * @return the length of the seed in seq1
	 */
	size_t
	decodeHelixSeedLength1( const size_t code ) const;

	/**
	 * Decodes the length of the seed within sequence 2 from an encoding
	 * generated with encodeSeedLength()
	 * @param code the lengths encoding
	 * @return the length of the seed in seq2
	 */
	size_t
	decodeHelixSeedLength2( const size_t code ) const;

protected:

	//! the used energy function
	const InteractionEnergy& energy;

	//! the helix constraint to be applied
	const HelixConstraint & helixConstraint;

	//! the helix mfe information for helix starting at (i1, i2)
	HelixHash helix;

	//! the helix mfe information for helix with seed starting at (i1, i2)
	HelixHash helixSeed;

	// seedHandler used in helixSeed computation
	SeedHandler * seedHandler;
};

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

inline
HelixHandlerNoBulgeMax::HelixHandlerNoBulgeMax(
		const InteractionEnergy & energy
		, const HelixConstraint & helixConstraint
		, SeedHandler * const seedHandler
)
		:
		energy(energy)
		, helixConstraint(helixConstraint)
		, helix()
		, helixSeed()
{
	if (seedHandler != NULL) {
		setSeedHandler( *seedHandler );
	}
}

////////////////////////////////////////////////////////////////////////////

inline
HelixHandlerNoBulgeMax::~HelixHandlerNoBulgeMax()
{
}

////////////////////////////////////////////////////////////////////////////

inline
const InteractionEnergy&
HelixHandlerNoBulgeMax::
getInteractionEnergy() const
{
	return energy;
}

////////////////////////////////////////////////////////////////////////////

inline
const HelixConstraint&
HelixHandlerNoBulgeMax::
getConstraint() const
{
	return helixConstraint;
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
HelixHandlerNoBulgeMax::
getHelixE(const size_t i1, const size_t i2) const
{
	// try to locate helix information
	HelixHash::const_iterator data = helix.find(BP(i1,i2));
	// return respective value
	return (data==helix.end()) ? E_INF : data->second.first;
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
HelixHandlerNoBulgeMax::
getHelixSeedE(const size_t i1, const size_t i2) const
{
	// try to locate helix information
	HelixHash::const_iterator data = helixSeed.find(BP(i1,i2));
	// return respective value
	return (data==helixSeed.end()) ? E_INF : data->second.first;
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerNoBulgeMax::
getHelixLength1(const size_t i1, const size_t i2) const
{
	// try to locate helix information
	HelixHash::const_iterator data = helix.find(BP(i1,i2));
	// return respective value
	return (data==helix.end()) ? 0 : data->second.second;
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerNoBulgeMax::
getHelixLength2(const size_t i1, const size_t i2) const
{
	// try to locate helix information
	HelixHash::const_iterator data = helix.find(BP(i1,i2));
	// return respective value
	return (data==helix.end()) ? 0 : data->second.second;
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerNoBulgeMax::
getHelixSeedLength1(const size_t i1, const size_t i2) const
{
	// try to locate helix information
	HelixHash::const_iterator data = helixSeed.find(BP(i1,i2));
	// return respective value
	return (data==helixSeed.end()) ? 0 : decodeHelixSeedLength1(data->second.second);
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerNoBulgeMax::
getHelixSeedLength2(const size_t i1, const size_t i2) const
{
	// try to locate helix information
	HelixHash::const_iterator data = helixSeed.find(BP(i1,i2));
	// return respective value
	return (data==helixSeed.end()) ? 0 : decodeHelixSeedLength2(data->second.second);
}

///////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerNoBulgeMax::
encodeHelixSeedLength( const size_t l1, const size_t l2 ) const
{
	return l1 + l2*(helixConstraint.getMaxLength1() + seedHandler->getConstraint().getMaxLength1()+1);
}

//////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerNoBulgeMax::
decodeHelixSeedLength1( const size_t code ) const
{
	return code % (helixConstraint.getMaxLength1() + seedHandler->getConstraint().getMaxLength1()+1);
}

//////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerNoBulgeMax::
decodeHelixSeedLength2( const size_t code ) const
{
	return code / (helixConstraint.getMaxLength1() + seedHandler->getConstraint().getMaxLength1()+1);
}

//////////////////////////////////////////////////////////////////////////

inline
void
HelixHandlerNoBulgeMax::setSeedHandler(SeedHandler & seedHandler) {
	this->seedHandler = &seedHandler;
}

//////////////////////////////////////////////////////////////////////////

} // namespace
#endif /* HELIXHANDLERNOBULGEMAX_H_ */
