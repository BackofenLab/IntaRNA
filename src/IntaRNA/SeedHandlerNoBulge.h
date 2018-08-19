
#ifndef INTARNA_SEEDHANDLERNOBULGE_H_
#define INTARNA_SEEDHANDLERNOBULGE_H_

#include "IntaRNA/InteractionEnergy.h"
#include "IntaRNA/SeedConstraint.h"
#include "IntaRNA/SeedHandler.h"

#include <boost/unordered_map.hpp>


namespace IntaRNA {

/**
 * Handler to provide seed interaction information for each intermolecular
 * index combination (= left end of seed) for perfectly stacked seeds only, ie.
 * seed contain no bulges.
 *
 */
class SeedHandlerNoBulge : public SeedHandler
{

protected:

	//! container to store stacking energies of seed base pairs
	typedef std::vector<E_type> StackingEnergyList;

	//! container type for sparse seed information
	typedef boost::unordered_map< Interaction::BasePair, E_type > SeedHash;

	//! container to store seeds' hybridization energies
	SeedHash seedForLeftEnd;


public:

	/**
	 * Construction
	 * @param energy the energy function to be used for seed prediction
	 * @param seedConstraint the seed constraint to be applied
	 */
	SeedHandlerNoBulge(
			const InteractionEnergy & energy
			, const SeedConstraint & seedConstraint
			);

	/**
	 * destruction
	 */
	virtual ~SeedHandlerNoBulge();

	/**
	 * Computes the seed information for the given interval boundaries
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
	traceBackSeed( Interaction & interaction, const size_t i1, const size_t i2);


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

	/**
	 * Checks and stores a given seed candidate
	 *
	 * @param j1 right-most seed position in seq1
	 * @param j2 right-most seed position in seq2
	 * @param bpE list of stacking energies for this seed
	 *
	 */
	void
	storeSeed( const size_t j1, const size_t j2, const StackingEnergyList & bpE );

};


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


inline
SeedHandlerNoBulge::SeedHandlerNoBulge(
		const InteractionEnergy & energy
		, const SeedConstraint & seedConstraint
		)
	:
		SeedHandler(energy,seedConstraint)
		, seedForLeftEnd()
{
#if INTARNA_IN_DEBUG_MODE
	if ( ! seedConstraint.getExplicitSeeds().empty()) {
		LOG(WARNING) <<"explicit seeds definitions not supported by noBulge-seed handler (and thus ignored)";
	}
	if ( seedConstraint.getMaxUnpaired1() > 0 ) { throw std::runtime_error("SeedHandlerNoBulge() : seedConstraint.maxUnpaired1 > 0"); }
	if ( seedConstraint.getMaxUnpaired2() > 0 ) { throw std::runtime_error("SeedHandlerNoBulge() : seedConstraint.maxUnpaired2 > 0"); }
#endif
}

////////////////////////////////////////////////////////////////////////////

inline
SeedHandlerNoBulge::~SeedHandlerNoBulge()
{
	seedForLeftEnd.clear();
}

//////////////////////////////////////////////////////////////////////////

inline
void
SeedHandlerNoBulge::
traceBackSeed( Interaction & interaction
		, const size_t i1
		, const size_t i2
		)
{
#if INTARNA_IN_DEBUG_MODE
	if ( E_isINF( getSeedE(i1,i2) ) ) throw std::runtime_error("SeedHandlerNoBulge::traceBackSeed(i1="+toString(i1)+",i2="+toString(i2)+") no seed known (E_INF)");
#endif

	// get number of base pairs within the seed
	const size_t seedBps = getConstraint().getBasePairs();

	// add seed base pairs
	// but exclude left- and right-most bp
	for (size_t k=1; (k+1)<seedBps; k++) {
		interaction.basePairs.push_back( energy.getBasePair(i1+k,i2+k) );
	}

}

//////////////////////////////////////////////////////////////////////////

inline
E_type
SeedHandlerNoBulge::
getSeedE( const size_t i1, const size_t i2 ) const
{
	// search for seed entry in hash
	SeedHash::const_iterator seedData = seedForLeftEnd.find( SeedHash::key_type(i1,i2) );
	// return according energy value
	return seedData == seedForLeftEnd.end() ? E_INF : seedData->second;
}

//////////////////////////////////////////////////////////////////////////

inline
size_t
SeedHandlerNoBulge::
getSeedLength1( const size_t i1, const size_t i2 ) const
{
#if INTARNA_IN_DEBUG_MODE
	if ( E_isINF( getSeedE(i1,i2) ) ) throw std::runtime_error("SeedHandlerNoBulge::getSeedLength1(i1="+toString(i1)+",i2="+toString(i2)+") no seed known (E_INF)");
#endif
	// check if
	return getConstraint().getBasePairs();
}

//////////////////////////////////////////////////////////////////////////

inline
size_t
SeedHandlerNoBulge::
getSeedLength2( const size_t i1, const size_t i2 ) const
{
#if INTARNA_IN_DEBUG_MODE
	if ( E_isINF( getSeedE(i1,i2) ) ) throw std::runtime_error("SeedHandlerNoBulge::getSeedLength2(i1="+toString(i1)+",i2="+toString(i2)+") no seed known (E_INF)");
#endif
	// check if
	return getConstraint().getBasePairs();
}

//////////////////////////////////////////////////////////////////////////

inline
void
SeedHandlerNoBulge::
storeSeed( const size_t j1, const size_t j2, const StackingEnergyList & bpE )
{
	const size_t seedBP = seedConstraint.getBasePairs();
	E_type seedEhybrid = 0;
//	bool stored = false;
//	std::stringstream tmp;
//	for(auto e=bpE.begin(); e!=bpE.end(); e++) { tmp <<" + "<<*e; }

	// check if EDs of full seed are within boundaries
	if (seedConstraint.getMaxED() >= energy.getED1(j1+1-seedBP, j1)
			&& seedConstraint.getMaxED() >= energy.getED2(j2+1-seedBP, j2) )
	{
		// compute seed energy
		for(auto e=bpE.begin(); e!=bpE.end(); e++) { seedEhybrid += *e; } // (left) stacking energies
		E_type seedEfull = energy.getE(j1+1-seedBP, j1, j2+1-seedBP, j2, seedEhybrid) + energy.getE_init(); // ED values etc.

		// store seed hybridization energy if overall E is below or equal to threshold
		if( seedEfull <= seedConstraint.getMaxE() ) {
//	if (E_isNotINF(getSeedE(j1+1-seedBP,j2+1-seedBP))) {
//		LOG(DEBUG)<<" seed duplicate for j "<<j1<<"-"<<j2<<" : ";
//		LOG(DEBUG)<<getSeedE(j1+1-seedBP,j2+1-seedBP) <<" <- "<<seedEhybrid<<" = "<<tmp.str();
//	}
//	stored = true;
			seedForLeftEnd[ SeedHash::key_type(j1+1-seedBP,j2+1-seedBP) ] = seedEhybrid ;
		}
	}

//	LOG(DEBUG)<<"storing "<<(j1+1-seedBP)<<"-"<<(j2+1-seedBP) <<" = "<<seedEhybrid <<" -> "<<(stored?"true":"false") <<" = "<<tmp.str();

}

//////////////////////////////////////////////////////////////////////////

} // namespace

#endif /* SEEDHANDLERNOBULGE_H_ */
