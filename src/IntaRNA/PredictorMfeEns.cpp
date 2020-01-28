
#include "IntaRNA/PredictorMfeEns.h"

#include <iostream>
#include <algorithm>

namespace IntaRNA {

////////////////////////////////////////////////////////////////////////////

PredictorMfeEns::PredictorMfeEns(
		const InteractionEnergy & energy
		, OutputHandler & output
		, PredictionTracker * predTracker
		)
	: PredictorMfe(energy,output,predTracker)
{
}

////////////////////////////////////////////////////////////////////////////

PredictorMfeEns::~PredictorMfeEns()
{
}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns::
initZ()
{
	// reset storage
	Z_partition.clear();
}

////////////////////////////////////////////////////////////////////////////

const PredictorMfeEns::Site2Z_hash &
PredictorMfeEns::
getZPartition() const {
	return Z_partition;
}

////////////////////////////////////////////////////////////////////////////

Z_type
PredictorMfeEns::
getHybridZ( const size_t i1, const size_t j1
	 , const size_t i2, const size_t j2)
{
	Interaction::Boundary key(i1, j1, i2, j2);
	if ( Z_partition.find(key) == Z_partition.end() ) {
		return Z_type(0);
	} else {
		return Z_partition[key];
	}
}

////////////////////////////////////////////////////////////////////////////

Z_type
PredictorMfeEns::
getZ( const size_t i1, const size_t j1
	 , const size_t i2, const size_t j2)
{
	Interaction::Boundary key(i1, j1, i2, j2);
	if ( Z_partition.find(key) == Z_partition.end() ) {
		return Z_type(0);
	} else {
		return Z_partition[key] * energy.getBoltzmannWeight(energy.getE(i1,j1,i2,j2, E_type(0)));
	}
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns::
reportZ( SeedHandler* seedHandler )
{
	if (predTracker != NULL) {
		predTracker->updateZ(this, seedHandler);
	}
}

//////////////////////////////////////////////////////////////////////////

E_type
PredictorMfeEns::
getNonOverlappingEnergy( const size_t si1, const size_t si2, const size_t si1p, const size_t si2p
	                    , const InteractionEnergy & energy, const SeedHandler & seedHandler ) {

#if INTARNA_IN_DEBUG_MODE
	// check indices
	if( !seedHandler.isSeedBound(si1,si2) )
		throw std::runtime_error("PredictorMfeEns2dSeedExtension::getNonOverlappingEnergy( si "+toString(si1)+","+toString(si2)+",..) is no seed bound");
	if( !seedHandler.isSeedBound(si1p,si2p) )
		throw std::runtime_error("PredictorMfeEns2dSeedExtension::getNonOverlappingEnergy( sip "+toString(si1p)+","+toString(si2p)+",..) is no seed bound");
	if( si1 > si1p )
		throw std::runtime_error("PredictorMfeEns2dSeedExtension::getNonOverlappingEnergy( si "+toString(si1)+","+toString(si2)+", sip "+toString(si1p)+","+toString(si2p)+",..) si1 > sj1 !");
	// check if loop-overlapping (i.e. share at least one loop)
	if( !seedHandler.areLoopOverlapping(si1,si2,si1p,si2p) ) {
		throw std::runtime_error("PredictorMfeEns2dSeedExtension::getNonOverlappingEnergy( si "+toString(si1)+","+toString(si2)+", sip "+toString(si1p)+","+toString(si2p)+",..) are not loop overlapping");
	}
#endif

	// identity check
	if( si1 == si1p ) {
		return E_type(0);
	}

	// trace seed at (si1,si2)
	Interaction interaction = Interaction(energy.getAccessibility1().getSequence(), energy.getAccessibility2().getAccessibilityOrigin().getSequence());
	interaction.basePairs.push_back( energy.getBasePair(si1, si2) );
	seedHandler.traceBackSeed( interaction, si1, si2 );

	E_type fullE = 0;
	size_t k1old = si1, k2old = si2;
	for (size_t i = 1; i < interaction.basePairs.size(); i++) {
		// get index of current base pair
		size_t k1 = energy.getIndex1(interaction.basePairs[i]);
		// check if overlap done
		if (k1 > si1p) break;
		size_t k2 = energy.getIndex2(interaction.basePairs[i]);
		// add hybridization energy
		fullE += energy.getE_interLeft(k1old,k1,k2old,k2);

		// store
		k1old = k1;
		k2old = k2;
	}
	return fullE;
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns::
fillHybridZ_left( const size_t si1, const size_t si2, const InteractionEnergy & energy
	             , const SeedHandler & seedHandler, const OutputConstraint & outConstraint
							 , Z2dMatrix & hybridZ_left )
{
#if INTARNA_IN_DEBUG_MODE
	// check indices
	if (!energy.areComplementary(si1,si2) )
		throw std::runtime_error("PredictorMfeEns2dSeedExtension::fillHybridZ_left("+toString(si1)+","+toString(si2)+",..) are not complementary");
#endif

	// global vars to avoid reallocation
	size_t i1,i2,k1,k2;

	// determine whether or not lonely base pairs are allowed or if we have to
	// ensure a stacking to the right of the left boundary (i1,i2)
	const size_t noLpShift = outConstraint.noLP ? 1 : 0;
	Z_type iStackZ = Z_type(1);

	// iterate over all window starts i1 (seq1) and i2 (seq2)
	for (size_t l1=0; l1 < hybridZ_left.size1(); l1++) {
		for (size_t l2=0; l2 < hybridZ_left.size2(); l2++) {
			i1 = si1-l1;
			i2 = si2-l2;

			// referencing cell access
			Z_type & curZ = hybridZ_left(si1-i1,si2-i2);

			// init current cell (0 if not just right-most (j1,j2) base pair)
			curZ = (i1==si1 && i2==si2) ? energy.getBoltzmannWeight(energy.getE_init()) : 0.0;

			// check if complementary (use global sequence indexing)
			if( i1<si1 && i2<si2 && energy.areComplementary(i1,i2) ) {

				// right-stacking of i if no-LP
				if (outConstraint.noLP) {
					// skip if no stacking possible
					if (!energy.areComplementary(i1+noLpShift,i2+noLpShift)) {
						continue;
					}
					// get stacking energy to avoid recomputation in recursion below
					iStackZ = energy.getBoltzmannWeight(energy.getE_interLeft(i1,i1+noLpShift,i2,i2+noLpShift));
					// check just stacked
					curZ += iStackZ + hybridZ_left(l1-noLpShift,l2-noLpShift);
				}

				// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
				for (k1=i1+noLpShift; k1++ < si1; ) {
					// ensure maximal loop length
					if (k1-i1-noLpShift > energy.getMaxInternalLoopSize1()+1) break;
					for (k2=i2+noLpShift; k2++ < si2; ) {
						// ensure maximal loop length
						if (k2-i2-noLpShift > energy.getMaxInternalLoopSize2()+1) break;
						// check if (k1,k2) are valid left boundary
						if ( ! Z_equal(hybridZ_left(si1-k1,si2-k2), 0.0) ) {
							curZ += (iStackZ
									* energy.getBoltzmannWeight(energy.getE_interLeft(i1+noLpShift,k1,i2+noLpShift,k2))
									* hybridZ_left(si1-k1,si2-k2));
						}
					} // k2
				} // k1

				// correction for left seeds
				if ( i1<si1 && i2<si2 && seedHandler.isSeedBound(i1, i2) ) {

					// check if seed is to be processed:
					bool substractThisSeed =
							// check if left of anchor seed
										( i1+seedHandler.getSeedLength1(i1,i2)-1 <= si1
										&& i2+seedHandler.getSeedLength2(i1,i2)-1 <= si2 )
							// check if overlapping with anchor seed
									||	seedHandler.areLoopOverlapping(i1,i2,si1,si2);

					if (substractThisSeed) {

						// iterate seeds in S region
						size_t sj1 = RnaSequence::lastPos, sj2 = RnaSequence::lastPos;
						size_t si1overlap = si1+1, si2overlap = si2+1;
						// find left-most loop-overlapping seed
						while( seedHandler.updateToNextSeed(sj1,sj2
								, i1, std::min(si1,i1+seedHandler.getSeedLength1(i1, i2)-2)
								, i2, std::min(si2,i2+seedHandler.getSeedLength2(i1, i2)-2)) )
						{
							// check if right of i1,i2 and overlapping
							if (sj1 > i1 && seedHandler.areLoopOverlapping(i1, i2, sj1, sj2)) {
								// update left-most loop-overlapping seed
								if (sj1 < si1overlap) {
									si1overlap = sj1;
									si2overlap = sj2;
								}
							}
						}

						// if we found an overlapping seed
						if (si1overlap <= si1) {
							// check if right side is non-empty
							if ( ! E_equal(hybridZ_left( si1-si1overlap, si2-si2overlap),0) ) {
								// compute Energy of loop S \ S'
								E_type nonOverlapE = getNonOverlappingEnergy(i1, i2, si1overlap, si2overlap, energy, seedHandler);
								// subtract energy.getBoltzmannWeight( nonOverlapE ) * hybridZ_left( si1overlap, si2overlap up to anchor seed [==1 if equal])
								Z_type correctionTerm = energy.getBoltzmannWeight( nonOverlapE )
														* hybridZ_left( si1-si1overlap, si2-si2overlap);
								curZ -= correctionTerm;
							// sanity insurance
								if (curZ < 0) {
									curZ = Z_type(0.0);
								}
							}
						} else {
							// get data for seed to be removed
							const Z_type seedZ_rm = energy.getBoltzmannWeight(seedHandler.getSeedE(i1, i2));
							const size_t sj1_rm = i1+seedHandler.getSeedLength1(i1,i2)-1;
							const size_t sj2_rm = i2+seedHandler.getSeedLength2(i1,i2)-1;
							// if no S'
							// substract seedZ * hybridZ_left(right end seed up to anchor seed)
							Z_type correctionTerm = seedZ_rm
									* hybridZ_left( si1-sj1_rm
												  , si2-sj2_rm );
							// if noLP : handle explicit loop right of current seed
							if (outConstraint.noLP) {
								for (k1=sj1_rm; k1++ < si1; ) {
									// ensure maximal loop length
									if (k1-sj1_rm > energy.getMaxInternalLoopSize1()+1) break;
									for (k2=sj2_rm; k2++ < si2; ) {
										// ensure at least one unpaired base in interior loop following the seed to be removed
										if (sj1_rm-k1 + sj2_rm-k2 == 2) {continue;}
										// ensure maximal loop length
										if (k2-sj2_rm > energy.getMaxInternalLoopSize2()+1) break;
										// check if (k1,k2) are valid left boundary
										if ( ! Z_equal(hybridZ_left(si1-k1,si2-k2), 0.0) ) {
											correctionTerm += (seedZ_rm
													* energy.getBoltzmannWeight(energy.getE_interLeft(sj1_rm,k1,sj2_rm,k2))
													* hybridZ_left(si1-k1,si2-k2) );
										}
									} // k2
								} // k1

							}
							curZ -= correctionTerm;
							// sanity ensurence
							if (curZ < 0) {
								curZ = Z_type(0.0);
							}
						}
					} // substractThisSeed
				}

			} // complementary

		} // i2
	} // i1

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns::
updateZ( const size_t i1, const size_t j1
		, const size_t i2, const size_t j2
		, const Z_type partZ
		, const bool isHybridZ )
{
	// check if something to be done
	if (Z_equal(partZ,0) || Z_isINF(Zall))
		return;
	// update overall partition function
	if (isHybridZ) {
#if INTARNA_IN_DEBUG_MODE
		if ( (std::numeric_limits<Z_type>::max() - (partZ*energy.getBoltzmannWeight(energy.getE(i1,j1,i2,j2, E_type(0))))) <= Zall) {
			LOG(WARNING) <<"PredictorMfeEns::updateZ() : partition function overflow! Recompile with larger partition function data type!";
		}
#endif
		// add ED penalties etc.
		Zall += partZ*energy.getBoltzmannWeight(energy.getE(i1,j1,i2,j2, E_type(0)));
	} else {
#if INTARNA_IN_DEBUG_MODE
		if ( (std::numeric_limits<Z_type>::max() - partZ) <= Zall) {
			LOG(WARNING) <<"PredictorMfeEns::updateZ() : partition function overflow! Recompile with larger partition function data type!";
		}
#endif
		// just increase
		Zall += partZ;
	}

	// store partial Z
	Interaction::Boundary key(i1,j1,i2,j2);
	auto keyEntry = Z_partition.find(key);
	if ( Z_partition.find(key) == Z_partition.end() ) {
		Z_partition[key] = partZ;
	} else {
		// update entry
		keyEntry->second += partZ;
	}

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns::
updateOptimaUsingZ()
{
	for (auto it = Z_partition.begin(); it != Z_partition.end(); ++it)
	{
		// if partition function is > 0
		if (Z_isNotINF(it->second) && it->second > 0) {
			PredictorMfe::updateOptima( it->first.i1, it->first.j1, it->first.i2, it->first.j2, energy.getE(it->second), true, false );
		}
	}
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns::
reportOptima()
{
	// update optima from Z information
	updateOptimaUsingZ();

	// call super-class function
	PredictorMfe::reportOptima();
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfeEns::
traceBack( Interaction & interaction )
{
	// temporary access
	const OutputConstraint & outConstraint = output.getOutputConstraint();
	// check if something to trace
	if (interaction.basePairs.size() < 2) {
		return;
	}

#if INTARNA_IN_DEBUG_MODE
	// sanity checks
	if ( interaction.basePairs.size() != 2 ) {
		throw std::runtime_error("PredictorMfeEns::traceBack() : given interaction does not contain boundaries only");
	}
#endif

	// ensure sorting
	interaction.sort();

	// check for single base pair interaction
	if (interaction.basePairs.at(0).first == interaction.basePairs.at(1).first) {
		// delete second boundary (identical to first)
		interaction.basePairs.resize(1);
		// update done
		return;
	}

#if INTARNA_IN_DEBUG_MODE
	// sanity checks
	if ( ! interaction.isValid() ) {
		throw std::runtime_error("PredictorMfeEns2d::traceBack() : given interaction is not valid");
	}
#endif

}

////////////////////////////////////////////////////////////////////////////


} // namespace
