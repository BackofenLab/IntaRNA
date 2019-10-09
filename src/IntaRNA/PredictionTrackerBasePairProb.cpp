
#include "IntaRNA/PredictionTrackerBasePairProb.h"

extern "C" {
	#include <ViennaRNA/vrna_config.h>
	#include <ViennaRNA/plotting/probabilities.h>
}

#include <boost/foreach.hpp>
#include <boost/regex.hpp>

namespace IntaRNA {

//////////////////////////////////////////////////////////////////////

PredictionTrackerBasePairProb::
PredictionTrackerBasePairProb(
		const InteractionEnergy & energy
		, const std::string & fileName
	)
 :	PredictionTracker()
	, energy(energy)
	, fileName(fileName)
	, probabilityThreshold(0.0001)
{
}

//////////////////////////////////////////////////////////////////////

PredictionTrackerBasePairProb::
~PredictionTrackerBasePairProb()
{
}

//////////////////////////////////////////////////////////////////////

void
PredictionTrackerBasePairProb::
updateOptimumCalled( const size_t i1, const size_t j1
					, const size_t i2, const size_t j2
					, const E_type curE
					)
{
}

//////////////////////////////////////////////////////////////////////

void
PredictionTrackerBasePairProb::
updateZ( PredictorMfeEns *predictor, SeedHandler *seedHandler )
{

	// sequence strings
	const std::string & rna1 = energy.getAccessibility1().getSequence().asString();
	const std::string & reverseRna2 = energy.getAccessibility2().getAccessibilityOrigin().getSequence().asString();

	size_t s1 = energy.size1();
	size_t s2 = energy.size2();
	size_t n1 = energy.getAccessibility1().getMaxLength();
	size_t n2 = energy.getAccessibility2().getMaxLength();

	Z_type maxZ = 0.0;
	Interaction::Boundary interactionBoundary;

	// initialize Z_partition
	const PredictorMfeEns::Site2Z_hash & Z_partition = predictor->getZPartition();
	// cleanup additional data
	Z_partitionMissing.clear();

	// create index of left/right boundaries
	for (auto it = Z_partition.begin(); it != Z_partition.end(); ++it) {
		// identify best interaction boundary
		Z_type Zstruct = it->second * energy.getBoltzmannWeight(energy.getE(it->first.i1, it->first.j1, it->first.i2, it->first.j2, E_type(0)));
		if (Zstruct > maxZ) {
			maxZ = Zstruct;
			interactionBoundary = it->first;
		}

		// create left and right index
		if (it->first.i1 != it->first.j1 && it->first.i2 != it->first.j2) {
			// encode left/right boundary
			Interaction::BasePair left(it->first.i1, it->first.i2);
			Interaction::BasePair right(it->first.j1, it->first.j2);
			// create left index
			rightExt[left].push_back(right);
			// create right index
			if (seedHandler != NULL) {
				leftExt[right].push_back(left);
			}
		}

	} // it (Z_partition)

	// compute all parts outside of seeds
	for (auto it = Z_partition.begin(); it != Z_partition.end(); ++it) {
		Interaction::BasePair leftEnd(it->first.i1, it->first.i2);
		for (auto k = rightExt[leftEnd].begin(); k != rightExt[leftEnd].end(); ++k) {
			if (k->first < it->first.j1 && k->second < it->first.j2) {
				computeMissingZ(	it->first.i1, k->first, it->first.i2, k->second
									, it->first.i1, it->first.j1, it->first.i2, it->first.j2
									, predictor, seedHandler);
			}
		}
		Interaction::BasePair rightEnd(it->first.j1, it->first.j2);
		for (auto k = leftExt[rightEnd].begin(); k != leftExt[rightEnd].end(); ++k) {
			if (k->first > it->first.i1 && k->second > it->first.i2) {
				computeMissingZ(	it->first.i1, k->first, it->first.i2, k->second
									, it->first.i1, it->first.j1, it->first.i2, it->first.j2
									, predictor, seedHandler);
			}
		}
	} // it (Z_partition)

	for (auto z = Z_partitionMissing.begin(); z != Z_partitionMissing.end(); z++) {
		LOG(DEBUG) <<" newZ( "<<z->first.i1<<":"<<z->first.i2<<", "<<z->first.j1<<":"<<z->first.j2 <<" ) = "<<z->second;
	}

	// compute missing Z values for seed-based predictions
	if (seedHandler != NULL
			&& ( ! seedHandler->getConstraint().getExplicitSeeds().empty()
					|| seedHandler->getConstraint().getBasePairs() >= 2) )
	{
		// iterate all seeds
		size_t si1 = RnaSequence::lastPos, si2 = RnaSequence::lastPos;
		while( seedHandler->updateToNextSeed(si1, si2) )
		{
			// compute missing Z values for seed base pairs

			const size_t sl1 = seedHandler->getSeedLength1(si1, si2);
			const size_t sl2 = seedHandler->getSeedLength2(si1, si2);
			const size_t sj1 = si1+sl1-1;
			const size_t sj2 = si2+sl2-1;
			const Interaction::BasePair siBP(si1, si2);
			const Interaction::BasePair sjBP(sj1, sj2);

			LOG(DEBUG)<<" internal of seed "<<si1<<":"<<si2;

			// trace all inner seed base pairs
			Interaction interaction = Interaction(energy.getAccessibility1().getSequence(), energy.getAccessibility2().getAccessibilityOrigin().getSequence());
			seedHandler->traceBackSeed( interaction, si1, si2 );

			// iterate all internal seed base pairs
			for (size_t i = 0; i < interaction.basePairs.size(); i++) {
				// get index of current base pair
				const size_t sk1 = energy.getIndex1(interaction.basePairs[i]);
				const size_t sk2 = energy.getIndex2(interaction.basePairs[i]);
				const Interaction::BasePair skBP(sk1, sk2);
			LOG(DEBUG)<<"   >> "<<sk1<<":"<<sk2<<"  Z "<<getHybridZ(sk1, sj1, sk2, sj2, predictor);

				// no extension
//				if (Z_equal(getHybridZ(sk1, sj1, sk2, sj2, predictor), 0)
//					|| Z_equal(getHybridZ(si1,sk1, si2, sk2, predictor), 0))
//				{
//			LOG(DEBUG)<<"  missing "<<sk1<<":"<<sk2<<" j "<<sj1<<":"<<sj2;
//				}
				computeMissingZ(sk1, sj1, sk2, sj2, si1, sj1, si2, sj2, predictor, seedHandler);
				// update right extension index
				if (std::find(rightExt[skBP].begin(), rightExt[skBP].end(), sjBP) == rightExt[skBP].end()) {
					rightExt[skBP].push_back(sjBP);
				}

				// right extensions
				// iterate all r=rightExt[si] to computeMissingZ(si,sk,r) and update rightExt[sk].push(r)
				for (auto seedRightExt = rightExt[siBP].begin(); seedRightExt != rightExt[siBP].end(); ++seedRightExt) {
					if (sk1 < seedRightExt->first && sk2 < seedRightExt->second
//						&& Z_equal(getHybridZ(sk1, seedRightExt->first, sk2, seedRightExt->second, predictor), 0)
					)
					{
						computeMissingZ(sk1, seedRightExt->first, sk2, seedRightExt->second, si1, seedRightExt->first, si2, seedRightExt->second, predictor, seedHandler);
						// update right extension index
						if (std::find(rightExt[skBP].begin(), rightExt[skBP].end(), *seedRightExt) == rightExt[skBP].end()) {
							rightExt[skBP].push_back(*seedRightExt);
						}
					}
				}

				// left extensions
				// iterate all l=leftEx[sj] to computeMissingZ(l,sk,sj)
				for (auto seedLeftExt = leftExt[sjBP].begin(); seedLeftExt != leftExt[sjBP].end(); ++seedLeftExt) {
					if (seedLeftExt->first < sk1 && seedLeftExt->second < sk2
//						&& Z_equal(getHybridZ(seedLeftExt->first, sk1, seedLeftExt->second, sk2, predictor), 0)
						)
					{
						computeMissingZ(seedLeftExt->first, sk1, seedLeftExt->second, sk2, seedLeftExt->first, sj1, seedLeftExt->second, sj2, predictor, seedHandler);
					}
				}

			}

		}

	}

	// Compute base-pair probabilities outside seeds
	computeBasePairProbs(predictor, Z_partition.begin(), Z_partition.end());
	// Compute base-pair probabilities within seeds
	computeBasePairProbs(predictor, Z_partitionMissing.begin(), Z_partitionMissing.end());

	// build plist
	struct vrna_elem_prob_s plist[structureProbs.size()+1];
	size_t i = 0;
	for (auto it = structureProbs.begin(); it != structureProbs.end(); ++it) {
		LOG(DEBUG) << "prob: " << it->first.first << ":" << it->first.second << ":" << it->first.first << ":" << it->first.second << " = " << it->second;
		if (it->second > probabilityThreshold) {
			plist[i].i = it->first.first + 1;
			plist[i].j = it->first.second + 1;
			plist[i].p = it->second;
			plist[i].type = 0; // base-pair prob
			i++;
		}
	}

	// create dot plot
	char *name = strdup(fileName.c_str());
	std::string comment =
	  "Intermolecular base-pair probabilities generated by "
	  INTARNA_PACKAGE_STRING
    " using Vienna RNA package "
    VRNA_VERSION;

	generateDotPlot(strdup(rna1.c_str()), strdup(reverseRna2.c_str()), name, plist, comment.c_str(), interactionBoundary);

}

////////////////////////////////////////////////////////////////////////////

void
PredictionTrackerBasePairProb::
computeBasePairProbs( PredictorMfeEns *predictor
										, const PredictorMfeEns::Site2Z_hash::const_iterator first
										, const PredictorMfeEns::Site2Z_hash::const_iterator last )
{
	const Z_type Zinit = energy.getBoltzmannWeight(energy.getE_init());

	for (auto it = first; it != last; ++it) {
		Z_type bpProb = it->second * energy.getBoltzmannWeight(energy.getE(it->first.i1, it->first.j1, it->first.i2, it->first.j2, E_type(0)));
//		LOG_IF((it->first.j1 == 2 && it->first.j2 == 2), DEBUG)
		LOG(DEBUG)
				<< " i: " << it->first.i1 << ":" << it->first.i2
				<< " j: " << it->first.j1 << ":" << it->first.j2
				;

		Interaction::BasePair iBP(it->first.i1, it->first.i2);
		Interaction::BasePair jBP(it->first.j1, it->first.j2);

		// init
		if (structureProbs.find(jBP) == structureProbs.end()) {
			structureProbs[jBP] = 0;
		}
		if (structureProbs.find(iBP) == structureProbs.end()) {
			structureProbs[iBP] = 0;
		}

		// no extension
		if (predictor->getZPartition().find( it->first ) != predictor->getZPartition().end()) {
			structureProbs[jBP] += (1 / predictor->getZall()) * bpProb;
			if (iBP != jBP) {
				structureProbs[iBP] += (1 / predictor->getZall()) * bpProb;
				LOG(DEBUG) <<"("<<iBP.first<<":"<<iBP.second<<")" << " - no ext: " << bpProb;
			}
//			if (it->first.j1 == 2 && it->first.j2 == 2)
				LOG(DEBUG) <<"("<<jBP.first<<":"<<jBP.second<<")" << " - no ext: " << bpProb;
		}

		assert( !Z_equal(it->second,0) );

		// extensions
		if (iBP != jBP) {
		for (auto it2 = rightExt[jBP].begin(); it2 != rightExt[jBP].end(); ++it2) {

			LOG(DEBUG) <<" rightExt("<<jBP.first<<":"<<jBP.second<<") = "<<it2->first<<":"<<it2->second;

			assert( !Z_equal(getHybridZ(it->first.j1, it2->first, it->first.j2, it2->second, predictor),0) );

			// ensure extension is valid (present in original Z data)
			if (predictor->getZPartition().find( Interaction::Boundary(it->first.i1, it2->first, it->first.i2, it2->second)) != predictor->getZPartition().end()) {

				assert( !Z_equal(it->second, 0) );
				bpProb = it->second * getHybridZ(it->first.j1, it2->first, it->first.j2, it2->second, predictor)
							// ED penalty
									 * energy.getBoltzmannWeight(energy.getE(it->first.i1, it2->first, it->first.i2, it2->second, E_type(0)))
									 / Zinit;

				structureProbs[jBP] += (1 / predictor->getZall()) * bpProb;
//				LOG_IF((it->first.j1 == 1 && it->first.j2 == 1), DEBUG)
				LOG(DEBUG)
				<<"("<<jBP.first<<":"<<jBP.second<<")" << " - ext at: " << it2->first << ":" << it2->second << "=" << bpProb << " | left: " << it->second << " | right = " << getHybridZ(it->first.j1, it2->first, it->first.j2, it2->second, predictor);
			}
		}
		} // if not single base pair
	}
}

////////////////////////////////////////////////////////////////////////////

void
PredictionTrackerBasePairProb::
computeMissingZ( const size_t i1, const size_t j1
							, const size_t i2, const size_t j2
							, const size_t l1, const size_t r1
							, const size_t l2, const size_t r2
							, PredictorMfeEns *predictor
							, SeedHandler* seedHandler )
{

	// ignore if invalid interation
	if ((i1 == j1 && i2 < j2) || (i2 == j2 && i1 < j1)) {
		LOG(DEBUG) << "ignoring invalid interaction";
		return;
	}

//	// check if unknown
//	if ( !Z_equal(getHybridZ(i1,j1,i2,j2,predictor),0) && !Z_equal(getHybridZ(i1,j1,i2,j2,predictor),0)) {
//		return;
//	}

	Z_type partZ = 0.0;
	const Z_type fullZ = getHybridZ(l1,r1,l2,r2,predictor);
	assert(!Z_equal(fullZ,0));

	const Z_type Zinit = energy.getBoltzmannWeight(energy.getE_init());

	// check if one side is existing -> easy computation
	if (j1==r1 && j2==r2) {
		partZ = getHybridZ(i1,r1,i2,r2,predictor);
		if (!Z_equal(partZ,0)) {
			updateHybridZ(l1, i1, l2, i2, fullZ/partZ*Zinit, predictor);
			return;
		}
		partZ = getHybridZ(l1, i1, l2, i2,predictor);
		if (!Z_equal(partZ,0)) {
			updateHybridZ(i1,r1,i2,r2, fullZ/partZ*Zinit, predictor);
			return;
		}
	}
	if (i1==l1 && i2==l2) {
		partZ = getHybridZ(l1,j1,l2,j2,predictor);
		if (!Z_equal(partZ,0)) {
			updateHybridZ(j1, r1, j2, r2, fullZ/partZ*Zinit, predictor);
			return;
		}
		partZ = getHybridZ(j1, r1, j2, r2,predictor);
		if (!Z_equal(partZ,0)) {
			updateHybridZ(l1,j1,l2,j2, fullZ/partZ*Zinit, predictor);
			return;
		}
	}

//	if (!isFullSeedinRegion(i1, j1, i2, j2, seedHandler)) {

//		// check if full seed outside of subregion (left or right)
//		if (isFullSeedinRegion(l1, i1, l2, i2, seedHandler)) {
//			// full seed left of subregion
//
//			// Case 2.1 (left)
//			// TODO: not working -> see test case 7
//			size_t seedCount = countNonOverlappingSeeds(l1, r1, l2, r2, seedHandler);
//
//			partZ = (
//				  getHybridZ(l1, r1, l2, r2, predictor)
//				/ getHybridZ(l1, i1, l2, i2, predictor)
//			) * energy.getBoltzmannWeight(energy.getE(i1,j1,i2,j2, E_type(0)))
//			* Zinit
//			/ seedCount;
//			updateHybridZ(i1, j1, i2, j2, partZ);
//
//		} else if (isFullSeedinRegion(j1, r1, j2, r2, seedHandler)) {
//			// full seed right of subregion
//
//			// Case 2.1 (right)
//			// TODO: not working -> see test case 7
//			size_t seedCount = countNonOverlappingSeeds(l1, r1, l2, r2, seedHandler);
//			if (i1 == 0 && i2 == 0 && j1 == 2 && j2 == 2) {
//				LOG(DEBUG) << "seeds: " << seedCount;
//				LOG(DEBUG) << getHybridZ(l1, r1, l2, r2, predictor);
//				LOG(DEBUG) << getHybridZ(j1, r1, j2, r2, predictor);
//			}
//
//
//			partZ = (
//				  getHybridZ(l1, r1, l2, r2, predictor)
//				/ getHybridZ(j1, r1, j2, r2, predictor)
//			) * energy.getBoltzmannWeight(energy.getE(i1,j1,i2,j2, E_type(0)))
//			* Zinit
//			/ seedCount;
//			updateHybridZ(i1, j1, i2, j2, partZ);
//
//		} else {

			// Case 2.2

	if( (l1==i1 && l2==i2) || (r1==j2 && r2==j2) ) {

			// check side of k
			size_t k1, k2;
//			bool kOnRight;
			if (i1 != l1) {
				// k at the left of region
//				kOnRight = false;
				k1 = i1;
				k2 = i2;
			} else {
//				kOnRight = true;
				k1 = j1;
				k2 = j2;
			}

			LOG(DEBUG) <<" split seed "<<l1<<":"<<l2<<" k "<<k1<<":"<<k2<<" "<<r1<<":"<<r2;
			// find leftmost overlapping seeds
			std::vector< std::pair <size_t, size_t> > overlappingSeeds = getLeftMostSeedsAtK(k1, k2, l1,l2, seedHandler);

			// loop over overlapping seeds
			for(std::vector< std::pair <size_t, size_t> >::iterator seedLeft = overlappingSeeds.begin(); seedLeft != overlappingSeeds.end(); ++seedLeft) {
				// calculate missing partition function
				const size_t sl1 = seedHandler->getSeedLength1(seedLeft->first, seedLeft->second);
				const size_t sl2 = seedHandler->getSeedLength2(seedLeft->first, seedLeft->second);

				LOG(DEBUG)<<"  >> seed "<<seedLeft->first<<":"<<seedLeft->second;

				// check if out of right boundary
				if ( seedLeft->first + sl1 -1 > r1 || seedLeft->first + sl2 -1 > r2 ) {
					continue;
				}


				// left part of seed until k
				partZ = energy.getBoltzmannWeight(getPartialSeedEnergy(seedLeft->first,seedLeft->second,seedLeft->first,k1,seedLeft->second,k2,seedHandler));
				// compute right part k-r
				assert( !Z_equal(getHybridZ(seedLeft->first,seedLeft->second,r1,r2,predictor),0));
				partZ = getHybridZ(seedLeft->first,seedLeft->second,r1,r2,predictor) / partZ;
				// store
				updateHybridZ(k1, r1, k2, r2, partZ, predictor);
				updateHybridZ(l1, k1, l2, k2, fullZ/partZ*Zinit, predictor);

//				if (i1==l1 && i2==l2) {
//					updateHybridZ(i1, k1, i2, k2, partZ, predictor);
//					updateHybridZ(l1, k1, l2, k2, partZ, predictor);
//					updateHybridZ(k1, j1, k2, j2, partZ, predictor);
//
//				} else
//				if (j1==r1 && j2==r2) {
//
//				}


//				if (kOnRight) {
//					if (seedLeft->first <= l1 && seedLeft->second <= l2) {
//						// region completely inside seed
//						partZ = energy.getBoltzmannWeight(getPartialSeedEnergy(seedLeft->first,seedLeft->second,l1,k1,l2,k2,seedHandler));
//					} else {
//						partZ += (
//							getHybridZ(l1, seedLeft->first+sl1-1, l2, seedLeft->second+sl2-1, predictor)
//							* energy.getBoltzmannWeight(energy.getED1(i1, k1) + energy.getED2(i2, k2))
//						) / energy.getBoltzmannWeight(getPartialSeedEnergy(seedLeft->first,seedLeft->second,k1,seedLeft->first+sl1-1,k2,seedLeft->second+sl2-1,seedHandler));
//					}
//				} else {
//					if (j1 <= seedLeft->first+sl1-1 && j2 <= seedLeft->second+sl2-1) {
//						// region completely inside seed
//						partZ = energy.getBoltzmannWeight(getPartialSeedEnergy(seedLeft->first,seedLeft->second,k1,j1,k2,j2,seedHandler));
//					} else {
//						partZ += (
//							getHybridZ(seedLeft->first, j1, seedLeft->second, j2, predictor)
//							* energy.getBoltzmannWeight(energy.getED1(k1, j1) + energy.getED2(k2, j2))
//						) / energy.getBoltzmannWeight(getPartialSeedEnergy(seedLeft->first,seedLeft->second,seedLeft->first,k1,seedLeft->second,k2,seedHandler));
//					}
//				}
			}

//			if (kOnRight) {
//				updateHybridZ(i1, k1, i2, k2, partZ);
//			} else {
//				updateHybridZ(k1, j1, k2, j2, partZ);
//			}

//		} else {
//			assert( false /* i-j not left/right bound */ );
		}

//	}

}

////////////////////////////////////////////////////////////////////////////

size_t
PredictionTrackerBasePairProb::
countNonOverlappingSeeds( const size_t i1, const size_t j1
				                , const size_t i2, const size_t j2
				                , SeedHandler *seedHandler )
{
	size_t count = 0;
	std::vector<Interaction::BasePair> seeds;
	size_t si1 = RnaSequence::lastPos, si2 = RnaSequence::lastPos;
	while( seedHandler->updateToNextSeed(si1,si2
			, i1, j1+1-seedHandler->getConstraint().getBasePairs()
			, i2, j2+1-seedHandler->getConstraint().getBasePairs()) )
	{
		const size_t sl1 = seedHandler->getSeedLength1(si1, si2);
		const size_t sl2 = seedHandler->getSeedLength2(si1, si2);

		// ignore seeds that are conflicting with region borders
		if ((si1 == i1 && si2 != i2) ||
		     (si1 != i1 && si2 == i2) ||
				 (si1+sl1-1 == j1 && si2+sl2-1 != j2) ||
				 (si1+sl1-1 != j1 && si2+sl2-1 == j2)) continue;

		Interaction::BasePair seed(si1, si2);
		bool found = false;
		for (auto it = seeds.begin(); it != seeds.end(); ++it) {
			if (seedHandler->areLoopOverlapping(si1, si2, it->first, it->second)) {
				found = true;
				break;
			}
		}

		if (!found) {
			seeds.push_back(seed);
			count++;
		}
	}
	return count;
}

////////////////////////////////////////////////////////////////////////////

bool
PredictionTrackerBasePairProb::
isFullSeedinRegion( const size_t i1, const size_t j1
				        	, const size_t i2, const size_t j2
				        	, SeedHandler* seedHandler )
{
	size_t si1 = RnaSequence::lastPos, si2 = RnaSequence::lastPos;
	while( seedHandler->updateToNextSeed(si1,si2
			, i1, j1+1-seedHandler->getConstraint().getBasePairs()
			, i2, j2+1-seedHandler->getConstraint().getBasePairs()) )
	{
		const size_t sl1 = seedHandler->getSeedLength1(si1, si2);
		const size_t sl2 = seedHandler->getSeedLength2(si1, si2);
		const size_t sj1 = si1+sl1-1;
		const size_t sj2 = si2+sl2-1;
		if (sj1 <= j1 && sj2 <= j2) {
			return true;
		}
	}
	return false;
}

////////////////////////////////////////////////////////////////////////////

std::vector< std::pair <size_t, size_t> >
PredictionTrackerBasePairProb::
getLeftMostSeedsAtK( const size_t k1, const size_t k2
				, const size_t i1min, const size_t i2min
									 , SeedHandler* seedHandler )
{
  std::vector< std::pair <size_t, size_t> > seeds;
	size_t maxSeedLength = seedHandler->getConstraint().getBasePairs();

	size_t i1 = std::max( i1min, k1 - std::min(k1, maxSeedLength-1));
	size_t i2 = std::max( i2min, k2 - std::min(k2, maxSeedLength-1));

	Interaction tmpInter(energy.getAccessibility1().getSequence(),energy.getAccessibility2().getAccessibilityOrigin().getSequence());
	// TODO works only for seedNoBulge : add "getMaxSeedLength1|2() to SeedHandler with implementations in subclasses"

	size_t si1 = RnaSequence::lastPos, si2 = RnaSequence::lastPos;
	while( seedHandler->updateToNextSeed(si1, si2, i1, k1, i2, k2) ) {
		// check if k is part of seed (inner base pair)
		if (   si1>=k1
			|| si2>=k2
			|| (k1 >= si1+seedHandler->getSeedLength1(si1,si2)-1)
			|| (k2 >= si2+seedHandler->getSeedLength2(si1,si2)-1) )
		{
			continue;
		}

		// check if k is truely bp of seed
		tmpInter.basePairs.clear();
		seedHandler->traceBackSeed(tmpInter,si1,si2);
		bool seedBpNotFound = true;
		for( auto bp = tmpInter.basePairs.begin(); seedBpNotFound && bp!=tmpInter.basePairs.end(); bp++) {
			seedBpNotFound = (bp->first != k1) || (bp->second != k2);
		}
		if (seedBpNotFound) {
			continue;
		}

		// check if a loop-overlapping seed already exists
		bool found = false;
		for (auto it = seeds.begin(); it != seeds.end(); ++it) {
			if (seedHandler->areLoopOverlapping(si1, si2, it->first, it->second)) {
				// replace if new seed is more left
				if (si1 < it->first && si2 < it->second) {
					*it = std::make_pair(si1, si2);
				}
				found = true;
				break;
			}
		}
		if (!found) {
			seeds.push_back(std::make_pair(si1, si2));
		}
	}
	return seeds;
}

////////////////////////////////////////////////////////////////////////////

E_type
PredictionTrackerBasePairProb::
getPartialSeedEnergy( const size_t si1, const size_t si2
										, const size_t i1, const size_t j1
										, const size_t i2, const size_t j2
										, SeedHandler* seedHandler )
{
	// ignore if invalid interation
	if ((i1 == j1 && i2 < j2) || (i2 == j2 && i1 < j1)) {
		return energy.getE_init();
	}

	// trace S
	Interaction interaction = Interaction(energy.getAccessibility1().getSequence(), energy.getAccessibility2().getAccessibilityOrigin().getSequence());
	const size_t sl1 = seedHandler->getSeedLength1(si1, si2);
	const size_t sl2 = seedHandler->getSeedLength2(si1, si2);
	interaction.basePairs.push_back( energy.getBasePair(si1, si2) );
	seedHandler->traceBackSeed( interaction, si1, si2 );
	interaction.basePairs.push_back( energy.getBasePair(si1+sl1-1, si2+sl2-1) );

	E_type partE = 0;
	size_t i1old = i1;
	size_t i2old = i2;

	for (size_t i = 0; i < interaction.basePairs.size(); i++) {
		// get index of current base pair
		size_t s1 = energy.getIndex1(interaction.basePairs[i]);
		size_t s2 = energy.getIndex2(interaction.basePairs[i]);

		if (s1 < i1 || s2 < i2) continue;
		if (s1 > j1 && s2 > j2) break;

		if (!E_isINF(energy.getE_interLeft(i1old,s1,i2old,s2))) {
			// add hybridization energy
			partE += energy.getE_interLeft(i1old,s1,i2old,s2);
			// store
			i1old = s1;
			i2old = s2;
		}
	}

	return partE + energy.getE_init();
}

////////////////////////////////////////////////////////////////////////////

Z_type
PredictionTrackerBasePairProb::
getHybridZ( const size_t i1, const size_t j1
					, const size_t i2, const size_t j2
					, PredictorMfeEns *predictor)
{
	Z_type partZ = predictor->getHybridZ(i1, j1, i2, j2);
	if (Z_equal(partZ, 0)) {
		Interaction::Boundary key(i1, j1, i2, j2);
		// check in original data
		if ( predictor->getZPartition().find(key) != predictor->getZPartition().end() ) {
			partZ = predictor->getZPartition().find(key)->second;
		} else
		// check in additional data
		if ( Z_partitionMissing.find(key) != Z_partitionMissing.end() ) {
			partZ = Z_partitionMissing[key];
		} else {
			// fall back
			partZ = Z_type(0);
		}
	}
	return partZ;
}

////////////////////////////////////////////////////////////////////////////

Z_type
PredictionTrackerBasePairProb::
getBasePairProb( const size_t i1, const size_t i2
					     , PredictorMfeEns *predictor)
{
	Interaction::BasePair key(i1, i2);
	if ( structureProbs.find(key) == structureProbs.end() ) {
		return Z_type(0);
	} else {
		return structureProbs[key];
	}
}

////////////////////////////////////////////////////////////////////////////

void
PredictionTrackerBasePairProb::
updateHybridZ( const size_t i1, const size_t j1
						 , const size_t i2, const size_t j2
						 , const Z_type partZ
						 , PredictorMfeEns *predictor )
{
	// store if unknown
	if (Z_equal(getHybridZ(i1,j1,i2,j2,predictor),0)) {
		Interaction::Boundary key(i1,j1,i2,j2);
		Z_partitionMissing[key] = partZ;
	}
}

////////////////////////////////////////////////////////////////////////////

bool
PredictionTrackerBasePairProb::
generateDotPlot( char *seq1, char *seq2, char *fileName
							 , plist *pl, const char *comment
						   , Interaction::Boundary interactionBoundary )
{
	FILE *file;
	file = fopen(fileName,"w");
	if (file == NULL) return false; /* failure */

	int bbox[4];
	bbox[0] = 0;
	bbox[1] = 0;
	bbox[2] = 72 * (strlen(seq1) + 3);
	bbox[3] = 72 * (strlen(seq2) + 3);

	fprintf(file,
          "%%!PS-Adobe-3.0 EPSF-3.0\n"
          "%%%%Creator: IntaRNA\n"
          "%%%%Title: RNA Dot Plot\n"
          "%%%%BoundingBox: %d %d %d %d\n"
          "%%%%DocumentFonts: Helvetica\n"
          "%%%%Pages: 1\n"
          "%%%%EndComments\n\n",
          bbox[0], bbox[1], bbox[2], bbox[3]);

  // comment

	if (comment) {
    fprintf(file, "%%%% %s\n", comment);
  }

	fprintf(file, "/DPdict 100 dict def\n");
	fprintf(file, "DPdict begin\n");

	// ps template

	fprintf(file, dotplotTemplate);
  fprintf(file, "end\n");
	fprintf(file, "DPdict begin\n");

	// sequences

	unsigned int i, length;
  length = strlen(seq1);
  fprintf(file, "/sequence1 { (\\\n");
  i = 0;
  while (i < length) {
    fprintf(file, "%.255s\\\n", seq1 + i);  /* no lines longer than 255 */
    i += 255;
  }
  fprintf(file, ") } def\n");
	fprintf(file, "/len { sequence1 length } bind def\n\n");
	length = strlen(seq2);
  fprintf(file, "/sequence2 { (\\\n");
  i = 0;
  while (i < length) {
    fprintf(file, "%.255s\\\n", seq2 + i);  /* no lines longer than 255 */
    i += 255;
  }
  fprintf(file, ") } def\n");
	fprintf(file, "/len2 { sequence2 length } bind def\n\n");

	fprintf(file, "72 72 translate\n"
	               "72 72 scale\n");

  fprintf(file, "/Helvetica findfont 0.95 scalefont setfont\n\n");

	fprintf(file,"drawseq1\n");
	fprintf(file,"drawseq2\n");

	// basepair data

	fprintf(file,"%%data starts here\n");

	fprintf(file, "\n%%draw the grid\ndrawgrid\n\n");
	fprintf(file,"%%start of base pair probability data\n");

  fprintf(file, "/coor [\n");

	for (plist *pl1 = pl; pl1->i > 0; pl1++) {
    if (pl1->type == 0) {
      fprintf(file, "%1.9f %d %d box\n", sqrt(pl1->p), pl1->i, pl1->j);
    }
  }

  fprintf(file, "] def\n");

	// print frame
	fprintf(file,
		      "0.03 setlinewidth\n\
           %1.1f %1.1f %zu %zu rectangle\n\
					 0 0 0 setrgbcolor\n\
           stroke\n", 0.5, 0.5, strlen(seq1), strlen(seq2));

	// print best interaction outline
	fprintf(file,
		      "0.03 setlinewidth\n\
           %1.1f %1.1f %zu %zu rectangle\n\
					 1 0 0 setrgbcolor\n\
           stroke\n", (float)interactionBoundary.i1 + 0.5, (float)interactionBoundary.i2 + 0.5, interactionBoundary.j1 - interactionBoundary.i1 + 1, interactionBoundary.j2 - interactionBoundary.i2 + 1);

	fprintf(file, "showpage\n"
	            "end\n"
              "%%%%EOF\n");

	fclose(file);
	return true; /* success */
}

////////////////////////////////////////////////////////////////////////////

} // namespace
