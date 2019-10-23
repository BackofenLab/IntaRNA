
#include "IntaRNA/PredictionTrackerBasePairProb.h"

extern "C" {
	#include <ViennaRNA/vrna_config.h>
	#include <ViennaRNA/plotting/probabilities.h>
}

#include <boost/foreach.hpp>
#include <boost/regex.hpp>

namespace IntaRNA {

const Interaction::BasePair bpToTrack(1,1);

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
	hasSeedhandler = (seedHandler != NULL);

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
	for (auto z = Z_partition.begin(); z != Z_partition.end(); ++z) {
		// identify best interaction boundary
		Z_type Zstruct = z->second * energy.getBoltzmannWeight(energy.getE(z->first.i1, z->first.j1, z->first.i2, z->first.j2, E_type(0)));
		if (Zstruct > maxZ) {
			maxZ = Zstruct;
			interactionBoundary = z->first;
		}

		Interaction::BasePair iBP(z->first.i1, z->first.i2);
		Interaction::BasePair jBP(z->first.j1, z->first.j2);

		// create left and jBP index
		if (z->first.i1 != z->first.j1 && z->first.i2 != z->first.j2) {
			// encode iBP/jBP boundary
			// create left index
			rightExt[iBP].insert(jBP);
			// create right index
			if (seedHandler != NULL) {
				leftExt[jBP].insert(iBP);
			}
		}

		// bp-prob without extension
		/*if (predictor->getZPartition().find( z->first ) != predictor->getZPartition().end()) {
			updateProb(jBP, z->second);
			LOG_IF(Interaction::BasePair(z->first.j1,z->first.j2)==bpToTrack, DEBUG) <<(z->first) << " - no ext at j: " <<  z->second;
			if (iBP != jBP) {
				updateProb(iBP, z->second);
				LOG_IF(Interaction::BasePair(z->first.i1,z->first.i2)==bpToTrack, DEBUG) <<(z->first) << " - no ext at i: " <<  z->second;
			}
		}*/


	} // it (Z_partition)

/*
	// compute all decompositions of Z not already known
	for (auto z1 = Z_partition.begin(); z1 != Z_partition.end(); ++z1) {
		Interaction::BasePair lBound(z1->first.i1, z1->first.i2);
		Interaction::BasePair rBound(z1->first.j1, z1->first.j2);
		// skip single base pairs
		if (lBound == rBound) {
			continue;
		}
		// decompose inside with left bound aligned
		auto rightExtofL = rightExt.find(lBound);
		if (rightExtofL != rightExt.end()) {
			for (auto k = rightExtofL->second.begin(); k != rightExtofL->second.end(); ++k) {
				// within boundary
				if (k->first < rBound.first && k->second < rBound.second) {
//					computeMissingZ( z1->first, *k, predictor, seedHandler);
					updateProb(*k, z1->second);
					// TODO : only needed if right bound of seed?!
					rightExt[*k].insert(rBound);
				}
			}
		}
		// decompose inside with right bound aligned
		auto leftExtofR = leftExt.find(rBound);
		if (leftExtofR != leftExt.end()) {
			for (auto l = leftExtofR->second.begin(); l != leftExtofR->second.end(); ++l) {
				// left out of boundary
				if (l->first < lBound.first && l->second < lBound.second) {
//					computeMissingZ( z1->first, *k, predictor, seedHandler);
					updateProb(lBound, getHybridZ( l->first, rBound.first, l->second, rBound.second, predictor));
//					rightExt[lBound].insert(rBound);
				}
			}
		}
	} // it (Z_partition)

	for (auto z = Z_partition.begin(); z != Z_partitionMissing.end(); z++) {
		LOG(DEBUG) <<" initZ: "<<z->first<<" = "<<z->second;
	}

	for (auto z = Z_partitionMissing.begin(); z != Z_partitionMissing.end(); z++) {
		LOG(DEBUG) <<" newZ: "<<z->first <<" = "<<z->second;
	}

	for (auto right = rightExt.begin(); right != rightExt.end(); right++) {
		for (auto bp = right->second.begin(); bp != right->second.end(); bp++) {
	  	LOG(DEBUG) <<" RightExt at: "<<right->first <<": "<<*bp;
		}
	}

	// compute missing Z values for seed-based predictions
	if (hasSeedhandler
			&& ( ! seedHandler->getConstraint().getExplicitSeeds().empty()
					|| seedHandler->getConstraint().getBasePairs() >= 2) )
	{

		std::set< Interaction::BasePair > unknownSeedBP;
		// TODO identify all inner-seed-basepairs that are not known in rightExt or leftExt

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

			//LOG(DEBUG)<<" internal of seed "<<si1<<":"<<si2;

			// trace all inner seed base pairs
			Interaction interaction = Interaction(energy.getAccessibility1().getSequence(), energy.getAccessibility2().getAccessibilityOrigin().getSequence());
			seedHandler->traceBackSeed( interaction, si1, si2 );

			// iterate all internal seed base pairs
			for (size_t i = 0; i < interaction.basePairs.size(); i++) {
				// get index of current base pair
				const size_t sk1 = energy.getIndex1(interaction.basePairs[i]);
				const size_t sk2 = energy.getIndex2(interaction.basePairs[i]);
				const Interaction::BasePair skBP(sk1, sk2);
			//LOG(DEBUG)<<"   >> "<<sk1<<":"<<sk2<<"  Z "<<getHybridZ(sk1, sj1, sk2, sj2, predictor);

				// no extension
//				if (Z_equal(getHybridZ(sk1, sj1, sk2, sj2, predictor), 0)
//					|| Z_equal(getHybridZ(si1,sk1, si2, sk2, predictor), 0))
//				{
//			LOG(DEBUG)<<"  missing "<<sk1<<":"<<sk2<<" j "<<sj1<<":"<<sj2;
//				}
				computeMissingZ(sk1, sj1, sk2, sj2, si1, sj1, si2, sj2, predictor, seedHandler);
				// update right extension index
				rightExt[skBP].insert(sjBP);

				// right extensions
				// iterate all r=rightExt[si] to computeMissingZ(si,sk,r) and update rightExt[sk].push(r)
				for (auto seedRightExt = rightExt[siBP].begin(); seedRightExt != rightExt[siBP].end(); ++seedRightExt) {
					if (sk1 < seedRightExt->first && sk2 < seedRightExt->second
//						&& Z_equal(getHybridZ(sk1, seedRightExt->first, sk2, seedRightExt->second, predictor), 0)
					)
					{
						computeMissingZ(sk1, seedRightExt->first, sk2, seedRightExt->second, si1, seedRightExt->first, si2, seedRightExt->second, predictor, seedHandler);
						// update right extension index
						rightExt[skBP].insert(*seedRightExt);
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


			// TODO compute their probability
			for (auto bp = unknownSeedBP.begin(); bp != unknownSeedBP.end(); bp++) {
				// TODO replace with sane left bound
				const size_t lMin1 = bp->first;
				const size_t lMin2 = bp->second;
				auto seedsWithBP = getLeftMostSeedsAtK(bp->first, bp->second, lMin1, lMin2, seedHandler);

				// compute respective probabilities
				for (auto seed = seedsWithBP.begin(); seed!=seedsWithBP.end(); seed++) {
					// TODO compute respective probability for this seed's decomposition and all its left/right extensions
				}
			}

		} // for all seeds

	} // if seedhandler
*/

	// Compute base-pair probabilities via combinations
	if (!hasSeedhandler) {
		computeBasePairProbsNoSeed(predictor);
	} else {
		// TODO: computeBasePairProbs
		throw std::runtime_error("PredictionTrackerBasePairProb not yet implemented for seeds");
	}

	// build plist
	struct vrna_elem_prob_s plist[structureProbs.size()+1];
	size_t i = 0;
	const Z_type Zall = predictor->getZall();
	for (auto sp = structureProbs.begin(); sp != structureProbs.end(); ++sp) {
		LOG_IF(Interaction::BasePair(sp->first.first,sp->first.second)==bpToTrack,DEBUG) << "prob: " << toString(bpToTrack) << " = " << sp->second <<"   Z = "<<sp->second*predictor->getZall();
		if ( (sp->second /Zall)  > probabilityThreshold) {
			plist[i].i = sp->first.first + 1;
			plist[i].j = sp->first.second + 1;
			plist[i].p = (sp->second /Zall);
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
computeBasePairProbsNoSeed( const PredictorMfeEns *predictor )
{
	for (auto z = predictor->getZPartition().begin(); z != predictor->getZPartition().end(); ++z) {
		assert( !Z_equal(z->second, 0) );

		Z_type bpProb = z->second * energy.getBoltzmannWeight(energy.getE(z->first.i1, z->first.j1, z->first.i2, z->first.j2, E_type(0)));

		Interaction::BasePair iBP(z->first.i1, z->first.i2);
		Interaction::BasePair jBP(z->first.j1, z->first.j2);

		// left end and single bp
		updateProb(iBP, bpProb);

		if (iBP < jBP) {
			// right end
			updateProb(jBP, bpProb);

			// inner (rightExt > jBP)
			for (auto right = rightExt[jBP].begin(); right != rightExt[jBP].end(); ++right) {
				Z_type innerProb = z->second * getHybridZ(z->first.j1, right->first, z->first.j2, right->second, predictor)
									 * energy.getBoltzmannWeight(energy.getE(z->first.i1, right->first, z->first.i2, right->second, -energy.getE_init()));
				updateProb(jBP, innerProb);
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////

void
PredictionTrackerBasePairProb::
computeBasePairProbs( const PredictorMfeEns *predictor
										, const PredictorMfeEns::Site2Z_hash::const_iterator first
										, const PredictorMfeEns::Site2Z_hash::const_iterator last )
{
	const Z_type Zinit = energy.getBoltzmannWeight(energy.getE_init());

	for (auto z = first; z != last; ++z) {
		Z_type bpProb = z->second * energy.getBoltzmannWeight(energy.getE(z->first.i1, z->first.j1, z->first.i2, z->first.j2, E_type(0)));

		Interaction::BasePair iBP(z->first.i1, z->first.i2);
		Interaction::BasePair jBP(z->first.j1, z->first.j2);

		LOG_IF(jBP==bpToTrack,DEBUG)<<"computeBpProb : rightExt of"
				<< " i: " << iBP
				<< " j: " << jBP
				;


		assert( !Z_equal(z->second,0) );

		Z_type probOfJ = 0;

		// extensions
		for (auto right = rightExt[jBP].begin(); right != rightExt[jBP].end(); ++right) {

			//LOG_IF((z->first.j1 == 2 && z->first.j2 == 2), DEBUG) <<" rightExt("<<jBP.first<<":"<<jBP.second<<") = "<<right->first<<":"<<right->second;

			assert( !Z_equal(getHybridZ(z->first.j1, right->first, z->first.j2, right->second, predictor),0) );

			// ensure extension is valid (present in original Z data)
			if (predictor->getZPartition().find( Interaction::Boundary(z->first.i1, right->first, z->first.i2, right->second)) != predictor->getZPartition().end()) {

				assert( !Z_equal(z->second, 0) );
				bpProb = z->second * getHybridZ(z->first.j1, right->first, z->first.j2, right->second, predictor)
							// ED penalty
									 * energy.getBoltzmannWeight(energy.getE(z->first.i1, right->first, z->first.i2, right->second, E_type(0)))
									 / Zinit;

				probOfJ += bpProb;
				LOG_IF(jBP==bpToTrack,DEBUG)
				<<"("<<jBP.first<<":"<<jBP.second<<")" << " - ext at: " << (*right) << "=" << bpProb << " | left: " << z->second << " | right = " << getHybridZ(z->first.j1, right->first, z->first.j2, right->second, predictor);
			}
		}
		// init if needed
		if (!Z_equal(probOfJ,0)) {
			updateProb(jBP, probOfJ);
		}
	}
}

////////////////////////////////////////////////////////////////////////////

std::pair< Z_type, Z_type >
PredictionTrackerBasePairProb::
computeMissingZseed( const size_t l1, const size_t k1, const size_t r1
				      	   , const size_t l2 , const size_t k2, const size_t r2
					         , const PredictorMfeEns *predictor
					         , const SeedHandler* seedHandler )
{
	Z_type leftZ = 0, rightZ = 0;

	const Z_type fullZ = getHybridZ(l1,r1,l2,r2,predictor);
	const Z_type Zinit = energy.getBoltzmannWeight(energy.getE_init());

	LOG_IF(Interaction::BasePair(k1,k2)==bpToTrack,DEBUG) <<" split seed "<<l1<<":"<<l2<<" k "<<k1<<":"<<k2<<" "<<r1<<":"<<r2;
	// find leftmost overlapping seeds
	std::vector< Interaction::BasePair > overlappingSeeds = getLeftMostSeedsAtK(k1, k2, l1,l2, seedHandler);
	LOG_IF(Interaction::BasePair(k1,k2)==bpToTrack,DEBUG) << " seedCount: " << overlappingSeeds.size();

	// loop over overlapping seeds
	for (auto seedLeft = overlappingSeeds.begin(); seedLeft != overlappingSeeds.end(); ++seedLeft) {
		// calculate missing partition function
		const size_t sl1 = seedHandler->getSeedLength1(seedLeft->first, seedLeft->second);
		const size_t sl2 = seedHandler->getSeedLength2(seedLeft->first, seedLeft->second);

		// check if out of right boundary
		if ( seedLeft->first + sl1 -1 > r1 || seedLeft->second + sl2 -1 > r2 ) {
			continue;
		}

		// find leftmost overlapping seeds
		std::vector< Interaction::BasePair > overlappingSeeds = getLeftMostSeedsAtK(k1, k2, l1,l2, seedHandler);

		// left part of seed until k
		Z_type partZ = energy.getBoltzmannWeight(getPartialSeedEnergy(seedLeft->first,seedLeft->second,seedLeft->first,k1,seedLeft->second,k2,seedHandler));
		// compute right part k-r
		assert( !Z_equal(getHybridZ(seedLeft->first,r1,seedLeft->second,r2,predictor),0));
		partZ = getHybridZ(seedLeft->first,r1,seedLeft->second,r2,predictor) / partZ;
		// store
		rightZ += partZ;
		leftZ += fullZ/partZ*Zinit;
	}
	return std::make_pair( leftZ, rightZ );

}

////////////////////////////////////////////////////////////////////////////

void
PredictionTrackerBasePairProb::
computeMissingZ( const Interaction::Boundary & boundFull
				       , const Interaction::BasePair & k
				       , const PredictorMfeEns *predictor
				       , const SeedHandler* seedHandler )
{

	// ignore if invalid interation
	if (!(boundFull.i1 < k.first && k.first < boundFull.j1)
		|| !(boundFull.i2 < k.second && k.second < boundFull.j2))
	{
		LOG(DEBUG) << "ignoring invalid interaction : "<<boundFull<<" k "<<k;
		return;
	}

	typedef Interaction::Boundary IB;

	const IB boundL = IB(boundFull.i1,k.first,boundFull.i2,k.second);
	const IB boundR = IB(k.first,boundFull.j1,k.second,boundFull.j2);

	Z_type ZL = getHybridZ(boundL,predictor);
	Z_type ZR = getHybridZ(boundR,predictor);

	if (Z_equal( ZL, 0 ) && Z_equal( ZR, 0 )) {
		LOG(DEBUG)<<"computeMissingZ("<<boundL<<","<<boundR<<"): both sides not known ";
	}

	auto fullZ = predictor->getZPartition().find(boundFull);
	assert(fullZ != predictor->getZPartition().end());
	const Z_type Zinit = energy.getBoltzmannWeight(energy.getE(k.first,k.first,k.second,k.second,energy.getE_init()));

	// store unknown part
	if ( Z_equal( ZL, 0 ) ) {
		updateHybridZ( boundL, fullZ->second / ZR * Zinit, *predictor );
	} else
	if ( Z_equal( ZR, 0 ) ) {
		updateHybridZ( boundR, fullZ->second / ZL * Zinit, *predictor );
	} else {
		// sanity check if both known
		LOG_IF( !Z_equal( (ZL * ZR / Zinit), fullZ->second ), DEBUG) <<" decomposition not correct: ZL("<<ZL<<" of "<<boundL<<") * ZR("<<ZR<<" of "<<boundR<<")/Zinit("<<Zinit<<") != Zfull("<<fullZ->second<<") for "<<boundFull;
	}

}

////////////////////////////////////////////////////////////////////////////

bool
PredictionTrackerBasePairProb::
isFullSeedinRegion( const size_t i1, const size_t j1
				        	, const size_t i2, const size_t j2
				        	, const SeedHandler* seedHandler )
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

std::vector< Interaction::BasePair >
PredictionTrackerBasePairProb::
getLeftMostSeedsAtK( const size_t k1, const size_t k2
				           , const size_t i1min, const size_t i2min
									 , const SeedHandler* seedHandler )
{
  std::vector<Interaction::BasePair> seeds;
	size_t maxSeedLength = seedHandler->getConstraint().getBasePairs();

	size_t i1 = std::max( i1min, k1 - std::min(k1, maxSeedLength-1));
	size_t i2 = std::max( i2min, k2 - std::min(k2, maxSeedLength-1));

	Interaction tmpInter(energy.getAccessibility1().getSequence(),energy.getAccessibility2().getAccessibilityOrigin().getSequence());
	// TODO works only for seedNoBulge : add "getMaxSeedLength1|2() to SeedHandler with implementations in subclasses"

	size_t si1 = RnaSequence::lastPos, si2 = RnaSequence::lastPos;
	while( seedHandler->updateToNextSeed(si1, si2, i1, k1, i2, k2) ) {
		// check if k is part of seed (inner base pair)
		if ( si1>=k1
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
			seedBpNotFound = (energy.getIndex1(*bp) != k1) || (energy.getIndex2(*bp) != k2);
		}
		if (seedBpNotFound) {
			continue;
		}

		// check if a loop-overlapping seed already exists
		bool found = false;
		for (auto seed = seeds.begin(); seed != seeds.end(); ++seed) {
			if (seedHandler->areLoopOverlapping(si1, si2, seed->first, seed->second)) {
				// replace if new seed is more left
				if (si1 < seed->first && si2 < seed->second) {
					*seed = Interaction::BasePair(si1, si2);
				}
				found = true;
				break;
			}
		}
		if (!found) {
			Interaction::BasePair seed(si1, si2);
			seeds.push_back(seed);
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
										, const SeedHandler* seedHandler )
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
					, const PredictorMfeEns *predictor)
{
	Interaction::Boundary boundary(i1, j1, i2, j2);
	return getHybridZ(boundary,predictor);
}

////////////////////////////////////////////////////////////////////////////

Z_type
PredictionTrackerBasePairProb::
getHybridZ( const Interaction::Boundary & boundary
					, const PredictorMfeEns *predictor)
{
	// check in original data
	if ( predictor->getZPartition().find(boundary) != predictor->getZPartition().end() ) {
		return predictor->getZPartition().find(boundary)->second;
	} else
	// check in additional data
	if ( Z_partitionMissing.find(boundary) != Z_partitionMissing.end() ) {
		return Z_partitionMissing.find(boundary)->second;
	} else {
		// fall back
		return Z_type(0);
	}
}

////////////////////////////////////////////////////////////////////////////

Z_type
PredictionTrackerBasePairProb::
getBasePairProb( const size_t i1, const size_t i2
					     , const PredictorMfeEns *predictor)
{
	Interaction::BasePair bp(i1, i2);
	if ( structureProbs.find(bp) == structureProbs.end() ) {
		return Z_type(0);
	} else {
		return structureProbs[bp] / predictor->getZall();
	}
}

////////////////////////////////////////////////////////////////////////////

void
PredictionTrackerBasePairProb::
updateHybridZ( const size_t i1, const size_t j1
						 , const size_t i2, const size_t j2
						 , const Z_type partZ
						 , const PredictorMfeEns & predictor )
{
	updateHybridZ(Interaction::Boundary(i1,j1,i2,j2),partZ,predictor);
}

////////////////////////////////////////////////////////////////////////////

void
PredictionTrackerBasePairProb::
updateHybridZ( const Interaction::Boundary & boundary
						 , const Z_type partZ
						 , const PredictorMfeEns & predictor )
{
	// store if unknown
	if (Z_partitionMissing.find( boundary) == Z_partitionMissing.end()) {
		Z_partitionMissing[boundary] = partZ;
	} else {
		LOG_IF( !Z_equal(Z_partitionMissing.find(boundary)->second,partZ), DEBUG )
				<<"updateHybridZ( "<<boundary<<" ) new val "<<partZ<<" != "<<Z_partitionMissing.find(boundary)->second<<" old val";
	}
}

////////////////////////////////////////////////////////////////////////////

bool
PredictionTrackerBasePairProb::
generateDotPlot( const char *seq1, const char *seq2, const char *fileName
							 , const plist *pl, const char *comment
						   , const Interaction::Boundary interactionBoundary )
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

	for (const plist *pl1 = pl; pl1->i > 0; pl1++) {
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
