
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
	const std::string & rna2 = energy.getAccessibility2().getSequence().asString();

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

	for (auto it = Z_partition.begin(); it != Z_partition.end(); ++it) {
		// identify best interaction boundary
		Z_type Zstruct = it->second * energy.getBoltzmannWeight(energy.getE(it->first.i1, it->first.j1, it->first.i2, it->first.j2, E_type(0)));
		if (Zstruct > maxZ) {
			maxZ = Zstruct;
			interactionBoundary = it->first;
		}

		// create left and right index
		if (it->first.i1 != it->first.j1 && it->first.i2 != it->first.j2) {
			Interaction::BasePair key(it->first.i1, it->first.i2);
			Interaction::BasePair entry(it->first.j1, it->first.j2);
			leftIndex[key].push_back(entry);

			if (seedHandler != NULL) {
	      // create right index
				Interaction::BasePair key(it->first.j1, it->first.j2);
				Interaction::BasePair entry(it->first.i1, it->first.i2);
				rightIndex[key].push_back(entry);
			}
		}

	} // it (Z_partition)

	// print indexes

	for (auto it = rightIndex.begin(); it != rightIndex.end(); ++it) {
		LOG(DEBUG) << "right: (" << it->first.first << ":" << it->first.second << ") " << it->second.size();
	}

	// ============

	if (seedHandler != NULL) {
		// iterate seeds
		size_t si1 = RnaSequence::lastPos, si2 = RnaSequence::lastPos;
		while( seedHandler->updateToNextSeed(si1, si2
				  , 0, s1-seedHandler->getConstraint().getBasePairs()
			  	, 0, s2-seedHandler->getConstraint().getBasePairs() ) )
		{
			// compute missing Z values for seed basepairs
			// and update left index
			size_t sl = seedHandler->getConstraint().getBasePairs();
			Interaction::BasePair seedKey(si1+sl-1, si2+sl-1);
			for (size_t i = 0; i < sl; i++ ) {

				// right extensions
				for (auto it = leftIndex[seedKey].begin(); it != leftIndex[seedKey].end(); ++it) {
					if (si1+i < it->first && si2+i < it->second && Z_equal(getHybridZ(si1+i, it->first, si2+i, it->second, predictor), 0)) {
						computeMissingZ(si1+i, it->first, si2+i, it->second, si1, it->first, si2, it->second, predictor, seedHandler);
						Interaction::BasePair key(si1+i, si2+i);
						Interaction::BasePair entry(it->first, it->second);
						leftIndex[key].push_back(entry);
					}
				}

				// left extensions
				for (auto it = rightIndex[seedKey].begin(); it != rightIndex[seedKey].end(); ++it) {
					if (si1+i > it->first && si2+i > it->second && Z_equal(getHybridZ(it->first, si1+i, it->second, si2+i, predictor), 0)) {
						computeMissingZ(it->first, si1+i, it->second, si2+i, it->first, si1+sl-1, it->second, si2+sl-1, predictor, seedHandler);
						Interaction::BasePair key(it->first, it->second);
						Interaction::BasePair entry(si1+i, si2+i);
						leftIndex[key].push_back(entry);
					}
				}

				// inside seed
				for (size_t j = i; j < seedHandler->getConstraint().getBasePairs(); j++ ) {
					if (Z_equal(getHybridZ(si1+i, si1+j, si2+i, si2+j, predictor), 0)) {
						computeMissingZ(si1+i, si1+j, si2+i, si2+j, si1, si1+sl-1, si2, si2+sl-1, predictor, seedHandler);
						if (j != i) {
							Interaction::BasePair key(si1+i, si2+i);
							Interaction::BasePair entry(si1+j, si2+j);
							leftIndex[key].push_back(entry);
						}
					}
				}

			}
		}

	}

	for (auto it = leftIndex.begin(); it != leftIndex.end(); ++it) {
		for (auto it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
		  LOG(DEBUG) << "left: (" << it->first.first << ":" << it->first.second << ") " << it2->first << ":" << it2->second;
		}
	}

	// print partitions

	for (auto it = Z_partition.begin(); it != Z_partition.end(); ++it) {
		LOG(DEBUG) << "Z: (" << it->first.i1 << ":" << it->first.j1 << ":" << it->first.i2 << ":" << it->first.j2 << ") " << it->second;
	}
	for (auto it = Z_partitionMissing.begin(); it != Z_partitionMissing.end(); ++it) {
		LOG(DEBUG) << "Zmissing: (" << it->first.i1 << ":" << it->first.j1 << ":" << it->first.i2 << ":" << it->first.j2 << ") " << it->second;
	}

	// ============

	// Compute basepair probabilities
	LOG(DEBUG) << "----------------------------------------------";
  computeBasePairProbs(predictor, Z_partition, Z_partition.begin(), Z_partition.end());
	LOG(DEBUG) << "----------------------------------------------";
	computeBasePairProbs(predictor, Z_partition, Z_partitionMissing.begin(), Z_partitionMissing.end());
	LOG(DEBUG) << "----------------------------------------------";

	// build plist
	struct vrna_elem_prob_s plist[structureProbs.size()+1];
	size_t i = 0;
	for (auto it = structureProbs.begin(); it != structureProbs.end(); ++it) {
		LOG(DEBUG) << "Z - prob: " << it->first.first << ":" << it->first.second << ":" << it->first.first << ":" << it->first.second << " - " << getHybridZ(it->first.first, it->first.second, it->first.first, it->first.second, predictor) << " = " << it->second;
		if (it->second > probabilityThreshold) {
			plist[i].i = it->first.first + 1;
			plist[i].j = it->first.second + 1;
			plist[i].p = it->second;
			plist[i].type = 0; // base-pair prob
			i++;
		}
	}

	// create dot plot
	std::string reverseRna2(rna2);
	std::reverse(reverseRna2.begin(), reverseRna2.end());
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
										, const PredictorMfeEns::Site2Z_hash & Z_partition
										, const PredictorMfeEns::Site2Z_hash::const_iterator first
										, const PredictorMfeEns::Site2Z_hash::const_iterator last )
{
	for (auto it = first; it != last; ++it) {
		// full Z of this site only
		Z_type bpProb = it->second * energy.getBoltzmannWeight(energy.getE(it->first.i1, it->first.j1, it->first.i2, it->first.j2, E_type(0)));
		Interaction::BasePair key(it->first.j1, it->first.j2);
		LOG(DEBUG) << " - key: " << it->first.j1 << ":" << it->first.j2;

		// no extension
		if (Z_partition.find( Interaction::Boundary(it->first.i1, it->first.j1, it->first.i2, it->first.j2)) != Z_partition.end()) {
			LOG(DEBUG) << "add " << bpProb;
			structureProbs[key] += (1 / predictor->getZall()) * bpProb;
		}

		// extensions
		for (auto it2 = leftIndex[key].begin(); it2 != leftIndex[key].end(); ++it2) {
			// ensure extension is valid (present in original Z data)
			if (Z_partition.find( Interaction::Boundary(it->first.i1, it2->first, it->first.i2, it2->second)) != Z_partition.end()) {
				if (Z_equal(it->second, 0)) {
					bpProb = getHybridZ(it->first.j1, it2->first, it->first.j2, it2->second, predictor)
								// ED penalty
										 * energy.getBoltzmannWeight(energy.getE(it->first.i1, it2->first, it->first.i2, it2->second, E_type(0)));
				} else {
					bpProb = it->second * getHybridZ(it->first.j1, it2->first, it->first.j2, it2->second, predictor)
								// ED penalty
										 * energy.getBoltzmannWeight(energy.getE(it->first.i1, it2->first, it->first.i2, it2->second, E_type(0)))
										 / energy.getBoltzmannWeight(energy.getE_init());
				}

				structureProbs[key] += (1 / predictor->getZall()) * bpProb;
				LOG(DEBUG) << "add " << bpProb;
			}
		}

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
	LOG(DEBUG) << "compute missing: " << i1 << ":" << j1 << ":" << i2 << ":" << j2;

	if (i1 == j1 && i2 == j2) {
		updateHybridZ(i1, j1, i2, j2, 0);
		return;
	}

	// ignore if invalid interation
	if ((i1 == j1 && i2 < j2) || (i2 == j2 && i1 < j1)) {
		LOG(DEBUG) << "ignoring invalid interaction";
		return;
	}

	// search for known partition
	/*bool foundOuter = false;
	size_t l1, r1, l2, r2;
	for (l1 = i1+1; !foundOuter && l1-- > 0;) {
		for (l2 = i2+1; !foundOuter && l2-- > 0;) {
			for (r1 = j1; !foundOuter && r1 < energy.size1(); r1++) {
				for (r2 = j2; !foundOuter && r2 < energy.size2(); r2++) {
					if (!Z_equal(getHybridZ(l1, r1, l2, r2, predictor), 0)) {
						foundOuter = true;
					}
				}
			}
		}
	}
	r1--;
	r2--;

	if (!foundOuter) {
		//throw std::runtime_error("Could not compute missing Z: no outer region found");
		return;
	}*/

	LOG(DEBUG) << l1 << ":" << r1 << ":" << l2 << ":" << r2;

	Z_type partZ = 0.0;

	if (isFullSeedinRegion(i1, j1, i2, j2, seedHandler)) {
		// Full seed in subregion

		// Case 1
		LOG(DEBUG) << "case1";

		// compute Z(l-i), Z(j-r)

		if (!Z_equal(getHybridZ(l1, i1, l2, i2, predictor), 0)) {
			partZ = (
					getHybridZ(l1, r1, l2, r2, predictor)
				/ getHybridZ(i1, r1, i2, r2, predictor)
			) * energy.getBoltzmannWeight(energy.getE(l1,i1,l2,i2, E_type(0)));
		}
		updateHybridZ(l1, i1, l2, i2, partZ);

		if (!Z_equal(getHybridZ(j1, r1, j2, r2, predictor), 0)) {
			partZ = (
					getHybridZ(l1, r1, l2, r2, predictor)
				/ getHybridZ(l1, j1, l2, j2, predictor)
			) * energy.getBoltzmannWeight(energy.getE(j1,r1,j2,r2, E_type(0)));
		}
		updateHybridZ(j1, r1, j2, r2, partZ);

	} else {

		// check if full seed outside of subregion (left or right)
		if (isFullSeedinRegion(l1, i1, l2, i2, seedHandler)) {
			// full seed left of subregion

			// Case 2.1 (left)
			LOG(DEBUG) << "case2.1 left";

			partZ = (
				  getHybridZ(l1, r1, l2, r2, predictor)
				/ getHybridZ(l1, i1, l2, i2, predictor)
			) * energy.getBoltzmannWeight(energy.getE(i1,j1,i2,j2, E_type(0)))
			* energy.getBoltzmannWeight(energy.getE_init());
			updateHybridZ(i1, j1, i2, j2, partZ);

		} else if (isFullSeedinRegion(j1, r1, j2, r2, seedHandler)) {
			// full seed right of subregion

			// Case 2.1 (right)
			LOG(DEBUG) << "case2.1 right";

			partZ = (
				  getHybridZ(l1, r1, l2, r2, predictor)
				/ getHybridZ(j1, r1, j2, r2, predictor)
			) * energy.getBoltzmannWeight(energy.getE(i1,j1,i2,j2, E_type(0)))
			* energy.getBoltzmannWeight(energy.getE_init());
			updateHybridZ(i1, j1, i2, j2, partZ);

		} else {

			// Case 2.2
			LOG(DEBUG) << "case2.2";

			// check side of k
			size_t k1, k2;
			bool kOnRight;
			if (i1 != l1 || i2 != l2) {
				// k at the left of region
				kOnRight = false;
				k1 = i1;
				k2 = i2;
			} else {
				kOnRight = true;
				k1 = j1;
				k2 = j2;
			}

			// find leftmost overlapping seeds
			std::vector< std::pair <size_t, size_t> > overlappingSeeds = getLeftMostSeedsAtK(k1, k2, seedHandler);

			// loop over overlapping seeds
			for(std::vector< std::pair <size_t, size_t> >::iterator it = overlappingSeeds.begin(); it != overlappingSeeds.end(); ++it) {
				// calculate missing partition function
				const size_t sl1 = seedHandler->getSeedLength1(it->first, it->second);
				const size_t sl2 = seedHandler->getSeedLength2(it->first, it->second);
				if (kOnRight && E_isNotINF(energy.getE_interLeft(k1, r1, k2, r2))) {
					LOG(DEBUG) << getHybridZ(l1, it->first+sl1-1, l2, it->second+sl2-1, predictor);
					LOG(DEBUG) << energy.getBoltzmannWeight(energy.getED1(i1, k1) + energy.getED2(i2, k2));
					LOG(DEBUG) << energy.getBoltzmannWeight(getPartialSeedEnergy(it->first,it->second,k1,r1,k2,r2,seedHandler));
					partZ += (
						getHybridZ(l1, it->first+sl1-1, l2, it->second+sl2-1, predictor)
						* energy.getBoltzmannWeight(energy.getED1(i1, k1) + energy.getED2(i2, k2))
					) / energy.getBoltzmannWeight(getPartialSeedEnergy(it->first,it->second,k1,r1,k2,r2,seedHandler));
				} else if (E_isNotINF(energy.getE_interLeft(l1, k1, l2, k2))) {
					partZ += (
						getHybridZ(it->first, r1, it->second, r2, predictor)
						* energy.getBoltzmannWeight(energy.getED1(k1, j1) + energy.getED2(k2, j2))
					) / energy.getBoltzmannWeight(getPartialSeedEnergy(it->first,it->second,l1,k1,l2,k2,seedHandler));
				}
	    }

			if (kOnRight) {
				updateHybridZ(i1, k1, i2, k2, partZ);
			} else {
				updateHybridZ(k1, j1, k2, j2, partZ);
			}

			// TODO: if no overlapping seeds
			// what now? :|

		}

	}

	LOG(DEBUG) << "Z: " << partZ;

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
									 , SeedHandler* seedHandler )
{
  std::vector< std::pair <size_t, size_t> > seeds;
	size_t maxSeedLength = seedHandler->getConstraint().getBasePairs();

	size_t i1 = (k1 < maxSeedLength) ? 0 : k1 - maxSeedLength;
	size_t i2 = (k2 < maxSeedLength) ? 0 : k2 - maxSeedLength;

	size_t si1 = RnaSequence::lastPos, si2 = RnaSequence::lastPos;
	while( seedHandler->updateToNextSeed(si1, si2, i1, k1, i2, k2) ) {
		// check if k is part of seed
		if (k1 - si1 != k2 - si2) {
			continue;
		}

		// check if a loop-overlapping seed already exists
		bool found = false;
		for(std::vector< std::pair <size_t, size_t> >::iterator it = seeds.begin(); it != seeds.end(); ++it) {
			if (seedHandler->areLoopOverlapping(si1, si2, it->first, it->second)) {
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

	LOG(DEBUG) << "search seed at: " << si1 << ":" << si2 << ":" << i1 << ":" << j1 << ":" << i2 << ":" << j2;

	for (size_t i = 0; i < interaction.basePairs.size(); i++) {
		// get index of current base pair
		size_t s1 = energy.getIndex1(interaction.basePairs[i]);
		size_t s2 = energy.getIndex2(interaction.basePairs[i]);

		if (s1 < i1 || s2 < i2) continue;
		if (s1 > j1 && s2 > j2) break;

		LOG(DEBUG) << "seed: " << i1old << ":" << s1 << ":" << i2old << ":" << s2;

		if (!E_isINF(energy.getE_interLeft(i1old,s1,i2old,s2))) {
			// add hybridization energy
			partE += energy.getE_interLeft(i1old,s1,i2old,s2);
			// store
			i1old = s1;
			i2old = s2;
		}
	}

	LOG(DEBUG) << "partE: " << energy.getBoltzmannWeight(partE);

	return partE;
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
						 , const Z_type partZ )
{
	LOG(DEBUG) << "update HybridZ: " << partZ;
	Interaction::Boundary key(i1,j1,i2,j2);
	Z_partitionMissing[key] = partZ;
}

////////////////////////////////////////////////////////////////////////////

bool
PredictionTrackerBasePairProb::
isSeedBp( const size_t i1, const size_t i2
        , SeedHandler* seedHandler )
{
	size_t si1 = RnaSequence::lastPos, si2 = RnaSequence::lastPos;
	size_t l1 = (i1 <= seedHandler->getConstraint().getBasePairs()) ? 0 : i1 - seedHandler->getConstraint().getBasePairs();
	size_t l2 = (i2 <= seedHandler->getConstraint().getBasePairs()) ? 0 : i2 - seedHandler->getConstraint().getBasePairs();
	while( seedHandler->updateToNextSeed(si1, si2, l1, i1, l2, i2) ) {
		// found seed overlapping basepair
		// check if basepair part of seed
		if (i1 - si1 == i2 - si2) {
			return true;
		}
	}
	return false;
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

	// print outline
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
