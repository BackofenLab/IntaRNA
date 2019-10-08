
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
			leftIndex[left].push_back(right);
			// create right index
			if (seedHandler != NULL) {
				rightIndex[right].push_back(left);
			}
		}

	} // it (Z_partition)

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

			// trace all seed base pairs
			Interaction interaction = Interaction(energy.getAccessibility1().getSequence(), energy.getAccessibility2().getAccessibilityOrigin().getSequence());
			const size_t sl1 = seedHandler->getSeedLength1(si1, si2);
			const size_t sl2 = seedHandler->getSeedLength2(si1, si2);
			interaction.basePairs.push_back( energy.getBasePair(si1, si2) );
			seedHandler->traceBackSeed( interaction, si1, si2 );
			interaction.basePairs.push_back( energy.getBasePair(si1+sl1-1, si2+sl2-1) );

			// iterate all internal seed base pairs
			for (size_t i = 0; i < interaction.basePairs.size(); i++) {
				// get index of current base pair
				size_t sk1 = energy.getIndex1(interaction.basePairs[i]);
				size_t sk2 = energy.getIndex2(interaction.basePairs[i]);

				// no extension
				if (Z_equal(getHybridZ(sk1, si1+sl1-1, sk2, si2+sl2-1, predictor), 0)) {
					computeMissingZ(sk1, si1+sl1-1, sk2, si2+sl2-1, si1, si1+sl1-1, si2, si2+sl2-1, predictor, seedHandler);
					Interaction::BasePair key(sk1, sk2);
					Interaction::BasePair entry(si1+sl1-1, si2+sl2-1);
					if (std::find(leftIndex[key].begin(), leftIndex[key].end(), entry) == leftIndex[key].end()) {
						leftIndex[key].push_back(entry);
					}
				}

				// right extensions
				// iterate all r=leftIndex[si] to computeMissingZ(si,sk,r) and update leftIndex[sk].push(r)
				Interaction::BasePair seedKey(si1, si2);
				for (auto it = leftIndex[seedKey].begin(); it != leftIndex[seedKey].end(); ++it) {
					if (sk1 < it->first && sk2 < it->second && Z_equal(getHybridZ(sk1, it->first, sk2, it->second, predictor), 0)) {
						computeMissingZ(sk1, it->first, sk2, it->second, si1, it->first, si2, it->second, predictor, seedHandler);
						Interaction::BasePair key(sk1, sk2);
						Interaction::BasePair entry(it->first, it->second);
						if (std::find(leftIndex[key].begin(), leftIndex[key].end(), entry) == leftIndex[key].end()) {
							leftIndex[key].push_back(entry);
						}
					}
				}

				// left extensions
				// iterate all l=rightIndex[sj] to computeMissingZ(l,sk,sj) and update leftIndex[l].push(sk)  ///////////rightIndex[sk].push(l)
				Interaction::BasePair seedKey2(si1+sl1-1, si2+sl2-1);
				for (auto it = rightIndex[seedKey2].begin(); it != rightIndex[seedKey2].end(); ++it) {
					if (it->first < sk1 && it->second < sk2 && Z_equal(getHybridZ(it->first, sk1, it->second, sk2, predictor), 0)) {
						computeMissingZ(it->first, sk1, it->second, sk2, it->first, si1+sl1-1, it->second, si2+sl2-1, predictor, seedHandler);
						Interaction::BasePair entry(sk1, sk2);
						Interaction::BasePair key(it->first, it->second);
						if (std::find(leftIndex[key].begin(), leftIndex[key].end(), entry) == leftIndex[key].end()) {
							leftIndex[key].push_back(entry);
						}
					}
				}

				// init basepair hybridZ
				updateHybridZ(sk1, sk1, sk2, sk2, 0);

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
	for (auto it = first; it != last; ++it) {
		Z_type bpProb = it->second * energy.getBoltzmannWeight(energy.getE(it->first.i1, it->first.j1, it->first.i2, it->first.j2, E_type(0)));
		Interaction::BasePair key(it->first.j1, it->first.j2);

		// no extension
		if (predictor->getZPartition().find( Interaction::Boundary(it->first.i1, it->first.j1, it->first.i2, it->first.j2)) != predictor->getZPartition().end()) {
			structureProbs[key] += (1 / predictor->getZall()) * bpProb;
		}

		// extensions
		for (auto it2 = leftIndex[key].begin(); it2 != leftIndex[key].end(); ++it2) {

			// ensure extension is valid (present in original Z data)
			if (predictor->getZPartition().find( Interaction::Boundary(it->first.i1, it2->first, it->first.i2, it2->second)) != predictor->getZPartition().end()) {

				// exclude wrong extensions of type /|/ given /|
				if (it2->first - it->first.i1 == it2->second - it->first.i2 &&
				  !(it2->first - it->first.j1 == it2->second - it->first.j2)) {
				  continue;
				}

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

	// ignore if invalid interation
	if ((i1 == j1 && i2 < j2) || (i2 == j2 && i1 < j1)) {
		LOG(DEBUG) << "ignoring invalid interaction";
		return;
	}

	Z_type partZ = 0.0;

	if (!isFullSeedinRegion(i1, j1, i2, j2, seedHandler)) {

		// check if full seed outside of subregion (left or right)
		if (isFullSeedinRegion(l1, i1, l2, i2, seedHandler)) {
			// full seed left of subregion

			// Case 2.1 (left)
			size_t seedCount = countNonOverlappingSeeds(l1, r1, l2, r2, seedHandler);

			partZ = (
				  getHybridZ(l1, r1, l2, r2, predictor)
				/ getHybridZ(l1, i1, l2, i2, predictor)
			) * energy.getBoltzmannWeight(energy.getE(i1,j1,i2,j2, E_type(0)))
			* energy.getBoltzmannWeight(energy.getE_init())
			/ seedCount;
			updateHybridZ(i1, j1, i2, j2, partZ);

		} else if (isFullSeedinRegion(j1, r1, j2, r2, seedHandler)) {
			// full seed right of subregion

			// Case 2.1 (right)
			size_t seedCount = countNonOverlappingSeeds(l1, r1, l2, r2, seedHandler);

			partZ = (
				  getHybridZ(l1, r1, l2, r2, predictor)
				/ getHybridZ(j1, r1, j2, r2, predictor)
			) * energy.getBoltzmannWeight(energy.getE(i1,j1,i2,j2, E_type(0)))
			* energy.getBoltzmannWeight(energy.getE_init())
			/ seedCount;
			updateHybridZ(i1, j1, i2, j2, partZ);

		} else {

			// Case 2.2

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
			std::vector< std::pair <size_t, size_t> > overlappingSeeds = getLeftMostSeedsAtK(j1, j2, seedHandler);

			// loop over overlapping seeds
			for(std::vector< std::pair <size_t, size_t> >::iterator it = overlappingSeeds.begin(); it != overlappingSeeds.end(); ++it) {
				// calculate missing partition function
				const size_t sl1 = seedHandler->getSeedLength1(it->first, it->second);
				const size_t sl2 = seedHandler->getSeedLength2(it->first, it->second);
				if (kOnRight) {
					if (it->first <= l1 && it->second <= l2) {
						// region completely inside seed
						partZ = energy.getBoltzmannWeight(getPartialSeedEnergy(it->first,it->second,l1,k1,l2,k2,seedHandler));
					} else {
						partZ += (
							getHybridZ(l1, it->first+sl1-1, l2, it->second+sl2-1, predictor)
							* energy.getBoltzmannWeight(energy.getED1(i1, k1) + energy.getED2(i2, k2))
						) / energy.getBoltzmannWeight(getPartialSeedEnergy(it->first,it->second,k1,it->first+sl1-1,k2,it->second+sl2-1,seedHandler));
					}
				} else {
					if (j1 <= it->first+sl1-1 && j2 <= it->second+sl2-1) {
						// region completely inside seed
						partZ = energy.getBoltzmannWeight(getPartialSeedEnergy(it->first,it->second,k1,j1,k2,j2,seedHandler));
					} else {
						partZ += (
							getHybridZ(it->first, j1, it->second, j2, predictor)
							* energy.getBoltzmannWeight(energy.getED1(k1, j1) + energy.getED2(k2, j2))
						) / energy.getBoltzmannWeight(getPartialSeedEnergy(it->first,it->second,it->first,k1,it->second,k2,seedHandler));
					}
				}
	    }

			if (kOnRight) {
				updateHybridZ(i1, k1, i2, k2, partZ);
			} else {
				updateHybridZ(k1, j1, k2, j2, partZ);
			}

		}

	}

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
									 , SeedHandler* seedHandler )
{
  std::vector< std::pair <size_t, size_t> > seeds;
	size_t maxSeedLength = seedHandler->getConstraint().getBasePairs();

	size_t i1 = (k1 < maxSeedLength) ? 0 : k1 - maxSeedLength + 1;
	size_t i2 = (k2 < maxSeedLength) ? 0 : k2 - maxSeedLength + 1;

	size_t si1 = RnaSequence::lastPos, si2 = RnaSequence::lastPos;
	while( seedHandler->updateToNextSeed(si1, si2, i1, k1, i2, k2) ) {
		// check if k is part of seed
		if (k1 - si1 != k2 - si2) {
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
						 , const Z_type partZ )
{
	Interaction::Boundary key(i1,j1,i2,j2);
	Z_partitionMissing[key] = partZ;
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
