
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

	// initialize Z_partition
	const PredictorMfeEns::Site2Z_hash & Z_partition = predictor->getZPartition();
	// cleanup additional data
	Z_partitionMissing.clear();

	// TODO fuse with leftIndex generation
	// TODO fuse with best interaction identification
	// compute missing Z values in case of a seed based prediction
	if (seedHandler != NULL) {
		// TODO replace with something iterating over elements of Z_partition
		// loop over window lengths
		for (size_t w1 = n1; w1 > 0; w1--) {
			for (size_t w2 = n2; w2 > 0; w2--) {
				// skip initial probabilities
				if (w1 == n1 && w2 == n2) continue;
				// shift window over sequence length
				for (size_t i1 = 0; i1 < s1-w1+1; i1++) {
					for (size_t i2 = 0; i2 < s2-w2+1; i2++) {
						size_t j1 = i1 + w1 - 1;
						size_t j2 = i2 + w2 - 1;

						// if seed-based prediction, compute missing Z values
						if (Z_equal(getHybridZ(i1, j1, i2, j2, predictor), 0)) {
							LOG(DEBUG) << "missing Z at " << i1 << ":" << j1 << ":" << i2 << ":" << j2;
							computeMissingZ(i1, j1, i2, j2, predictor, seedHandler);
						} else {
							LOG(DEBUG) << "found Z at " << i1 << ":" << j1 << ":" << i2 << ":" << j2 << " with Z = " << getHybridZ(i1, j1, i2, j2, predictor);
						}
					} // i1
				} // i2
			} // w1
		} // w2
	}


	// build left side index of all true
	// TODO : iterate o
	for (auto it = Z_partition.begin(); it != Z_partition.end(); ++it) {
		size_t i1 = it->first.i1;
		size_t i2 = it->first.i2;
		Interaction::BasePair key(i1, i2);
		if ( leftIndex.find(key) == leftIndex.end() ) {
			std::vector<Interaction::BasePair> basePairs;
			for (auto it2 = Z_partition.begin(); it2 != Z_partition.end(); ++it2) {
				if (it2->first.i1 == i1 && it2->first.i2 == i2
					&& it2->first.j1 != i1 && it2->first.j2 != i2)
				{
					Interaction::BasePair entry(it2->first.j1, it2->first.j2);
					basePairs.push_back(entry);
				}
			}
			leftIndex[key] = basePairs;
		}
	}

	// TODO iterate all seeds : handle seed-BP computation if needed (update leftindex)

			// check if basepair is part of seed
			/*if (seedHandler != NULL) {
				auto seeds = getLeftMostSeedsAtK(k1, k2, seedHandler);
				if (seeds.size() > 0) {
					Interaction::Boundary probKey(k1, k1, k2, k2);
					const size_t sl1 = seedHandler->getSeedLength1(seeds[0].first, seeds[0].second);
					const size_t sl2 = seedHandler->getSeedLength2(seeds[0].first, seeds[0].second);
					const size_t sj1 = seeds[0].first+sl1-1;
					const size_t sj2 = seeds[0].second+sl2-1;
					structureProbs[probKey] = (getHybridZ(seeds[0].first, sj1, seeds[0].second, sj2, predictor) / predictor->getZall());
					continue;
				}
			}*/

	// TODO iterate over both
	for (auto it = Z_partition.begin(); it != Z_partition.end(); ++it) {
		// full Z of this site only
		Z_type bpProb = it->second * energy.getBoltzmannWeight(energy.getE(it->first.i1, it->first.j1, it->first.i2, it->first.j2, E_type(0)));;
		// TODO : refactor to function + iterate over both Z_partition and Z_partitionMissing
		Z_type extendedProb = 0.0;
		Interaction::BasePair key(it->first.j1, it->first.j2);
		for (auto it2 = leftIndex[key].begin(); it2 != leftIndex[key].end(); ++it2) {
			// ensure extension is a valid (present in original Z data)
			if (Z_partition.find( Interaction::Boundary(it->first.i1, it2->first, it->first.i2, it2->second) != Z_partition.end())) {
				extendedProb
						+= getHybridZ(it->first.j1, it2->first, it->first.j2, it2->second, predictor)
							// ED penalty
			             * energy.getBoltzmannWeight(energy.getE(it->first.i1, it2->first, it->first.i2, it2->second, E_type(0)));
			}
		}

		if (!Z_equal(extendedProb, 0)) {
			// all right extensions
			bpProb += it->second * extendedProb / energy.getBoltzmannWeight(energy.getE_init());
		}

		// final partition function
		// todo reuse "key" as key
		Interaction::Boundary probKey(it->first.j1, it->first.j1, it->first.j2, it->first.j2);
		structureProbs[probKey] += (1 / predictor->getZall()) * bpProb;
	}

	// create plist
	struct vrna_elem_prob_s plist[structureProbs.size()+1];
	size_t i = 0;
	Z_type maxZ = 0.0;
	Interaction::Boundary interactionBoundary;

	// search best interaction boundary
  for (auto it = Z_partition.begin(); it != Z_partition.end(); ++it) {
		Z_type Zstruct = it->second * energy.getBoltzmannWeight(energy.getE(it->first.i1, it->first.j1, it->first.i2, it->first.j2, E_type(0)));
		if (Zstruct > maxZ) {
			maxZ = Zstruct;
			interactionBoundary = it->first;
		}
	}

	// build plist
	for (auto it = structureProbs.begin(); it != structureProbs.end(); ++it) {
		LOG(DEBUG) << "Z - prob: " << it->first.i1 << ":" << it->first.j1 << ":" << it->first.i2 << ":" << it->first.j2 << " - " << getHybridZ(it->first.i1, it->first.j1, it->first.i2, it->first.j2, predictor) << " = " << it->second;
		if (it->first.i1 == it->first.j1 && it->first.i2 == it->first.j2 && it->second > probabilityThreshold) {
			plist[i].i = it->first.i1 + 1;
			plist[i].j = it->first.i2 + 1;
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
computeMissingZ( const size_t i1, const size_t j1
							, const size_t i2, const size_t j2
							, PredictorMfeEns *predictor
							, SeedHandler* seedHandler )
{
	// ignore if invalid interation
	if ((i1 == j1 && i2 < j2) || (i2 == j2 && i1 < j1)) {
		LOG(DEBUG) << "ignoring invalid interaction";
		return;
	}

	// search for known partition
	bool foundOuter = false;
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
	}

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
		if (s1 > j1 && s2 < j2) break;

		LOG(DEBUG) << "seed: " << i1old << ":" << s1 << ":" << i2old << ":" << s2;

		if (!E_isINF(energy.getE_interLeft(i1old,s1,i2old,s2))) {
			// add hybridization energy
			partE += energy.getE_interLeft(i1old,s1,i2old,s2);
			// store
			i1old = s1;
			i2old = s2;
		}
	}

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
getBasePairProb( const size_t i1, const size_t j1
				       , const size_t i2, const size_t j2
					     , PredictorMfeEns *predictor)
{
	Interaction::Boundary key(i1, j1, i2, j2);
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
