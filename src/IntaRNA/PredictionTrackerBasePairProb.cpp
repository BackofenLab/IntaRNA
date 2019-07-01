
#include "IntaRNA/PredictionTrackerBasePairProb.h"

extern "C" {
	#include <ViennaRNA/plotting/probabilities.h>
}

#include <boost/foreach.hpp>
#include <boost/regex.hpp>

namespace IntaRNA {

//////////////////////////////////////////////////////////////////////

PredictionTrackerBasePairProb::
PredictionTrackerBasePairProb(
		const InteractionEnergy & energy
		, const std::string & outStreamName
	)
 :	PredictionTracker()
	, energy(energy)
	, outStream(NULL)
	, deleteOutStream(true)
	, overallZ(0.0)
{
	// open stream
	if (!outStreamName.empty()) {
		outStream = newOutputStream( outStreamName );
		if (outStream == NULL) {
			throw std::runtime_error("PredictionTrackerBasePairProb() : could not open output stream '"+outStreamName+"' for writing");
		}
	}

}

//////////////////////////////////////////////////////////////////////

PredictionTrackerBasePairProb::
PredictionTrackerBasePairProb(
		const InteractionEnergy & energy
		, std::ostream & outStream
	)
 :	PredictionTracker()
	, energy(energy)
	, outStream(&outStream)
	, deleteOutStream(false)
	, overallZ(0.0)
{
}

//////////////////////////////////////////////////////////////////////

PredictionTrackerBasePairProb::
~PredictionTrackerBasePairProb()
{
	if (deleteOutStream) {
		// clean up if file pointers were created in constructor
		deleteOutputStream( outStream );
	}
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
updateZ( PredictorMfeEns *predictor )
{

	// sequence strings
	const std::string & rna1 = energy.getAccessibility1().getSequence().asString();
	const std::string & rna2 = energy.getAccessibility2().getSequence().asString();

	// write probabilities to streams
	size_t s1 = energy.size1();
	size_t s2 = energy.size2();
	size_t n1 = energy.getAccessibility1().getMaxLength();
	size_t n2 = energy.getAccessibility2().getMaxLength();

	// calculate initial probabilities (for max window length)
	for (size_t i1 = 0; i1 < s1-n1+1; i1++) {
		for (size_t i2 = 0; i2 < s2-n2+1; i2++) {
			StructureProb sProb;
			sProb.i1 = i1;
			sProb.j1 = i1 + n1 - 1;
			sProb.i2 = i2;
			sProb.j2 = i2 + n2 - 1;
			sProb.prob = predictor->getZ(i1, i1 + n1 - 1, i2, i2 + n2 - 1) / predictor->getOverallHybridZ();
			structureProbs[generateMapKey(i1, i1 + n1 -1, i2, i2 + n2 - 1)] = sProb;
		}
	}

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
					float prob = 0.0;
					// loop over combinations of (0..i), (j..|s|)
					for (size_t l1 = 0; l1 <= i1; l1++) {
						for (size_t l2 = 0; l2 <= i2; l2++) {
							for (size_t r1 = j1; r1 < s1; r1++) {
								// ensure maximal loop length
								if (r1-l1 > energy.getMaxInternalLoopSize1()+1) break;
								for (size_t r2 = j2; r2 < s2; r2++) {
									// ensure maximal loop length
									if (r2-l2 > energy.getMaxInternalLoopSize2()+1) break;

									if (l1 == i1 && l2 == i2 && j1 == r1 && j2 == r2) {
										prob += predictor->getZ(i1, j1, i2, j2) / predictor->getOverallHybridZ();
									} else {
										// get outer probability
										size_t key = generateMapKey(l1, r1, l2, r2);
										if ( structureProbs.find(key) == structureProbs.end() ) {
											throw std::runtime_error("PredictionTrackerBasePairProb() : could not find outer probability for interaction " + toString(l1) + ":" + toString(r1) + ":" + toString(l2) + ":" + toString(r2));
										}

										if (!Z_equal(predictor->getZ(l1, r1, l2, r2), 0)) {
											prob += (structureProbs[key].prob
													     * energy.getBoltzmannWeight(energy.getE_interLeft(l1,i1,l2,i2))
															 * predictor->getZ(i1, j1, i2, j2)
															 * energy.getBoltzmannWeight(energy.getE_interLeft(j1,r1,j2,r2))
												     ) / predictor->getZ(l1, r1, l2, r2);
										}
									}

								} // r2
							} // r1
						} // l2
					} // l1

					// store probability
					StructureProb sProb;
					sProb.i1 = i1;
					sProb.j1 = j1;
					sProb.i2 = i2;
					sProb.j2 = j2;
					sProb.prob = prob;
					structureProbs[generateMapKey(i1, j1, i2, j2)] = sProb;

				} // i2
			} // i1
		} // w2
	} // w1

	// create plist
	struct vrna_elem_prob_s plist1[s1*s2+1];
	struct vrna_elem_prob_s plist2[s1*s2+1];
	size_t i = 0;

	for (std::unordered_map<size_t, StructureProb >::const_iterator it = structureProbs.begin(); it != structureProbs.end(); ++it)
	{
		if (it->second.i1 == it->second.j1 && it->second.i2 == it->second.j2 && it->second.prob > 0.0001) {
			LOG(DEBUG) << it->second.i1 << ":" << it->second.j1 << ":" << it->second.i2 << ":" << it->second.j2 << " = " << it->second.prob;
			plist1[i].i = it->second.i1;
			plist1[i].j = s1 + it->second.i2;
			plist1[i].p = it->second.prob;
			plist1[i].type = 0; // base-pair prob
			plist2[i].i = it->second.i1;
			plist2[i].j = s1 + it->second.i2;
			plist2[i].p = it->second.prob;
			plist2[i].type = 0; // base-pair prob
			i++;
		}
	}

	// create dot plot
	std::string reverseRna2(rna2);
	std::reverse(reverseRna2.begin(), reverseRna2.end());
	char *rna = strdup((rna1 + reverseRna2).c_str());
	std::string outName = "plot";
	std::string comment = "";
	PS_dot_plot_list(rna, &outName[0u], plist1, plist2, &comment[0u]);

}

////////////////////////////////////////////////////////////////////////////

size_t
PredictionTrackerBasePairProb::
generateMapKey( const size_t i1, const size_t j1
					, const size_t i2, const size_t j2 ) const
{
	size_t maxLength = std::max(energy.getAccessibility1().getMaxLength(), energy.getAccessibility2().getMaxLength());
	size_t key = 0;
	key += i1;
	key += j1 * pow(maxLength, 1);
	key += i2 * pow(maxLength, 2);
	key += j2 * pow(maxLength, 3);
	return key;
}

////////////////////////////////////////////////////////////////////////////

} // namespace
