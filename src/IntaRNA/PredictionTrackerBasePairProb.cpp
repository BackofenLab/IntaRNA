
#include "IntaRNA/PredictionTrackerBasePairProb.h"

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
	, mfeEnsPredictor(NULL)
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
	, mfeEnsPredictor(NULL)
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
updateZ( PredictorMfeEns *predictor )
{
	mfeEnsPredictor = predictor;

	// write probabilities to streams
	size_t n1 = energy.getAccessibility1().getMaxLength();
	size_t n2 = energy.getAccessibility2().getMaxLength();

	// recursively compute structure probabilities
	float initProb = (mfeEnsPredictor->getZ(0, n1-1, 0, n2-1) * energy.getBoltzmannWeight( energy.getE(0,n1-1,0,n2-1,E_type(0)) )) / mfeEnsPredictor->getOverallZ();
	float prob = computeStructureProbRecursive(0, n1-1, 0, n2-1, initProb);

	// print base-pair probabilities
	for (std::unordered_map<size_t, StructureProb >::const_iterator it = structureProbs.begin(); it != structureProbs.end(); ++it)
	{
		if (it->second.i1 == it->second.j1 && it->second.i2 == it->second.j2) {
			LOG(DEBUG) << it->second.i1 << ":" << it->second.j1 << ":" << it->second.i2 << ":" << it->second.j2 << " = " << it->second.prob;
		}
	}

	(*outStream) << "stuff" << "\n";

	outStream->flush();

	if (deleteOutStream) {
		// clean up if file pointers were created in constructor
		deleteOutputStream( outStream );
	}
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

float
PredictionTrackerBasePairProb::
computeStructureProbRecursive( const size_t k1, const size_t l1
					, const size_t k2, const size_t l2, const float structProb )
{
	size_t n1 = energy.getAccessibility1().getMaxLength();
	size_t n2 = energy.getAccessibility2().getMaxLength();
	size_t w1 = l1 - k1;
	size_t w2 = l2 - k2;

	if (w1 == 0 || w2 == 0) {
		// store base-pair prob
		size_t key = generateMapKey(k1, l1, k2, l2);
		float prob = computeSubStructureProb(k1, l1, k2, l2, k1, l1, k2, l2, structProb);
		if ( structureProbs.find(key) == structureProbs.end() ) {
			// create new entry
			StructureProb sProb;
			sProb.i1 = k1;
			sProb.j1 = l1;
			sProb.i2 = k2;
			sProb.j2 = l2;
			sProb.prob = prob;
			structureProbs[key] = sProb;
		} else {
			// update entry
			StructureProb & sProb = structureProbs[key];
			sProb.prob += prob;
		}
	} else {
		// calculate sub probs
		float prob = 0.0;

		for (size_t i1 = 0; i1 <= n1-w1; i1++) {
			for (size_t i2 = 0; i2 <= n2-w2; i2++) {
				size_t j1 = i1 + w1 - 1;
				size_t j2 = i2 + w2 - 1;

				float subProb = computeSubStructureProb(i1, j1, i2, j2, k1, l1, k2, l2, structProb);
				prob += subProb;

				computeStructureProbRecursive(i1, j1, i2, j2, subProb);
			}
		}
	}



}

////////////////////////////////////////////////////////////////////////////

float
PredictionTrackerBasePairProb::
computeSubStructureProb( const size_t i1, const size_t j1
					, const size_t i2, const size_t j2
					, const size_t k1, const size_t l1
					, const size_t k2, const size_t l2, const float structProb )
{
	float prob = 0;
	Z_type partZ = mfeEnsPredictor->getZ(i1, j1, i2, j2);
	Z_type partZprime = mfeEnsPredictor->getZ(k1, l1, k2, l2);

	// case 1 (i=k, j=l)
	if (i1 == k1 && i2 == k2 && j1 == l1 && j2 == l2) {
		prob = (partZ * energy.getBoltzmannWeight( energy.getE(i1,j1,i2,j2,E_type(0)) )) / mfeEnsPredictor->getOverallZ();
	}

	// case 2 (i>k, j=l)
	if ((i1 > k1 || i2 > k2) && j1 == l1 && j2 == l2) {
		if (!Z_equal(partZprime, 0)) {
			prob = (structProb * partZ * energy.getBoltzmannWeight(energy.getE_interLeft(k1,i1,k2,i2))) / partZprime;
		}
	}

	// case 3 (i=k, j<l)
	if (i1 == k1 && i2 == k2 && (j1 < l1 || j2 < l2)) {
		if (!Z_equal(partZprime, 0)) {
			prob = (structProb * partZ * energy.getBoltzmannWeight(energy.getE_interLeft(j1,l1,j2,l2))) / partZprime;
		}
	}

	// case 4 (i<k, j<l)
	if ((i1 < k1 || i2 < k2) && (j1 < l1 || j2 < l2)) {
		if (!Z_equal(partZprime, 0)) {
			prob = (structProb * partZ * energy.getBoltzmannWeight(energy.getE_interLeft(k1,i1,k2,i2)) * energy.getBoltzmannWeight(energy.getE_interLeft(j1,l1,j2,l2))) / partZprime;
		}
	}

	return prob;
}

////////////////////////////////////////////////////////////////////////////


} // namespace
