
#include "IntaRNA/PredictionTrackerBasePairProb.h"
#include "IntaRNA/PredictorMfeEns.h"

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

	LOG(DEBUG) << "Hello Tracker";
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
	// write probabilities to streams

  /*

	// overall structure probability
	size_t key = generateMapKey(i1, j1, i2, j2);
	StructureProb sProb;
	sProb.i1 = i1;
	sProb.j1 = j1;
	sProb.i2 = i2;
	sProb.j2 = j2;
	sProb.prob = (Z_partitions[key].partZ * energy.getBoltzmannWeight( energy.getE(i1,j1,i2,j2,E_type(0)) )) / overallZ;
	structureProbs[key] = sProb;

	// last seen partition
	size_t i1last = i1;
	size_t j1last = j1;
	size_t i2last = i2;
	size_t j2last = j2;
	size_t lastKey = key;

	// compute structure probabilities
	for (size_t l1 = j1+1; l1 --> i1;) {
		for (size_t l2 = j2+1; l2 --> i2;) {
			for (size_t k1 = i1; k1 <= l1; k1++) {
				for (size_t k2 = i2; k2 <= l2; k2++) {
					// check if hybrid partition available
					key = generateMapKey(k1, l1, k2, l2);

					// compute probability and update last seen
					computeStructureProb(k1, l1, k2, l2, i1last, j1last, i2last, j2last, key, lastKey);
					i1last = k1;
					j1last = l1;
					i2last = k2;
					j2last = l2;
					lastKey = key;
				}
			}
		}
	}

	for (std::unordered_map<size_t, StructureProb >::const_iterator it = structureProbs.begin(); it != structureProbs.end(); ++it)
	{
		if (it->second.i1 == it->second.j1 && it->second.i2 == it->second.j2) {
			LOG(DEBUG) << it->second.i1 << ":" << it->second.j1 << ":" << it->second.i2 << ":" << it->second.j2 << " = " << it->second.prob;
		}
	}
	*/

	(*outStream) << "stuff" << "\n";

	outStream->flush();

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
updateZ( const PredictorMfeEns * predictor )
{
	LOG(DEBUG) << predictor->getOverallZ();
	/*
	// update overallZ
	overallZ += partZ * energy.getBoltzmannWeight( energy.getE(i1,j1,i2,j2,E_type(0)) );

	// store partial hybridZ
	size_t key = generateMapKey(i1, j1, i2, j2);
	if ( Z_partitions.find(key) == Z_partitions.end() ) {
		// create new entry
		ZPartition zPartition;
		zPartition.i1 = i1;
		zPartition.j1 = j1;
		zPartition.i2 = i2;
		zPartition.j2 = j2;
		zPartition.partZ = partZ;
		Z_partitions[key] = zPartition;
	} else {
		// update entry
		ZPartition & zPartition = Z_partitions[key];
		zPartition.partZ += partZ;
	}
	*/
}

////////////////////////////////////////////////////////////////////////////

size_t
PredictionTrackerBasePairProb::
generateMapKey( const size_t i1, const size_t j1
					, const size_t i2, const size_t j2 )
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

void
PredictionTrackerBasePairProb::
computeStructureProb( const size_t i1, const size_t j1
					, const size_t i2, const size_t j2
					, const size_t i1last, const size_t j1last
					, const size_t i2last, const size_t j2last, const size_t key, const size_t lastKey )
{
	Z_type partZ = 0;
	Z_type partZprime = 0;
	if ( Z_partitions.find(key) != Z_partitions.end() ) {
		partZ = Z_partitions[key].partZ;
	}
	if ( Z_partitions.find(lastKey) != Z_partitions.end() ) {
		partZprime = Z_partitions[lastKey].partZ;
	}

	float ps = (partZ * energy.getBoltzmannWeight( energy.getE(i1,j1,i2,j2,E_type(0)) )) / overallZ;

	float psPrime = (partZprime * energy.getBoltzmannWeight( energy.getE(i1last,j1last,i2last,j2last,E_type(0)) )) / overallZ;

	float b = 0;
	if (!Z_equal(partZprime, 0)) {
		b = (psPrime * partZ * energy.getBoltzmannWeight(energy.getE_interLeft(i1last,j1last,i1,j1)) * energy.getBoltzmannWeight(energy.getE_interLeft(i2,j2,i2last,j2last))) / partZprime;
	}

	StructureProb sProb;
	sProb.i1 = i1;
	sProb.j1 = j1;
	sProb.i2 = i2;
	sProb.j2 = j2;
	sProb.prob = ps + b;
	structureProbs[key] = sProb;
}

////////////////////////////////////////////////////////////////////////////


} // namespace
