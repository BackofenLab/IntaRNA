
#ifndef INTARNA_PREDICTORMFE2DSEEDEXTENSION_H_
#define INTARNA_PREDICTORMFE2DSEEDEXTENSION_H_

#include "IntaRNA/PredictorMfe2d.h"
#include "IntaRNA/SeedHandlerIdxOffset.h"
#include <vector>

namespace IntaRNA {

/**
 * Implements seed-based, space-efficient interaction prediction.
 *
 * Note, for each seed start (i1,i2) only the mfe seed is considered for the
 * overall interaction computation instead of considering all possible seeds
 * starting at (i1,i2).
 *
 * Note, strongly underestimates Zall since some recursion cases (that are
 * excluded from Zall updating) might produce interaction duplicates (which
 * would lead to an overestiamte of Zall).
 *
 * @author Frank Gelhausen
 * @author Martin Raden
 *
 */
class PredictorMfe2dSeedExtension: public PredictorMfe2d {

protected:

	//! matrix type to hold the mfe energies for interaction site starts
	typedef PredictorMfe2d::E2dMatrix E2dMatrix;

public:


	/**
	 * Constructs a predictor and stores the energy and output handler
	 *
	 * @param energy the interaction energy handler
	 * @param output the output handler to report mfe interactions to
	 * @param predTracker the prediction tracker to be used or NULL if no
	 *         tracking is to be done; if non-NULL, the tracker gets deleted
	 *         on this->destruction.
	 * @param seedHandler the seed handler to be used
	 */
	PredictorMfe2dSeedExtension(
			const InteractionEnergy & energy
			, OutputHandler & output
			, PredictionTracker * predTracker
			, SeedHandler * seedHandler );


	/**
	 * data cleanup
	 */
	virtual ~PredictorMfe2dSeedExtension();


	/**
	 * Computes the mfe for the given sequence ranges (i1-j1) in the first
	 * sequence and (i2-j2) in the second sequence and reports it to the output
	 * handler.
	 *
	 * Each considered interaction contains a seed according to the seed handler
	 * constraints.
	 *
	 * @param r1 the index range of the first sequence interacting with r2
	 * @param r2 the index range of the second sequence interacting with r1
	 *
	 */
	virtual
	void
	predict( const IndexRange & r1 = IndexRange(0,RnaSequence::lastPos)
			, const IndexRange & r2 = IndexRange(0,RnaSequence::lastPos) );


protected:


	//! access to the interaction energy handler of the super class
	using PredictorMfe2d::energy;

	//! access to the output handler of the super class
	using PredictorMfe2d::output;

	//! the seed handler (with idx offset)
	SeedHandlerIdxOffset seedHandler;

	//! energy of all interaction hybrids that start on the left side of the seed
	E2dMatrix hybridE_left;

	//! energy of all interaction hybrids that start on the right side of the seed
	E2dMatrix hybridE_right;

protected:

	//! One-based complementary-pair indices; zero denotes an unavailable pair.
	Matrix<size_t> extensionPairIndex;
	//! Seed-independent loop energies and flags; E_INF is a cacheable result.
	std::vector<E_type> extensionLoopEnergy;
	std::vector<unsigned char> extensionLoopKnown;
	//! Indexed positive loop displacements on the two increasing local strands.
	size_t extensionLoopWidth1 = 0, extensionLoopWidth2 = 0;

	/** Discard cached values before a prediction changes the local offsets. */
	void resetExtensionCache();
	/** Prepare reusable seed-independent pair/loop
	 * information for this prediction. Cache storage is limited to 16 MiB;
	 * domains exceeding that limit retain direct energy evaluation.
	 * @param n number of active positions on sequence 1
	 * @param m number of active positions on reversed sequence 2
	 */
	void prepareExtensionWorkspace(size_t n, size_t m);
	/** Check pair feasibility in the current prediction's local coordinates.
	 * @param i1 position in sequence 1
	 * @param i2 position in reversed sequence 2
	 * @return whether the pair is complementary and accessible
	 */
	bool extensionPairIsValid(size_t i1, size_t i2) const;
	/** Reuse a seed-independent loop energy, with direct evaluation outside
	 * the bounded cache. Both pairs use increasing local strand coordinates.
	 * @param i1 left pair's position in sequence 1
	 * @param j1 right pair's position in sequence 1
	 * @param i2 left pair's position in reversed sequence 2
	 * @param j2 right pair's position in reversed sequence 2
	 * @return loop energy in hundredths of kcal/mol, or E_INF
	 */
	E_type extensionLoopE(size_t i1, size_t j1, size_t i2, size_t j2);

	/**
	 * Computes all entries of the hybridE matrix for interactions starting in
	 * i1 and i2 and report all valid interactions to updateOptima()
	 *
	 * @param j1 start of the interaction within seq 1
	 * @param j2 start of the interaction within seq 2
	 *
	 */
	void
	fillHybridE_left( const size_t j1, const size_t j2 );

	/**
	 * Computes all entries of the hybridE matrix for interactions starting in
	 * i1 and i2 and report all valid interactions to updateOptima()
	 *
	 * @param i1 start of the interaction within seq 1
	 * @param i2 start of the interaction within seq 2
	 *
	 */
	void
	fillHybridE_right( const size_t i1, const size_t i2 );

	/**
	 * Fills a given interaction (boundaries given) with the according
	 * hybridizing base pairs using hybridE_seed.
	 * @param interaction IN/OUT the interaction to fill
	 */
	virtual
	void
	traceBack( Interaction & interaction );

	void
	printMatrix( const E2dMatrix & matrix );

};

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


inline
bool
PredictorMfe2dSeedExtension::extensionPairIsValid(size_t i1, size_t i2) const
{
	return extensionPairIndex.size1() ? extensionPairIndex(i1,i2) != 0
			: energy.areComplementary(i1,i2);
}

inline
void
PredictorMfe2dSeedExtension::
printMatrix( const E2dMatrix & matrix )
{
	for (int i = 0; i < matrix.size1(); i++) {
		std::cout << "| ";
		for (int j = 0; j < matrix.size2(); j++) {
			std::cout << matrix(i, j) << " | ";
		}
		std::cout << std::endl;
	}
}

} // namespace

#endif /* INTARNA_PREDICTORMFE2DSEEDEXTENSION_H_ */
