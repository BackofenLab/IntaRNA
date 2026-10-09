
#ifndef INTARNA_PREDICTORMFEENS2D_H_
#define INTARNA_PREDICTORMFEENS2D_H_

#include "IntaRNA/PredictorMfeEns.h"
#include "IntaRNA/BasePairProbabilities.h"
#include "IntaRNA/Interaction.h"

#include "IntaRNA/Matrix.h"

namespace IntaRNA {

/**
 * Memory efficient ensemble predictor for RNAup-like computation, i.e. full
 * DP-implementation without seed-heuristic, using 2D matrices
 *
 * @author Martin Raden
 * @author Frank Gelhausen
 *
 */
class PredictorMfeEns2d: public PredictorMfeEns {

protected:

	//! matrix type to hold the partition functions for interaction site starts
	typedef Matrix<Z_type> Z2dMatrix;

public:

	/**
	 * Constructs a predictor and stores the energy and output handler
	 *
	 * @param energy the interaction energy handler
	 * @param output the output handler to report mfe interactions to
	 * @param predTracker the prediction tracker to be used or NULL if no
	 *         tracking is to be done; if non-NULL, the tracker gets deleted
	 *         on this->destruction.
	 * @param pairProbabilities optional non-owning sink; requires needZall and
	 *         must remain pending until all disjoint regions have succeeded.
	 *         Collection reverses this class's fillHybridZ recurrence; subclasses
	 *         changing that recurrence must also supply their own collection.
	 */
	PredictorMfeEns2d( const InteractionEnergy & energy
					, OutputHandler & output
					, PredictionTracker * predTracker
					, BasePairProbabilities * pairProbabilities = nullptr );

	virtual ~PredictorMfeEns2d();

	/**
	 * Computes the mfe for the given sequence ranges (i1-j1) in the first
	 * sequence and (i2-j2) in the second sequence and reports it to the output
	 * handler.
	 *
	 * @param r1 the index range of the first sequence interacting with r2
	 * @param r2 the index range of the second sequence interacting with r1
	 *
	 */
	virtual
	void
	predict( const IndexRange & r1 = IndexRange(0,RnaSequence::lastPos)
			, const IndexRange & r2 = IndexRange(0,RnaSequence::lastPos)
			);

protected:

	//! optional sequence-pair owner, never deleted or finalized by this predictor
	BasePairProbabilities * pairProbabilities;

	//! access to the interaction energy handler of the super class
	using PredictorMfeEns::energy;

	//! access to the output handler of the super class
	using PredictorMfeEns::output;

	//! energy of all interaction hybrids that end in position p (seq1) and
	//! q (seq2)
	Z2dMatrix hybridZ;

protected:

	/** Compute one region, allowing predict() to fail the sink on any exception. */
	void predictRegion(const IndexRange & r1,const IndexRange & r2);

	/** Reverse the current fixed-right hybrid recurrence. Adjoint seeds obey the
	 * same site filters and numerical-zero tests as updateZ(). Implicit noLP
	 * stacking partners are counted in addition to the explicit DP-state pairs.
	 * @param j1 right boundary in local target coordinates
	 * @param j2 right boundary in local reversed-query coordinates
	 * @param outside scratch adjoints, same dimensions as hybridZ
	 * @param masses accumulated numerators in original regional orientation
	 * @param denominator accumulated boundary contributions, checked against Zall
	 */
	void accumulatePairMasses(size_t j1,size_t j2,Z2dMatrix & outside,
			Z2dMatrix & masses,Z_type & denominator) const;

	/**
	 * Computes all entries of the hybridE matrix for interactions ending in
	 * p=j1 and q=j2 and reports every complete boundary via updateCompleteZ(),
	 * which invokes the virtual updateZ() hook.
	 *
	 * @param j1 end of the interaction within seq 1
	 * @param j2 end of the interaction within seq 2
	 * @param i1init smallest value for i1
	 * @param i2init smallest value for i2
	 * @param callUpdateZ whether or not updateCompleteZ() is to be called
	 *
	 */
	virtual
	void
	fillHybridZ( const size_t j1, const size_t j2
				, const size_t i1init, const size_t i2init
				, const bool callUpdateZ
				);

};

} // namespace

#endif /* INTARNA_PREDICTORMFEENS2D_H_ */
