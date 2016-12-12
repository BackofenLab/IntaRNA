
#ifndef PREDICTOR_H_
#define PREDICTOR_H_

#include "general.h"
#include "InteractionEnergyIdxOffset.h"

#include "OutputHandler.h"

/**
 * RNA-RNA interaction prediction handler interface.
 *
 * @author Martin Mann 2014
 *
 */
class Predictor {

public:

	/**
	 * Constructs an RNA-RNA interaction prediction handler and sets the
	 * central data members.
	 *
	 * @param energy the handler for interaction energy computation
	 * @param output the output handler for identified interactions
	 *
	 */
	Predictor( const InteractionEnergy & energy, OutputHandler & output );

	/**
	 * destruction
	 */
	virtual ~Predictor();

	/**
	 * Computes the predictors optimization target for the given sequence
	 * ranges in the first sequence and second sequence.
	 * The according optimal interaction is given to the output handler.
	 *
	 * @param r1 the index range of the first sequence interacting with r2
	 * @param r2 the index range of the second sequence interacting with r1
	 * @param reportMax the maximal number of (sub)optimal interactions to be
	 *            reported to the output handler
	 * @param reportNonOverlapping whether or not the reported interactions
	 *            should be non-overlapping or not
	 *
	 */
	virtual
	void
	predict( const IndexRange & r1 = IndexRange(0,RnaSequence::lastPos)
			, const IndexRange & r2 = IndexRange(0,RnaSequence::lastPos)
			, const size_t reportMax = 1
			, const bool reportNonOverlapping = true ) = 0;

	/**
	 * Computes the maximal width of an interaction for a given site width and
	 * maximal size of interaction loops.
	 *
	 * @param w the width to compute the maximal interacting width for
	 * @param maxLoopSize the maximal size of loops within interactions
	 *
	 * @return 1 + (w-1)*(maxLoopSize+1): if w>0; 0 otherwise
	 */
	static
	size_t
	getMaxInteractionWidth( const size_t w, const size_t maxLoopSize );

protected:

	//! energy computation handler
	InteractionEnergyIdxOffset energy;

	//! interaction output handler
	OutputHandler & output;

};

#endif /* PREDICTOR_H_ */
