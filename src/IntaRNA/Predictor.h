
#ifndef PREDICTOR_H_
#define PREDICTOR_H_

#include "IntaRNA/general.h"
#include "IntaRNA/InteractionEnergyIdxOffset.h"

#include "IntaRNA/OutputConstraint.h"
#include "IntaRNA/OutputHandler.h"

#include "IntaRNA/PredictionTracker.h"

namespace IntaRNA {

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
	 * @param predTracker the prediction tracker to be used or NULL if no
	 *         tracking is to be done; if non-NULL, the tracker gets deleted
	 *         on this->destruction.
	 *
	 */
	Predictor( const InteractionEnergy & energy
				, OutputHandler & output
				, PredictionTracker * predTracker
				);

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
	 * @param outConstraint constrains the interactions reported to the output handler
	 *
	 */
	virtual
	void
	predict( const IndexRange & r1 = IndexRange(0,RnaSequence::lastPos)
			, const IndexRange & r2 = IndexRange(0,RnaSequence::lastPos)
			, const OutputConstraint & outConstraint = OutputConstraint() ) = 0;

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

	//! prediction tracker to be used
	PredictionTracker * predTracker;


	/**
	 * Initializes the list of best solutions to be filled by updateOptima()
	 *
	 * @param outConstraint constrains the interactions reported to the output handler
	 */
	virtual
	void
	initOptima( const OutputConstraint & outConstraint ) = 0;

	/**
	 * Updates the global the list of best solutions found so far
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param j1 the index of the first sequence interacting with j2
	 * @param i2 the index of the second sequence interacting with i1
	 * @param j2 the index of the second sequence interacting with j1
	 * @param energy the energy of the interaction site
	 * @param isHybridE whether or not the given energy is only the
	 *        hybridization energy (init+loops) or the total interaction energy
	 */
	virtual
	void
	updateOptima( const size_t i1, const size_t j1
				, const size_t i2, const size_t j2
				, const E_type energy
				, const bool isHybridE ) = 0;




	/**
	 * Pushes the optimal and suboptimal solutions to the output handler.
	 *
	 * @param outConstraint constrains the interactions reported to the output handler
	 */
	virtual
	void
	reportOptima( const OutputConstraint & outConstraint ) = 0;


};


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

inline
Predictor::Predictor( const InteractionEnergy & energy
					, OutputHandler & output
					, PredictionTracker * predTracker )
:
	energy(energy)
	, output(output)
	, predTracker(predTracker)
{
}

////////////////////////////////////////////////////////////////////////////

inline
Predictor::~Predictor()
{
	CLEANUP(predTracker);
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
Predictor::
getMaxInteractionWidth( const size_t w, const size_t maxLoopSize )
{
	// check if no window at all
	if (w==0)
		return 0;
	// compute maximal interaction width, ie each position pairs with max loop size
	return 1 + ((w-1)*(maxLoopSize+1));
}

////////////////////////////////////////////////////////////////////////////

} // namespace

#endif /* PREDICTOR_H_ */
