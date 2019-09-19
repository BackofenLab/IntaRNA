
#ifndef INTARNA_PREDICTOR_H_
#define INTARNA_PREDICTOR_H_

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
	 *
	 */
	virtual
	void
	predict( const IndexRange & r1 = IndexRange(0,RnaSequence::lastPos)
			, const IndexRange & r2 = IndexRange(0,RnaSequence::lastPos) ) = 0;

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

	/**
	 * Access to the current overall partition function covering
	 * all interactions of the last predict() call.
	 *
	 * @return the overall hybridization partition function
	 */
	Z_type
	getZall() const;

protected:

	//! energy computation handler
	InteractionEnergyIdxOffset energy;

	//! interaction output handler
	OutputHandler & output;

	//! prediction tracker to be used
	PredictionTracker * predTracker;

	//! the overall partition function since initZ() was last called.
	//! its value is updated by updateZ()
	Z_type Zall;


	/**
	 * Initializes the list of best solutions to be filled by updateOptima()
	 */
	virtual
	void
	initOptima() = 0;

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
	 * @param incrementZall whether or not Zall is to be incremented (if needed)
	 */
	virtual
	void
	updateOptima( const size_t i1, const size_t j1
				, const size_t i2, const size_t j2
				, const E_type energy
				, const bool isHybridE
				, const bool incrementZall ) = 0;

	/**
	 * Updates the global Zall partition function.
	 *
	 * Note: should not be called for interactions for which
	 * updateOptima(..,incrementZall=true) is called. Otherwise, these
	 * interactions are considered multiple times.
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param j1 the index of the first sequence interacting with j2
	 * @param i2 the index of the second sequence interacting with i1
	 * @param j2 the index of the second sequence interacting with j1
	 * @param energy the energy of the interaction site
	 * @param isHybridE whether or not the given energy is only the
	 *        hybridization energy (init+loops) or the total interaction energy
	 * @param incrementZall whether or not Zall is to be incremented (if needed)
	 */
	virtual
	void
	updateZall( const size_t i1, const size_t j1
				, const size_t i2, const size_t j2
				, const E_type energy
				, const bool isHybridE );



	/**
	 * Pushes the optimal and suboptimal solutions to the output handler.
	 */
	virtual
	void
	reportOptima() = 0;

	/**
	 * Increments the overall partition function Zall with the given value
	 * @param partZ the increment for Zall
	 */
	void
	incrementZall( const Z_type partZ );

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
	, Zall(0)
{
}

////////////////////////////////////////////////////////////////////////////

inline
Predictor::~Predictor()
{
	 INTARNA_CLEANUP(predTracker);
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

inline
Z_type
Predictor::
getZall() const
{
	return Zall;
}

////////////////////////////////////////////////////////////////////////////

inline
void
Predictor::
incrementZall( const Z_type partZ )
{
#if INTARNA_IN_DEBUG_MODE
	if ( (std::numeric_limits<Z_type>::max() - partZ) <= Zall) {
		LOG(WARNING) <<"PredictorMfeEns::incrementZall() : partition function overflow! Recompile with larger partition function data type!";
	}
#endif
	// increment overall partition function
	Zall += partZ;
}

////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////

inline
void
Predictor::
updateZall( const size_t i1, const size_t j1
		, const size_t i2, const size_t j2
		, const E_type interE
		, const bool isHybridE )
{
	// ignore if not needed
	if (!output.getOutputConstraint().needZall) {
		return;
	}

	// ignore invalid reports
	if (E_isINF(interE) || interE >= E_MAX) {
		return;
	}

	// check GU ends if needed
	if (output.getOutputConstraint().noGUend && (energy.isGU(i1,i2) || energy.isGU(j1,j2)) ) {
		return;
	}

	// increment Zall with BW of overall energy
	incrementZall(
			energy.getBoltzmannWeight(
					isHybridE ?
							energy.getE( i1,j1, i2,j2, interE )
							: interE
			));
}

////////////////////////////////////////////////////////////////////////////

} // namespace

#endif /* PREDICTOR_H_ */
