
#ifndef INTARNA_PREDICTORMFEENS2DHEURISTIC_H_
#define INTARNA_PREDICTORMFEENS2DHEURISTIC_H_

#include "IntaRNA/PredictorMfeEns2d.h"
#include "IntaRNA/Interaction.h"

#include <boost/numeric/ublas/matrix.hpp>

namespace IntaRNA {

/**
 * Memory efficient interaction predictor that uses a heuristic to
 * find the mfe or a close-to-mfe interaction.
 *
 * To this end, for each interaction start i1,i2 only the optimal right side
 * interaction with boundaries j1,j2 is considered in the recursion instead of
 * all possible interaction ranges.
 *
 * This yields a quadratic time and space complexity.
 *
 * @author Martin Mann
 *
 */
class PredictorMfeEns2dHeuristic: public PredictorMfeEns2d {

protected:

	/**
	 * Describes the currently best interaction found for a left interaction
	 * boundary i1,i2
	 */
	class BestInteraction {
	public:

		//! init data
		BestInteraction( const Z_type Z=Z_INF, const size_t j1=RnaSequence::lastPos, const size_t j2=RnaSequence::lastPos )
			: Z(Z), j1(j1), j2(j2)
		{}

	public:
		//! energy of the interaction
		Z_type Z;
		//! right end of the interaction in seq1
		size_t j1;
		//! right end of the interaction in seq2
		size_t j2;
	};

	//! matrix type to hold the mfe energies and boundaries for interaction site starts
	typedef boost::numeric::ublas::matrix<BestInteraction> Z2dMatrix;

public:

	/**
	 * Constructs a predictor and stores the energy and output handler
	 *
	 * @param energy the interaction energy handler
	 * @param output the output handler to report mfe interactions to
	 * @param predTracker the prediction tracker to be used or NULL if no
	 *         tracking is to be done; if non-NULL, the tracker gets deleted
	 *         on this->destruction.
	 */
	PredictorMfeEns2dHeuristic( const InteractionEnergy & energy
							, OutputHandler & output
							, PredictionTracker * predTracker );

	virtual ~PredictorMfeEns2dHeuristic();

	/**
	 * Computes the mfe for the given sequence ranges (i1-j1) in the first
	 * sequence and (i2-j2) in the second sequence and reports it to the output
	 * handler.
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
			, const OutputConstraint & outConstraint = OutputConstraint() );

protected:

	//! access to the interaction energy handler of the super class
	using PredictorMfeEns2d::energy;

	//! access to the output handler of the super class
	using PredictorMfeEns2d::output;

	//! access to the list of reported interaction ranges of the super class
	using PredictorMfeEns2d::reportedInteractions;

	//! energy of all interaction hybrids starting in i1,i2
	Z2dMatrix hybridZ;

protected:

	/**
	 * Computes all entries of the hybridE matrix
	 * and reports all valid interactions via updateOptima()
	 * @param outConstraint constrains the interactions reported to the output handler
	 */
	virtual
	void
	fillHybridZ( const OutputConstraint & outConstraint );

	/**
	 * Identifies the next best interaction (containing a seed)
	 * with an energy equal to or higher
	 * than the given interaction. The new interaction will not overlap any
	 * index range stored in reportedInteractions.
	 *
	 * @param curBest IN/OUT the current best interaction to be replaced with one
	 *        of equal or higher energy not overlapping with any reported
	 *        interaction so far; an interaction with energy E_INF is set, if
	 *        there is no better interaction left
	 */
	virtual
	void
	getNextBest( Interaction & curBest );

};

} // namespace

#endif /* INTARNA_PREDICTORMFEENS2DHEURISTIC_H_ */
