
#ifndef PREDICTORMAXPROB_H_
#define PREDICTORMAXPROB_H_

#include "Predictor.h"
#include "Interaction.h"

#include <boost/numeric/ublas/matrix.hpp>

/**
 * Computes the interaction site with maximal probability among all interaction
 * sites
 *
 * @author Martin Mann
 *
 */
class PredictorMaxProb: public Predictor {

protected:

	//! matrix type to cover the energies for different interaction site widths
	typedef boost::numeric::ublas::matrix<E_type> E2dMatrix;

	//! full 4D DP-matrix for computation to hold all start position combinations
	//! first index = start positions (i1,i2) of (seq1,seq2)
	//! second index = interaction window sizes (w1,w2) or NULL if (i1,i2) not complementary
	typedef boost::numeric::ublas::matrix< E2dMatrix* > E4dMatrix;


public:

	/**
	 * Constructs a predictor and stores the energy and output handler
	 *
	 * @param energy the interaction energy handler
	 * @param output the output handler to report optimal interactions to
	 */
	PredictorMaxProb( const InteractionEnergy & energy, OutputHandler & output );

	virtual ~PredictorMaxProb();

	/**
	 * Computes the interaction site with maximal probability
	 * for the given sequence ranges (i1-j1) in the first
	 * sequence and (i2-j2) in the second sequence and reports it to the output
	 * handler.
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param j1 the index of the first sequence interacting with j2
	 * @param i2 the index of the second sequence interacting with i1
	 * @param j2 the index of the second sequence interacting with j1
	 */
	void
	predict( const size_t i1 = 0, const size_t j1 = RnaSequence::lastPos
				, const size_t i2 = 0, const size_t j2 = RnaSequence::lastPos);

protected:

	//! access to the interaction energy handler of the super class
	using Predictor::energy;

	//! access to the output handler of the super class
	using Predictor::output;

	//! partition function of all interaction hybrids computed by the recursion with indices
	//! hybridZ(i1,i2)->(w1,w2), with interaction start i1 (seq1) and i2 (seq2) and
	//! ineraction end j1=i1+w1 and j2=j2+w2
	//! NOTE: hybridZ(i1,i2)==NULL if not complementary(seq1[i1],seq2[i2])
	E4dMatrix hybridZ;

	//! the overall partition function = sum or all hybridZ entries
	double Z;

	//! interaction boundaries with maximal probability
	Interaction maxProbInteraction;

	//! offset for indices in sequence 1 for current computation
	size_t i1offset;
	//! offset for indices in sequence 2 for current computation
	size_t i2offset;

protected:

	/**
	 * Removes all temporary data structures and resets the predictor
	 */
	void
	clear();

	/**
	 * computes all entries of the hybridZ matrix
	 */
	void
	fillHybridZ( const InteractionEnergy & energy );

	/**
	 * Initializes the interaction site with maximal probability
	 */
	virtual
	void
	initMaxProbInteraction();

	/**
	 * updates the global optimum if needed
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param j1 the index of the first sequence interacting with j2
	 * @param i2 the index of the second sequence interacting with i1
	 * @param j2 the index of the second sequence interacting with j1
	 * @param curZ the partition function for the interaction site
	 */
	virtual
	void
	updateMaxProbInteraction( const size_t i1, const size_t j1
			, const size_t i2, const size_t j2
			, const E_type curZ );

	/**
	 * Provides the Boltzmann weight for a given energy.
	 * @param energ the energy the Boltzmann weight is to be computed for
	 * @return the Boltzmann weight, i.e. exp( - energy / RT );
	 */
	E_type
	getBoltzmannWeight( const E_type energy );

};

#endif /* PREDICTORMAXPROB_H_ */
