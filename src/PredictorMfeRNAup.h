
#ifndef PREDICTORMFERNAUP_H_
#define PREDICTORMFERNAUP_H_

#include "Predictor.h"
#include "Interaction.h"

#include <boost/numeric/ublas/matrix.hpp>

/**
 * Predictor for RNAup-like computation, i.e. full DP-implementation without
 * seed-heuristic
 *
 * @author Martin Mann
 *
 */
class PredictorMfeRNAup: public Predictor {

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
	 * @param output the output handler to fill with mfe interactions
	 */
	PredictorMfeRNAup( const InteractionEnergy & energy, OutputHandler & output );

	virtual ~PredictorMfeRNAup();

	/**
	 * Computes the mfe for the given sequence ranges (i1-j1) in the first
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

	//! energy of all interaction hybrids computed by the recursion with indices
	//! hybridE(i1,i2)->(w1,w2), with interaction start i1 (seq1) and i2 (seq2) and
	//! ineraction end j1=i1+w1 and j2=j2+w2
	//! NOTE: hybridE(i1,i2)==NULL if not complementary(seq1[i1],seq2[i2])
	E4dMatrix hybridE;

	//! mfe interaction boundaries
	Interaction mfeInteraction;

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
	 * computes all entries of the hybridE matrix
	 */
	void
	fillHybridE( const InteractionEnergy & energy );

	/**
	 * Fills a given interaction (boundaries given) with the according
	 * hybridizing base pairs.
	 * @param interaction IN/OUT the interaction to fill
	 */
	void
	traceBack( Interaction & interaction ) const;

	/**
	 * Initializes the global energy minimum
	 */
	virtual
	void
	initMfe();

	/**
	 * updates the global optimum to be the mfe interaction if needed
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param j1 the index of the first sequence interacting with j2
	 * @param i2 the index of the second sequence interacting with i1
	 * @param j2 the index of the second sequence interacting with j1
	 * @param eH the energy of the hybridization only (just base pairs etc)
	 * @param eE the energy of the hybridization's dangling ends (left+right)
	 * @param eD the energy needed to make it accessible (seq1[i-j]+seq2[i-j])
	 */
	virtual
	void
	updateMfe( const size_t i1, const size_t j1
			, const size_t i2, const size_t j2
			, const E_type eH, const E_type eE, const E_type eD );

};

#endif /* PREDICTORMFERNAUP_H_ */
