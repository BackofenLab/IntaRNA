/*
 * PredictorMfe2dSeed.h
 *
 *  Created on: 13.10.2016
 *      Author: Mmann
 */

#ifndef PREDICTORMFE2DSEED_H_
#define PREDICTORMFE2DSEED_H_

#include "PredictorMfe2d.h"
#include "SeedConstraint.h"

//#define BOOST_DISABLE_ASSERTS // define to disable dimension checks
#include <boost/multi_array.hpp>

/**
 * Implements seed-based space-efficient interaction prediction.
 *
 * Note, for each seed start (i1,i2) only the mfe seed is considered for the
 * overall interaction computation instead of considering all possible seeds
 * starting at (i1,i2).
 *
 * @author Martin Mann
 *
 */
class PredictorMfe2dSeed: public PredictorMfe2d {

protected:

	//! matrix type to hold the mfe energies for interaction site starts
	typedef PredictorMfe2d::E2dMatrix E2dMatrix;

	//! 5D matrix type to hold the mfe energies for seed interactions
	//! of the ranges i1..(i1+bp+u1-1) with i2..(i2+bp+u2-1), with
	//! i1,i2 = the start index of the seed in seq1/2
	//! bp = the number of base pairs within the seed
	//! bpInbetween = the number of base pairs enclosed by left and right base pair, ie. == (bp-2)
	//! u1/u2 = the number of unpaired positions within the seed,
	//! using the index [i1][i2][bpInbetween][u1][u2] or a SeedIndex object
	typedef boost::multi_array<E_type,5> SeedRecMatrix;

	//! defines the seed data {{ i1, i2, bpInbetween, u1, u2 }} to access elements of
	//! the SeedRecMatrix
	typedef boost::array<SeedRecMatrix::index, 5> SeedIndex;

	//! matrix to store the seed information for each seed left side (i1,i2);
	//! it holds both the energy (first) as well as the length of the seed using
	//! the length combination using encodeSeedLength()
	typedef boost::numeric::ublas::matrix< std::pair<E_type, size_t> > SeedMatrix;

public:


	/**
	 * Constructs a predictor and stores the energy and output handler
	 *
	 * @param energy the interaction energy handler
	 * @param output the output handler to report mfe interactions to
	 * @param seedConstraint the seed constraints to be applied
	 */
	PredictorMfe2dSeed(
			const InteractionEnergy & energy
			, OutputHandler & output
			, const SeedConstraint & seedConstraint );


	/**
	 * data cleanup
	 */
	virtual ~PredictorMfe2dSeed();


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

	//! access to the mfe interaction of the super class
	using PredictorMfe2d::mfeInteraction;

	//! access to the index offset in seq1 of the super class
	using PredictorMfe2d::i1offset;

	//! access to the index offset in seq2 of the super class
	using PredictorMfe2d::i2offset;

	// TODO provide all data structures as arguments to make predict() call threadsafe

	//! energy of all interaction hybrids that end in position p (seq1) and
	//! q (seq2)
	using PredictorMfe2d::hybridE_pq;

	//! the current range of computed entries within hybridE_pq set by initHybridE()
	using PredictorMfe2d::hybridErange;

	//! the seed constraints
	const SeedConstraint & seedConstraint;

	//! for fixed interaction end p=j1,q=j2: each cell (i1,i2) provides the mfe
	//! for the interaction i1..j1 with i2..j2 given that the range contains
	//! a valid seed interaction
	E2dMatrix hybridE_pq_seed;

	//! the recursion data for the computation of a seed interaction
	//! i1..(i1+bpInbetween+u1-1) with i2..(i2+bpInbetween+u2-1)
	//! using the indexing [i1][i2][bpInbetween][u1][u2]
	SeedRecMatrix seedE_rec;

	//! the seed mfe information for seeds starting at (i1,i2)
	//! TODO replace with sparse data structure
	SeedMatrix seed;

protected:

	/**
	 * does nothing but to ignore the calls from fillHybridE()
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param j1 the index of the first sequence interacting with j2
	 * @param i2 the index of the second sequence interacting with i1
	 * @param j2 the index of the second sequence interacting with j1
	 * @param hybridE ignored
	 */
	virtual
	void
	updateMfe( const size_t i1, const size_t j1
			, const size_t i2, const size_t j2
			, const E_type energy );

	/**
	 * Computes the seed matrix for the given interval boundaries
	 * @param i1 the index of the first sequence interacting with i2
	 * @param j1 the index of the first sequence interacting with j2
	 * @param i2 the index of the second sequence interacting with i1
	 * @param j2 the index of the second sequence interacting with j1
	 */
	void
	fillSeed(const size_t i1, const size_t j1, const size_t i2, const size_t j2);

	/**
	 * Computes all entries of the hybridE_seed matrix for interactions ending in
	 * p=j1 and q=j2 and report all valid interactions to updateMfe()
	 *
	 * @param energy the energy function to use
	 * @param j1 end of the interaction within seq 1
	 * @param j2 end of the interaction within seq 2
	 * @param i1min smallest value for i1
	 * @param i2min smallest value for i2
	 *
	 */
	void
	fillHybridE_seed( const size_t j1, const size_t j2, const size_t i1min=0, const size_t i2min=0  );

	/**
	 * Fills a given interaction (boundaries given) with the according
	 * hybridizing base pairs using hybridE_seed.
	 * @param interaction IN/OUT the interaction to fill
	 */
	void
	traceBack( Interaction & interaction );

	/**
	 * Fills a given interaction with the according
	 * hybridizing base pairs of the provided seed interaction (excluding
	 * the right-most seed base pair)
	 *
	 * @param interaction IN/OUT the interaction to fill
	 * @param i1 the seed left end in seq 1
	 * @param i2 the seed left end in seq 2
	 * @param bp the number of base pairs (bp+2) within the seed so far
	 * @param u1 the number of unpaired bases within seq 1
	 * @param u2 the number of unpaired bases within seq 2
	 */
	void
	traceBackSeed( Interaction & interaction
			, const size_t i1, const size_t i2, const size_t bp
			, const size_t u1, const size_t u2 );

	/**
	 * Provides the seed energy during recursion.
	 *
	 * NOTE: internally a ring-list data structure is used which reuses memory
	 * instead of allocating mem for all possible parameter combinations. Thus,
	 * you have to call the method in appropriate order depending on your seed
	 * recursion.
	 *
	 * @param i1 the seed left end in seq 1
	 * @param i2 the seed left end in seq 2
	 * @param bpInbetween the number of seed base pairs enclosed by the left-
	 *        and right-most base pair, ie. bpSeed-2
	 * @param u1 the number of unpaired bases within seq 1
	 * @param u2 the number of unpaired bases within seq 2
	 *
	 * @return the energy of the according (sub)seed
	 */
	E_type
	getSeedE( const size_t i1, const size_t i2, const size_t bpInbetween, const size_t u1, const size_t u2 );

	/**
	 * Fills the seed energy during recursion.
	 *
	 * NOTE: internally a ring-list data structure is used which reuses memory
	 * instead of allocating mem for all possible parameter combinations. Thus,
	 * you have to call the method in appropriate order depending on your seed
	 * recursion.
	 *
	 * @param i1 the seed left end in seq 1
	 * @param i2 the seed left end in seq 2
	 * @param bpInbetween the number of seed base pairs enclosed by the left-
	 *        and right-most base pair, ie. bpSeed-2
	 * @param u1 the number of unpaired bases within seq 1
	 * @param u2 the number of unpaired bases within seq 2
	 * @param E the energy value to be set
	 */
	void
	setSeedE( const size_t i1, const size_t i2, const size_t bpInbetween, const size_t u1, const size_t u2, const E_type E );

	/**
	 * Encodes the seed lengths into one number
	 * @param l1 the length of the seed in seq1
	 * @param l2 the length of the seed in seq2
	 * @return the combined encoding = (l1 + l2*(max_l1+1))
	 */
	size_t
	encodeSeedLength( const size_t l1, const size_t l2 ) const;

	/**
	 * Decodes the length of the seed within sequence 1 from an encoding
	 * generated with encodeSeedLength()
	 * @param code the lengths encoding
	 * @return the length of the seed in seq1
	 */
	size_t
	decodeSeedLength1( const size_t code ) const;

	/**
	 * Decodes the length of the seed within sequence 2 from an encoding
	 * generated with encodeSeedLength()
	 * @param code the lengths encoding
	 * @return the length of the seed in seq2
	 */
	size_t
	decodeSeedLength2( const size_t code ) const;

};

#endif /* PREDICTORMFE2DSEED_H_ */
