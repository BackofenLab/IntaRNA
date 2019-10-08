
#ifndef INTARNA_PREDICTIONTRACKERBASEPAIRPROB_H_
#define INTARNA_PREDICTIONTRACKERBASEPAIRPROB_H_

#include "IntaRNA/PredictionTracker.h"
#include "IntaRNA/InteractionEnergy.h"
#include "IntaRNA/Interaction.h"
#include "IntaRNA/PredictorMfeEns.h"

#include <iostream>

#include <boost/algorithm/string.hpp>

namespace IntaRNA {

/**
 * Collects partition function parts Z(i,j), computes base-pair probabilities
 * and prints them as a dot plot
 *
 * The information is written to stream on destruction.
 *
 * @author Frank Gelhausen
 *
 */
class PredictionTrackerBasePairProb: public PredictionTracker
{

public:

	typedef std::unordered_map<Interaction::BasePair, Z_type, Interaction::BasePairHash> BasePair2Prob_hash;

	typedef std::unordered_map<Interaction::BasePair, std::list<Interaction::BasePair>, Interaction::BasePairHash> BasePairIndex;

public:

	/**
	 * Constructs a PredictionTracker that collects probability information
	 * for an interaction by computing the Boltzmann probabilities and
	 * generates a basepair-probabilities dotplot
	 *
	 * @param energy the energy function used for energy calculation
	 * @param fileName the name of the generated postscript file containing
	 * the dotplot
	 */
	PredictionTrackerBasePairProb(
				const InteractionEnergy & energy
				, const std::string & fileName
			);

	/**
	 * destruction: write the probabilities to stream.
	 */
	virtual ~PredictionTrackerBasePairProb();


	/**
	 * Updates the probability information for each Predictor.updateOptima() call.
	 *
	 * @param i1 the index of the first sequence interacting with i2
	 * @param j1 the index of the first sequence interacting with j2
	 * @param i2 the index of the second sequence interacting with i1
	 * @param j2 the index of the second sequence interacting with j1
	 * @param energy the overall energy of the interaction site
	 */
	virtual
	void
	updateOptimumCalled( const size_t i1, const size_t j1
						, const size_t i2, const size_t j2
						, const E_type energy
						);

	/**
	 * Updates the probability information.
	 *
	 * @param predictor the predictor providing the probability information
	 * @param seedHandler the seedHandler of the predictor (NULL if noseed predictior)
	 */
	virtual
	void
	updateZ( PredictorMfeEns *predictor, SeedHandler* seedHandler ) override;

	/**
	 * Compute partition function of the
	 * the subregion (i1, j1, i2, j2)
	 * @param i1 region index
	 * @param j1 region index
	 * @param i2 region index
	 * @param j2 region index
	 * @param predictor the predictor providing the probability information
	 * @param seedHandler the seedHandler of the predictor
	 */
	void
	computeMissingZ( const size_t i1, const size_t j1
								, const size_t i2, const size_t j2
								, const size_t l1, const size_t r1
								, const size_t l2, const size_t r2
								, PredictorMfeEns *predictor
								, SeedHandler* seedHandler );

	/**
	 * Check if full seed exists in the region (i1, j1, i2, j2)
	 * @param i1 region index
	 * @param j1 region index
	 * @param i2 region index
	 * @param j2 region index
	 * @param seedHandler the seedHandler of the predictor
	 *
	 * @return true if region contains a full seed
	 */
	bool
	isFullSeedinRegion( const size_t i1, const size_t j1
					        	, const size_t i2, const size_t j2
					        	, SeedHandler* seedHandler );

	/**
	 * Get leftmost seeds containing basepair (k1, k2)
	 * @param k1 base pair index
	 * @param k2 base pair index
	 * @param seedHandler the seedHandler of the predictor
	 *
	 * @return vector of pairs containing the left-most seed of each
	 *         loop-overlapping cluster of seeds containing base pair k
	 */
	std::vector< std::pair <size_t, size_t> >
	getLeftMostSeedsAtK( const size_t k1, const size_t k2
					           , SeedHandler* seedHandler );

	/**
	 * Get partial energy of seed (si1, si2) at region (i1, j1, i2, j2)
	 * @param i1 region index
	 * @param j1 region index
	 * @param i2 region index
	 * @param j2 region index
	 * @param seedHandler the seedHandler of the predictor
	 *
	 * @return energy of seed at given region
	 */
	E_type
	getPartialSeedEnergy( const size_t si1, const size_t si2
		                  , const size_t i1, const size_t j1
											, const size_t i2, const size_t j2
											, SeedHandler* seedHandler );
	/**
	 * Access to the current partition function covering
	 * the interaction at region (i1, j1, i2, j2).
	 * @param i1 region index
	 * @param j1 region index
	 * @param i2 region index
	 * @param j2 region index
	 * @param predictor the predictor providing the probability information
	 *
	 * @return the hybridization partition function at given region
	 */
	Z_type
	getHybridZ( const size_t i1, const size_t j1
	          , const size_t i2, const size_t j2
	          , PredictorMfeEns *predictor);

	/**
	 * Access to the base pair probability
	 * of basepair (i1, i2)
	 * @param i1 index in first sequence
	 * @param i2 index in second sequence
	 * @param predictor the predictor providing the probability information
	 *
	 * @return the base pair probability of given basepair
	 */
	Z_type
	getBasePairProb( const size_t i1, const size_t i2
	               , PredictorMfeEns *predictor);

	/**
	 * Set the current partition function covering
	 * the interaction at region (i1, j1, i2, j2).
	 * @param i1 region index
	 * @param j1 region index
	 * @param i2 region index
	 * @param j2 region index
	 * @param partZ the new partition function for givent region
	 */
	void
	updateHybridZ( const size_t i1, const size_t j1
						 	 , const size_t i2, const size_t j2
							 , const Z_type partZ );

  /**
	 * Generates a dotplot of the given base pair probabilities
	 * @param seq1 first RNA sequence
	 * @param seq2 second RNA sequence
	 * @param fileName name of the output file
	 * @param pl plist containing base pair probabilities
	 * @param comment comment to include in postscript
	 * @param interactionBoundary boundary of the predicted interaction
	 *
	 * @return false in case of failure
	 */
	bool
	generateDotPlot( char *seq1, char *seq2, char *fileName
	               , plist *pl, const char *comment
								 , Interaction::Boundary interactionBoundary );

	/**
	 * Compute basepair probabilities and store in structureProbs
	 * @param predictor the predictor providing the probability information
	 * @param iterator start of partition function
	 * @param iterator end of partition function
	 */
	void
	computeBasePairProbs( PredictorMfeEns *predictor
		                  , const PredictorMfeEns::Site2Z_hash::const_iterator first
	                    , const PredictorMfeEns::Site2Z_hash::const_iterator last );

protected:

	//! energy handler used for predictions
	const InteractionEnergy & energy;

	//! filename of the generated dotplot
	const std::string fileName;

	//! threshold used to draw probabilities in dotplot
	const Z_type probabilityThreshold;

	//! map storing structure probabilities
	BasePair2Prob_hash structureProbs;

	//! map storing missing Z partitions for a given interaction
	PredictorMfeEns::Site2Z_hash Z_partitionMissing;

	//! left side index
	BasePairIndex leftIndex;

	//! right side index
	BasePairIndex rightIndex;

	//! postscript template for dotplots
	const char* dotplotTemplate =
		"/box { %%size x y box - draws box centered on x,y\n\
		   2 index 0.5 mul sub            %% x -= 0.5\n\
		   exch 2 index 0.5 mul sub exch  %% y -= 0.5\n\
		   3 -1 roll dup rectfill\n\
		} bind def\n\
		\n\
		/drawseq1 { %% print sequence1\n\
		[ [0.7 -0.3 ]\n\
		  [0.7 0.7 len2 add]\n\
		] {\n\
		   gsave\n\
		    aload pop translate\n\
		    0 1 len 1 sub {\n\
		     dup 0 moveto\n\
		     sequence1 exch 1 getinterval\n\
		     show\n\
		    } for\n\
		   grestore\n\
		  } forall\n\
		} bind def\n\
    \n\
		/drawseq2 { %% print sequence2\n\
		[ [-0.3 len2 sub -0.4 -90]\n\
		  [-0.3 len2 sub 0.7 len add -90]\n\
		] {\n\
		   gsave\n\
		    aload pop rotate translate\n\
		    0 1 len2 1 sub {\n\
		     dup 0 moveto\n\
		     sequence2 exch 1 getinterval\n\
		     show\n\
		    } for\n\
		   grestore\n\
		  } forall\n\
		} bind def\n\
    \n\
		/rectangle {%% x y w h RT -\n\
			   %% draw a rectangle size w h at x y\n\
				 4 -2 roll moveto %% lower left corner\n\
				 dup 0 exch rlineto %% to upper left\n\
				 exch 0 rlineto %% to upper right\n\
				 neg 0 exch rlineto %% to lower right\n\
				 closepath\n\
		} def\n\
		/drawgrid{\n\
		  gsave\n\
		  0.5 dup translate\n\
		  1 %% len log 0.9 sub cvi 10 exch exp %% grid spacing\n\
		  0 exch len {\n\
		     dup\n\
				 dup cvi 10 mod 0 eq {\n\
				   0.01 setlinewidth\n\
					 [1 0] 0 setdash\n\
				 } {\n\
				   0.01 setlinewidth\n\
					 [0.3 0.7] 0.15 setdash\n\
				 } ifelse\n\
		     0 moveto\n\
		     len2 lineto %% vertical\n\
		     stroke\n\
		  } for\n\
      \n\
		  1 %% len log 0.9 sub cvi 10 exch exp %% grid spacing\n\
		  0 exch len2 {\n\
		     dup\n\
				 dup cvi 10 mod 0 eq {\n\
				   0.01 setlinewidth\n\
					 [1 0] 0 setdash\n\
				 } {\n\
				   0.01 setlinewidth\n\
					 [0.3 0.7] 0.15 setdash\n\
				 } ifelse\n\
		     len2 exch sub 0 exch moveto\n\
		     len exch len2 exch sub lineto %% horizontal\n\
		     stroke\n\
		  } for\n\
      \n\
		  grestore\n\
		} bind def\n";

};

//////////////////////////////////////////////////////////////////////

} // namespace

#endif /* INTARNA_PREDICTIONTRACKERBASEPAIRPROB_H_ */
