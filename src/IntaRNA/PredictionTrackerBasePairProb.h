
#ifndef INTARNA_PREDICTIONTRACKERBASEPAIRPROB_H_
#define INTARNA_PREDICTIONTRACKERBASEPAIRPROB_H_

#include "IntaRNA/PredictionTracker.h"
#include "IntaRNA/InteractionEnergy.h"
#include "IntaRNA/Interaction.h"
#include "IntaRNA/PredictorMfeEns.h"
#include "IntaRNA/PredictorMfeEns2dSeedExtension.h"

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

	typedef std::unordered_map<Interaction::Boundary, Z_type, Interaction::Boundary::Hash, Interaction::Boundary::Equal> Site2Z_hash;
	typedef std::unordered_map<Interaction::BasePair, Z_type, Interaction::BasePair::Hash, Interaction::BasePair::Equal> BasePair2Prob_hash;
	typedef std::unordered_map<Interaction::BasePair, std::set<Interaction::BasePair>, Interaction::BasePair::Hash, Interaction::BasePair::Equal> BasePairIndex;
	typedef boost::numeric::ublas::matrix<Z_type> Z2dMatrix;

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
	               , const PredictorMfeEns *predictor);

protected:

	//! energy handler used for predictions
	const InteractionEnergy & energy;

	//! filename of the generated dotplot
	const std::string fileName;

	//! threshold used to draw probabilities in dotplot
	const Z_type probabilityThreshold;

	//! map storing structure probabilities
	BasePair2Prob_hash structureProbs;

	//! left side index
	BasePairIndex rightExt;

	//! right side index
	BasePairIndex leftExt;

	//! flag for seed-based predictors
	bool isSeedPredictor;

	//! maximum postscript width/height in ps units
	const size_t maxDotPlotSize;

	//! partition function of all interaction hybrids that start on the right side of the seed including E_init
	Z2dMatrix hybridZ;

	//! map storing the missing partitions of ZH for all considered interaction sites
	Site2Z_hash ZH_partition_missing;

	/**
	 * Access to the given partition function covering
	 * the given boundary
	 * @param Site2Z_hash partition function hash
	 * @param Boundary boundary
	 * @param addZInit whether or not to add Z_init
	 *
	 * @return the partition function at given boundary
	 */
	Z_type
	getZPartitionValue( const Site2Z_hash *Zpartition, const Interaction::Boundary & boundary, const bool addZInit );

	/**
	 * Compute ZH partition function for given region
	 * @param predictor the predictor providing the probability information
	 * @param seedHandler the seedHandler of the predictor
	 * @param i1 region index
	 * @param j1 region index
	 * @param i2 region index
	 * @param j2 region index
	 *
	 * @return the ZH partition function at given region
	 */
	Z_type
	getZHPartition( const PredictorMfeEns2dSeedExtension *predictor, const SeedHandler* seedHandler
	              , const size_t i1, const size_t j1
	              , const size_t i2, const size_t j2 );

	/**
	 * Compute ZR partition function for given region
	 * @param predictor the predictor providing the probability information
	 * @param seedHandler the seedHandler of the predictor
	 * @param i1 region index
	 * @param j1 region index
	 * @param i2 region index
	 * @param j2 region index
	 * @param si1 index of seed bordering the left side of ZR
	 * @param si2 index of seed bordering the left side of ZR
	 *
	 * @return the ZR partition function at given region
	 */
	Z_type
	getZRPartition( const PredictorMfeEns2dSeedExtension *predictor, const SeedHandler* seedHandler
	              , const size_t i1, const size_t j1
	              , const size_t i2, const size_t j2
								, const size_t si1, const size_t si2 );

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
	          , const PredictorMfeEns *predictor);
	/**
	 * Access to the current partition function covering
	 * the interaction at region (i1, j1, i2, j2).
	 * @param boundary the region of interest
	 * @param predictor the predictor providing the probability information
	 *
	 * @return the hybridization partition function at given region
	 */
	Z_type
	getHybridZ(  const Interaction::Boundary & boundary
	           , const PredictorMfeEns *predictor);

	void
	updateProb( const Interaction::BasePair & bp, const Z_type prob ) {
		if (structureProbs.find(bp)==structureProbs.end()) {
			structureProbs[bp] = prob;
		} else {
			structureProbs[bp] += prob;
		}
	}

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
	generateDotPlot( const char *seq1, const char *seq2, const char *fileName
	               , const plist *pl, const char *comment
								 , const Interaction::Boundary interactionBoundary );

	/**
	 * Compute basepair probabilities and store in structureProbs
	 * @param predictor the predictor providing the probability information
	 * @param seedHandler the seedHandler of the predictor
	 */
	void
	computeBasePairProbs( const PredictorMfeEns2dSeedExtension *predictor, const SeedHandler* seedHandler );

	/**
	 * Compute basepair probabilities for no-seed predictions and store in structureProbs
	 * @param predictor the predictor providing the probability information
	 */
	void
	computeBasePairProbsNoSeed( const PredictorMfeEns *predictor );

	/**
	 * Computes hybridZ
	 *
	 * Note: (i1,i2) have to be complementary (right-most base pair of seed)
	 *
	 * @param l1 start of the interaction within seq 1
	 * @param si1 start of anchor seed in seq 1
	 * @param l2 start of the interaction within seq 2
	 * @param si2 start of anchor seed in seq 2
	 * @param seedHandler the seedHandler of the predictor
	 * @return Z of region
	 */
	Z_type
	fillHybridZ( const size_t l1, const size_t si1, const size_t l2, const size_t si2, const SeedHandler* seedHandler );

	//! postscript template for dotplots
	const char* dotplotTemplate =
		"/box { %%size x y box - draws box centered on x,y\n\
		   2 index 0.5 mul sub            %% x -= 0.5\n\
		   exch 2 index 0.5 mul sub exch  %% y -= 0.5\n\
		   3 -1 roll dup rectfill\n\
		} bind def\n\
		/boxgray { %%size x y box - draws box centered on x,y\n\
			 0 0 1 5 index sub sethsbcolor  %% grayscale\n\
		   1 index 0.5 sub  %% x -= 0.5   s x y x'\n\
		   1 index 0.5 sub  %% y -= 0.5   s x y x' y'\n\
		   5 2 roll  %% x' y' s x y\n\
		   pop pop pop\n\
		   1 dup rectfill\n\
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
