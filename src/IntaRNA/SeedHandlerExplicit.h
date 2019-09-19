
#ifndef INTARNA_SEEDHANDLEREXPLICIT_H_
#define INTARNA_SEEDHANDLEREXPLICIT_H_


#include "IntaRNA/SeedHandler.h"

#include <boost/unordered_map.hpp>

namespace IntaRNA
{

/**
 * Enables the explicit definition of possible seed interactions.
 *
 * Only the provided seeds will be considered for interaction prediction.
 *
 */
class SeedHandlerExplicit : public SeedHandler
{
public:

	/**
	 * Construction
	 * @param energy the energy function to be used for seed prediction
	 * @param seedConstraint the seed constraint to be applied that contains the
	 *        explicit seed encodings
	 */
	SeedHandlerExplicit(
			const InteractionEnergy & energy
			, const SeedConstraint & seedConstraint
			);

	/**
	 * destruction
	 */
	virtual ~SeedHandlerExplicit();


	/**
	 * Computes the seed information for the given interval boundaries
	 * @param i1 the first index of seq1 that might interact
	 * @param j1 the last index of seq1 that might interact
	 * @param i2 the first index of seq2 that might interact
	 * @param j2 the last index of seq2 that might interact
	 * @return the number of potential seed interactions
	 */
	virtual
	size_t
	fillSeed(const size_t i1, const size_t j1, const size_t i2, const size_t j2);

	/**
	 * Identifies the base pairs of the seed interaction starting at i1,i2
	 * and writes them to the provided container
	 *
	 * NOTE: the left- and right-most base pair is excluded!
	 *
	 * @param interaction the container to add the base pairs too
	 * @param i1 the start of the seed in seq1
	 * @param i2 the start of the seed in seq2
	 */
	virtual
	void
	traceBackSeed( Interaction & interaction, const size_t i1, const size_t i2) const;


	/**
	 * Access to the energy of a seed with left-most base pair (i1,i2)
	 * @param i1 the left most interacting base of seq1
	 * @param i2 the left most interacting base of seq2
	 * @return the energy of a seed starting at (i1,i2) or E_INF if none possible
	 */
	virtual
	E_type
	getSeedE( const size_t i1, const size_t i2 ) const;

	/**
	 * Checks whether or not a given base pair is the left-most base pair of
	 * any seed
	 * @param i1 the interacting base of seq1
	 * @param i2 the interacting base of seq2
	 * @return true if (i1,i2) is the left most base pair of some seed; false
	 *         otherwise
	 */
	virtual
	bool
	isSeedBound( const size_t i1, const size_t i2 ) const;

	/**
	 * Access to the length in seq1 of the seed with left-most base pair (i1,i2)
	 * @param i1 the left most interacting base of seq1
	 * @param i2 the left most interacting base of seq2
	 * @return the length in seq1 of the mfe seed starting at (i1,i2) or 0 if none possible
	 */
	virtual
	size_t
	getSeedLength1( const size_t i1, const size_t i2 ) const;

	/**
	 * Access to the length in seq2 of the seed with left-most base pair (i1,i2)
	 * @param i1 the left most interacting base of seq1
	 * @param i2 the left most interacting base of seq2
	 * @return the length in seq2 of the seed starting at (i1,i2) or 0 if none possible
	 */
	virtual
	size_t
	getSeedLength2( const size_t i1, const size_t i2 ) const;



	/**
	 * checks whether or not a given explicit seed encoding is valid or not
	 *
	 * @param seed the seed encoding to check
	 * @return an empty string for valid encodings; or an error description
	 *         otherwise
	 */
	static
	std::string
	checkSeedEncoding( const std::string & seed );

	/**
	 * parses the given explicit seed encoding for the maximal number of
	 * base pairs within a seed
	 * @param seedEncoding the explicit seed encoding to parse
	 * @return the maximal number of base pairs among all encoded seeds
	 */
	static
	size_t
	getSeedMaxBP( const std::string & seedEncoding );

	/**
	 * Replace the input variables i1 and i2 to values to within the given range
	 * that correspond to
	 *
	 * - the first seed (if the given index pair is no valid seed start or one
	 *   of the indices is out of range bounds)
	 * - the next seed according to some seed order
	 *
	 * The indices are not updated if the last seed within the range is given
	 * or no seed within the range could be found.
	 * It returns whether or not the input variables have been updated.
	 *
	 * Note, if changed, only the seed left-most base pair is within the range
	 * but the full seed indices might exceed i1max or i2max.
	 *
	 * @param i1 seq1 seed index to be changed
	 * @param i2 seq2 seed index to be changed
	 * @param i1min first position within seq1 (inclusive)
	 * @param i1max last position within seq1 (inclusive)
	 * @param i2min first position within seq2 (inclusive)
	 * @param i2max last position within seq2 (inclusive)
	 * @return true if the input variables have been changed; false otherwise
	 */
	virtual
	bool
	updateToNextSeed( size_t & i1, size_t & i2
			, const size_t i1min = 0, const size_t i1max = RnaSequence::lastPos
			, const size_t i2min = 0, const size_t i2max = RnaSequence::lastPos
			) const;

public:

	/**
	 * data structure to store the information of an explicit seed
	 */
	class SeedData {
	public:

		/**
		 * default construction
		 */
		SeedData();

		/**
		 * construction based on a given seed encoding
		 *
		 * Invalid encodings will produce an invalid Seed object (isValid()==false)
		 *
		 * @param seedEncoding the explicit seed to be parsed (has to be valid!)
		 * @param energyFunction the energy function to be used for index and energy setup
		 */
		SeedData( const std::string & seedEncoding, const InteractionEnergy & energyFunction );

		//! 5' start position in sequence 1
		size_t start1;
		//! 5' start position in sequence 2
		size_t start2;
		//! dot-bar encoding for sequence 1
		std::string dotBar1;
		//! dot-bar encoding for sequence 2
		std::string dotBar2;
		// flag about bp complementarity
		bool allBpCompl;
		//! energy of the seed interaction excluding the right base pair
		E_type energy;

		/**
		 * checks whether or not the based data is encoding a valid interaction
		 */
		bool
		isValid() const;

	};

protected:

	//! container to store
	boost::unordered_map< Interaction::BasePair, SeedData > seedForLeftEnd;

};

} /* namespace IntaRNA */

#endif /* INTARNA_SEEDHANDLEREXPLICIT_H_ */
