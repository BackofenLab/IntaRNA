
#ifndef INTERACTION_H_
#define INTERACTION_H_

#include <vector>
#include <utility>

#include "general.h"
#include "RnaSequence.h"

// dummy declaration to avoid inclusion loop
class InteractionRange;


/**
 * Data structure to describe and store an Interaction
 */
class Interaction {

public:

	//! type of a base pair index encoding
	typedef std::pair<size_t,size_t> BasePair;

	//! type of a vector encoding base pair indices that are interacting
	typedef std::vector<BasePair> PairingVec;

public:

	//! the first interaction partner
	const RnaSequence * s1;

	//! the second interaction partner
	const RnaSequence * s2;

	//! interacting indices
	PairingVec basePairs;

	//! energy of the interaction (can be NaN)
	E_type energy;

	//! optional: range of seed interaction and full seed energy
	InteractionRange * seedRange;

	/**
	 * construction
	 * @param s1 the sequence of the first interaction partner
	 * @param s2 the sequence of the second interaction partner
	 */
	Interaction( const RnaSequence & s1, const RnaSequence & s2 );

	/**
	 * construction from interaction range, ie. forming base pairs at the
	 * beginning and end of the interaction range.
	 *
	 * NOTE: the ends have to be complementary nucleotides
	 *
	 * @param range the interaction range to be used for initialization
	 */
	Interaction( const InteractionRange & range );

	/**
	 * destruction
	 */
	virtual ~Interaction();

	/**
	 * checks whether or not the current interaction is non-empty, nested and
	 * sorted.
	 * @return true if the interaction encoding is non-empty, nested and sorted;
	 *         false otherwise
	 */
	bool
	isValid() const;

	/**
	 * Checks whether or not the interaction contains no base pairs.
	 * @return true if the interaction encodes no base pair; false otherwise
	 */
	bool
	isEmpty() const;

	/**
	 * Sorts the stored interacting base pairs.
	 */
	void
	sort();

	/**
	 * Clears all interaction information.
	 */
	void
	clear();

	/**
	 * Sets the seedRange member according to the given data
	 * @param ij1 left most base pair in seed
	 * @param ij2 right most base pair in seed
	 * @param energy full energy of the seed interaction
	 */
	void
	setSeedRange( const BasePair ij1, const BasePair ij2, const E_type energy );


	/**
	 * Creates an interaction with one base pair for each interaction range
	 * boundary.
	 *
	 * NOTE: the ends have to be complementary nucleotides
	 *
	 * NOTE: the interaction energy is reset too.
	 *
	 * @param range the interaction range to get the boundaries from
	 * @return the altered range object (*this)
	 */
	Interaction &
	operator= ( const InteractionRange & range );

	/**
	 * Prints the interacting base pairs to stream
	 * @param out the ostream to write to
	 * @param i the Interaction object to add
	 * @return the altered stream out
	 */
	friend std::ostream& operator<<(std::ostream& out, const Interaction& i);

	/**
	 * Prints the interacting base pair to stream
	 * @param out the ostream to write to
	 * @param bp the Interaction base pair object to add
	 * @return the altered stream out
	 */
	friend std::ostream& operator<<(std::ostream& out, const Interaction::BasePair& bp);

};

#endif /* INTERACTION_H_ */
