
#ifndef INTERACTION_H_
#define INTERACTION_H_

#include <vector>
#include <utility>

#include "RnaSequence.h"

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
	const RnaSequence & s1;

	//! the second interaction partner
	const RnaSequence & s2;


	/**
	 * construction
	 * @param s1 the sequence of the first interaction partner
	 * @param s2 the sequence of the second interaction partner
	 */
	Interaction( const RnaSequence & s1, const RnaSequence & s2 );

	/**
	 * destruction
	 */
	virtual ~Interaction();

	/**
	 * Adds an interaction base pair to the container
	 * @param i1 index in sequence 1
	 * @param i2 index in sequence 2
	 */
	void
	addInteraction( const size_t i1, const size_t i2 );

	/**
	 * checks whether or not the current interaction is non-empty, nested and
	 * sorted.
	 * @return true if the interaction encoding is non-empty, nested and sorted;
	 *         false otherwise
	 */
	bool
	isValid() const;

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
	 * Access to the interacting base pairs between s1 and s2.
	 * @return the list of interacting base pairs (idx s1, idx s2)
	 */
	const PairingVec &
	getBasePairs() const;

protected:

	//! interacting indices
	PairingVec interaction;

};

#endif /* INTERACTION_H_ */
