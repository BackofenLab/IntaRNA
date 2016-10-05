/*
 * RnaSequence.h
 *
 *  Created on: 24.06.2014
 *      Author: Mmann
 */

#ifndef RNASEQUENCE_H_
#define RNASEQUENCE_H_

#include <locale>
#include <string>
#include <vector>

#include "general.h"




/**
 * Represents an RNA sequence within IntaRNA and holds all necessary utility
 * functions.
 *
 * The integer sequence encoding uses the Vienna RNA package encoding functions.
 *
 * @author Martin Mann 2014
 */
class RnaSequence {

public:

	//! type for sequence string representation
	typedef std::string String_type;

	//! type for integer encoding of a single sequence letter
	//! (currently same as in Vienna package)
	typedef short Code_type;

	//! type for sequence integer encoded representation
	typedef std::vector<Code_type> CodeSeq_type;


	/**
	 * Allowed nucleotide single letter character alphabet and its lower case
	 * variants.
	 */
	const static std::string SequenceAlphabet;

	//! constant that serves as a placeholder for the last position within a
	//! sequence independent of its real size
	static const size_t lastPos;

public:

	/**
	 * Creates a sequence.
	 * @param id the name of the sequence
	 * @param seqString the nucleotide string encoding the sequence
	 *
	 */
	RnaSequence(const std::string& id
			, const std::string & seqString );

	/**
	 * Destruction and garbage collection
	 */
	virtual ~RnaSequence();


	/////////////////////  DATA ACCESS  //////////////////////////////

	/**
	 * Access to the sequence's ID.
	 * @return the ID of the sequence.
	 */
	const std::string&
	getId() const;

	/**
	 * Access to the length of the sequence.
	 * @return the length of the sequence.
	 */
	size_t
	size() const;

	/**
	 * Access to the sequence in character encoding.
	 * @return the character encoding of the sequence
	 */
	const String_type&
	asString() const;

	/**
	 * Access to the sequence in integer coding using the Vienna RNA package
	 * encoding.
	 * @return the integer encoding of the sequence
	 */
	const CodeSeq_type&
	asCodes() const;

	/**
	 * Whether or not the sequence contains ambiguous nucleotide encodings.
	 * @return true if the sequence contains ambiguous nucleotide encodings;
	 *         false otherwise
	 */
	bool
	isAmbiguous() const;

	/**
	 * Whether or not a specific sequence position shows an ambiguous nucleotide
	 * encoding.
	 * @param i the sequence position of interest
	 * @return true if the sequence contains an ambiguous nucleotide encoding
	 *         at position i;
	 *         false otherwise
	 */
	bool
	isAmbiguous( const size_t i ) const;

	/**
	 * prints the sequence id and the sequence to stream
	 * @out the ostream to write to
	 * @rna the RnaSequence object to add
	 * @return the altered stream out
	 */
	friend std::ostream& operator<<(std::ostream& out, const RnaSequence& rna);

public:

	/////////////////////  STATIC UTILITY  //////////////////////////////

	/**
	 * Utility function that converts a string into the internal sequence's
	 * string representation. This covers an upper case conversion.
	 * @param seqString sequence string
	 * @return the internally used sequence string representation
	 */
	static
	String_type
	getUpperCase( const std::string & seqString );


	/**
	 * Utility function that converts a sequence's string representation into
	 * its integer encoding
	 *
	 * NOTE: encoding covers only upper case characters! Otherwise an exception
	 * is raised.
	 *
	 * @param seqString the string to encode
	 * @param return the integer encoding
	 * @throw std::runtime_error if an unsupported nucleotide character is given
	 */
	static
	CodeSeq_type
	getCodeForString( const String_type& seqString );

	/**
	 * Utility function that converts a nucleotide's char representation into
	 * its integer encoding.
	 *
	 * NOTE: encoding covers only upper case characters! Otherwise an exception
	 * is raised.
	 *
	 * @param nucleotide the char to encode
	 * @param return the integer encoding
	 * @throw std::runtime_error if an unsupported nucleotide character is given
	 */
	static
	Code_type
	getCodeForChar( const char nucleotide );


	/**
	 * Utility function that tests whether or not a given sequence is a valid
	 * RNA sequence
	 */
	static
	bool
	isValidSequence( const std::string& sequence );


	/**
	 * Whether or not two positions of two RNAs are complementary, ie. can
	 * form a base pair
	 * @param s1 first RNA
	 * @param s2 second RNA
	 * @param p1 position within s1
	 * @param p2 position within s2
	 * @return true if s1[p1] can form a base pair with s2[p2]; false otherwise
	 */
	static
	bool
	areComplementary( const RnaSequence & s1, const RnaSequence & s2,
					const size_t p1, const size_t p2 );

protected:

	/////////////////////  DATA MEMBERS  //////////////////////////////


	//! the locale to use for integer encoding
	static std::locale codeLocale;

	//! ID of this sequence
	std::string id;

	//! The sequence's string representation.
	String_type seqString;

	//! Integer encoding of the sequence.
	CodeSeq_type seqCode;

	//! Whether or not the sequence contains ambiguous nucleotide encodings
	bool ambiguous;

};

#endif /* RNASEQUENCE_H_ */
