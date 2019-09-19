
#ifndef INTARNA_RNASEQUENCE_H_
#define INTARNA_RNASEQUENCE_H_

#include <locale>
#include <string>
#include <vector>

#include "IntaRNA/general.h"


#ifndef VIENNA_RNA_PAIR_MAT_H
#define VIENNA_RNA_PAIR_MAT_H
extern "C" {
	#include <ViennaRNA/pair_mat.h>
}
#endif

namespace IntaRNA {

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
	 * Allowed nucleotide single letter character alphabet according to IUPAC
	 * codes and their lower case variants.
	 */
	const static std::string SequenceAlphabetIUPAC;

	/**
	 * Allowed nucleotide single letter character alphabet (uppercase only)
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
	 * @param idxPos0 input/output index of the first sequence position
	 */
	RnaSequence(const std::string& id
			, const std::string & seqString
			, const long idxPos0 = 1 );

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
	 * Provides the input/output index of a given sequence position.
	 * @param i the position of interest
	 * @return the input/output index shifted by idxPos0
	 */
	long
	getInOutIndex( const size_t i ) const;

	/**
	 * Provides the internal sequence position of an input/output index given
	 * that the first position represents idxPos0.
	 * @param i the input/output index of interest
	 * @return the internal position, i.e. i shifted by idxPos0
	 */
	size_t
	getIndex( const long i ) const;

	/**
	 * Provides the reverse index of a given sequence position.
	 * @param i the index of interest
	 * @return the reverse index, i.e. (seq.size()-i-1)
	 */
	size_t
	getReversedIndex( const size_t i ) const;

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
	 * Checks for equality
	 * @param rna2 the RnaSequence to compare to
	 * @return true if both are describing the same sequence with the same identifier.
	 */
	const bool operator == ( const RnaSequence &rna2 ) const;

	/**
	 * prints the sequence id and the sequence to stream
	 * @param out the ostream to write to
	 * @param rna the RnaSequence object to add
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
	 * @return the integer encoding
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
	 * @return the integer encoding
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
	 * Utility function that tests whether or not a given sequence is a valid
	 * RNA sequence according to IUPAC encoding
	 */
	static
	bool
	isValidSequenceIUPAC( const std::string& sequence );


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

	/**
	 * Whether or not two positions of two RNAs are forming a GU base pair
	 * @param s1 first RNA
	 * @param s2 second RNA
	 * @param p1 position within s1
	 * @param p2 position within s2
	 * @return true if (s1[p1]=G and s2[p2]==U) or (s1[p1]=U and s2[p2]==G);
	 *         false otherwise
	 */
	static
	bool
	isGU( const RnaSequence & s1, const RnaSequence & s2,
					const size_t p1, const size_t p2 );

protected:

	/////////////////////  DATA MEMBERS  //////////////////////////////


	//! the locale to use for integer encoding
	static std::locale codeLocale;

	//! codes for G and U to check for GU base pairs
	static int bpGUcodes[];

	//! ID of this sequence
	std::string id;

	//! The sequence's string representation.
	String_type seqString;

	//! Integer encoding of the sequence.
	CodeSeq_type seqCode;

	//! Whether or not the sequence contains ambiguous nucleotide encodings
	bool ambiguous;

	//! Input/output index of the first sequence position
	long idxPos0;

};


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////


inline
RnaSequence::RnaSequence(
		const std::string & id
		, const std::string & seqString
		, const long idxPos0 )
 :
	id(id)
	, seqString(getUpperCase(seqString))
	, seqCode(getCodeForString(this->seqString))
	, ambiguous(this->seqString.find('N')!=std::string::npos)
	, idxPos0(idxPos0)
{
#if INTARNA_IN_DEBUG_MODE
	if (id.size() == 0) {
		throw std::runtime_error("RnaSequence::RnaSequence : id empty");
	}
	if (seqString.size() == 0) {
		throw std::runtime_error("RnaSequence::RnaSequence : seqString empty");
	}
#endif

}

/////////////////////////////////////////////////////////////////////////////

inline
RnaSequence::~RnaSequence()
{
}

/////////////////////////////////////////////////////////////////////////////

inline
const std::string&
RnaSequence::
getId() const
{
	return id;
}

/////////////////////////////////////////////////////////////////////////////

inline
size_t
RnaSequence::
size() const
{
	return seqString.size();
}

/////////////////////////////////////////////////////////////////////////////

inline
long
RnaSequence::
getInOutIndex( const size_t i ) const
{
#if INTARNA_IN_DEBUG_MODE
	if (i >= size()) {
		throw std::runtime_error("RnaSequence::getInOutIndex : index "+toString(i)+" >= length "+toString(size()));
	}
#endif
	// get in/out index
	long p = idxPos0 + (long)i;
	// check for -+1 index transition
	if (idxPos0<0 && p>=0) {
		p++;
	}
	// final in/out index
	return p;
}

/////////////////////////////////////////////////////////////////////////////

inline
size_t
RnaSequence::
getIndex( const long i ) const
{
	// check lower bounds
	if (i < idxPos0) {
		throw std::runtime_error("RnaSequence::getIndex : index "+toString(i)+" < idxPos0 "+toString(idxPos0));
	}
	// shift to internal indexing
	size_t p = (size_t)(i - idxPos0);
	// check for -+1 index transition
	if (idxPos0<0 && i>=0) {
		assert(i>0);
		p--;
	}
	// check upper bound
	if (p >= size()) {
		throw std::runtime_error("RnaSequence::getIndex : index "+toString(i)+" relates to "+toString(p)+" >= length "+toString(size()));
	}
	// return internal index
	return p;
}

/////////////////////////////////////////////////////////////////////////////

inline
size_t
RnaSequence::
getReversedIndex( const size_t i ) const
{
#if INTARNA_IN_DEBUG_MODE
	if (i >= size()) {
		throw std::runtime_error("RnaSequence::getReversedIndex : index "+toString(i)+" >= seq.length "+toString(size()));
	}
#endif

	return this->size() -i -1;
}

/////////////////////////////////////////////////////////////////////////////

inline
const
RnaSequence::
String_type&
RnaSequence::
asString() const
{
	return seqString;
}

/////////////////////////////////////////////////////////////////////////////

inline
const
RnaSequence::
CodeSeq_type&
RnaSequence::
asCodes() const
{
	return seqCode;
}

/////////////////////////////////////////////////////////////////////////////

inline
bool
RnaSequence::
isAmbiguous() const
{
	return ambiguous;
}

/////////////////////////////////////////////////////////////////////////////

inline
bool
RnaSequence::
isAmbiguous( const size_t i ) const
{
	// check for ambiguous nucleotide encoding
	return this->seqString.at(i) == 'N';
}

/////////////////////////////////////////////////////////////////////////////

inline
RnaSequence::
String_type
RnaSequence::
getUpperCase( const std::string & seqString )
{
#if INTARNA_IN_DEBUG_MODE
	if (!isValidSequenceIUPAC(seqString)) throw std::runtime_error("RnaSequence::getUpperCase() : the given sequence contains non-IUPAC codes : '"+seqString+"'");
#endif

	// create container to fill
	String_type seqRet(seqString.size(),'_');

	for (size_t i=0; i<seqString.size(); ++i)
	{
		// get upper case characters
		seqRet[i] = std::toupper(seqString.at(i),codeLocale);
		if (seqRet[i]=='T') {
			seqRet[i] = 'U';
		}
		// overwrite non-ACGU characters with N = ambiguous (no distinction needed)
		if (SequenceAlphabet.find(seqRet[i]) == -1) {
			seqRet[i] = 'N';
		}
	}

	return seqRet;
}

/////////////////////////////////////////////////////////////////////////////

inline
RnaSequence::
CodeSeq_type
RnaSequence::
getCodeForString( const String_type& seqString )
{
	// create container to fill and init with 'X'
	CodeSeq_type seqCode(seqString.size());

	for (size_t i=0; i<seqString.size(); ++i)
	{
		seqCode[i] = getCodeForChar( seqString.at(i) );
	}

	return seqCode;
}

/////////////////////////////////////////////////////////////////////////////

inline
bool
RnaSequence::
isValidSequence( const std::string & sequence )
{
	// check whether or not the string contains unsupported characters
	return sequence.find_first_not_of(SequenceAlphabet) == std::string::npos;
}

/////////////////////////////////////////////////////////////////////////////

inline
bool
RnaSequence::
isValidSequenceIUPAC( const std::string & sequence )
{
	// check whether or not the string contains unsupported characters
	return sequence.find_first_not_of(SequenceAlphabetIUPAC) == std::string::npos;
}

/////////////////////////////////////////////////////////////////////////////

inline
std::ostream&
operator<<(std::ostream& out, const RnaSequence& rna)
{
	// add ID(SEQUENCE) to stream
	out <<rna.id <<'(' <<rna.asString() <<')';
	return out;
}

/////////////////////////////////////////////////////////////////////////////

inline
RnaSequence::
Code_type
RnaSequence::
getCodeForChar( const char nucleotide )
{
#if INTARNA_IN_DEBUG_MODE
	// check if nucleotide character is NOT supported
	if (SequenceAlphabet.find(nucleotide) == std::string::npos)
		throw std::runtime_error("RnaSequence::getCodeForChar() : unsupported nucleotide character '"+toString(nucleotide)+"' in sequence");
#endif

	// otherwise get encoding:
	// use nucleotide character encoding from Vienna RNA package
	return (Code_type)encode_char(nucleotide);
}

/////////////////////////////////////////////////////////////////////////////

inline
bool
RnaSequence::
areComplementary( const RnaSequence & s1, const RnaSequence & s2,
					const size_t p1, const size_t p2 )
{
#if INTARNA_IN_DEBUG_MODE
	// check if valid positions
	if (p1>=s1.size() || p2>=s2.size())
		throw std::runtime_error("RnaSequence::areComplementary : index positions p1/p2 ("
				+ toString(p1)+"/"+toString(p2)
				+ ") are out of bounds s1/s2 ("
				+ toString(s1.size())+"/"+toString(s2.size())
				+")"
				);
#endif

	// check via VRNA util
	return BP_pair[s1.seqCode.at(p1)][s2.seqCode.at(p2)] > 0;
}

/////////////////////////////////////////////////////////////////////////////

inline
bool
RnaSequence::
isGU( const RnaSequence & s1, const RnaSequence & s2,
					const size_t p1, const size_t p2 )
{
#if INTARNA_IN_DEBUG_MODE
	// check if valid positions
	if (p1>=s1.size() || p2>=s2.size())
		throw std::runtime_error("RnaSequence::areComplementary : index positions p1/p2 ("
				+ toString(p1)+"/"+toString(p2)
				+ ") are out of bounds s1/s2 ("
				+ toString(s1.size())+"/"+toString(s2.size())
				+")"
				);
#endif

	// get bp code
	int bpCode = BP_pair[s1.seqCode.at(p1)][s2.seqCode.at(p2)];

	// check if GU pair
	return bpCode == bpGUcodes[0] || bpCode == bpGUcodes[1] ;
}

/////////////////////////////////////////////////////////////////////////////

inline
const bool
RnaSequence::
operator == ( const RnaSequence &rna2 ) const
{
	return	// check if same object (pointer)
			this == &rna2
			|| (
				// ensure same lengths
				this->size() == rna2.size()
				// ids are identical (most likely shorter, i.e. faster check)
				&& this->id == rna2.id
				// sequences identical
				&& this->seqString == rna2.seqString
			);
}

/////////////////////////////////////////////////////////////////////////////

} // namespace


#endif /* RNASEQUENCE_H_ */
