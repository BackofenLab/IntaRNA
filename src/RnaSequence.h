
#ifndef RNASEQUENCE_H_
#define RNASEQUENCE_H_

#include <locale>
#include <string>
#include <vector>

#include "general.h"


#ifndef VIENNA_RNA_PAIR_MAT_H
#define VIENNA_RNA_PAIR_MAT_H
extern "C" {
	#include "ViennaRNA/pair_mat.h"
}
#endif



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


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////


inline
RnaSequence::RnaSequence(
		const std::string & id
		, const std::string & seqString )
 :
	id(id)
	, seqString(getUpperCase(seqString))
	, seqCode(getCodeForString(this->seqString))
	, ambiguous(this->seqString.find('N')!=std::string::npos)
{
#if IN_DEBUG_MODE
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
#if IN_DEBUG_MODE
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
#if IN_DEBUG_MODE
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
#if IN_DEBUG_MODE
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




#endif /* RNASEQUENCE_H_ */
