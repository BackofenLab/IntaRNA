
#ifndef INTARNA_OUTPUTHANDLERCSV_H_
#define INTARNA_OUTPUTHANDLERCSV_H_

#include "IntaRNA/general.h"

#include "IntaRNA/OutputHandler.h"
#include "IntaRNA/InteractionEnergy.h"

#include <list>
#include <map>

#include <boost/algorithm/string.hpp>

namespace IntaRNA {

/**
 * OutputHandler that stores interactions in CSV stream output.
 */
class OutputHandlerCsv : public OutputHandler
{

public:

	//! the column types supported
	enum ColType {
		id1 = 0, //!< id of first sequence
		id2, //!< id of second sequence
		seq1, //!< full first sequence
		seq2, //!< full second sequence
		subseq1, //!< interacting subsequence of first sequence
		subseq2, //!< interacting subsequence of second sequence
		subseqDP, //!< hybrid subsequences compatible with hybridDP
		subseqDB, //!< hybrid subsequences compatible with hybridDB
		start1, //!< start index of hybrid in seq1
		end1, //!< end index of hybrid in seq1
		start2, //!< start index of hybrid in seq2
		end2, //!< end index of hybrid in seq2
		hybridDP, //!< hybrid in VRNA dot-bracket notation
		hybridDB, //!< hybrid in dot-bar notation
		hybridDPfull, //!< hybrid in VRNA dot-bracket notation for full sequence lengths
		hybridDBfull, //!< hybrid in dot-bar notation for full sequence lengths
		E, //!< overall hybridization energy
		ED1, //!< ED value of seq1
		ED2, //!< ED value of seq2
		Pu1, //!< probability to be accessible for seq1
		Pu2, //!< probability to be accessible for seq2
		E_init, //!< initiation energy
		E_loops, //!< sum of loop energies (excluding E_init)
		E_dangleL, //!< dangling end contribution of base pair (start1,end2)
		E_dangleR, //!< dangling end contribution of base pair (end1,start2)
		E_endL, //!< penalty of closing base pair (start1,end2)
		E_endR, //!< penalty of closing base pair (end1,start2)
		E_hybrid, //!< energy of hybridization only = E - ED1 - ED2
		E_norm, //!< length normalized energy = E/ln(length(seq1)*length(seq2))
		E_hybridNorm, //!< length normalized energy of hybridization only = E_hybrid / ln(length(seq1)*length(seq2))
		seedStart1, //!< start index of the seed in seq1
		seedEnd1, //!< end index of the seed in seq1
		seedStart2, //!< start index of the seed in seq2
		seedEnd2, //!< end index of the seed in seq2
		seedE, //!< overall hybridization energy of the seed only (excluding rest)
		seedED1, //!< ED value of seq1 of the seed only (excluding rest)
		seedED2, //!< ED value of seq2 of the seed only (excluding rest)
		seedPu1, //!< probability of seed region to be accessible for seq1
		seedPu2, //!< probability of seed region to be accessible for seq2
		ColTypeNumber //!< number of column types
	};

	//! list of ColTypes
	typedef std::list<ColType> ColTypeList;

	/**
	 * Access to the mapping of ColTypes to according strings
	 */
	static
	const std::map<ColType,std::string> &
	getColType2string();

protected:

	//! mapping of ColTypes to according strings
	static std::map<ColType,std::string> colType2string;

	/**
	 * initializes the mapping of ColType elements to respective strings
	 */
	static
	void
	initColType2string() {
		if (colType2string.empty()) {
			// fill container : ensure strings are non-redundant!
			colType2string[id1] = "id1";
			colType2string[id2] = "id2";
			colType2string[seq1] = "seq1";
			colType2string[seq2] = "seq2";
			colType2string[start1] = "start1";
			colType2string[end1] = "end1";
			colType2string[start2] = "start2";
			colType2string[end2] = "end2";
			colType2string[subseq1] = "subseq1";
			colType2string[subseq2] = "subseq2";
			colType2string[subseqDP] = "subseqDP";
			colType2string[subseqDB] = "subseqDB";
			colType2string[hybridDP] = "hybridDP";
			colType2string[hybridDB] = "hybridDB";
			colType2string[hybridDPfull] = "hybridDPfull";
			colType2string[hybridDBfull] = "hybridDBfull";
			colType2string[E] = "E";
			colType2string[ED1] = "ED1";
			colType2string[ED2] = "ED2";
			colType2string[Pu1] = "Pu1";
			colType2string[Pu2] = "Pu2";
			colType2string[E_init] = "E_init";
			colType2string[E_loops] = "E_loops";
			colType2string[E_dangleL] = "E_dangleL";
			colType2string[E_dangleR] = "E_dangleR";
			colType2string[E_endL] = "E_endL";
			colType2string[E_endR] = "E_endR";
			colType2string[E_hybrid] = "E_hybrid";
			colType2string[E_norm] = "E_norm";
			colType2string[E_hybridNorm] = "E_hybridNorm";
			colType2string[seedStart1] = "seedStart1";
			colType2string[seedEnd1] = "seedEnd1";
			colType2string[seedStart2] = "seedStart2";
			colType2string[seedEnd2] = "seedEnd2";
			colType2string[seedE] = "seedE";
			colType2string[seedED1] = "seedED1";
			colType2string[seedED2] = "seedED2";
			colType2string[seedPu1] = "seedPu1";
			colType2string[seedPu2] = "seedPu2";
			// ensure filling is complete
			for (size_t i=0; i<ColTypeNumber; i++) {
				if ( colType2string.find( static_cast<ColType>(i) ) == colType2string.end() ) {
					throw std::runtime_error("OutputHandlerCsv::initColType2string() : ColType "+toString(i)+" without string representative");
				}
			}
		}
	}

public:

	/**
	 * Construct a CSV table output handler for interaction reporting.
	 *
	 * @param out the stream to write to
	 * @param energy the interaction energy object used for computation
	 * @param colOrder the order and list of columns to be printed
	 * @param colSep the column separator to be used in CSV output
	 * @param printHeader whether or not to print header information = col names
	 */
	OutputHandlerCsv( std::ostream & out
						, const InteractionEnergy & energy
						, const ColTypeList colOrder
						, const std::string& colSep = ";"
						, const bool printHeader = false
						);

	/**
	 * destruction
	 */
	virtual ~OutputHandlerCsv();

	/**
	 * Write a given RNA-RNA interaction in simple text format to the output
	 * stream.
	 *
	 * @param interaction the interaction to output
	 */
	virtual
	void
	add( const Interaction & interaction );

	/**
	 * Handles a given RNA-RNA interaction range as a
	 * RNA-RNA interaction with two base pairs and writes it in simple
	 * text format to the output stream.
	 *
	 * @param range the interaction range to add
	 */
	virtual
	void
	add( const InteractionRange & range );

	/**
	 * Converts a list of Coltypes to their string representation.
	 * @param colTypes the list to convert
	 * @return the string representation of the list
	 */
	static
	std::string
	list2string( const ColTypeList & colTypes );

	/**
	 * Converts string representation of a list of Coltypes to the respective
	 * list.
	 *
	 * @param listString the string representation of the list
	 * @return the parsed list
	 *
	 * @throws std::runtime_error in case of parsing problems
	 */
	static
	ColTypeList
	string2list( const std::string & listString );

	/**
	 * Generates the header line for a given list of columns
	 * @param colTypes the list of column types to consider
	 * @param colSep the column separator to be used
	 * @return the header line of the CSV output
	 */
	static
	std::string
	getHeader( const ColTypeList & colTypes, const std::string& colSep = ";" );

protected:

	//! the output stream to write to
	std::ostream & out;

	//! the interaction energy function used for interaction computation
	const InteractionEnergy & energy;

	//! the sequence of columns to be reported
	const std::list< ColType > colOrder;

	//! the column separator to be used
	std::string colSep;



};


//////////////////////////////////////////////////////////////////////////

inline
const std::map<OutputHandlerCsv::ColType,std::string> &
OutputHandlerCsv::
getColType2string()
{
	return colType2string;
}

//////////////////////////////////////////////////////////////////////////

inline
void
OutputHandlerCsv::
add( const InteractionRange & range )
{
	// forward to interaction reporting
	add( Interaction(range) );
}

//////////////////////////////////////////////////////////////////////////

inline
std::string
OutputHandlerCsv::
list2string( const ColTypeList & colTypes )
{
	// init string encodings
	initColType2string();
	// convert list
	std::stringstream ret;
	for (auto it = colTypes.begin(); it!=colTypes.end(); it++) {
		if (it!=colTypes.begin())
			ret <<',';
		ret <<colType2string[*it];
	}
	// return compiled string
	return ret.str();
}

//////////////////////////////////////////////////////////////////////////

inline
OutputHandlerCsv::
ColTypeList
OutputHandlerCsv::
string2list( const std::string & stringEncoding )
{
	// init string encodings
	initColType2string();

	ColTypeList list;
	// empty string : give full list
	if (stringEncoding.empty()) {
		// generate list of all types
		for (size_t c = 0; c < (size_t)ColTypeNumber; c++) {
			list.push_back( static_cast<ColType>(c) );
		}
	}
	else {
		// find split position
		size_t startPos = 0, splitPos = std::string::npos;
		while (startPos != splitPos) {
			splitPos = stringEncoding.find(',',startPos);
			// get current type string
			std::string curTypeString = stringEncoding.substr(startPos,splitPos-(splitPos==std::string::npos?0:startPos));
			// trim leading/trailing whitespaces
			boost::trim(curTypeString);
			// try to find type encoding
			bool notFound = true;
			for (auto it = colType2string.begin(); notFound && it != colType2string.end(); it++ ) {
				// check if type found (case insensitive for being user friendly
				if ( boost::iequals( it->second, curTypeString) ) {
					list.push_back( it->first );
					notFound = false;
				}
			}
			// check if error to be raised
			if (notFound) {
				throw std::runtime_error("OutputHandlerCsv::string2list("+stringEncoding+") contains unsupported ColType encoding '"+curTypeString+"'");
			}
			// update start of next interval encoding to parse
			startPos = splitPos + (splitPos != std::string::npos ? 1 : 0);
		}
	}
	return list;
}

//////////////////////////////////////////////////////////////////////////

inline
std::string
OutputHandlerCsv::
getHeader( const OutputHandlerCsv::ColTypeList & colList, const std::string& colSep )
{
	// init string encodings
	initColType2string();

	// generate header via stream
	std::stringstream header;
	for (auto col = colList.begin(); col != colList.end(); col++) {
		header	<<(col==colList.begin()?"":colSep)
				<<colType2string[*col];
	}
	header <<'\n';
	// return string
	return header.str();
}

//////////////////////////////////////////////////////////////////////////

} // namespace

#endif /* OUTPUTHANDLERCSV_H_ */
