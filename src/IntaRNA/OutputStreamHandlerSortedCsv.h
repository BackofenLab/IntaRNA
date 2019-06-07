
#ifndef OUTPUTSTREAMHANDLERSORTEDCSV_H_
#define OUTPUTSTREAMHANDLERSORTEDCSV_H_

#include "IntaRNA/OutputStreamHandler.h"

#include <sstream>


namespace IntaRNA
{

/**
 * Provides access to a temporal output stream for interaction reporting, which
 * is assumed to be populated by an OutputHandlerCsv.
 * On every flush(), the temporary stream is parsed and the sorted output is
 * reported to the final output stream.
 *
 * @author: Martin Raden 2019
 */
class OutputStreamHandlerSortedCsv : public OutputStreamHandler
{
public:

	/**
	 * construction
	 * @param outStreamHandler the OutputStreamHandler to write output to (non-NULL) and which is to destroy on destruction
	 * @param colToSort the CSV column to sort
	 * @param sortLexOrder whether or not lexicographic (or numerical) sorting is to be done
	 * @param colSep the column separator used within the output
	 * @param csvWithHeader whether or not the CSV output will contain a header line
	 * @param listSep the separator used within a single column to represent multiple values
	 */
	OutputStreamHandlerSortedCsv( OutputStreamHandler * outStreamHandler
								, const size_t colToSort
								, const bool sortLexOrder
								, const std::string colSep
								, const bool csvWithHeader
								, const std::string listSep = ""
								);

	/**
	 * destruction:
	 * - parses the unsorted CSV output
	 * - writes the sorted CSV output to the final output stream
	 * - calls flush() on the final output
	 * - destruction of the underlying OutputStreamHandler
	 */
	virtual ~OutputStreamHandlerSortedCsv();

	/**
	 * Access to the ostream to which unsorted CSV output is to be written.
	 * @return the output stream to write unsorted CSV output to
	 */
	virtual
	std::ostream& getOutStream();


protected:

	//! the underlying OutputStreamHandler to which the final output is reported
	OutputStreamHandler* outStreamHandler;

	//! the stream to write output to
	std::stringstream outStreamUnsorted;

	//! the CSV column to sort
	const size_t colToSort;

	//! whether or not lexicographic (or numerical) sorting is to be done
	const bool sortLexOrder;

	//! the column separator used within the output
	const std::string colSep;

	//! whether or not the CSV output will contain a header line
	const bool csvWithHeader;

	//! the parsed CSV header (empty if not parsed yet)
	std::string csvHeader;

	//! the separator used within a single column to represent multiple values
	const std::string listSep;
};

/////////////////////////////////////////////////////////////////////////////

inline
OutputStreamHandlerSortedCsv::
OutputStreamHandlerSortedCsv( OutputStreamHandler * osh
			, const size_t colToSort
			, const bool sortLexOrder
			, const std::string colSep
			, const bool csvWithHeader
			, const std::string listSep
		)
 :	OutputStreamHandler( osh == NULL ? &std::cout : &(osh->getOutStream()) )
	, outStreamHandler(osh)
	, outStreamUnsorted()
	, colToSort(colToSort)
	, sortLexOrder(sortLexOrder)
	, colSep(colSep)
	, csvWithHeader(csvWithHeader)
	, csvHeader()
	, listSep(listSep)
{
	if (osh == NULL) {
		throw std::runtime_error("OutputStreamHandlerSortedCsv() : outStreamHandler == NULL");
	}
	if (colSep.empty()) {
		throw std::runtime_error("OutputStreamHandlerSortedCsv() : colSep is empty");
	}
}

/////////////////////////////////////////////////////////////////////////////

inline
std::ostream&
OutputStreamHandlerSortedCsv::
getOutStream()
{
	return outStreamUnsorted;
}

/////////////////////////////////////////////////////////////////////////////



} /* namespace IntaRNA */

#endif /* OUTPUTSTREAMHANDLERSORTEDCSV_H_ */
