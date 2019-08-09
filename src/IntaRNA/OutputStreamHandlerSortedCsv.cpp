
#include "IntaRNA/OutputStreamHandlerSortedCsv.h"

#include <vector>
#include <utility>

#include <boost/algorithm/string/iter_find.hpp>
#include <boost/algorithm/string/finder.hpp>

namespace IntaRNA {

/////////////////////////////////////////////////////////////////////////////

OutputStreamHandlerSortedCsv::
~OutputStreamHandlerSortedCsv()
{

	//! type to capture the row's value to sort on (first) and the overall row data (second)
	typedef std::pair< std::string, std::string > RowData;

	//! list of all reported CSV rows to be sorted for final output
	std::vector< RowData > rows;

	// parse unsorted stream
	std::string curRow;
	while( outStreamUnsorted.rdbuf()->in_avail() > 0 ) {
		// parse next row
		std::getline( outStreamUnsorted, curRow, '\n' );
		if (csvWithHeader && csvHeader.empty()) {
			csvHeader = curRow;
			continue;
		}
		// extract value for sorting
		std::vector<std::string> values;
		boost::iter_split(values, curRow, boost::algorithm::first_finder(colSep));
		// store data
		rows.push_back( RowData( (values.size() > colToSort) ? values.at(colToSort) : "", curRow ) );
	}

	// sort all rows
	if (sortLexOrder) {
		// string-based lex-order sorting
		std::sort(rows.begin(), rows.end());
	} else {
		// numerical sort
		if (listSep.empty()) {
			// single value sorting
			std::sort(rows.begin(), rows.end(), [&](const RowData& a, const RowData& b) {
				return std::stod(a.first) < std::stod(b.first);
			});
		} else {
			// first-value sorting of possible list of values
			std::sort(rows.begin(), rows.end(), [&](const RowData& a, const RowData& b) {
				// check on first element of list
				return
					std::stod(a.first.substr( 0, a.first.find(listSep) ))
					< std::stod(b.first.substr( 0, b.first.find(listSep) ));
			});
		}
	}

#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_outputStreamUpdate)
#endif
	{
		// write header if needed
		if (csvWithHeader) {
			(*outStream) <<csvHeader <<'\n';
		}
		// write sorted output to stream
		std::for_each( rows.begin(), rows.end(), [&](const RowData & d) { (*outStream) << d.second <<"\n"; });
	}

	// flush output
	outStream->flush();

	// disconnect outstream to avoid double deletion
	outStream = NULL;
	// delete wrapped handler
	INTARNA_CLEANUP(outStreamHandler);
}

/////////////////////////////////////////////////////////////////////////////



} /* namespace IntaRNA */

