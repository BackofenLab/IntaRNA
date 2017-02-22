
#include "general.h"

#include <fstream>

#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>


/////////////////////////////////////////////////////////////////////

std::ostream *
newOutputStream( const std::string & out )
{
	// check if empty or whitespace string
	if (boost::regex_match( out, boost::regex("^\\s*$"), boost::match_perl)) {
		return NULL;
	}
	// open according stream
	if (boost::iequals(out,"STDOUT")) {
		return & std::cout;
	} else
	if (boost::iequals(out,"STDERR")) {
		return & std::cerr;
	} else {
		// open file stream
		std::fstream * outFileStream = new std::fstream();
		outFileStream->open( out.c_str(), std::ios_base::out );
		if (!outFileStream->is_open()) {
			delete outFileStream;
			return NULL;
		} else {
			// set output stream
			return outFileStream;
		}
	}
}

/////////////////////////////////////////////////////////////////////

void
deleteOutputStream( std::ostream * outStream )
{
	// check if something to be done
	if (outStream == NULL) {
		return;
	}
	// flush content
	outStream->flush();

	// check if to be closed and deleted
	std::fstream * outFile = dynamic_cast<std::fstream *>(outStream);
	if (outFile != NULL) {
		// close and delete file handle
		outFile->close();
		CLEANUP(outFile);
#if IN_DEBUG_MODE
	} else {
		// sanity check
		if (outStream != &std::cout && outStream != &std::cerr) {
			throw std::runtime_error("deleteOutputStream() : no file, nor STDOUT/STDERR");
		}
#endif
	}
}

/////////////////////////////////////////////////////////////////////



