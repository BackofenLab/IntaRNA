
#include "IntaRNA/general.h"

#include <fstream>

#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>

namespace IntaRNA {

/////////////////////////////////////////////////////////////////////

std::ostream *
newOutputStream( const std::string & out )
{
	// check if empty or whitespace string
	if (boost::regex_match( out, boost::regex(R"(^\s*$)"), boost::match_perl)) {
		return NULL;
	}
	// open according stream
	if (boost::iequals(out,"STDOUT")) {
		return & std::cout;
	} else
	if (boost::iequals(out,"STDERR")) {
		return & std::cerr;
	} else {
		// file output
		namespace bio = boost::iostreams;
		bio::filtering_ostream* fstream = new bio::filtering_ostream();
		BOOST_IOS::openmode fopenmode = BOOST_IOS::out;

		// gzipped output file stream
		if (out.size()>3 && boost::iequals(out.substr(out.size()-3,3),".gz")) {
			// gzip compression
			fstream->push( bio::gzip_compressor() );
			// binary output
			fopenmode = BOOST_IOS::out | BOOST_IOS::binary;
		}

		// register final file
		fstream->push( bio::file_descriptor_sink( out, fopenmode ) );

		// check if all went fine so far
		if (fstream->is_complete()) {
			return fstream;
		} else {
			INTARNA_CLEANUP(fstream);
			return NULL;
		}
	}
}

/////////////////////////////////////////////////////////////////////

void
deleteOutputStream( std::ostream *& outStream )
{
	// check if something to be done
	if (outStream == NULL) {
		return;
	}
	// flush content
	outStream->flush();

	// handle file output
	namespace bio = boost::iostreams;
	bio::filtering_ostream * outFile = dynamic_cast<bio::filtering_ostream *>(outStream);
	if (outFile != NULL) {
		// ensure devices are closed on destruction
		outFile->set_auto_close(true);
		// close all file handles
		outFile->clear();
		// delete stream
		INTARNA_CLEANUP(outStream);
	}

	// ensure NULL setting
	outStream = NULL;
}

/////////////////////////////////////////////////////////////////////

std::istream *
newInputStream( const std::string & in )
{
	// check if empty or whitespace string
	if (boost::regex_match( in, boost::regex(R"(^\s*$)"), boost::match_perl)) {
		return NULL;
	}
	// open according stream
	if (boost::iequals(in,"STDIN")) {
		return & std::cin;
	} else {
		// file input
		namespace bio = boost::iostreams;
		bio::filtering_istream* fstream = new bio::filtering_istream();
		BOOST_IOS::openmode fopenmode = BOOST_IOS::in;

		// gzipped input file stream
		if (in.size()>3 && boost::iequals(in.substr(in.size()-3,3),".gz")) {
			// gzip compression
			fstream->push( bio::gzip_decompressor() );
			// binary input
			fopenmode = BOOST_IOS::in | BOOST_IOS::binary;
		}

		// register final file
		fstream->push( bio::file_descriptor_source( in, fopenmode ) );

		// check if all went fine so far
		if (fstream->is_complete()) {
			return fstream;
		} else {
			INTARNA_CLEANUP(fstream);
			return NULL;
		}
	}
}

/////////////////////////////////////////////////////////////////////

void
deleteInputStream( std::istream *& inStream )
{
	// check if something to be done
	if (inStream == NULL) {
		return;
	}

	// handle file input
	namespace bio = boost::iostreams;
	bio::filtering_ostream * inFile = dynamic_cast<bio::filtering_ostream *>(inStream);
	if (inFile != NULL) {
		// ensure devices are closed on destruction
		inFile->set_auto_close(true);
		// close all file handles
		inFile->clear();
		// delete stream
		INTARNA_CLEANUP(inStream);
	}

	// ensure NULL setting
	inStream = NULL;
}

/////////////////////////////////////////////////////////////////////


} // namespace
