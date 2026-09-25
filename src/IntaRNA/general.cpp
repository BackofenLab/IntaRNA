
#include "IntaRNA/general.h"

#include <string>
#include <memory>
#include <algorithm>
#include <fstream>
#include <iostream>

#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/algorithm/string/predicate.hpp>

namespace IntaRNA {

/////////////////////////////////////////////////////////////////////

std::ostream* newOutputStream(const std::string& out)
{
    // 1. Whitespace-Check
    if (out.empty() || out.find_first_not_of(" \t\n\r\f\v") == std::string::npos) {
        return NULL;
    }

    // 2. Standard-Streams
    if (boost::iequals(out, "STDOUT")) {
        return &std::cout;
    } 
    if (boost::iequals(out, "STDERR")) {
        return &std::cerr;
    }

    // 3. File-Streams (with optional gzip compression)
    namespace bio = boost::iostreams;
    
    // Puffergröße für hohe Performance & reduziertes File-Handle-Overhead (z.B. 512 KB)
    constexpr std::streamsize BUFFER_SIZE = 512 * 1024;

    // Verwende std::unique_ptr für Exception-Safety (Leak-Schutz)
    auto fstream = std::make_unique<bio::filtering_ostream>();
    BOOST_IOS::openmode fopenmode = BOOST_IOS::out;

	// Gzip-Erkennung
	if (boost::iends_with(out, ".gz")) {
		// Gzip-Kompressor mit explizit vergrößertem Puffer hinzufügen
		fstream->push(bio::gzip_compressor(), BUFFER_SIZE);
		fopenmode |= BOOST_IOS::binary;
	}

	// File Sink öffnen
	bio::file_descriptor_sink sink(out, fopenmode);

	// Prüfen, ob die Datei tatsächlich geöffnet werden konnte
	if (!sink.is_open()) {
		return NULL;
	}

	// Sink mit großem Puffer zur Pipeline hinzufügen
	fstream->push(sink, BUFFER_SIZE);

	// Prüfen, ob die Pipeline bereit ist und der Stream sich in einem korrekten Zustand befindet
	if (fstream->is_complete() && fstream->good()) {
		return fstream.release(); // Ownership an den Aufrufer übergeben
	}

    return NULL;
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
	bio::filtering_istream * inFile = dynamic_cast<bio::filtering_istream *>(inStream);
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
