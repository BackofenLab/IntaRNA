
#ifndef OUTPUTSTREAMHANDLER_H_
#define OUTPUTSTREAMHANDLER_H_

#include "IntaRNA/general.h"

namespace IntaRNA
{

/**
 * Provides access to the output stream for interaction reporting.
 *
 * @author: Martin Raden 2019
 */
class OutputStreamHandler
{
public:

	/**
	 * construction
	 * @param outStream the stream to write output to (non-NULL)
	 */
	OutputStreamHandler( std::ostream * outStream );

	/**
	 * destruction: if the provided outStream is a file stream, it is closed and
	 * destroyed.
	 */
	virtual ~OutputStreamHandler();

	/**
	 * Access to the ostream to which output is to be written.
	 * @return the output stream to write to
	 */
	virtual
	std::ostream& getOutStream();

protected:

	//! the stream to write output to
	std::ostream * outStream;
};

/////////////////////////////////////////////////////////////////////////////

inline
OutputStreamHandler::
OutputStreamHandler( std::ostream * os )
 :	outStream(os)
{
	if (outStream == NULL) {
		throw std::runtime_error("OutputStreamHandler : provided outStream is NULL");
	}
}

/////////////////////////////////////////////////////////////////////////////

inline
OutputStreamHandler::
~OutputStreamHandler()
{
	// cleanup output stream
	deleteOutputStream(outStream);
}

/////////////////////////////////////////////////////////////////////////////

inline
std::ostream&
OutputStreamHandler::
getOutStream()
{
	return *outStream;
}

/////////////////////////////////////////////////////////////////////////////



} /* namespace IntaRNA */

#endif /* OUTPUTSTREAMHANDLER_H_ */
