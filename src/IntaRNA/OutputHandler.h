
#ifndef INTARNA_OUTPUTHANDLER_H_
#define INTARNA_OUTPUTHANDLER_H_

#include "IntaRNA/general.h"
#include "IntaRNA/Interaction.h"
#include "IntaRNA/InteractionRange.h"
#include "IntaRNA/OutputConstraint.h"
#include <string>

namespace IntaRNA {

/**
 * Defines an output and storage interface to enable different output formats
 * for all predicted RNA-RNA interactions.
 *
 * @author Martin Mann 2014
 *
 */
class OutputHandler {

public:

	/**
	 * Construction
	 * @param outConstraint the output constraints to be heeded
	 */
	OutputHandler( const OutputConstraint & outConstraint );

	/**
	 * Destruction
	 */
	virtual ~OutputHandler();


	/**
	 * Adds a given RNA-RNA interaction to the storage/output.
	 *
	 * NOTE: the given interaction object and its source object might be deleted
	 * later on. Thus, a deep copy is needed in order to make them completely
	 * independent.
	 *
	 * @param interaction the interaction to add
	 * @param outConstraint the output constraint applied to find the reported
	 *        interaction
	 */
	virtual
	void
	add( const Interaction & interaction ) = 0;

	/**
	 * Returns the number of reported interactions.
	 * @return the number of reported interactions
	 */
	virtual
	size_t
	reported() const;


	/**
	 * returns the reversed string
	 * @param str the string to reverse
	 * @return the reversed string
	 */
	static
	std::string
	reverse( const std::string & str );

	/**
	 * Increments the overall partition function with the given value.
	 *
	 * Note, take care that multiple increments have to represent disjoint subsets
	 * of all interactions for the same pair of sequences, otherwise the
	 * aggregated overall partition function will be wrong!
	 *
	 * @param subZ increment to be added to the overall partition function
	 */
	virtual
	void
	incrementZ( const Z_type subZ );

	/**
	 * Access to the partition function (initialized with 0) aggregated via
	 * incrementZ().
	 *
	 * @return the aggregated partition function Z
	 */
	virtual
	Z_type
	getZ() const;

	/**
	 * Access to the output constraints to be applied
	 * @return the OutputConstraint object to be heeded
	 */
	const OutputConstraint &
	getOutputConstraint() const;


protected:

	//! the output constraints to be applied
	const OutputConstraint outConstraint;

	//! number of reported interactions
	size_t reportedInteractions;

	//! overall partition function for the sequences provided in Z_energy
	Z_type Z;

};


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

inline
OutputHandler::OutputHandler( const OutputConstraint & outConstraint )
	: outConstraint(outConstraint)
	, reportedInteractions(0)
	, Z(0)
{
}

////////////////////////////////////////////////////////////////////////////

inline
OutputHandler::~OutputHandler() {
}

////////////////////////////////////////////////////////////////////////////

inline
std::string
OutputHandler::
reverse( const std::string & str )
{
	return std::string(str.rbegin(), str.rend());
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
OutputHandler::
reported() const
{
	return reportedInteractions;
}

////////////////////////////////////////////////////////////////////////////

inline
Z_type
OutputHandler::
getZ() const
{
	return Z;
}

////////////////////////////////////////////////////////////////////////////

inline
const OutputConstraint &
OutputHandler::
getOutputConstraint() const
{
	return outConstraint;
}

////////////////////////////////////////////////////////////////////////////

inline
void
OutputHandler::
incrementZ( const Z_type subZ )
{
#if INTARNA_IN_DEBUG_MODE
	if (std::numeric_limits<Z_type>::max() - subZ <= Z) {
		LOG(WARNING) <<"OutputHandler::incrementZ() : partition function overflow! Recompile with larger partition function data type!";
	}
#endif
	// increment partition function
#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_outputHandlerIncrementZ)
#endif
	{Z += subZ;}
}

////////////////////////////////////////////////////////////////////////////

} // namespace

#endif /* OUTPUTHANDLER_H_ */
