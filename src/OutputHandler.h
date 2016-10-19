
#ifndef OUTPUTHANDLER_H_
#define OUTPUTHANDLER_H_

#include "general.h"
#include "Interaction.h"
#include "InteractionRange.h"
#include <string>

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
	 */
	OutputHandler();

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
	 */
	virtual
	void
	add( const Interaction & interaction ) = 0;

	/**
	 * Adds a given RNA-RNA interaction range to the storage/output.
	 *
	 * NOTE: the given interaction object and its source object might be deleted
	 * later on. Thus, a deep copy is needed in order to make them completely
	 * independent.
	 *
	 * @param range the interaction range to add
	 */
	virtual
	void
	add( const InteractionRange & range ) = 0;

	/**
	 * returns the reversed string
	 * @param str the string to reverse
	 * @return the reversed string
	 */
	static
	std::string
	reverse( const std::string & str );

};

#endif /* OUTPUTHANDLER_H_ */
