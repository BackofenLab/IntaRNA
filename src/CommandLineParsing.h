/*
 * CommandLineParsing.h
 *
 *  Created on: 18.06.2014
 *      Author: Mmann
 */

#ifndef COMMANDLINEPARSING_H_
#define COMMANDLINEPARSING_H_

#include "config.h"
#include "general.h"
#include "RnaSequence.h"


#include <boost/program_options.hpp>

	namespace po = boost::program_options;

#include <iostream>
#include <cstdarg>

#include "Accessibility.h"
#include "Energy.h"
#include "Predictor.h"
#include "OutputHandler.h"

/**
 * Central handler for all command line arguments etc.
 *
 */
class CommandLineParsing {
public:

	//! type for a set of sequences
	typedef std::vector< RnaSequence > RnaSequenceVec;

public:

	/**
	 * Constructs a commandline argument parser for IntaRNA.
	 *
	 * @param logStream the stream to write validation log messages to
	 */
	CommandLineParsing( std::ostream& logStream );
	virtual ~CommandLineParsing();

	/**
	 * Parses the commandline arguments as passed to the 'main' method
	 * @param argc the number of arguments
	 * @param argv the argument array
	 *
	 * @return the parsing code: 1 (>0) if all went fine, otherwise the exit code to return
	 *
	 */
	int
	parse( int argc, char ** argv );

public:

	/////////  GETTERS  ///////////////////////////////////////////////////

	/**
	 * Access to the minimal number of inter-molecular base pairs in the seed.
	 * @return the minimal number of base pairs in the seed region
	 */
	int getSeedMinBp() const;

	/**
	 * Parses the query parameter and returns all parsed sequences.
	 * @return the set of parsed query sequences
	 */
	const RnaSequenceVec& getQuerySequences() const;

	/**
	 * Parses the target parameter and returns all parsed sequences.
	 * @return the set of parsed target sequences
	 */
	const RnaSequenceVec& getTargetSequences() const;

	/**
	 * Returns a newly allocated Accessibility object for the given query
	 * sequence according to the user defined parameters.
	 * @param sequenceNumber the number of the sequence within the vector
	 *         returned by getQuerySequences()
	 * @return a newly allocated Accessibility object or NULL in error case
	 */
	Accessibility* getQueryAccessibility( const size_t sequenceNumber ) const;

	/**
	 * Returns a newly allocated Accessibility object for the given target
	 * sequence according to the user defined parameters.
	 * @param sequenceNumber the number of the sequence within the vector
	 *         returned by getTargetSequences()
	 * @return a newly allocated Accessibility object or NULL in error case
	 */
	Accessibility* getTargetAccessibility( const size_t sequenceNumber ) const;

	/**
	 * Returns a newly allocated Energy object according to the user defined
	 * parameters.
	 * @param accQuery the accessibility object of the query sequence
	 * @param accTarget the (reversed) accessibility object of the target sequence
	 * @return the newly allocated Energy object or NULL in error case
	 */
	Energy* getEnergyHandler( const Accessibility& accQuery, const ReverseAccessibility& accTarget ) const;

protected:

	/////////  PRIVATE STUFF  ////////////////////////////////////////////////

	/**
	 * Limits and values for a number parameter
	 */
	template <typename T>
	class NumberParameter {
	public:
		//! the value of the parameter
	  T val;
		//! the minimally allowed value
	  const T min;
		//! the maximally allowed value
	  const T max;
		//! the default value
	  const T def;
		/**
		 * construction feeding the members
		 * @param min the minimally allowed value
		 * @param max the maximally allowed value
		 * @param def the default value
		 */
	  NumberParameter( const T min, const T max, const T def )
	   : val(def), min(min), max(max), def(def)
	  {}
	    //! checks if the given value is in the allowed range [min,max]
	    //! @param value the value to check
	    //! @return true if in range; false otherwise
	  bool isInRange(const T& value) const {
		  return value >= min && value <= max;
	  }
	    //! checks whether or not val is in the allowed range [min,max]
	    //! @return true if in range; false otherwise
	  bool isInRange() const {
		  return isInRange(val);
	  }
	};

	/**
	 * Allowed alphabet for a single char parameter
	 */
	class CharParameter {
	public:
		  //! the value of the parameter
		char val;
		  //! the set of allowed values for this parameter as a string
		const std::string alphabet;
		  //! the default value of the parameter
		const char def;
		/**
		 * Construction and member setup
		 * @param alphabet the allowed set of character values
		 * @param def the default value (has to be part of the alphabet)
		 */
		CharParameter( const std::string& alphabet, const char def )
		  : val(def), alphabet(alphabet), def(def)
		{
			if (alphabet.find(def) == std::string::npos) {
				throw std::runtime_error("CharParameter() : default value '"+toString(def)+"' is not within alphabet '"+alphabet+"'");
			}
		}
		//! checks if the given value is in the allowed alphabet
		//! @param value the value to check
		//! @return true if in alphabet; false otherwise
		bool isInAlphabet(const char value) const {
		  return alphabet.find(value) != std::string::npos;
		}
		//! checks whether or not val is in the allowed alphabet
		//! @return true if in alphabet; false otherwise
		bool isInAlphabet() const {
		  return isInAlphabet(val);
		}

	};


	//! query specific options
	po::options_description opts_query;
	//! target specific options
	po::options_description opts_target;
	//! seed specific options
	po::options_description opts_seed;
	//! interaction/energy specific options
	po::options_description opts_inter;
	//! general options
	po::options_description opts_general;

	//! overall option list
	po::options_description opts_cmdline_all;

	//! central result code to be set by validate_* functions in error case
	int parsingCode;
	//! value to mark that parsingCode was not set yet, ie. parse() not called
	const static int parsingCodeNotSet;

	//! stream to be used for log messages in validate_* functions
	std::ostream & logStream;

	//! the query command line argument
	std::string queryArg;
	//! the container holding all query sequences
	RnaSequenceVec query;
	//! accessibility computation mode for query sequences
	CharParameter qAcc;
	//! constraint for accessibility computation for query sequences
	std::string qAccConstr;
	//! window length to be considered accessible/interacting within query
	NumberParameter<int> qIntLenMax;
	//! maximal internal loop length to be considered accessible/interacting within query
	NumberParameter<int> qIntLoopMax;

	//! the target command line argument
	std::string targetArg;
	//! the container holding all target sequences
	RnaSequenceVec target;
	//! accessibility computation mode for target sequences
	CharParameter tAcc;
	//! constraint for accessibility computation for target sequences
	std::string tAccConstr;
	//! window length to be considered accessible/interacting within target
	NumberParameter<int> tIntLenMax;
	//! maximal internal loop length to be considered accessible/interacting within target
	NumberParameter<int> tIntLoopMax;

	//! minimal number of base pairs in seed
	NumberParameter<int> seedMinBP;

	//! the temperature to be used for energy computations
	NumberParameter<T_type> temperature;

	//! the selected energy model
	CharParameter energy;

protected:

	/**
	 * Validates the query sequence argument.
	 * @param value the argument value to validate
	 */
	void validate_query(const std::string & value);

	/**
	 * Validates the query accessibility argument.
	 * @param value the argument value to validate
	 */
	void validate_qAcc(const char & value);

	/**
	 * Validates the query accessibility constraint argument.
	 * @param value the argument value to validate
	 */
	void validate_qAccConstr(const std::string & value);

	/**
	 * Validates the query's maximal accessibility argument.
	 * @param value the argument value to validate
	 */
	void validate_qIntLenMax(const int & value);

	/**
	 * Validates the query's maximal internal loop length argument.
	 * @param value the argument value to validate
	 */
	void validate_qIntLoopMax(const int & value);

	/**
	 * Validates the target sequence argument.
	 * @param value the argument value to validate
	 */
	void validate_target(const std::string & value);

	/**
	 * Validates the target accessibility argument.
	 * @param value the argument value to validate
	 */
	void validate_tAcc(const char & value);

	/**
	 * Validates the target accessibility constraint argument.
	 * @param value the argument value to validate
	 */
	void validate_tAccConstr(const std::string & value);

	/**
	 * Validates the target's maximal accessibility argument.
	 * @param value the argument value to validate
	 */
	void validate_tIntLenMax(const int & value);

	/**
	 * Validates the target's maximal internal loop length argument.
	 * @param value the argument value to validate
	 */
	void validate_tIntLoopMax(const int & value);

	/**
	 * Validates the seedMinBP arguments.
	 * @param value the argument value to validate
	 */
	void validate_seedMinBP(const int & value);

	/**
	 * Validates the temperature arguments.
	 * @param value the argument value to validate
	 */
	void validate_temperature(const T_type & value);

	/**
	 * Validates a CharParameter.
	 * @param argName the name of the parameter (for exception handling)
	 * @param param the parameter object
	 * @param value the value of the parameter to validate
	 */
	void validate_charArgument(const std::string & argName, const CharParameter& param, const char & value);


	/**
	 * Validates a NumberParameter.
	 * @param argName the name of the parameter (for exception handling)
	 * @param param the parameter object
	 * @param value the value of the parameter to validate
	 */
	template <typename T>
	void validate_numberArgument(const std::string & name, const NumberParameter<T> & param, const T& value)
	{
		// alphabet check
		if ( ! param.isInRange(value) ) {
			logStream <<"\n "<<name<<" = " <<value <<" : has to be in the range [" <<param.min <<","<<param.max<<"]\n";
			parsingCode = std::min(-1,parsingCode);
		}
	}


	/**
	 * Validates a sequence arguments.
	 * @param argument the argument name of this parameter
	 * @param value the argument value to validate
	 */
	void validate_sequenceArgument(const std::string& argument, const std::string & value);


	/**
	 * Validates a structure constraint argument.
	 * @param argument the argument name of this parameter
	 * @param value the argument value to validate
	 */

	void validate_structureConstraintArgument(const std::string& argument, const std::string & value);

	/**
	 * Parses the target parameter and returns all parsed sequences.
	 * @param paramName the name of the parameter (for exception handling)
	 * @param paramArg the given argument for the parameter
	 * @param sequences the container to fill
	 * TODO exception handling
	 */
	void parseSequences(const std::string & paramName,
					const std::string& paramArg,
					RnaSequenceVec& sequences );

	/**
	 * Checks whether or not a sequence container holds a specific number of
	 * sequences.
	 * @param paramName the name of the parameter that enforces this constraint
	 *        (for exception handling)
	 * @param sequences the sequence container to check
	 * @param min the minimal number of sequences allowed
	 * @param max the maximal number of sequences allowed
	 * @return true if the constraint is met; false otherwise
	 * */
	bool validateSequenceNumber(const std::string & paramName
						, const RnaSequenceVec & sequences
						, const size_t min
						, const size_t max );

	/**
	 * Checks whether or not the sequences in the container conform to the
	 * allowed sequence alphabet.
	 *
	 * @param paramName the name of the parameter that enforces this constraint
	 *        (for exception handling)
	 * @param sequences the sequence container to check
	 * @return true if the constraint is met; false otherwise
	 * */
	bool validateSequenceAlphabet(const std::string & paramName
						, const RnaSequenceVec & sequences );

	/**
	 * Checks whether or not any command line argument were parsed. Throws a
	 * std::runtime_error if not.
	 */
	void checkIfParsed() const;

};

inline
void
CommandLineParsing::checkIfParsed() const {
	if (parsingCode == parsingCodeNotSet) {
		throw std::runtime_error("CommandLineParsing::getQuerySequences() : parse() function was not called yet");
	}
}


#endif /* COMMANDLINEPARSING_H_ */
