
#ifndef COMMANDLINEPARSING_H_
#define COMMANDLINEPARSING_H_

#include "general.h"
#include "RnaSequence.h"

#include <boost/regex.hpp>
#include <boost/program_options.hpp>

#include <iostream>
#include <cstdarg>

#include "Accessibility.h"
#include "InteractionEnergy.h"
#include "OutputHandler.h"
#include "Predictor.h"
#include "SeedConstraint.h"
#include "VrnaHandler.h"

/**
 * Central handler for all command line arguments etc.
 *
 */
class CommandLineParsing {
public:

	//! type for a list of sequences
	typedef std::vector< RnaSequence > RnaSequenceVec;
	//! type for a list of ranges for the sequences
	typedef std::vector<IndexRangeList> IndexRangeListVec;

	//! different exit codes for parsing
	enum ReturnCode {
		KEEP_GOING = -1,
		STOP_ALL_FINE = 0,
		STOP_PARSING_ERROR = 1,
		NOT_PARSED_YET = 999
	};


	//! regular expression used to parse range list string encodings
	static const boost::regex regexRangeEncoding;

public:

	/**
	 * Constructs a commandline argument parser for IntaRNA.
	 *
	 * @param logStream the stream to write validation log messages to
	 */
	CommandLineParsing();
	virtual ~CommandLineParsing();

	/**
	 * Parses the commandline arguments as passed to the 'main' method
	 * @param argc the number of arguments
	 * @param argv the argument array
	 *
	 * @return the parsing code: KEEP_GOING if all went fine, otherwise the exit code to return
	 *
	 */
	ReturnCode
	parse( int argc, char ** argv );

public:

	/////////  GETTERS  ///////////////////////////////////////////////////

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
	 * Access to the ranges to screen for interactions for the query with the
	 * according sequence number.
	 * @param sequenceNumber the number of the sequence within the vector
	 *         returned by getQuerySequences()
	 * @return the range list for the according sequence.
	 */
	const IndexRangeList& getQueryRanges( const size_t sequenceNumber ) const;

	/**
	 * Access to the ranges to screen for interactions for the target with the
	 * according sequence number.
	 * @param sequenceNumber the number of the sequence within the vector
	 *         returned by getTargetSequences()
	 * @return the range list for the according sequence.
	 */
	const IndexRangeList& getTargetRanges( const size_t sequenceNumber ) const;

	/**
	 * Returns a newly allocated Energy object according to the user defined
	 * parameters.
	 * @param accQuery the accessibility object of the query sequence
	 * @param accTarget the (reversed) accessibility object of the target sequence
	 * @return the newly allocated Energy object to be deleted by the calling
	 * function or NULL in error case
	 */
	InteractionEnergy* getEnergyHandler( const Accessibility& accQuery, const ReverseAccessibility& accTarget ) const;

	/**
	 * Provides a newly allocated output handler according to the user request.
	 *
	 * @param energy the energy handler used for interaction computation
	 *
	 * @return the newly allocated OutputHandler object to be deleted by the
	 * calling function
	 */
	OutputHandler* getOutputHandler(const InteractionEnergy & energy) const;

	/**
	 * Provides a newly allocated predictor according to the user defined
	 * parameters
	 * @param energy the interaction energy handler to be used
	 * @param output the output handler to be used
	 * @return the newly allocated Predictor object to be deleted by the calling
	 * function
	 */
	Predictor* getPredictor( const InteractionEnergy & energy
			, OutputHandler & output ) const;


	/**
	 * Provides the seed constraint according to the user settings
	 * @return the user defined seed constraints
	 */
	const SeedConstraint & getSeedConstraint() const;

	/**
	 * Access to the set folding temperature in Celsius.
	 * @return the chosen temperature in Celsisus
	 */
	T_type getTemperature() const;

	/**
	 * The constraints to be applied to the interaction output generation
	 * @return the output constraints to be applied
	 */
	OutputConstraint getOutputConstraint() const;


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
	boost::program_options::options_description opts_query;
	//! target specific options
	boost::program_options::options_description opts_target;
	//! seed specific options
	boost::program_options::options_description opts_seed;
	//! interaction/energy specific options
	boost::program_options::options_description opts_inter;
	//! general options
	boost::program_options::options_description opts_general;
	//! output options
	boost::program_options::options_description opts_output;

	//! overall option list
	boost::program_options::options_description opts_cmdline_all;

	//! central result code to be set by validate_* functions in error case
	ReturnCode parsingCode;

	//! the query command line argument
	std::string queryArg;
	//! the container holding all query sequences
	RnaSequenceVec query;
	//! accessibility computation mode for query sequences
	CharParameter qAcc;
	//! window length for query accessibility computation (plFold)
	NumberParameter<int> qAccW;
	//! maximal base pair span for query accessibility computation (plFold)
	NumberParameter<int> qAccL;
	//! constraint for accessibility computation for query sequences
	std::string qAccConstr;
	//! window length to be considered accessible/interacting within query
	NumberParameter<int> qIntLenMax;
	//! maximal internal loop length to be considered accessible/interacting within query
	NumberParameter<int> qIntLoopMax;
	//! the string encoding of the interaction intervals for the query(s)
	std::string qRegionString;
	//! the list of interaction intervals for each query sequence
	IndexRangeListVec qRegion;

	//! the target command line argument
	std::string targetArg;
	//! the container holding all target sequences
	RnaSequenceVec target;
	//! accessibility computation mode for target sequences
	CharParameter tAcc;
	//! window length for target accessibility computation (plFold)
	NumberParameter<int> tAccW;
	//! maximal base pair span for target accessibility computation (plFold)
	NumberParameter<int> tAccL;
	//! constraint for accessibility computation for target sequences
	std::string tAccConstr;
	//! window length to be considered accessible/interacting within target
	NumberParameter<int> tIntLenMax;
	//! maximal internal loop length to be considered accessible/interacting within target
	NumberParameter<int> tIntLoopMax;
	//! the string encoding of the interaction intervals for the target(s)
	std::string tRegionString;
	//! the list of interaction intervals for each target sequence
	IndexRangeListVec tRegion;

	//! whether or not a seed is to be required for an interaction or not
	bool noSeedRequired;
	//! number of base pairs in seed
	NumberParameter<int> seedBP;
	//! max overall unpaired in seed
	NumberParameter<int> seedMaxUP;
	//! max unpaired in query's seed
	NumberParameter<int> seedMaxUPq;
	//! max unpaired in target's seed
	NumberParameter<int> seedMaxUPt;
	//! max energy of a seed to be considered
	NumberParameter<E_type> seedMaxE;
	//! the final seed constraint to be used
	mutable SeedConstraint * seedConstraint;

	//! the temperature to be used for energy computations
	NumberParameter<T_type> temperature;

	//! the prediction mode
	NumberParameter<int> predMode;

	//! the selected energy model
	CharParameter energy;
	//! the provided energy parameter file of the VRNA package
	std::string energyFile;

	//! number of (sub)optimal interactions to report
	NumberParameter<int> oNumber;
	//! whether or not reported interactions can to be overlapping
	NumberParameter<int> oOverlap;

	//! the vienna energy parameter handler initialized by #parse()
	mutable VrnaHandler vrnaHandler;

protected:

	////////////  INDIVIDUAL TESTS  //////////////////

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
	 * Validates the query accessibility sliding window size argument.
	 * @param value the argument value to validate
	 */
	void validate_qAccW(const int & value);

	/**
	 * Validates the query accessibility maximal loop length argument.
	 * @param value the argument value to validate
	 */
	void validate_qAccL(const int & value);

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
	 * Validates the query's region argument.
	 * @param value the argument value to validate
	 */
	void validate_qRegion(const std::string & value);

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
	 * Validates the target accessibility sliding window size argument.
	 * @param value the argument value to validate
	 */
	void validate_tAccW(const int & value);

	/**
	 * Validates the target accessibility maximal loop length argument.
	 * @param value the argument value to validate
	 */
	void validate_tAccL(const int & value);

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
	 * Validates the target's region argument.
	 * @param value the argument value to validate
	 */
	void validate_tRegion(const std::string & value);

	/**
	 * Validates the seedBP argument.
	 * @param value the argument value to validate
	 */
	void validate_seedBP(const int & value);

	/**
	 * Validates the seedMaxUP argument.
	 * @param value the argument value to validate
	 */
	void validate_seedMaxUP(const int & value);

	/**
	 * Validates the seedMaxUPq argument.
	 * @param value the argument value to validate
	 */
	void validate_seedMaxUPq(const int & value);

	/**
	 * Validates the seedMaxUPt argument.
	 * @param value the argument value to validate
	 */
	void validate_seedMaxUPt(const int & value);

	/**
	 * Validates the seedMaxE argument.
	 * @param value the argument value to validate
	 */
	void validate_seedMaxE(const E_type & value);

	/**
	 * Validates the temperature argument.
	 * @param value the argument value to validate
	 */
	void validate_temperature(const T_type & value);

	/**
	 * Validates the prediction mode argument.
	 * @param value the argument value to validate
	 */
	void validate_predMode(const int & value);

	/**
	 * Validates the temperature argument.
	 * @param value the argument value to validate
	 */
	void validate_energy(const char & value);

	/**
	 * Validates the energy parameter file argument.
	 * @param value the argument value to validate
	 */
	void validate_energyFile(const std::string & value);

	/**
	 * Validates the oNumber argument.
	 * @param value the argument value to validate
	 */
	void validate_oNumber(const int & value);

	/**
	 * Validates the oOverlap argument.
	 * @param value the argument value to validate
	 */
	void validate_oOverlap(const int & value);


	////////////  GENERIC TESTS  /////////////////

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
			LOG(ERROR) <<name<<" = " <<value <<" : has to be in the range [" <<param.min <<","<<param.max<<"]";
			parsingCode = std::max(ReturnCode::STOP_PARSING_ERROR,parsingCode);
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
	 * Parses the parameter value and returns all parsed sequences.
	 * @param paramName the name of the parameter (for exception handling)
	 * @param paramArg the given argument for the parameter
	 * @param sequences the container to fill
	 */
	void parseSequences(const std::string & paramName,
					const std::string& paramArg,
					RnaSequenceVec& sequences );

	/**
	 * Parses the parameter input stream from FASTA format and returns all
	 * parsed sequences.
	 * @param paramName the name of the parameter (for exception handling)
	 * @param input the input stream from where to read the FASTA data
	 * @param sequences the container to fill
	 */
	void parseSequencesFasta( const std::string & paramName,
					std::istream& input,
					RnaSequenceVec& sequences);

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
	 * Checks whether or not a fil with the given file name exists and is
	 * readable.
	 * @param filename the file name to check
	 * @return true if the file can be accessed and is readable; false otherwise
	 */
	bool validateFile( const std::string & filename );

	/**
	 * Checks whether or not a range input is either a valid string encoding
	 * or an accessible file name.
	 * @param argName the name of the argument the @p value is for
	 * @param value the argument's value to be parsed
	 */
	bool
	validateRegion( const std::string & argName
					, const std::string & value );

	/**
	 * Parses a given range input and pushes the parsed ranges to the given
	 * container.
	 * @param argName the name of the argument the @p value is for
	 * @param value the argument's value to be parsed
	 * @param sequences the parsed sequences the regions are for
	 * @param rangeList the container to feed the parsed ranges to (individually for each file if read from BED input)
	 */
	void
	parseRegion( const std::string & argName
				, const std::string & value
				, const RnaSequenceVec & sequences
				, IndexRangeListVec & rangeList );

	/**
	 * Checks whether or not any command line argument were parsed. Throws a
	 * std::runtime_error if not.
	 */
	void checkIfParsed() const;

	/**
	 * Updates the overall parsing code according to the given one
	 * @param currentParsingCode the current parsing code to be used for the update
	 */
	void updateParsingCode( const ReturnCode currentParsingCode );

};

inline
void
CommandLineParsing::checkIfParsed() const {
	if (parsingCode == ReturnCode::NOT_PARSED_YET) {
		throw std::runtime_error("CommandLineParsing::checkIfParsed() : parse() function was not called yet");
	}
}


#endif /* COMMANDLINEPARSING_H_ */
