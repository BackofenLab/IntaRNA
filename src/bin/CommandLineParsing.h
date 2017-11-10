
#ifndef COMMANDLINEPARSING_H_
#define COMMANDLINEPARSING_H_

#include "IntaRNA/general.h"
#include "IntaRNA/RnaSequence.h"

#include <boost/regex.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

#include <iostream>
#include <cstdarg>

#include "IntaRNA/Accessibility.h"
#include "IntaRNA/InteractionEnergy.h"
#include "IntaRNA/OutputHandler.h"
#include "IntaRNA/Predictor.h"
#include "IntaRNA/SeedConstraint.h"
#include "IntaRNA/SeedHandler.h"
#include "IntaRNA/SeedHandlerExplicit.h"
#include "IntaRNA/VrnaHandler.h"

using namespace IntaRNA;

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
	 * @param accTarget the accessibility object of the target sequence
	 * @param accQuery the (reversed) accessibility object of the query sequence
	 * @return the newly allocated Energy object to be deleted by the calling
	 * function or NULL in error case
	 */
	InteractionEnergy* getEnergyHandler( const Accessibility& accTarget, const ReverseAccessibility& accQuery ) const;

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
	 * @param energy the interaction energy handler to be used
	 * @return the user defined seed constraints
	 */
	const SeedConstraint & getSeedConstraint( const InteractionEnergy & energy ) const;

	/**
	 * Provides a newly allocated seed handler object according to the user settings
	 *
	 * NOTE: the calling function has to remove the returned object!
	 *
	 * @param energy the interaction energy handler to be used
	 * @return a newly allocated seed handler respective the user defined seed constraints
	 */
	SeedHandler * getSeedHandler( const InteractionEnergy & energy ) const;

	/**
	 * Access to the set folding temperature in Celsius.
	 * @return the chosen temperature in Celsius
	 */
	T_type getTemperature() const;

	/**
	 * The constraints to be applied to the interaction output generation
	 * @return the output constraints to be applied
	 */
	OutputConstraint getOutputConstraint() const;

	/**
	 * The stream to write the interaction output to
	 * @return the output stream to write interaction output to
	 */
	std::ostream & getOutputStream() const;

	/**
	 * Writes the query accessibility to file/stream if requested
	 */
	void
	writeQueryAccessibility( const Accessibility & acc ) const;

	/**
	 * Writes the query accessibility to file/stream if requested
	 */
	void
	writeTargetAccessibility( const Accessibility & acc ) const;

#if INTARNA_MULITHREADING
	/**
	 * Number of threads to be used for parallel processing of
	 * query-target-combinations.
	 * @return number of threads to be used (>0)
	 */
	size_t
	getThreads() const;
#endif

protected:

	/////////  PRIVATE STUFF  ////////////////////////////////////////////////

	/**
	 * Defines codes for all registered prefixes for the --out argument
	 */
	enum OutPrefixCode {
		OP_EMPTY,
		OP_qMinE,
		OP_tMinE,
		OP_pMinE,
		OP_qAcc,
		OP_tAcc,
		OP_qPu,
		OP_tPu,
		OP_UNKNOWN
	};

	/**
	 * Provides the OutPrefixCode for a given prefix string or UNKNOWN if not
	 * registered. The mapping is done case-insensitive.
	 *
	 * @param outPrefix the prefix to encode
	 * @return the according OutPrefixCode or UNKNOWN if the prefix is not registered.
	 */
	static
	OutPrefixCode
	getCodeForOutPrefix( const std::string & outPrefix )
	{
		if (outPrefix.empty())
			return OutPrefixCode::OP_EMPTY;
		std::string prefLC = boost::to_lower_copy( outPrefix );
		if (prefLC == "qmine")	{ return OutPrefixCode::OP_qMinE; } else
		if (prefLC == "tmine")	{ return OutPrefixCode::OP_tMinE; } else
		if (prefLC == "pmine")	{ return OutPrefixCode::OP_pMinE; } else
		if (prefLC == "qacc")	{ return OutPrefixCode::OP_qAcc; } else
		if (prefLC == "tacc")	{ return OutPrefixCode::OP_tAcc; } else
		if (prefLC == "qpu")	{ return OutPrefixCode::OP_qPu; } else
		if (prefLC == "tpu")	{ return OutPrefixCode::OP_tPu; } else
		// not known
		return OutPrefixCode::OP_UNKNOWN;
	}

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

	//! whether or not STDIN was already requested by one of the following
	//! arguments
	bool stdinUsed;

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

	//! short option list
	boost::program_options::options_description opts_cmdline_short;

	//! central result code to be set by validate_* functions in error case
	ReturnCode parsingCode;

	//! the query command line argument
	std::string queryArg;
	//! the container holding all query sequences
	RnaSequenceVec query;
	//! subset of query sequence indices to be processed
	IndexRangeList qSet;
	//! string encoding of qSet
	std::string qSetString;
	//! accessibility computation mode for query sequences
	CharParameter qAcc;
	//! window length for query accessibility computation (plFold)
	NumberParameter<int> qAccW;
	//! maximal base pair span for query accessibility computation (plFold)
	NumberParameter<int> qAccL;
	//! constraint for accessibility computation for query sequences
	std::string qAccConstr;
	//! the file/stream to read the query's accessibility data from
	std::string qAccFile;
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
	//! subset of target sequence indices to be processed
	IndexRangeList tSet;
	//! string encoding of tSet
	std::string tSetString;
	//! accessibility computation mode for target sequences
	CharParameter tAcc;
	//! window length for target accessibility computation (plFold)
	NumberParameter<int> tAccW;
	//! maximal base pair span for target accessibility computation (plFold)
	NumberParameter<int> tAccL;
	//! constraint for accessibility computation for target sequences
	std::string tAccConstr;
	//! the file/stream to read the query's accessibility data from
	std::string tAccFile;
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
	//! explicit seed encodings (optional)
	std::string seedTQ;
	//! number of base pairs in seed
	NumberParameter<int> seedBP;
	//! max overall unpaired in seed
	NumberParameter<int> seedMaxUP;
	//! max unpaired in query's seed
	NumberParameter<int> seedQMaxUP;
	//! max unpaired in target's seed
	NumberParameter<int> seedTMaxUP;
	//! max energy of a seed to be considered
	NumberParameter<E_type> seedMaxE;
	//! minimal unpaired probability (per sequence) of a seed to be considered
	NumberParameter<E_type> seedMinPu;
	//! intervals in query for seed search
	std::string seedQRange;
	//! intervals in target for seed search
	std::string seedTRange;
	//! the final seed constraint to be used
	mutable SeedConstraint * seedConstraint;

	//! the temperature to be used for energy computations
	NumberParameter<T_type> temperature;

	//! the prediction target (mfe-single-site, max-prob-site, ..)
	CharParameter pred;
	//! the prediction mode (heuristic, space-efficient, exact)
	CharParameter predMode;
#if INTARNA_MULITHREADING
	//! number of threads = number of parallel predictors running
	NumberParameter<int> threads;
#endif

	//! the selected energy model
	CharParameter energy;
	//! the provided energy parameter file of the VRNA package
	std::string energyFile;

	//! where to write the output to and for each in what format
//	std::string out;
	std::vector< std::string > out;
	//! provides the parsed stream name for each out prefix
	std::map<OutPrefixCode,std::string> outPrefix2streamName;
	//! output stream
	std::ostream * outStream;
	//! output mode
	CharParameter outMode;
	//! number of (sub)optimal interactions to report
	NumberParameter<int> outNumber;
	//! whether or not reported interactions can to be overlapping
	CharParameter outOverlap;
	//! deltaE to mfe allowed to report an interaction
	NumberParameter<double> outDeltaE;
	//! max E allowed to report an interaction
	NumberParameter<double> outMaxE;
	//! the CSV column selection
	std::string outCsvCols;
	//! the CSV column selection
	static const std::string outCsvCols_default;

	//! (optional) file name for log output
	std::string logFileName;

	//! the vienna energy parameter handler initialized by #parse()
	mutable VrnaHandler vrnaHandler;

protected:

	/**
	 * sets the stdinUsed member to true if so far false or raises an exception
	 * if it is already true.
	 * @return true if stdinUsed was false so far; false otherwise (error logged)
	 */
	bool setStdinUsed();

	////////////  INDIVIDUAL TESTS  //////////////////

	/**
	 * Validates the query sequence argument.
	 * @param value the argument value to validate
	 */
	void validate_query(const std::string & value);

	/**
	 * Validates the query's qSet argument.
	 * @param value the argument value to validate
	 */
	void validate_qSet(const std::string & value);

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
	 * Validates the qAccFile argument.
	 * @param value the argument value to validate
	 */
	void validate_qAccFile(const std::string & value);

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
	 * Validates the target's tSet argument.
	 * @param value the argument value to validate
	 */
	void validate_tSet(const std::string & value);

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
	 * Validates the tAccFile argument.
	 * @param value the argument value to validate
	 */
	void validate_tAccFile(const std::string & value);

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
	 * Validates the explicit seed argument.
	 * @param value the argument value to validate
	 */
	void validate_seedTQ(const std::string & value);

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
	 * Validates the seedQMaxUP argument.
	 * @param value the argument value to validate
	 */
	void validate_seedQMaxUP(const int & value);

	/**
	 * Validates the seedTMaxUP argument.
	 * @param value the argument value to validate
	 */
	void validate_seedTMaxUP(const int & value);

	/**
	 * Validates the seedMaxE argument.
	 * @param value the argument value to validate
	 */
	void validate_seedMaxE(const E_type & value);

	/**
	 * Validates the seedMinPu argument.
	 * @param value the argument value to validate
	 */
	void validate_seedMinPu(const E_type & value);

	/**
	 * Validates the seedQRange argument.
	 * @param value the argument value to validate
	 */
	void validate_seedQRange(const std::string & value);

	/**
	 * Validates the seedTRange argument.
	 * @param value the argument value to validate
	 */
	void validate_seedTRange(const std::string & value);

	/**
	 * Validates the temperature argument.
	 * @param value the argument value to validate
	 */
	void validate_temperature(const T_type & value);

	/**
	 * Validates the prediction target argument.
	 * @param value the argument value to validate
	 */
	void validate_pred(const char & value);

	/**
	 * Validates the prediction mode argument.
	 * @param value the argument value to validate
	 */
	void validate_predMode(const char & value);

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
	 * Validates the out argument.
	 *
	 * Furthermore, stores all valid out mappings in outprefix2streamName
	 *
	 * @param list the list of argument values to validate
	 */
	void validate_out(const std::vector<std::string> & list);

	/**
	 * Validates the outMode argument.
	 * @param value the argument value to validate
	 */
	void validate_outMode(const char & value);

	/**
	 * Validates the outNumber argument.
	 * @param value the argument value to validate
	 */
	void validate_outNumber(const int & value);

	/**
	 * Validates the outOverlap argument.
	 * @param value the argument value to validate
	 */
	void validate_outOverlap(const char & value);

	/**
	 * Validates the outDeltaE argument.
	 * @param value the argument value to validate
	 */
	void validate_outDeltaE(const double & value);

	/**
	 * Validates the outMaxE argument.
	 * @param value the argument value to validate
	 */
	void validate_outMaxE(const double & value);

	/**
	 * Validates the outCsvCols argument.
	 * @param value the argument value to validate
	 */
	void validate_outCsvCols(const std::string & value);

#if INTARNA_MULITHREADING
	/**
	 * Validates the threads argument.
	 * @param value the argument value to validate
	 */
	void validate_threads( const int & value);
#endif

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
	 * Validates the string encoding of an IndexRangeList encoding parameter
	 * @param argName the name of the parameter (for exception handling)
	 * @param value the value of the parameter to validate
	 * @param indexMin the minimal value of an index allowed
	 * @param indexMax the maximal value of an index allowed
	 */
	void validate_indexRangeList(const std::string & argName
								, const std::string & value
								, const size_t indexMin
								, const size_t indexMax );

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
	 * @param seqSubset the indices of the input sequences to store (all other
	 *                  ignored)
	 */
	void parseSequences(const std::string & paramName,
					const std::string& paramArg,
					RnaSequenceVec& sequences,
					const IndexRangeList & seqSubset );

	/**
	 * Parses the parameter input stream from FASTA format and returns all
	 * parsed sequences.
	 * @param paramName the name of the parameter (for exception handling)
	 * @param input the input stream from where to read the FASTA data
	 * @param sequences the container to fill
	 * @param seqSubset the indices of the input sequences to store (all other
	 *                  ignored)
	 */
	void parseSequencesFasta( const std::string & paramName,
					std::istream& input,
					RnaSequenceVec& sequences,
					const IndexRangeList & seqSubset );

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
	 * Validates whether or not a given output target value is either STDOUT,
	 * STDERR, or non-empty (= file name)
	 */
	void
	validate_outputTarget( const std::string &argName , const std::string & value);

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

	/**
	 * Triggers an initial output of the output handler if needed.
	 */
	void initOutputHandler();

	/**
	 * Writes the accessibility to file or stream if requested by the user
	 * @param acc the accessibility data assigned
	 * @param fileOrStream the name of file/stream to write to
	 * @param (true) writes ED values, (false) writes Pu values
	 */
	void writeAccessibility( const Accessibility& acc, const std::string & fileOrStream, const bool writeED ) const;

	/**
	 * Adds a generic file prefix for input/output files for the given query
	 * and/or target sequence. Empty strings as well as STDOUT/STDERR are
	 * ignored, i.e. directly returned.
	 *
	 * @param fileName the file name that might have to be extended
	 * @param target the target sequence or NULL if no target is involved in naming
	 * @param query the query sequence or NULL if no query is involved in naming
	 *
	 * @return the prefixed file name to be used
	 */
	std::string
	getFullFilename( const std::string & fileName
					, const RnaSequence * target
					, const RnaSequence * query ) const;

};



////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


inline
bool
CommandLineParsing::
setStdinUsed() {
	if (stdinUsed) {
		LOG(ERROR) <<"STDIN named more than once as input source, which is not supported";
		updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
		return false;
	}
	stdinUsed = true;
	return true;
}

////////////////////////////////////////////////////////////////////////////

inline
void
CommandLineParsing::checkIfParsed() const {
	if (parsingCode == ReturnCode::NOT_PARSED_YET) {
		throw std::runtime_error("CommandLineParsing::checkIfParsed() : parse() function was not called yet");
	}
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_query(const std::string & value)
{
	validate_sequenceArgument("query",value);
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_qSet(const std::string & value) {
	// clear current qSet data
	qSet.clear();
	// parse input value
	if (!value.empty()) {
		// check regex
		if (!boost::regex_match(value, IndexRangeList::regex, boost::match_perl) ) {
			LOG(ERROR) <<"qSet"<<" = " <<value <<" : is not in the format 'from1-to1,from2-to2,..'";
			parsingCode = std::max(ReturnCode::STOP_PARSING_ERROR,parsingCode);
		} else {
			// parse and store subset definitions
			qSet = IndexRangeList(value);
		}
	}
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_qIntLenMax(const int & value)
{
	// forward check to general method
	validate_numberArgument("qIntW", qIntLenMax, value);
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_qIntLoopMax(const int & value)
{
	// forward check to general method
	validate_numberArgument("qIntL", qIntLoopMax, value);
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_qAcc(const char & value)
{
	// forward check to general method
	validate_charArgument("qAcc", qAcc, value);
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_qAccW(const int & value)
{
	// forward check to general method
	validate_numberArgument("qAccW", qAccW, value);
	// check lower bound
	if (qAccW.val > 0 && qAccW.val < 3) {
		LOG(ERROR) <<"\n qAccW = " <<value <<" : has to be 0 or > 3";
		updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
	}
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_qAccL(const int & value)
{
	// forward check to general method
	validate_numberArgument("qAccL", qAccL, value);
	// check lower bound
	if (qAccL.val > 0 && qAccL.val < 3) {
		LOG(ERROR) <<"qAccL = " <<value <<" : has to be 0 or > 3";
		updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
	}
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_qAccConstr(const std::string & value)
{
	// forward check to general method
	validate_structureConstraintArgument("qAccConstr", value);
	// check if no sliding window computation requested
	if (qAccW.val > 0 || qAccL.val > 0) {
		LOG(ERROR) <<"query accessibility constraint not possible for sliding window computation (qAccL/W > 0)";
		updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
	}
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_qAccFile(const std::string & value)
{
	// if not empty
	if (!value.empty()) {
		// if not STDIN
		if ( boost::iequals(value,"STDIN") ) {
			if (getTargetSequences().size()>1) {
				LOG(ERROR) <<"reading quary accessibilities for multiple sequences from '"<<value<<"' is not supported";
				updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
			} else {
				setStdinUsed();
			}
		} else {
			// should be file
			for (auto rna = getQuerySequences().begin(); rna != getQuerySequences().end(); rna++) {
				if ( ! validateFile( getFullFilename(value, NULL, &(*rna) ) ) ) {
					LOG(ERROR) <<"query accessibility file '"<<value<<"' could not be found";
					updateParsingCode(ReturnCode::STOP_PARSING_ERROR);

				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_qRegion(const std::string & value) {
	// check and store region information
	validateRegion( "qRegion", value );
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_target(const std::string & value)
{
	validate_sequenceArgument("target",value);
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_tSet(const std::string & value) {
	// clear current qSet data
	tSet.clear();
	// parse input value
	if (!value.empty()) {
		// check regex
		if (!boost::regex_match(value, IndexRangeList::regex, boost::match_perl) ) {
			LOG(ERROR) <<"tSet"<<" = " <<value <<" : is not in the format 'from1-to1,from2-to2,..'";
			parsingCode = std::max(ReturnCode::STOP_PARSING_ERROR,parsingCode);
		} else {
			// parse and store subset definitions
			tSet = IndexRangeList(value);
		}
	}
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_tIntLenMax(const int & value)
{
	// forward check to general method
	validate_numberArgument("tIntLenMax", tIntLenMax, value);
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_tIntLoopMax(const int & value)
{
	// forward check to general method
	validate_numberArgument("tIntLoopMax", tIntLoopMax, value);
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_tAcc(const char & value)
{
	// forward check to general method
	validate_charArgument("tAcc", tAcc, value);
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_tAccW(const int & value)
{
	// forward check to general method
	validate_numberArgument("tAccW", tAccW, value);
	// check lower bound
	if (tAccW.val > 0 && tAccW.val < 3) {
		LOG(ERROR) <<"tAccW = " <<value <<" : has to be 0 or > 3";
		updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
	}
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_tAccL(const int & value)
{
	// forward check to general method
	validate_numberArgument("tAccL", tAccL, value);
	// check lower bound
	if (tAccL.val > 0 && tAccL.val < 3) {
		LOG(ERROR) <<"tAccL = " <<value <<" : has to be 0 or > 3";
		updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
	}
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_tAccConstr(const std::string & value)
{
	// forward check to general method
	validate_structureConstraintArgument("tAccConstr", value);
	// check if no sliding window computation requested
	if (tAccW.val > 0 || tAccL.val > 0) {
		LOG(ERROR) <<"query accessibility constraint not possible for sliding window computation (tAccL/W > 0)";
		updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
	}
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_tAccFile(const std::string & value)
{
	// if not empty
	if (!value.empty()) {
		// if not STDIN
		if ( boost::iequals(value,"STDIN") ) {
			if (getTargetSequences().size()>1) {
				LOG(ERROR) <<"reading target accessibilities for multiple sequences from '"<<value<<"' is not supported";
				updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
			} else {
				setStdinUsed();
			}
		} else {
			// should be file
			for (auto rna = getTargetSequences().begin(); rna != getTargetSequences().end(); rna++) {
				if ( ! validateFile( getFullFilename(value, &(*rna), NULL ) ) ) {
					LOG(ERROR) <<"target accessibility file '"<<value<<"' could not be found";
					updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_tRegion(const std::string & value) {
	// check and store region information
	validateRegion( "tRegion", value );
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_seedTQ(const std::string & value) {
	if (!value.empty()) {
		// split by commas
		size_t p2 = 0;
		size_t p1 = value.find_first_not_of(',',p2);
		size_t non_empty_seeds = 0;
		// do for all substrings
		while( p1 != std::string::npos ) {
			// find end of current substring
			p2 = value.find_first_of(',', p1 + 1);
			const size_t length = (p2 == std::string::npos ? value.size() : p2) - p1;
			// check substring
			std::string errMsg = SeedHandlerExplicit::checkSeedEncoding( value.substr(p1, length));
			non_empty_seeds++;
			if (!errMsg.empty()) {
				LOG(ERROR) <<"explicit seed encoding '" <<value.substr(p1, length) <<"' : "<<errMsg;
				parsingCode = std::max(ReturnCode::STOP_PARSING_ERROR,parsingCode);
			}
			// update p1
			p1 = value.find_first_not_of(',',p2);
		}
		if (non_empty_seeds == 0) {
			LOG(ERROR) <<"explicit seed encoding '" <<value <<"' does not contain any seed information";
			parsingCode = std::max(ReturnCode::STOP_PARSING_ERROR,parsingCode);
		}
	}
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_seedBP(const int & value) {
	// forward check to general method
	validate_numberArgument("seedBP", seedBP, value);
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_seedMaxUP(const int & value) {
	// forward check to general method
	validate_numberArgument("seedMaxUP", seedMaxUP, value);
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_seedQMaxUP(const int & value) {
	// forward check to general method
	validate_numberArgument("seedQMaxUP", seedQMaxUP, value);
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_seedTMaxUP(const int & value) {
	// forward check to general method
	validate_numberArgument("seedTMaxUP", seedTMaxUP, value);
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_seedMaxE(const E_type & value) {
	// forward check to general method
	validate_numberArgument("seedMaxE", seedMaxE, value);
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_seedMinPu(const E_type & value) {
	// forward check to general method
	validate_numberArgument("seedMinPu", seedMinPu, value);
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_seedQRange(const std::string & value) {
	if (!value.empty()) {
		// check regex
		if (!boost::regex_match(value, IndexRangeList::regex, boost::match_perl) ) {
			LOG(ERROR) <<"seedQRange"<<" = " <<value <<" : is not in the format 'from1-to1,from2-to2,..'";
			parsingCode = std::max(ReturnCode::STOP_PARSING_ERROR,parsingCode);
		}
	}
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_seedTRange(const std::string & value) {
	if (!value.empty()) {
		// check regex
		if (!boost::regex_match(value, IndexRangeList::regex, boost::match_perl) ) {
			LOG(ERROR) <<"seedTRange"<<" = " <<value <<" : is not in the format 'from1-to1,from2-to2,..'";
			parsingCode = std::max(ReturnCode::STOP_PARSING_ERROR,parsingCode);
		}
	}
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_temperature(const T_type & value) {
	// forward check to general method
	validate_numberArgument("temperature", temperature, value);
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_pred(const char & value)
{
	// forward check to general method
	validate_charArgument("mode", pred, value);
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_predMode(const char & value)
{
	// forward check to general method
	validate_charArgument("mode", predMode, value);
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_energy(const char & value)
{
	// forward check to general method
	validate_charArgument("energy", energy, value);
}

////////////////////////////////////////////////////////////////////////////

inline
void
CommandLineParsing::
validate_energyFile(const std::string & value)
{
	// check if file exists and is readable
	if (!validateFile( value )) {
		LOG(ERROR) <<"provided VRNA energy parameter file '" <<value <<"' could not be processed.";
		updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
	}
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_outMode(const char & value) {
	// forward check to general method
	validate_charArgument("outMode", outMode, value);
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_outputTarget(
		const std::string &argName
		, const std::string & value)
{
	// check if standard stream
	if (boost::iequals(value,"STDOUT") || boost::iequals(value,"STDERR")) {
		return;
	}
	// check if empty filename
	if ( value.empty() ) {
		LOG(ERROR) <<argName<<" : no output defined (empty)";
		updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
	}
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_out(const std::vector<std::string> & list) {

	// check for uniqueness of argument prefixes
	bool emptyPrefixSeen = false;
	for (auto v=list.begin(); parsingCode != ReturnCode::STOP_PARSING_ERROR && v!=list.end(); v++) {
		// get prefix
		std::string curPref = (v->find(':')==std::string::npos ? "" : v->substr(0,v->find(':')));
		// get code for prefix
		OutPrefixCode curPrefCode = getCodeForOutPrefix( curPref );

		switch (curPrefCode) {
		// handle unknown prefix
		case OutPrefixCode::OP_UNKNOWN : {
			LOG(ERROR) <<"--out : prefix '"<<curPref<<"' is not supported.. maybe misspelled?";
			updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
			break;
		}
		// handle empty prefix
		case OutPrefixCode::OP_EMPTY : {
			// check if already seen
			if (emptyPrefixSeen) {
				LOG(ERROR) <<"--out : specified more than once without prefix";
				updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
				break;
			}
			// keep track that empty prefix was already seen
			emptyPrefixSeen = true;
			// cut of leading ':' if needed
			std::string streamName = (v->find(':')==std::string::npos ? *v : v->substr(v->find(':')+1));
			// forward check to general method
			validate_outputTarget( "--out", streamName );
			// store stream name
			outPrefix2streamName[curPrefCode] = streamName;
			// proceed to next argument
			break;
		}
		// handle all other prefixes
		default : {
			// check if prefix was already seen
			if ( !outPrefix2streamName.at(curPrefCode).empty()) {
				LOG(ERROR) <<"--out : specified more than once with prefix '"<<curPref<<"'";
				updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
				break;
			}
			// store prefix to identify another existence
			std::string streamName = v->substr(v->find(':')+1);
			// forward check to general method
			validate_outputTarget( "--out="+curPref+":", streamName );
			// store stream name
			outPrefix2streamName[curPrefCode] = streamName;
			break;
		}
		} // switch curPrefCode
	} // for all arguments
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_outNumber(const int & value) {
	// forward check to general method
	validate_numberArgument("oNumber", outNumber, value);
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_outOverlap(const char & value) {
	// forward check to general method
	validate_charArgument("outOverlap", outOverlap, value);
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_outDeltaE(const double & value) {
	// forward check to general method
	validate_numberArgument("outDeltaE", outDeltaE, value);
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_outMaxE(const double & value) {
	// forward check to general method
	validate_numberArgument("outMaxE", outMaxE, value);
}

////////////////////////////////////////////////////////////////////////////

#if INTARNA_MULITHREADING
inline
void CommandLineParsing::validate_threads(const int & value)
{
	// forward check to general method
	validate_numberArgument("threads", threads, value);
}
#endif

////////////////////////////////////////////////////////////////////////////

inline
void
CommandLineParsing::
writeQueryAccessibility( const Accessibility & acc ) const
{
	// forward to generic function
	if (!outPrefix2streamName.at(OutPrefixCode::OP_qAcc).empty()) {
		VLOG(2) <<"writing ED values for query '"<<acc.getSequence().getId()<<"' to "<<outPrefix2streamName.at(OutPrefixCode::OP_qAcc);
		writeAccessibility( acc
				// get file name prefixed with sequence number if needed
				, getFullFilename(outPrefix2streamName.at(OutPrefixCode::OP_qAcc), NULL, &(acc.getSequence()))
				, true );
	}
	if (!outPrefix2streamName.at(OutPrefixCode::OP_qPu).empty()) {
		VLOG(2) <<"writing unpaired probabilities for query '"<<acc.getSequence().getId()<<"' to "<<outPrefix2streamName.at(OutPrefixCode::OP_qPu);
		writeAccessibility( acc
				// get file name prefixed with sequence number if needed
				, getFullFilename(outPrefix2streamName.at(OutPrefixCode::OP_qPu), NULL, &(acc.getSequence()))
				, false );
	}
}

////////////////////////////////////////////////////////////////////////////

inline
void
CommandLineParsing::
writeTargetAccessibility( const Accessibility & acc ) const
{
	// forward to generic function
	if (!outPrefix2streamName.at(OutPrefixCode::OP_tAcc).empty()) {
		VLOG(2) <<"writing ED values for target '"<<acc.getSequence().getId()<<"' to "<<outPrefix2streamName.at(OutPrefixCode::OP_tAcc);
		writeAccessibility( acc
				// get file name prefixed with sequence number if needed
				, getFullFilename(outPrefix2streamName.at(OutPrefixCode::OP_tAcc), &(acc.getSequence()), NULL)
				, true );
	}
	if (!outPrefix2streamName.at(OutPrefixCode::OP_tPu).empty()) {
		VLOG(2) <<"writing unpaired probabilities for target '"<<acc.getSequence().getId()<<"' to "<<outPrefix2streamName.at(OutPrefixCode::OP_tPu);
		writeAccessibility( acc
				// get file name prefixed with sequence number if needed
				, getFullFilename(outPrefix2streamName.at(OutPrefixCode::OP_tPu), &(acc.getSequence()), NULL)
				, false );
	}
}

////////////////////////////////////////////////////////////////////////////

#if INTARNA_MULITHREADING
inline
size_t
CommandLineParsing::
getThreads() const
{
	return threads.val;
}
#endif

////////////////////////////////////////////////////////////////////////////

inline
std::string
CommandLineParsing::
getFullFilename( const std::string & fileNamePath, const RnaSequence * target, const RnaSequence * query ) const
{
	// do nothing for empty file names
	if (fileNamePath.empty()) {
		return fileNamePath;
	}
	// exclude stream names from prefixing
	if (boost::iequals(fileNamePath,"STDOUT") || boost::iequals(fileNamePath,"STDERR")) {
		return fileNamePath;
	}

	// generate file ID if necessary
	std::string fileID = "";

	// generate target only
	if (target != NULL && query == NULL) {
		if (getTargetSequences().size() > 1) {
			fileID += "s";
//			prefix += "t";
			// search for index of the target sequence
			for (size_t t = 0; t < getTargetSequences().size(); t++) {
				if (getTargetSequences().at(t) == *target) {
					// indexing starts with 1
					fileID += toString(t+1);
					break;
				}
			}
		}
	} else
	// generate query only
	if (query != NULL && target == NULL) {
		if (getQuerySequences().size() > 1) {
			fileID += "s";
//			prefix += "q";
			// search for index of the query sequence
			for (size_t q = 0; q < getQuerySequences().size(); q++) {
				if (getQuerySequences().at(q) == *query) {
					// indexing starts with 1
					fileID += toString(q+1);
					break;
				}
			}
		}
	} else
	// generate combined part
	{
		if (getQuerySequences().size() > 1 || getTargetSequences().size() > 1) {
			fileID += "t";
			// search for index of the target sequence
			for (size_t t = 0; t < getTargetSequences().size(); t++) {
				if (getTargetSequences().at(t) == *target) {
					// indexing starts with 1
					fileID += toString(t+1);
					break;
				}
			}
			fileID += "q";
			// search for index of the query sequence
			for (size_t q = 0; q < getQuerySequences().size(); q++) {
				if (getQuerySequences().at(q) == *query) {
					// indexing starts with 1
					fileID += toString(q+1);
					break;
				}
			}
		}
	}

	if (fileID.empty()) {
		return fileNamePath;
	} else {

		// get position in fileNamePath where the
		const size_t startOfExtension = fileNamePath.size() - boost::filesystem::path(fileNamePath).extension().string().size();

		// return compiled file name including file ID
		return fileNamePath.substr(0,startOfExtension)
				+ "-" + fileID
				+ fileNamePath.substr(startOfExtension);
	}
}


////////////////////////////////////////////////////////////////////////////




#endif /* COMMANDLINEPARSING_H_ */
