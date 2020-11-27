
#ifndef COMMANDLINEPARSING_H_
#define COMMANDLINEPARSING_H_

#include "IntaRNA/general.h"
#include "IntaRNA/RnaSequence.h"

#include <boost/regex.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/spirit/include/qi.hpp>

#include <iostream>
#include <cstdarg>

#include "IntaRNA/Accessibility.h"
#include "IntaRNA/InteractionEnergy.h"
#include "IntaRNA/HelixConstraint.h"
#include "IntaRNA/HelixHandler.h"
#include "IntaRNA/OutputStreamHandler.h"
#include "IntaRNA/OutputHandler.h"
#include "IntaRNA/Predictor.h"
#include "IntaRNA/SeedConstraint.h"
#include "IntaRNA/SeedHandler.h"
#include "IntaRNA/SeedHandlerExplicit.h"
#include "IntaRNA/PredictionTrackerSpotProb.h"
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

	// TODO if extended: also extend both getPersonality() and getPersonalityName()
	enum Personality {
		IntaRNA,		// default
		IntaRNA1,		// IntaRNA v1 like setup
		IntaRNA2,		// IntaRNA v2 like setup
		IntaRNA3,		// default IntaRNA v3 setup
		IntaRNAens,		// ensemble-based prediction
		IntaRNAsTar,	// sRNA-target prediction (optimized parameter)
		IntaRNAseed,  	// seed-only predictions
		IntaRNAhelix,  	// helix-block-based predictions
		IntaRNAduplex,	// RNAhybrid/RNAduplex-like
		IntaRNAexact	// RNAup-like exact predictions
	};


	/**
	 * Identifies the requested personality from the given call arguments
	 */
	static
	std::string
	getPersonalityName( Personality p ) {
		switch(p) {
		// the following case list is parsed by ROOT/Makefile.am to generate
		// respective binary links for each personality
		// thus: keep the one-line format!
		case IntaRNA : return "IntaRNA";
		case IntaRNA1 : return "IntaRNA1";
		case IntaRNA2 : return "IntaRNA2";
		case IntaRNA3 : return "IntaRNA3";
		case IntaRNAens : return "IntaRNAens";
		case IntaRNAsTar : return "IntaRNAsTar";
		case IntaRNAseed : return "IntaRNAseed";
		case IntaRNAhelix : return "IntaRNAhelix";
		case IntaRNAduplex : return "IntaRNAduplex";
		case IntaRNAexact : return "IntaRNAexact";
		default : return "unknown";
		}
	}

	template <typename InputIterator = std::string::const_iterator>
	struct unescaped_string
			: boost::spirit::qi::grammar<InputIterator, std::string(char const*)>
			{
				boost::spirit::qi::rule<InputIterator, std::string(char const*)> unesc_str;
				boost::spirit::qi::symbols<char const, char const> unesc_char;

				unescaped_string()
	    		  : unescaped_string::base_type(unesc_str)
				{
					using namespace boost::spirit;
					unesc_char.add(R"(\a)", '\a')(R"(\b)", '\b')(R"(\f)", '\f')(R"(\n)", '\n')
	                    		  (R"(\r)", '\r')(R"(\t)", '\t')(R"(\v)", '\v')
	                    		  (R"(\\)", '\\')(R"(\')", '\'')(R"(\")", '\"')
	                    		  ;

					unesc_str = qi::lit(qi::_r1)
								>> *(unesc_char | qi::graph | qi::space | R"(\x)" >> qi::hex)
								>>  qi::lit(qi::_r1)
								;
				}

				static
				std::string
				getUnescaped( const std::string & rawString )
				{
					using namespace boost::spirit;
					std::string parsed;
					char const* quote = ""; // leading and trailing quote string to parse
					unescaped_string myGrammar;
					qi::parse( rawString.begin(), rawString.end(), myGrammar(quote), parsed);
					return parsed;
				}
			};

	/**
	 * Identifies the requested personality from the given call arguments
	 */
	static
	Personality
	getPersonality( int argc, char** argv );

public:

	/**
	 * Constructs a commandline argument parser for IntaRNA.
	 *
	 * @param personality the Personality for which to provide and parse parameters
	 */
	CommandLineParsing( const Personality personality );
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
	 * @param energy the energy handler later used for prediction (needed for region postprocessing)
	 * @param sequenceNumber the number of the sequence within the vector
	 *         returned by getQuerySequences()
	 * @param acc the accessibility object for the sequence of the given number
	 * @return the range list for the according sequence.
	 */
	const IndexRangeList& getQueryRanges( const InteractionEnergy & energy, const size_t sequenceNumber, const Accessibility & acc ) const;

	/**
	 * Access to the ranges to screen for interactions for the target with the
	 * according sequence number.
	 * @param energy the energy handler later used for prediction (needed for region postprocessing)
	 * @param sequenceNumber the number of the sequence within the vector
	 *         returned by getTargetSequences()
	 * @param acc the accessibility object for the sequence of the given number
	 * @return the range list for the according sequence.
	 */
	const IndexRangeList& getTargetRanges( const InteractionEnergy & energy, const size_t sequenceNumber, const Accessibility & acc ) const;

	/**
	 * Access to the maximal window width of a query/target sequence range to
	 * be used for prediction using overlapping windows to save memory.
	 *
	 * @return the maximal window width to be used for windows-based
	 *         computations
	 */
	const size_t getWindowWidth() const;

	/**
	 * Access to the window overlap of a query/target sequence range to
	 * be used for prediction using overlapping windows to save memory.
	 *
	 * @return the window overlap to be used for window-based computations
	 */
	const size_t getWindowOverlap() const;

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
	 * Furthermore, it checks whether Zall has to be computed to generate the
	 * output.
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
	const HelixConstraint & getHelixConstraint( const InteractionEnergy & energy ) const;

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
	Z_type getTemperature() const;

	/**
	 * The constraints to be applied to the interaction output generation
	 * @param energy the interaction energy to be used for computation
	 * @return the output constraints to be applied
	 */
	OutputConstraint getOutputConstraint( const InteractionEnergy & energy ) const;

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

	/**
	 * Whether or not output is to be written for each region combination
	 * @return true if output is to be written for each region combination;
	 *         false otherwise (best for each query-target combination only)
	 */
	bool
	reportBestPerRegion() const;

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
		OP_qSpotProb,
		OP_tSpotProb,
		OP_pMinE,
		OP_qAcc,
		OP_tAcc,
		OP_qPu,
		OP_tPu,
		OP_spotProb,
		OP_spotProbAll,
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
		if (prefLC == "qspotprob")	{ return OutPrefixCode::OP_qSpotProb; } else
		if (prefLC == "tspotprob")	{ return OutPrefixCode::OP_tSpotProb; } else
		if (prefLC == "pmine")	{ return OutPrefixCode::OP_pMinE; } else
		if (prefLC == "qacc")	{ return OutPrefixCode::OP_qAcc; } else
		if (prefLC == "tacc")	{ return OutPrefixCode::OP_tAcc; } else
		if (prefLC == "qpu")	{ return OutPrefixCode::OP_qPu; } else
		if (prefLC == "tpu")	{ return OutPrefixCode::OP_tPu; } else
		if (prefLC == "spotprob")	{ return OutPrefixCode::OP_spotProb; } else
		// not known
		return OutPrefixCode::OP_UNKNOWN;
	}

	/**
	 * Limits and values for a number parameter
	 */
	template <typename T>
	class NumberParameter {
	public:
		//! name of the parameter
	  std::string name;
		//! the value of the parameter
	  T val;
		//! the minimally allowed value
	  const T min;
		//! the maximally allowed value
	  const T max;
		//! the default value
	  T def;
		/**
		 * construction feeding the members
		 * @param name the name of the parameter within the CLI
		 * @param min the minimally allowed value
		 * @param max the maximally allowed value
		 * @param def the default value
		 */
	  NumberParameter( const std::string& name, const T min, const T max, const T def )
	   : name(name), val(def), min(min), max(max), def(def)
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
		//! whether or not the value was set
		//! @return val != def
		bool isSet() const {
			return val != def;
		}

	};

	/**
	 * Allowed alphabet for a single char parameter
	 */
	class CharParameter {
	public:
		  //! the name of the parameter in the CLI
		std::string name;
		  //! the value of the parameter
		char val;
		  //! the set of allowed values for this parameter as a string
		const std::string alphabet;
		  //! the default value of the parameter
		char def;
		/**
		 * Construction and member setup
		 * @param name name of the parameter in the CLI
		 * @param alphabet the allowed set of character values
		 * @param def the default value (has to be part of the alphabet)
		 */
		CharParameter( const std::string& name, const std::string& alphabet, const char def )
		  : name(name), val(def), alphabet(alphabet), def(def)
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
		//! whether or not the value was set
		//! @return val != def
		bool isSet() const {
			return val != def;
		}

	};

	//! what is the requested personality for which we parse the parameters
	Personality personality;
	//! might hold the personality string after parsing if given via parameter
	std::string personalityParamValue;

	//! whether or not STDIN was already requested by one of the following
	//! arguments
	bool stdinUsed;

	//! query specific options
	boost::program_options::options_description opts_query;
	//! target specific options
	boost::program_options::options_description opts_target;
	//! helix specific options
	boost::program_options::options_description opts_helix;
	//! seed specific options
	boost::program_options::options_description opts_seed;
	//! SHAPE reactivity data specific options
	boost::program_options::options_description opts_shape;
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
	//! the id (prefix) to be used for query naming
	std::string qId;
	//! in/output index of pos 0 (of all queries)
	NumberParameter<long> qIdxPos0;
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
	mutable IndexRangeListVec qRegion;
	//! maximal length of automatically detected highly accessible regions for
	//! for query sequences; if 0, no automatic detection is done
	NumberParameter<int> qRegionLenMax;
	//! optional file name that contains structure probing reactivity data for
	//! the query sequence (e.g. SHAPE data) to guide accessibility prediction
	std::string qShape;
	//! optional encoding what method is to be used to convert the data from
	//! qShape to pseudo energies for according accessibility prediction
	std::string qShapeMethod;
	//! optional encoding how data from qShape is converted into pairing
	//! probabilities for according accessibility prediction
	std::string qShapeConversion;

	//! the target command line argument
	std::string targetArg;
	//! the container holding all target sequences
	RnaSequenceVec target;
	//! the id (prefix) to be used for target naming
	std::string tId;
	//! in/output index of pos 0 (of all targets)
	NumberParameter<long> tIdxPos0;
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
	mutable IndexRangeListVec tRegion;
	//! maximal length of automatically detected highly accessible regions for
	//! for target sequences; if 0, no automatic detection is done
	NumberParameter<int> tRegionLenMax;
	//! optional file name that contains structure probing reactivity data for
	//! the target sequence (e.g. SHAPE data) to guide accessibility prediction
	std::string tShape;
	//! optional encoding what method is to be used to convert the data from
	//! tShape to pseudo energies for according accessibility prediction
	std::string tShapeMethod;
	//! optional encoding how data from tShape is converted into pairing
	//! probabilities for according accessibility prediction
	std::string tShapeConversion;

	// META PARAMETER applied to both query and target
	//! accessibility computation mode
	CharParameter acc;
	//! window length for accessibility computation (plFold)
	NumberParameter<int> accW;
	//! maximal base pair span for accessibility computation (plFold)
	NumberParameter<int> accL;
	//! whether or not lonely base pairs are allowed in accessibility computation
	bool accNoLP;
	//! whether or not GU base pairs are allowed to close loops in accessibility computation
	bool accNoGUend;
	//! window length to be considered accessible/interacting
	NumberParameter<int> intLenMax;
	//! maximal internal loop length to be considered accessible/interacting
	NumberParameter<int> intLoopMax;


	//! the minimal number of base pairs allowed in the helix (>2)
	NumberParameter<int> helixMinBP;
	//! the maximal number of base pairs allowed in the helix (>helixMinBP)
	NumberParameter<int> helixMaxBP;
	//! maximal internal loop size in the helix computation (0-2)
	NumberParameter<int> helixMaxIL;
	//! minimal unpaired probability (per sequence) of a helix to be considered
	NumberParameter<Z_type> helixMinPu;
	//! maximal energy of a helix to be considered
	NumberParameter<E_kcal_type> helixMaxE;
	//! when set, full helix energy is to be used for energy checks
	bool helixFullE;
	//! the final helix constraint to be used
	mutable HelixConstraint * helixConstraint;

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
	NumberParameter<E_kcal_type> seedMaxE;
	//! minimal unpaired probability (per sequence) of a seed to be considered
	NumberParameter<Z_type> seedMinPu;
	//! max hybridization energy of a seed to be considered
	NumberParameter<E_kcal_type> seedMaxEhybrid;
	//! whether or not GU base pairs are allowed within seeds
	bool seedNoGU;
	//! whether or not GU base pairs are allowed at seed ends
	bool seedNoGUend;
	//! intervals in query for seed search
	std::string seedQRange;
	//! intervals in target for seed search
	std::string seedTRange;
	//! the final seed constraint to be used
	mutable SeedConstraint * seedConstraint;

	//! the temperature to be used for energy computations
	NumberParameter<Z_type> temperature;

	//! the interaction model to predict in (mfe-single-site, helix-based-single-site, max-prob-site, ..)
	CharParameter model;
	//! the prediction mode (heuristic, space-efficient, exact)
	CharParameter mode;
#if INTARNA_MULITHREADING
	//! number of threads = number of parallel predictors running
	NumberParameter<int> threads;
#endif
	//! the window width to be used for window-based computations
	NumberParameter<int> windowWidth;
	//! the window overlap to be used for window-based computations
	NumberParameter<int> windowOverlap;

	//! the selected energy model
	CharParameter energy;
	//! the provided energy parameter file of the VRNA package
	std::string energyFile;

	//! the energy shift to be applied
	NumberParameter<E_kcal_type> energyAdd;
	//! whether or not the overall energy covers dangling end contributions
	bool energyNoDangles;


	//! where to write the output to and for each in what format
//	std::string out;
	std::vector< std::string > out;
	//! provides the parsed stream name for each out prefix
	std::map<OutPrefixCode,std::string> outPrefix2streamName;
	//! output mode
	CharParameter outMode;
	//! number of (sub)optimal interactions to report
	NumberParameter<int> outNumber;
	//! whether or not reported interactions can to be overlapping
	CharParameter outOverlap;
	//! deltaE to mfe allowed to report an interaction
	NumberParameter<E_kcal_type> outDeltaE;
	//! max E allowed to report an interaction
	NumberParameter<E_kcal_type> outMaxE;
	//! min unpaired prob of an interacting subsequence allowed
	NumberParameter<Z_type> outMinPu;
	//! whether or not only the best seed is to be reported
	bool outBestSeedOnly;
	//! whether or not only lonely (non-stacked) inter-molecular base pairs are allowed
	bool outNoLP;
	//! whether or not GU base pairs are allowed at interaction and helix ends
	bool outNoGUend;
	//! the column separator to be used for tabular output
	std::string outSep;
	//! the list separator within individual columns of tabular output
	static const std::string outCsvLstSep;
	//! the CSV column selection
	std::string outCsvCols;
	//! the column ID from outCsvCols to be used for sorting the output
	std::string outCsvSort;
	//! the CSV column selection
	static const std::string outCsvCols_default;
	//! whether or not best interaction output should be provided independently
	//! for all region combinations or only the best for each query-target
	//! combination
	bool outPerRegion;
	//! for SpotProb output : spots to be tracked
	std::string outSpotProbSpots;
	//! whether or not Zall is needed for output generation
	mutable bool outNeedsZall;
	//! whether or not interaction BPs are needed for output generation
	mutable bool outNeedsBPs;

	//! (optional) file name for log output
	std::string logFileName;
	//! (optional) file name for input parameter configuration file
	std::string configFileName;

	//! the vienna energy parameter handler initialized by #parse()
	mutable VrnaHandler vrnaHandler;

	//! the handler of the final output stream
	OutputStreamHandler * outStreamHandler;

protected:

	/**
	 * sets the stdinUsed member to true if so far false or raises an exception
	 * if it is already true.
	 * @return true if stdinUsed was false so far; false otherwise (error logged)
	 */
	bool setStdinUsed();


	/**
	 * resets the default value of a parameter member to the given value.
	 * if the value is changed, a respective VLOG output is produced.
	 * @param param the parameter member to update
	 * @param value the new default value to set
	 */
	template <typename Param, typename Value>
	void
	resetParamDefault( Param & param, Value value ) {
		if (param.def != value) {
			param.def = value;
			param.val = value;
			VLOG(1) <<"  "<<param.name<<"=" <<value;
		}
	}

	/**
	 * resets the default value of a parameter member to the given value.
	 * if the value is changed, a respective VLOG output is produced.
	 * @param param the parameter member to update
	 * @param value the new default value to set
	 * @param paramName the name of the parameter as used in the CLI
	 */
	template <typename ParamType>
	void
	resetParamDefault( ParamType & param, ParamType value, const std::string paramName ) {
		if (param != value) {
			param = value;
			VLOG(1) <<"  "<<paramName<<"=" <<value;
		}
	}


	////////////  INDIVIDUAL TESTS  //////////////////

	/**
	 * Validates the query's qSet argument.
	 * @param value the argument value to validate
	 */
	void validate_qSet(const std::string & value);

	/**
	 * Validates the qAccFile argument.
	 * @param value the argument value to validate
	 */
	void validate_qAccFile(const std::string & value);

	/**
	 * Validates the SHAPE reactivity data file.
	 * @param paramName the name of the argument
	 * @param value the filename of the query's SHAPE reactivity data
	 */
	void validate_shape(const std::string& paramName, const std::string & value );

	/**
	 * Validates the method to transform SHAPE reactivity data to
	 * pseudo energies.
	 * @param paramName the name of the argument
	 * @param value the query's SHAPE method encoding
	 */
	void validate_shapeMethod(const std::string& paramName,  const std::string & value );

	/**
	 * Validates the SHAPE reactivity data conversion method encoding.
	 * @param paramName the name of the argument
	 * @param value the query's SHAPE conversion method encoding
	 */
	void validate_shapeConversion(const std::string& paramName,  const std::string & value );

	/**
	 * Validates the target's tSet argument.
	 * @param value the argument value to validate
	 */
	void validate_tSet(const std::string & value);

	/**
	 * Validates the tAccFile argument.
	 * @param value the argument value to validate
	 */
	void validate_tAccFile(const std::string & value);

	/**
	 * Validates the explicit seed argument.
	 * @param value the argument value to validate
	 */
	void validate_seedTQ(const std::string & value);

	/**
	 * Validates the seed?Range argument.
	 * @param paramName the name of the argument
	 * @param value the argument value to validate
	 */
	void validate_seedRange(const std::string& paramName, const std::string & value);

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
	 * Validates the tabular column separator argument.
	 *
	 * @param sep the argument value to validate
	 */
	void validate_outSep(const std::string & sep);

	/**
	 * Validates the outCsvCols argument.
	 * @param value the argument value to validate
	 */
	void validate_outCsvCols(const std::string & value);

	/**
	 * Validates the outCsvSort argument.
	 * @param value the argument value to validate
	 */
	void validate_outCsvSort(const std::string & value);

	/**
	 * Validates the configFileName argument.
	 * @param value the argument value to validate
	 */
	void validate_configFileName( const std::string & value);

	////////////  GENERIC TESTS  /////////////////

	/**
	 * Validates a CharParameter.
	 * @param argName the name of the parameter (for exception handling)
	 * @param param the parameter object
	 * @param value the value of the parameter to validate
	 */
	void validate_charArgument(const CharParameter& param, const char & value);


	/**
	 * Validates a NumberParameter.
	 * @param param the parameter object
	 * @param value the value of the parameter to validate
	 */
	template <typename T>
	void validate_numberArgument(const NumberParameter<T> & param, const T& value)
	{
		// range check
		if ( ! param.isInRange(value) ) {
			LOG(ERROR) <<param.name<<" = " <<value <<" : has to be in the range [" <<param.min <<","<<param.max<<"]";
			updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
		}
	}

	/**
	 * Validates a NumberParameter and excludes a certain value range.
	 * @param param the parameter object
	 * @param value the value of the parameter to validate
	 * @param minExcl the minimum of the range to be excluded
	 * @param minExcl the maximum of the range to be excluded
	 */
	template <typename T>
	void validate_numberArgumentExcludeRange(const NumberParameter<T> & param, const T& value, const T& minExcl, const T& maxExcl)
	{
		// standard check
		validate_numberArgument( param, value );

		// check excluded range
		if (param.val >= minExcl && param.val <= maxExcl) {
			LOG(ERROR) <<"\n "<<param.name<<" = " <<value <<" : has to be < "<<minExcl <<" or > "<<maxExcl;
			updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
		}
	}


	/**
	 * Validates the string encoding of an IndexRangeList encoding parameter
	 * @param argName the name of the parameter (for exception handling)
	 * @param value the value of the parameter to validate
	 * @param seq the RnaSequence to be used for mapping in/output indices to
	 *          internal index positions
	 */
	void validate_indexRangeList(const std::string & argName
								, const std::string & value
								, const RnaSequence & seq );

	/**
	 * Validates a sequence arguments.
	 * @param argument the argument name of this parameter
	 * @param value the argument value to validate
	 */
	void validate_sequenceArgument(const std::string& argument, const std::string & value);

	/**
	 * Validates a sequence id arguments.
	 * @param argument the argument name of this parameter
	 * @param value the argument value to validate
	 */
	void validate_id(const std::string& argument, const std::string & value);


	/**
	 * Validates a structure constraint argument.
	 * @param argument the argument name of this parameter
	 * @param value the argument value to validate
	 */

	void validate_structureConstraintArgument(const std::string& argument, const std::string & value);


	/**
	 * Parses the parameter value and returns all parsed sequences.
	 * @param paramName the name of the parameter (for exception handling)
	 * @param idPrefix the id (FASTA prefix) to be used for sequence id setup
	 * @param paramArg the given argument for the parameter
	 * @param sequences the container to fill
	 * @param seqSubset the indices of the input sequences to store (all other
	 *                  ignored)
	 * @param idxPos0 input/output index of first sequence position to be used
	 */
	void parseSequences(const std::string & paramName,
					const std::string& idPrefix,
					const std::string& paramArg,
					RnaSequenceVec& sequences,
					const IndexRangeList & seqSubset,
					const long idxPos0 );

	/**
	 * Parses the parameter input stream from FASTA format and returns all
	 * parsed sequences.
	 * @param paramName the name of the parameter (for exception handling)
	 * @param idPrefix the FASTA-id prefix to be used for sequence id setup
	 * @param input the input stream from where to read the FASTA data
	 * @param sequences the container to fill
	 * @param seqSubset the indices of the input sequences to store (all other
	 *                  ignored)
	 * @param idxPos0 input/output index of first sequence position to be used
	 */
	void parseSequencesFasta( const std::string & paramName,
					const std::string& idPrefix,
					std::istream& input,
					RnaSequenceVec& sequences,
					const IndexRangeList & seqSubset,
					const long idxPos0 );

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
void CommandLineParsing::validate_id(const std::string& paramName, const std::string & value)
{
	std::string invalidChars = "\n\r";
	for (auto c : invalidChars) {
		if (value.find(c) != std::string::npos) {
			LOG(ERROR) <<paramName<<" = '" <<value <<"' contains a line break";
			updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
		}
	}
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
			updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
		} else {
			// parse and store subset definitions
			qSet = IndexRangeList(value);
		}
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
				LOG(ERROR) <<"reading query accessibilities for multiple sequences from '"<<value<<"' is not supported";
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
void CommandLineParsing::validate_shape(const std::string& paramName, const std::string & value )
{
	if (!validateFile( value )) {
		LOG(ERROR) <<"Can not access/read SHAPE reactivity file for "<<paramName<<"='" <<value <<"'";
		updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
	}
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_shapeMethod( const std::string& paramName, const std::string & value )
{

	if (!boost::regex_match(value, AccessibilityConstraint::regexShapeMethod, boost::match_perl) ) {
		LOG(ERROR) <<"SHAPE method encoding "<<paramName<<"='" <<value <<"' is not valid.";
		updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
	}
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_shapeConversion( const std::string& paramName, const std::string & value )
{
	if (!boost::regex_match(value, AccessibilityConstraint::regexShapeConversion, boost::match_perl) ) {
		LOG(ERROR) <<"SHAPE conversion method encoding "<<paramName<<"'" <<value <<"' is not valid.";
		updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
	}
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
			updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
		} else {
			// parse and store subset definitions
			tSet = IndexRangeList(value);
		}
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
				updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
			}
			// update p1
			p1 = value.find_first_not_of(',',p2);
		}
		if (non_empty_seeds == 0) {
			LOG(ERROR) <<"explicit seed encoding '" <<value <<"' does not contain any seed information";
			updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
		}
	}
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_seedRange(const std::string& paramName, const std::string & value) {
	if (!value.empty()) {
		// check regex
		if (!boost::regex_match(value, IndexRangeList::regex, boost::match_perl) ) {
			LOG(ERROR) <<paramName<<" = " <<value <<" : is not in the format 'from1-to1,from2-to2,..'";
			updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
		}
	}
}

////////////////////////////////////////////////////////////////////////////

inline
void
CommandLineParsing::
validate_energyFile(const std::string & value)
{
	// check for supported default models
	if (value == std::string(VrnaHandler::Turner99)
		|| value == std::string(VrnaHandler::Turner04)
		|| value == std::string(VrnaHandler::Andronescu07)
	)
	{
		return;
	}
	// check if file exists and is readable
	if (!validateFile( value )) {
		LOG(ERROR) <<"provided VRNA energy parameter file '" <<value <<"' could not be processed.";
		updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
	}
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

	///////////  ASSUME IT IS A FILENAME/-PATH  ///////////////

	// check if empty filename
	if ( value.empty() ) {
		LOG(ERROR) <<argName<<" : no output defined (empty)";
		updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
		return;
	}

	// check if parent directory exists
	boost::filesystem::path p(value);
	if ( !p.parent_path().empty() && ! boost::filesystem::exists(p.parent_path()) ) {
		LOG(ERROR) <<argName<<" : parent directory '"<<p.parent_path().string()<<"' of output file '"<<value<<"' does not exist";
		updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
		return;
	}

	// check if output file exists and can not be overwritten
	if ( boost::filesystem::exists(p) ) {
		// start overwriting and catch if needed
		bool fileCanNotBeOverwritten = false;
		try {
			// open dummy file stream to check if writeable
			std::ofstream file(value);
			if (!file) {
				fileCanNotBeOverwritten = true;
			}
			file.close();
		} catch (std::exception & ex) {
			fileCanNotBeOverwritten = true;
		}
		if (fileCanNotBeOverwritten) {
			LOG(ERROR) <<argName<<" : of output file '"<<value<<"' exists and can not be overwritten";
			updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
			return;
		}
	}
}

////////////////////////////////////////////////////////////////////////////

inline
void CommandLineParsing::validate_outSep(const std::string & outSep) {

	// check if empty
	if ( outSep.empty() ) {
		LOG(ERROR) <<"no column separator for tabular CSV output defined (empty)";
		updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
		return;
	}

	// check if distinct from value list separator
	if ( outSep == outCsvLstSep ) {
		LOG(ERROR) <<"the selected column separator '"<<outSep<<"' for tabular CSV output equals the value separator used to multi-value output. Please change!";
		updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
		return;
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

			// handle SpotProb setup
			if (curPrefCode == OP_spotProb) {
				// if spots are defined
				if (streamName.find(':') != std::string::npos) {
					// get spots
					outSpotProbSpots = streamName.substr(0,streamName.find(':'));
					// check if valid spot encoding
					if (!boost::regex_match(outSpotProbSpots, PredictionTrackerSpotProb::regexSpotString, boost::match_perl)){
						// check sanity of spot encodings (indexing starts with 1)
						LOG(ERROR) <<"--out : spot encoding '"<<curPref<<"' is not valid.";
						updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
						break;
					}
					// get stream name (remove spot encoding prefix)
					streamName = streamName.substr(outSpotProbSpots.size()+1);
				} else {
					// change PrefCode to reflect exhaustive spot prob computation
					curPrefCode = OP_spotProbAll;
				}
			}

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
void CommandLineParsing::validate_configFileName( const std::string & value )
{
	if (!validateFile( value )) {
		LOG(ERROR) <<"Can not access/read parameter file '" <<value <<"'";
	}
}

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
	return threads.val == 0 ? threads.max : threads.val;
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

inline
bool
CommandLineParsing::
reportBestPerRegion() const
{
	return outPerRegion;
}

////////////////////////////////////////////////////////////////////////////




#endif /* COMMANDLINEPARSING_H_ */
