
#include "CommandLineParsing.h"

#include "general.h"

#include <cmath>
#include <stdexcept>

#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

#include "AccessibilityConstraint.h"

#include "AccessibilityDisabled.h"
#include "AccessibilityVrna.h"

#include "InteractionEnergyBasePair.h"
#include "InteractionEnergyVrna.h"

#include "PredictorMfe2d.h"
#include "PredictorMfe4d.h"
#include "PredictorMaxProb.h"
#include "PredictorMfe2dSeed.h"

#include "OutputHandlerText.h"





////////////////////////////////////////////////////////////////////////////

const boost::regex CommandLineParsing::regexRangeEncoding("^([123456789]\\d*-[123456789]\\d*,)*[123456789]\\d*-[123456789]\\d*$");

////////////////////////////////////////////////////////////////////////////

CommandLineParsing::CommandLineParsing()
	:
	opts_query("Query"),
	opts_target("Target"),
	opts_seed("Seed"),
	opts_inter("Interaction"),
	opts_general("General"),
	opts_cmdline_all(),

	parsingCode(NOT_PARSED_YET),

	queryArg(""),
	query(),
	qAcc("NF",'F'),
	qAccW( 0, 99999, 150),
	qAccL( 0, 99999, 100),
	qAccConstr(""),
	qIntLenMax( 0, 99999, 0),
	qIntLoopMax( 0, 100, 16),
	qRegionString(""),
	qRegion(),

	targetArg(""),
	target(),
	tAcc("NF",'F'),
	tAccW( 0, 99999, 150),
	tAccL( 0, 99999, 100),
	tAccConstr(""),
	tIntLenMax( 0, 99999, 0),
	tIntLoopMax( 0, 100, 16),
	tRegionString(""),
	tRegion(),

	noSeedRequired(false),
	seedBP(2,20,7),
	seedMaxUP(0,20,0),
	seedMaxUPq(-1,20,-1),
	seedMaxUPt(-1,20,-1),
	seedMaxE(-999,+999,0),
	seedConstraint(NULL),

	temperature(0,100,37),

	predMode(1,3,1),

	energy("BF",'F'),
	energyFile(""),

	vrnaHandler()

{
	using namespace boost::program_options;


	////  SEQUENCE OPTIONS  ////////////////////////////////////

	opts_query.add_options()
		("query,q"
			, value<std::string>(&queryArg)
				->required()
				->notifier(boost::bind(&CommandLineParsing::validate_query,this,_1))
			, "either an RNA sequence or the stream/file name from where to read the query sequences; use 'STDIN' to read from standard input stream")
		("qAcc"
			, value<char>(&(qAcc.val))
				->default_value(qAcc.def)
				->notifier(boost::bind(&CommandLineParsing::validate_qAcc,this,_1))
			, std::string("accessibility computation : 'N'o accessibility contributions, or 'F'ull accessibility computation").c_str())
		("qAccW"
			, value<int>(&(qAccW.val))
				->default_value(qAccW.def)
				->notifier(boost::bind(&CommandLineParsing::validate_qAccW,this,_1))
			, std::string("accessibility computation : sliding window size for query accessibility computation (arg in range ["+toString(qAccW.min)+","+toString(qAccW.max)+"]; 0 defaults to the full sequence length)").c_str())
		("qAccL"
			, value<int>(&(qAccL.val))
				->default_value(qAccL.def)
				->notifier(boost::bind(&CommandLineParsing::validate_qAccL,this,_1))
			, std::string("accessibility computation : sliding window size for query accessibility computation (arg in range ["+toString(qAccL.min)+","+toString(qAccL.max)+"]; 0 defaults to sliding window size 'qAccW')").c_str())
		("qAccConstr"
			, value<std::string>(&(qAccConstr))
				->notifier(boost::bind(&CommandLineParsing::validate_qAccConstr,this,_1))
			, std::string("accessibility computation : structure constraint for each sequence position: '.' no constraint, '"+toString(AccessibilityConstraint::dotBracket_accessible)+"' unpaired, '"+toString(AccessibilityConstraint::dotBracket_blocked)+"' blocked. Note, blocked positions are excluded from interaction prediction and considered unpaired!").c_str())
		("qIntW"
			, value<int>(&(qIntLenMax.val))
				->default_value(qIntLenMax.def)
				->notifier(boost::bind(&CommandLineParsing::validate_qIntLenMax,this,_1))
			, std::string("interaction site : maximal window size to be considered for interaction (and thus accessible) within the query (arg in range ["+toString(qIntLenMax.min)+","+toString(qIntLenMax.max)+"]; 0 defaults to the full sequence length)").c_str())
		("qIntL"
			, value<int>(&(qIntLoopMax.val))
				->default_value(qIntLoopMax.def)
				->notifier(boost::bind(&CommandLineParsing::validate_qIntLoopMax,this,_1))
			, std::string("interaction site : maximal number of unpaired bases between neighbored interacting bases to be considered in interactions within the query (arg in range ["+toString(qIntLoopMax.min)+","+toString(qIntLoopMax.max)+"]; 0 enforces stackings only)").c_str())
		("qRegion"
			, value<std::string>(&(qRegionString))
				->notifier(boost::bind(&CommandLineParsing::validate_qRegion,this,_1))
			, std::string("interaction site : query regions to be considered for interaction prediction. Either given as BED file (for multi-sequence FASTA input) or in the format 'from1-to1,from2-to2,..' assuming indexing starts with 1").c_str())
		;

	opts_target.add_options()
		("target,t"
			, value<std::string>(&targetArg)
				->required()
				->notifier(boost::bind(&CommandLineParsing::validate_target,this,_1))
				, "either an RNA sequence or the stream/file name from where to read the target sequences; use 'STDIN' to read from standard input stream")
		("tAcc"
			, value<char>(&(tAcc.val))
				->default_value(tAcc.def)
				->notifier(boost::bind(&CommandLineParsing::validate_tAcc,this,_1))
			, std::string("accessibility computation : 'N'o accessibility contributions, or 'F'ull accessibility computation").c_str())
		("tAccW"
			, value<int>(&(tAccW.val))
				->default_value(tAccW.def)
				->notifier(boost::bind(&CommandLineParsing::validate_tAccW,this,_1))
			, std::string("accessibility computation : sliding window size for query accessibility computation (arg in range ["+toString(tAccW.min)+","+toString(tAccW.max)+"]; 0 defaults to the full sequence length)").c_str())
		("tAccL"
			, value<int>(&(tAccL.val))
				->default_value(tAccL.def)
				->notifier(boost::bind(&CommandLineParsing::validate_tAccL,this,_1))
			, std::string("accessibility computation : sliding window size for query accessibility computation (arg in range ["+toString(tAccL.min)+","+toString(tAccL.max)+"]; 0 defaults to sliding window size 'tAccW')").c_str())
		("tAccConstr"
			, value<std::string>(&(tAccConstr))
				->notifier(boost::bind(&CommandLineParsing::validate_tAccConstr,this,_1))
			, std::string("accessibility computation : structure constraint for each sequence position: '.' no constraint, '"+toString(AccessibilityConstraint::dotBracket_accessible)+"' unpaired, '"+toString(AccessibilityConstraint::dotBracket_blocked)+"' blocked. Note, blocked positions are excluded from interaction prediction and considered unpaired!").c_str())
		("tIntLenMax"
			, value<int>(&(tIntLenMax.val))
				->default_value(tIntLenMax.def)
				->notifier(boost::bind(&CommandLineParsing::validate_tIntLenMax,this,_1))
			, std::string("interaction site : maximal window size to be considered for interaction (and thus accessible) within the target (arg in range ["+toString(tIntLenMax.min)+","+toString(tIntLenMax.max)+"]; 0 defaults to the full sequence length)").c_str())
		("tIntLoopMax"
			, value<int>(&(tIntLoopMax.val))
				->default_value(tIntLoopMax.def)
				->notifier(boost::bind(&CommandLineParsing::validate_tIntLoopMax,this,_1))
			, std::string("interaction site : maximal number of unpaired bases between neighbored interacting bases to be considered in interactions within the target (arg in range ["+toString(tIntLoopMax.min)+","+toString(tIntLoopMax.max)+"]; 0 enforces stackings only)").c_str())
		("tRegion"
			, value<std::string>(&(tRegionString))
				->notifier(boost::bind(&CommandLineParsing::validate_tRegion,this,_1))
			, std::string("interaction site : target regions to be considered for interaction prediction. Either given as BED file (for multi-sequence FASTA input) or in the format 'from1-to1,from2-to2,..' assuming indexing starts with 1").c_str())
		;

	////  SEED OPTIONS  ////////////////////////////////////


	opts_seed.add_options()
	    ("noSeed", "if present, no seed is enforced within the predicted interactions")
	    ;

	opts_seed.add_options()
		("seedBP"
			, value<int>(&(seedBP.val))
				->default_value(seedBP.def)
				->notifier(boost::bind(&CommandLineParsing::validate_seedBP,this,_1))
			, std::string("number of inter-molecular base pairs within the seed region (arg in range ["+toString(seedBP.min)+","+toString(seedBP.max)+"])").c_str())
		;

	opts_seed.add_options()
		("seedMaxUP"
			, value<int>(&(seedMaxUP.val))
				->default_value(seedMaxUP.def)
				->notifier(boost::bind(&CommandLineParsing::validate_seedMaxUP,this,_1))
			, std::string("maximal overall number (query+target) of unpaired bases within the seed region (arg in range ["+toString(seedMaxUP.min)+","+toString(seedMaxUP.max)+"])").c_str())
		;

	opts_seed.add_options()
		("seedMaxUPq"
			, value<int>(&(seedMaxUPq.val))
				->default_value(seedMaxUPq.def)
				->notifier(boost::bind(&CommandLineParsing::validate_seedMaxUPq,this,_1))
			, std::string("maximal number of unpaired bases within the query's seed region (arg in range ["+toString(seedMaxUPq.min)+","+toString(seedMaxUPq.max)+"]); if -1 the value of seedMaxUP is used.").c_str())
		;

	opts_seed.add_options()
		("seedMaxUPt"
			, value<int>(&(seedMaxUPt.val))
				->default_value(seedMaxUPt.def)
				->notifier(boost::bind(&CommandLineParsing::validate_seedMaxUPt,this,_1))
			, std::string("maximal number of unpaired bases within the target's seed region (arg in range ["+toString(seedMaxUPt.min)+","+toString(seedMaxUPt.max)+"]); if -1 the value of seedMaxUP is used.").c_str())
		;

	opts_seed.add_options()
		("seedMaxE"
			, value<E_type>(&(seedMaxE.val))
				->default_value(seedMaxE.def)
				->notifier(boost::bind(&CommandLineParsing::validate_seedMaxE,this,_1))
			, std::string("maximal energy a seed region may have (arg in range ["+toString(seedMaxE.min)+","+toString(seedMaxE.max)+"]).").c_str())
		;

	////  INTERACTION/ENERGY OPTIONS  ////////////////////////

	opts_inter.add_options()
		("mode,m"
			, value<int>(&(predMode.val))
				->default_value(predMode.def)
				->notifier(boost::bind(&CommandLineParsing::validate_predMode,this,_1))
			, std::string("prediction mode : 1 (default) : MFE (2D), 2 : MFE (4D), 3 : MaxProb (4D)").c_str())
		("energy,e"
			, value<char>(&(energy.val))
				->default_value(energy.def)
				->notifier(boost::bind(&CommandLineParsing::validate_energy,this,_1))
			, std::string("energy computation : 'B'ase pair == -1, or 'F'ull VRNA-based computation").c_str())
		("energyVRNA"
			, value<std::string>(&energyFile)
				->notifier(boost::bind(&CommandLineParsing::validate_energyFile,this,_1))
			, std::string("energy parameter file of VRNA package to be used").c_str())
		("temperature"
			, value<T_type>(&(temperature.val))
				->default_value(temperature.def)
				->notifier(boost::bind(&CommandLineParsing::validate_temperature,this,_1))
			, std::string("temperature in Celsius to setup the VRNA energy parameters (arg in range ["+toString(temperature.min)+","+toString(temperature.max)+"])").c_str())
		;

	// TODO parse energy function selection

	////  GENERAL OPTIONS  ////////////////////////////////////

	opts_general.add_options()
	    ("version", "print version")
	    ("verbose,v", "verbose output")
	    ("default-log-file", "name of log file to be used for output")
	    ("help,h", "show the help page with all available parameters")
	    ;

	////  GENERAL OPTIONS  ////////////////////////////////////

	opts_cmdline_all.add(opts_query).add(opts_target).add(opts_seed).add(opts_inter).add(opts_general);


}

////////////////////////////////////////////////////////////////////////////

CommandLineParsing::~CommandLineParsing() {

	CLEANUP(seedConstraint);

}

////////////////////////////////////////////////////////////////////////////

CommandLineParsing::ReturnCode
CommandLineParsing::
parse(int argc, char** argv)
{

	// init: nothing parsed yet
	parsingCode = ReturnCode::NOT_PARSED_YET;

	using namespace boost::program_options;

	variables_map vm;
	try {
		store( parse_command_line(argc, argv, opts_cmdline_all), vm);
		// parsing fine so far
		parsingCode = ReturnCode::KEEP_GOING;
	} catch (error& e) {
		LOG(ERROR) <<e.what() << " : run with '--help' for allowed arguments";
		updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
	}

	// if parsing was successful, check for help request
	if (parsingCode == ReturnCode::KEEP_GOING) {
		if (vm.count("help")) {
			std::cout
				<<"\nIntaRNA predicts RNA-RNA interactions.\n"
				<<"\nThe following program arguments are supported:\n"
				<< opts_cmdline_all
				<< "\n";
			parsingCode = ReturnCode::STOP_ALL_FINE;
			return parsingCode;
		}
		if (vm.count("version")) {
			std::cout <<PACKAGE_STRING << "\n";
			parsingCode = ReturnCode::STOP_ALL_FINE;
			return parsingCode;
		}
	}

	// if parsing was successful, continue with additional checks
	if (parsingCode == ReturnCode::KEEP_GOING) {
		try {
			// run all notifier checks
			notify(vm);
		} catch (required_option& e) {
			LOG(ERROR) <<"mandatory option '"<<e.get_option_name() << "' not provided";
			updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
		} catch (error& e) {
			LOG(ERROR) <<e.what();
			updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
		}
	}

	// if parsing was successful, continue with final parsing
	if (parsingCode == ReturnCode::KEEP_GOING) {
		try {

			// parse the sequences
			parseSequences("query",queryArg,query);
			// parse the sequences
			parseSequences("target",targetArg,target);

			// check for minimal sequence length (>=seedBP)
			for( size_t i=0; i<query.size(); i++) {
				if (query.at(i).size() < seedBP.val) {
					throw error("length of query sequence "+toString(i+1)+" is below minimal number of seed base pairs (seedBP="+toString(seedBP.val)+")");
				}
			}
			for( size_t i=0; i<target.size(); i++) {
				if (target.at(i).size() < seedBP.val) {
					throw error("length of target sequence "+toString(i+1)+" is below minimal number of seed base pairs (seedBP="+toString(seedBP.val)+")");
				}
			}

			// parse regions to be used for interaction prediction
			parseRegion( "qRegion", qRegionString, query, qRegion );
			parseRegion( "tRegion", tRegionString, target, tRegion );

			// check qAccConstr - query sequence compatibility
			if (vm.count("qAccConstr") > 0) {
				// only for single sequence input supported
				if (validateSequenceNumber("qAccConstr",query,1,1)) {
					// check length
					if (qAccConstr.size() != query.at(0).size()) {
						throw error("qAccConstr and query sequence differ in size");
					}
				} else {
					// TODO report error
					NOTIMPLEMENTED("check not implemented");
				}
			} else {
				// generate empty constraint
				qAccConstr = std::string(query.at(0).size(),'.');
			}
			// check tAccConstr - target sequence compatibility
			if (vm.count("tAccConstr") > 0) {
				// only for single sequence input supported
				if (validateSequenceNumber("tAccConstr",target,1,1)) {
					// check length
					if (tAccConstr.size() != target.at(0).size()) {
						throw error("tAccConstr and target sequence differ in size");
					}
				} else {
					// TODO report error
					NOTIMPLEMENTED("check not implemented");
				}
			} else {
				// generate empty constraint
				tAccConstr = std::string(target.at(0).size(),'.');
			}

			// check energy setup
			if (vm.count("energyVRNA") > 0 && energy.val != 'F') {
				throw error("energyVRNA provided but no (F)ull energy computation requested");
			}

			// check seed setup
			noSeedRequired = vm.count("noSeed") > 0;
			if (noSeedRequired) {
				// input sanity check : maybe seed constraints defined -> warn
				if (seedBP.val != seedBP.def) LOG(INFO) <<"no seed constraint wanted, but seedBP provided (will be ignored)";
				if (seedMaxUP.val != seedMaxUP.def) LOG(INFO) <<"no seed constraint wanted, but seedMaxUP provided (will be ignored)";
				if (seedMaxUPq.val != seedMaxUPq.def) LOG(INFO) <<"no seed constraint wanted, but seedMaxUPq provided (will be ignored)";
				if (seedMaxUPt.val != seedMaxUPt.def) LOG(INFO) <<"no seed constraint wanted, but seedMaxUPt provided (will be ignored)";
				if (seedMaxE.val != seedMaxE.def) LOG(INFO) <<"no seed constraint wanted, but seedMaxE provided (will be ignored)";
			}

		} catch (error& e) {
			LOG(ERROR) <<e.what();
			updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
		} catch (std::exception& e) {
			LOG(ERROR) <<e.what();
			updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
		}
	}

	// setup new VRNA handler with the given arguments
	if ( energy.val == 'F') {
		vrnaHandler = VrnaHandler( temperature.val, (energyFile.size() > 0 ? & energyFile : NULL) );
	}


	// return validate_* dependent parsing code
	return parsingCode;
}

////////////////////////////////////////////////////////////////////////////

void CommandLineParsing::validate_query(const std::string & value)
{
	validate_sequenceArgument("query",value);
}

////////////////////////////////////////////////////////////////////////////

void CommandLineParsing::validate_qIntLenMax(const int & value)
{
	// forward check to general method
	validate_numberArgument("qIntW", qIntLenMax, value);
}

////////////////////////////////////////////////////////////////////////////

void CommandLineParsing::validate_qIntLoopMax(const int & value)
{
	// forward check to general method
	validate_numberArgument("qIntL", qIntLoopMax, value);
}

////////////////////////////////////////////////////////////////////////////

void CommandLineParsing::validate_qAcc(const char & value)
{
	// forward check to general method
	validate_charArgument("qAcc", qAcc, value);
}

////////////////////////////////////////////////////////////////////////////

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

void CommandLineParsing::validate_qAccL(const int & value)
{
	// forward check to general method
	validate_numberArgument("qAccL", qAccL, value);
	// check upper bound
	if (qAccL.val > qAccW.val) {
		LOG(ERROR) <<"qAccL = " <<value <<" : has to be <= qAccW (=" <<qAccW.val<<")";
		updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
	}
	// check lower bound
	if (qAccL.val > 0 && qAccL.val < 3) {
		LOG(ERROR) <<"qAccL = " <<value <<" : has to be 0 or > 3";
		updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
	}
}

////////////////////////////////////////////////////////////////////////////

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

void CommandLineParsing::validate_qRegion(const std::string & value) {
	// check and store region information
	validateRegion( "qRegion", value );
}

////////////////////////////////////////////////////////////////////////////

void CommandLineParsing::validate_target(const std::string & value)
{
	validate_sequenceArgument("target",value);
}


////////////////////////////////////////////////////////////////////////////

void CommandLineParsing::validate_tIntLenMax(const int & value)
{
	// forward check to general method
	validate_numberArgument("tIntLenMax", tIntLenMax, value);
}

////////////////////////////////////////////////////////////////////////////

void CommandLineParsing::validate_tIntLoopMax(const int & value)
{
	// forward check to general method
	validate_numberArgument("tIntLoopMax", tIntLoopMax, value);
}

////////////////////////////////////////////////////////////////////////////

void CommandLineParsing::validate_tAcc(const char & value)
{
	// forward check to general method
	validate_charArgument("tAcc", tAcc, value);
}

////////////////////////////////////////////////////////////////////////////

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

void CommandLineParsing::validate_tAccL(const int & value)
{
	// forward check to general method
	validate_numberArgument("tAccL", tAccL, value);
	// check upper bound
	if (tAccL.val > tAccW.val) {
		LOG(ERROR) <<"tAccL = " <<value <<" : has to be <= tAccW (=" <<tAccW.val<<")";
		updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
	}
	// check lower bound
	if (tAccL.val > 0 && tAccL.val < 3) {
		LOG(ERROR) <<"tAccL = " <<value <<" : has to be 0 or > 3";
		updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
	}
}

////////////////////////////////////////////////////////////////////////////

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

void CommandLineParsing::validate_tRegion(const std::string & value) {
	// check and store region information
	validateRegion( "tRegion", value );
}

////////////////////////////////////////////////////////////////////////////

void CommandLineParsing::validate_seedBP(const int & value) {
	// forward check to general method
	validate_numberArgument("seedBP", seedBP, value);
}

////////////////////////////////////////////////////////////////////////////

void CommandLineParsing::validate_seedMaxUP(const int & value) {
	// forward check to general method
	validate_numberArgument("seedMaxUP", seedMaxUP, value);
}

////////////////////////////////////////////////////////////////////////////

void CommandLineParsing::validate_seedMaxUPq(const int & value) {
	// forward check to general method
	validate_numberArgument("seedMaxUPq", seedMaxUPq, value);
}

////////////////////////////////////////////////////////////////////////////

void CommandLineParsing::validate_seedMaxUPt(const int & value) {
	// forward check to general method
	validate_numberArgument("seedMaxUPt", seedMaxUPt, value);
}

////////////////////////////////////////////////////////////////////////////

void CommandLineParsing::validate_seedMaxE(const E_type & value) {
	// forward check to general method
	validate_numberArgument("seedMaxE", seedMaxE, value);
}

////////////////////////////////////////////////////////////////////////////

void CommandLineParsing::validate_temperature(const T_type & value) {
	// forward check to general method
	validate_numberArgument("temperature", temperature, value);
}

////////////////////////////////////////////////////////////////////////////

void CommandLineParsing::validate_predMode(const int & value)
{
	// forward check to general method
	validate_numberArgument("mode", predMode, value);

	// TODO check if mode with seed and "no-seed" requested

}

////////////////////////////////////////////////////////////////////////////

void CommandLineParsing::validate_energy(const char & value)
{
	// forward check to general method
	validate_charArgument("energy", energy, value);
}

////////////////////////////////////////////////////////////////////////////

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

void
CommandLineParsing::
validate_charArgument(const std::string & name, const CommandLineParsing::CharParameter& param, const char & value)
{
	// alphabet check
	if ( ! param.isInAlphabet(value) ) {
		LOG(ERROR) <<""<<name<<" = " <<value <<" : has to be one of '" <<param.alphabet <<"'";
		updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
	}
}

////////////////////////////////////////////////////////////////////////////

void
CommandLineParsing::
validate_sequenceArgument(const std::string & name, const std::string & value)
{
	if (value.compare("STDIN") != 0) {

		// check if it is a sequence
		if ( RnaSequence::isValidSequence(value) ) {

			// do nothing, all fine

		} else { // it is NOT a sequence! can be a file name!

			// check if a file of this name exists and is readable
			if ( ! validateFile( value ) ) {
				LOG(ERROR) <<""<<name<<" '"<<value <<"' : is neither STDIN, a file name, or a sequence!";
				updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
			}
		}

	}

	// TODO : check for ambiguous nucleotides 'N' --> give warning that these are ignored for pairing in interactions
}


////////////////////////////////////////////////////////////////////////////

void
CommandLineParsing::
validate_structureConstraintArgument(const std::string & name, const std::string & value)
{
	// check if valid alphabet
	if (value.find_first_not_of(AccessibilityConstraint::dotBracketAlphabet) != std::string::npos) {
		LOG(ERROR) <<""<<name<<" '"<<value <<"' : contains invalid characters!";
		updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
	} else {
		// check for base pair balance / nestedness
		int bpStackLvl = 0;
		for (std::string::const_iterator c = value.begin(); c!=value.end(); ++c) {
			switch (*c) {
			case '(' : ++bpStackLvl; break;
			case ')' : --bpStackLvl; break;
			default: break;
			}
			if (bpStackLvl<0) {
				LOG(ERROR) <<""<<name<<" '"<<value <<"' : unbalanced base pairs!";
				updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
			}
		}

	}
}

////////////////////////////////////////////////////////////////////////////

bool
CommandLineParsing::
validateFile( const std::string & filename )
{
	// check if file accessible
	if ( !boost::filesystem::exists( filename ) )
	{
		LOG(ERROR) <<"Can't find the file '"<<filename<<"'!";
		updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
		return false;
	}
	// check if it is a file
	if ( !boost::filesystem::is_regular_file( filename ) )
	{
		LOG(ERROR) <<"'"<<filename<<"' : Is no file!";
		updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
		return false;
	}
	return true;
}

////////////////////////////////////////////////////////////////////////////

bool
CommandLineParsing::
validateRegion( const std::string & argName, const std::string & value )
{
	// check if nothing given
	if (value.empty()) {
		return true;
	} else
	// check if direct range input
	if (boost::regex_match( value, regexRangeEncoding, boost::match_perl )) {
		return true;
	} else
	// might be BED file input
	if ( validateFile( value ) ) {
		return true;
	} else
	// no valid input
	{
		LOG(ERROR) <<"the argument for "<<argName<<" is neither a valid range string encoding nor a file that can be found";
		updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
	}
	return false;
}

////////////////////////////////////////////////////////////////////////////

void
CommandLineParsing::
parseRegion( const std::string & argName, const std::string & value, const RnaSequenceVec & sequences, IndexRangeListVec & rangeList )
{
	// check if nothing given
	if (value.empty()) {
		// ensure range list size sufficient
		rangeList.resize( sequences.size() );
		size_t s=0;
		BOOST_FOREACH( IndexRangeList & r, rangeList ) {
			// clear old data if any
			r.clear();
			assert(sequences.at(s).size()>0);
			// push full range
			r.push_back( IndexRange(0,sequences.at(s++).size()-1) );
		}
		return;
	} else
	// check direct range input
	if (boost::regex_match( value, regexRangeEncoding, boost::match_perl )) {
		// ensure single sequence input
		if(sequences.size() != 1) {
			LOG(ERROR) <<argName <<" : string range list encoding provided but more than one sequence present.";
			updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
			return;
		}
		// ensure range list size sufficient
		rangeList.resize(1);
		// clear first range list
		rangeList[0].clear();
		// decompose string into ranges
		boost::char_separator<char> commaSeparator(",");
		boost::tokenizer< boost::char_separator<char> > rangeStrings(value, commaSeparator);
		IndexRange toAdd(0,0);
		size_t dashPos = 0;
		BOOST_FOREACH(const std::string & rangeEncoding, rangeStrings)
		{
			// decompose range, check and store
			try {
				dashPos = rangeEncoding.find('-');
				// parse indices
				toAdd.from = boost::lexical_cast<size_t>( rangeEncoding.substr(0, dashPos) ) -1;
				toAdd.to = boost::lexical_cast<size_t>( rangeEncoding.substr( dashPos+1 ) ) -1;
				if (toAdd.to>=sequences.at(0).size()) {
					LOG(ERROR) <<argName <<" : upper boundary of range '"<<rangeEncoding<<"' exceeds the sequence's length "<<sequences.at(0).size();
					updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
				}
				// check sanity
				if (toAdd.isAscending()) {
					// insert range
					rangeList[0].insert( toAdd );
				} else {
					LOG(ERROR) <<argName <<" : error while parsing '"<<rangeEncoding<<"' : range not ascending";
					updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
				}
			}
			catch(const boost::bad_lexical_cast & e) {
				LOG(ERROR) <<"CommandLineParsing::parseRegion() : error while parsing '"<<rangeEncoding<<"' for argument "<<argName<<" : "<<e.what();
				updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
			}
		}
		return;
	}
	// might be BED file input
	if ( validateFile( value ) ) {
		NOTIMPLEMENTED("BED file input for index range list not implemented");
		return;
	}
	assert(false) /*should never happen*/;
}

////////////////////////////////////////////////////////////////////////////

const CommandLineParsing::RnaSequenceVec &
CommandLineParsing::
getQuerySequences() const
{
	checkIfParsed();
	return query;
}

////////////////////////////////////////////////////////////////////////////

const CommandLineParsing::RnaSequenceVec &
CommandLineParsing::
getTargetSequences() const
{
	checkIfParsed();
	return target;
}

////////////////////////////////////////////////////////////////////////////

Accessibility*
CommandLineParsing::
getQueryAccessibility( const size_t sequenceNumber ) const
{
	checkIfParsed();
	// input check
	if (sequenceNumber >= getQuerySequences().size()) {
		throw std::runtime_error("CommandLineParsing::getQueryAccessibility : sequence number "+toString(sequenceNumber)+" is out of range (<"+toString(getQuerySequences().size())+")");
	}
	const RnaSequence& seq = getQuerySequences().at(sequenceNumber);

	// construct selected accessibility object
	switch(qAcc.val) {
	case 'N' : return new AccessibilityDisabled( seq, qIntLenMax.val );
	case 'F' : {
		// create temporary constraint object (will be copied)
		AccessibilityConstraint accConstraint(qAccConstr);
		return new AccessibilityVrna( seq, vrnaHandler, qIntLenMax.val, qAccW.val, qAccL.val, &accConstraint);
	}
	default :
		NOTIMPLEMENTED("CommandLineParsing::getQueryAccessibility : qAcc = '"+toString(qAcc.val)+"' is not supported");
	}
	return NULL;
}

////////////////////////////////////////////////////////////////////////////

Accessibility*
CommandLineParsing::
getTargetAccessibility( const size_t sequenceNumber ) const
{
	checkIfParsed();
	// input check
	if (sequenceNumber >= getTargetSequences().size()) {
		throw std::runtime_error("CommandLineParsing::getTargetAccessibility : sequence number "+toString(sequenceNumber)+" is out of range (<"+toString(getTargetSequences().size())+")");
	}
	const RnaSequence& seq = getTargetSequences().at(sequenceNumber);
	switch(tAcc.val) {
	case 'N' : return new AccessibilityDisabled( seq, tIntLenMax.val );
	case 'F' : {
			// create temporary constraint object (will be copied)
			AccessibilityConstraint accConstraint(tAccConstr);
			return new AccessibilityVrna( seq, vrnaHandler, tIntLenMax.val, tAccW.val, tAccL.val, &accConstraint);
		}
	default :
		NOTIMPLEMENTED("CommandLineParsing::getTargetAccessibility : tAcc = '"+toString(tAcc.val)+"' is not supported");
	}
	return NULL;
}

////////////////////////////////////////////////////////////////////////////

InteractionEnergy*
CommandLineParsing::
getEnergyHandler( const Accessibility& accQuery, const ReverseAccessibility& accTarget ) const
{
	checkIfParsed();
	switch( energy.val ) {
	case 'B' : return new InteractionEnergyBasePair( accQuery, accTarget, qIntLoopMax.val, tIntLoopMax.val );
	case 'F' : return new InteractionEnergyVrna( accQuery, accTarget, vrnaHandler, qIntLoopMax.val, tIntLoopMax.val );
	default :
		NOTIMPLEMENTED("CommandLineParsing::getEnergyHandler : energy = '"+toString(energy.val)+"' is not supported");
	}
	return NULL;

}

////////////////////////////////////////////////////////////////////////////

void
CommandLineParsing::
parseSequences(const std::string & paramName,
					const std::string& paramArg,
					RnaSequenceVec& sequences )
{
	// clear sequence container
	sequences.clear();

	// read FASTA from STDIN stream
	if (paramArg.compare("STDIN") == 0) {
		parseSequencesFasta(paramName, std::cin, sequences);
	} else
	if (RnaSequence::isValidSequence(paramArg)) {
		// direct sequence input
		sequences.push_back(RnaSequence(paramName,paramArg));
	} else
	{
		// open file handle
		std::ifstream infile(paramArg);
		try {
			if(!infile.good()){
				LOG(ERROR) <<"FASTA parsing of "<<paramName<<" : could not open FASTA file  '"<<paramArg<<"'";
				updateParsingCode( ReturnCode::STOP_PARSING_ERROR );
			} else {
				parseSequencesFasta(paramName, infile, sequences);
			}
		} catch (std::exception & ex) {
			LOG(ERROR) <<"error while FASTA parsing of "<<paramName<<" : "<<ex.what();
			updateParsingCode( ReturnCode::STOP_PARSING_ERROR );
		}
		// close stream
		infile.close();
	}

	// holds current validation status to supress checks once a validation failed
	bool valid = true;

	// ensure at least one sequence was parsed
	valid = valid && validateSequenceNumber(paramName, sequences, 1, 9999);
	// validate alphabet
	valid = valid && validateSequenceAlphabet(paramName, sequences);

}


////////////////////////////////////////////////////////////////////////////

void
CommandLineParsing::
parseSequencesFasta( const std::string & paramName,
					std::istream& input,
					RnaSequenceVec& sequences)
{
	// temporary variables
	std::string line, name, sequence;
	int trimStart = 0;

	// read linewise
	while( ! std::getline( input, line ).eof() ) {
		// ignore empty lines
		if( line.empty() ) {
			continue;
		}
		// check of ID marker
		if ( line[0] == '>' ){ // Identifier marker
			if( !name.empty() ){ // we had read a name before
				// store last sequence
				// check if data complete
				if (sequence.empty())  {
					LOG(ERROR) <<"FASTA parsing of "<<paramName<<" : no sequence for ID '"<<name<<"'";
					updateParsingCode( ReturnCode::STOP_PARSING_ERROR );
				} else {
					// store sequence
					sequences.push_back( RnaSequence( name, sequence ) );
				}
				// clear name data
				name.clear();
			}
			// read new name
			if( !line.empty() ){
				// trim leading '>' plus successive and trailing whitespaces
				trimStart = line.find_first_not_of(" \t",1);
				line = line.substr( trimStart, std::max(0,(int)line.find_last_not_of(" \t\n\r")+1-trimStart) );
				name = line;
			}
			// clear sequence data
			sequence.clear();
		} else
		// has to be a sequence
		if( !name.empty() ){
			// trim leading/trailing whitespaces
			trimStart = line.find_first_not_of(" \t");
			line = line.substr( trimStart, std::max(0,(int)line.find_last_not_of(" \t\n\r")+1-trimStart) );
			// check for enclosed whitespaces
			if( line.find(' ') != std::string::npos ){ // Invalid sequence--no spaces allowed
				LOG(ERROR) <<"FASTA parsing of "<<paramName<<" : sequence for ID '"<<name<<"' contains spaces";
				updateParsingCode( ReturnCode::STOP_PARSING_ERROR );
				name.clear();
				sequence.clear();
			} else {
				// store sequence line
				sequence += line;
			}
		} else {
			LOG(ERROR) <<"FASTA parsing of "<<paramName<<" : found sequence data without leading ID";
			updateParsingCode( ReturnCode::STOP_PARSING_ERROR );
		}
	}
	// store last sequence
	// check if data complete
	if (sequence.empty())  {
		LOG(ERROR) <<"FASTA parsing of "<<paramName<<" : no sequence for ID '"<<name<<"'";
		updateParsingCode( ReturnCode::STOP_PARSING_ERROR );
	} else {
		// store sequence
		sequences.push_back( RnaSequence( name, sequence ) );
	}
}

////////////////////////////////////////////////////////////////////////////

bool
CommandLineParsing::
validateSequenceNumber( const std::string& paramName,
				const RnaSequenceVec& sequences,
				const size_t min,
				const size_t max)
{
	// check sequence number boundaries
	if (sequences.size() < min) {
		LOG(ERROR) <<""<<paramName<<" requires at least "<<min <<" sequences!";
		updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
		return false;
	} else
	if (sequences.size() > max) {
		LOG(ERROR) <<""<<paramName<<" allows at most "<<max <<" sequences!";
		updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
		return false;
	}
	// boundaries are met
	return true;
}

////////////////////////////////////////////////////////////////////////////

bool
CommandLineParsing::
validateSequenceAlphabet( const std::string& paramName,
				const RnaSequenceVec& sequences)
{
	bool allValid = true;
	// check each sequence
	for (int i=0; i<sequences.size(); i++) {
		// check if valid
		if (! RnaSequence::isValidSequence(sequences.at(i).asString())) {
			LOG(ERROR) <<"sequence " <<(i+1)<<" for parameter "<<paramName<<" is not valid!";
			updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
			allValid = false;
		}
	}
	return allValid;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

T_type
CommandLineParsing::
getTemperature() const {
	return temperature.val;
}

////////////////////////////////////////////////////////////////////////////

Predictor*
CommandLineParsing::
getPredictor( const InteractionEnergy & energy, OutputHandler & output ) const
{
	if (noSeedRequired) {
		switch ( predMode.val ) {
		case 1 :  return new PredictorMfe2d( energy, output );
		case 2 :  return new PredictorMfe4d( energy, output );
		case 3 :  return new PredictorMaxProb( energy, output );
		default: throw std::runtime_error("CommandLineParsing::getPredictor() : unknown predMode value "+toString(predMode.val));
		}
	} else {
		switch ( predMode.val ) {
		case 1 :  return new PredictorMfe2dSeed( energy, output, getSeedConstraint() );
		case 2 :  NOTIMPLEMENTED("mode "+toString(predMode.val)+" not implemented for seed constraint (try --noSeed)"); return NULL;
		case 3 :  NOTIMPLEMENTED("mode "+toString(predMode.val)+" not implemented for seed constraint (try --noSeed)"); return NULL;
		default: throw std::runtime_error("CommandLineParsing::getPredictor() : unknown predMode value "+toString(predMode.val));
		}
	}
}

////////////////////////////////////////////////////////////////////////////

OutputHandler*
CommandLineParsing::
getOutputHandler( const InteractionEnergy & energy ) const
{
	// TODO add according arguments and parsing
	return new OutputHandlerText(std::cout, energy );
}

////////////////////////////////////////////////////////////////////////////

void
CommandLineParsing::
updateParsingCode( const ReturnCode currentParsingCode )
{
	parsingCode = std::max( parsingCode, currentParsingCode );
}

////////////////////////////////////////////////////////////////////////////

const SeedConstraint &
CommandLineParsing::
getSeedConstraint() const
{
	if (seedConstraint == NULL) {
		// setup according to user data
		seedConstraint = new SeedConstraint(
							  seedBP.val
							, seedMaxUP.val
							, seedMaxUPq.val<0 ? seedMaxUP.val : seedMaxUPq.val
							, seedMaxUPt.val<0 ? seedMaxUP.val : seedMaxUPt.val
							, seedMaxE.val
						);
	}
	return *seedConstraint;
}

////////////////////////////////////////////////////////////////////////////

const IndexRangeList&
CommandLineParsing::
getQueryRanges( const size_t sequenceNumber ) const
{
#if IN_DEBUG_MODE
	if (sequenceNumber>=qRegion.size())
		throw std::runtime_error("CommandLineParsing::getQueryRanges("+toString(sequenceNumber)+") out of bounds");
	if (qRegion.at(sequenceNumber).empty())
		throw std::runtime_error("CommandLineParsing::getQueryRanges("+toString(sequenceNumber)+") is empty");
#endif
	return qRegion.at(sequenceNumber);
}

////////////////////////////////////////////////////////////////////////////

const IndexRangeList&
CommandLineParsing::
getTargetRanges( const size_t sequenceNumber ) const
{
#if IN_DEBUG_MODE
	if (sequenceNumber>=tRegion.size())
		throw std::runtime_error("CommandLineParsing::getTargetRanges("+toString(sequenceNumber)+") out of bounds");
	if (tRegion.at(sequenceNumber).empty())
		throw std::runtime_error("CommandLineParsing::getTargetRanges("+toString(sequenceNumber)+") is empty");
#endif
	return tRegion.at(sequenceNumber);
}

////////////////////////////////////////////////////////////////////////////



