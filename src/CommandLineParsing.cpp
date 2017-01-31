
#include "CommandLineParsing.h"

#include "general.h"

#include <cmath>
#include <stdexcept>
#include <fstream>

#include <omp.h>

#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

#include "AccessibilityConstraint.h"

#include "AccessibilityDisabled.h"
#include "AccessibilityFromStream.h"
#include "AccessibilityVrna.h"

#include "InteractionEnergyBasePair.h"
#include "InteractionEnergyVrna.h"

#include "PredictorMfe2dHeuristic.h"
#include "PredictorMfe2d.h"
#include "PredictorMfe4d.h"
#include "PredictorMaxProb.h"

#include "PredictorMfe2dHeuristicSeed.h"
#include "PredictorMfe2dSeed.h"
#include "PredictorMfe4dSeed.h"

// TODO maxprob seed

#include "OutputHandlerText.h"
#include "OutputHandlerCsv.h"
#include "OutputHandlerIntaRNA1.h"





////////////////////////////////////////////////////////////////////////////

const std::string CommandLineParsing::outCsvCols_default = "id1,start1,end1,id2,start2,end2,subseqDP,hybridDP,E";

////////////////////////////////////////////////////////////////////////////

CommandLineParsing::CommandLineParsing()
	:
	stdinUsed(false),
	opts_query("Query"),
	opts_target("Target"),
	opts_seed("Seed"),
	opts_inter("Interaction"),
	opts_general("General"),
	opts_output("Output"),
	opts_cmdline_all(),
	opts_cmdline_short(),

	parsingCode(NOT_PARSED_YET),

	queryArg(""),
	query(),
	qAcc("NCPE",'C'),
	qAccW( 0, 99999, 150),
	qAccL( 0, 99999, 100),
	qAccConstr(""),
	qAccFile(""),
	qIntLenMax( 0, 99999, 0),
	qIntLoopMax( 0, 30, 16),
	qRegionString(""),
	qRegion(),

	targetArg(""),
	target(),
	tAcc("NCPE",'C'),
	tAccW( 0, 99999, 150),
	tAccL( 0, 99999, 100),
	tAccConstr(""),
	tAccFile(""),
	tIntLenMax( 0, 99999, 0),
	tIntLoopMax( 0, 30, 16),
	tRegionString(""),
	tRegion(),

	noSeedRequired(false),
	seedBP(2,20,7),
	seedMaxUP(0,20,0),
	seedMaxUPq(-1,20,-1),
	seedMaxUPt(-1,20,-1),
	seedMaxE(-999,+999,0),
	seedMinPu(0,1,0),
	seedRangeq(""),
	seedRanget(""),
	seedConstraint(NULL),

	temperature(0,100,37),

	predMode( PredictionMode_min, PredictionMode_max, HEURISTIC),
	threads( 1, omp_get_max_threads(), 1),

	energy("BF",'F'),
	energyFile(""),

	out("STDOUT"),
	outStream(&(std::cout)),
	outMode( OutputMode_min, OutputMode_max, DETAILED ),
	outNumber( 0, 1000, 1),
	outOverlap( OutputConstraint::ReportOverlap::OVERLAP_NONE, OutputConstraint::ReportOverlap::OVERLAP_BOTH, OutputConstraint::ReportOverlap::OVERLAP_SEQ2 ),
	outDeltaE( 0.0, 100.0, 100.0),
	outMaxE( -999.0, 0.0, 0.0),
	outCsvCols(outCsvCols_default),
	outAccFileq(""),
	outAccFilet(""),
	outPuFileq(""),
	outPuFilet(""),

	vrnaHandler()

{
	using namespace boost::program_options;


	////  SEQUENCE OPTIONS  ////////////////////////////////////

	opts_query.add_options()
		("query,q"
			, value<std::string>(&queryArg)
				->required()
				->notifier(boost::bind(&CommandLineParsing::validate_query,this,_1))
			, "either an RNA sequence or the stream/file name from where to read the query sequences (should be the shorter sequences to increase efficiency); use 'STDIN' to read from standard input stream; sequences have to use IUPAC nucleotide encoding")
		("qAcc"
			, value<char>(&(qAcc.val))
				->default_value(qAcc.def)
				->notifier(boost::bind(&CommandLineParsing::validate_qAcc,this,_1))
			, std::string("accessibility computation : 'N'o accessibility contributions"
					", 'C' computation of accessibilities"
					", 'P' RNAplfold unpaired probability output from --qAccFile"
					", 'E' ED values in RNAplfold Pu-like format from --qAccFile"
					).c_str())
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
		;
	opts_cmdline_short.add(opts_query);
	opts_query.add_options()
		("qAccConstr"
			, value<std::string>(&(qAccConstr))
				->notifier(boost::bind(&CommandLineParsing::validate_qAccConstr,this,_1))
			, std::string("accessibility computation : structure constraint for each sequence position: '.' no constraint, '"+toString(AccessibilityConstraint::dotBracket_accessible)+"' unpaired, '"+toString(AccessibilityConstraint::dotBracket_blocked)+"' blocked. Note, blocked positions are excluded from interaction prediction and considered unpaired!").c_str())
		("qAccFile"
			, value<std::string>(&(qAccFile))
				->notifier(boost::bind(&CommandLineParsing::validate_qAccFile,this,_1))
			, std::string("accessibility computation : if --qAcc is to be read from file, the file/stream to be parsed. Used 'STDIN' if to read from standard input stream.").c_str())
		("qIntLenMax"
			, value<int>(&(qIntLenMax.val))
				->default_value(qIntLenMax.def)
				->notifier(boost::bind(&CommandLineParsing::validate_qIntLenMax,this,_1))
			, std::string("interaction site : maximal window size to be considered"
					" for interaction (and thus accessible) within the query"
					" (arg in range ["+toString(qIntLenMax.min)+","+toString(qIntLenMax.max)+"];"
					" 0 defaults to the full sequence length)"
					" If --qAccW is provided, the smaller window size of both is used."
					).c_str())
		("qIntLoopMax"
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
				, "either an RNA sequence or the stream/file name from where to read the target sequences (should be the longer sequences to increase efficiency); use 'STDIN' to read from standard input stream; sequences have to use IUPAC nucleotide encoding")
		("tAcc"
			, value<char>(&(tAcc.val))
				->default_value(tAcc.def)
				->notifier(boost::bind(&CommandLineParsing::validate_tAcc,this,_1))
			, std::string("accessibility computation : 'N'o accessibility contributions"
					", 'C' computation of accessibilities"
					", 'P' RNAplfold unpaired probability output from --qAccFile"
					", 'E' ED values in RNAplfold Pu-like format from --qAccFile"
					).c_str())
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
		;
	opts_cmdline_short.add(opts_target);
	opts_target.add_options()
		("tAccConstr"
			, value<std::string>(&(tAccConstr))
				->notifier(boost::bind(&CommandLineParsing::validate_tAccConstr,this,_1))
			, std::string("accessibility computation : structure constraint for each sequence position: '.' no constraint, '"+toString(AccessibilityConstraint::dotBracket_accessible)+"' unpaired, '"+toString(AccessibilityConstraint::dotBracket_blocked)+"' blocked. Note, blocked positions are excluded from interaction prediction and considered unpaired!").c_str())
		("tAccFile"
			, value<std::string>(&(tAccFile))
				->notifier(boost::bind(&CommandLineParsing::validate_tAccFile,this,_1))
			, std::string("accessibility computation : if --tAcc is to be read from file, the file/stream to be parsed. Used 'STDIN' if to read from standard input stream.").c_str())
		("tIntLenMax"
			, value<int>(&(tIntLenMax.val))
				->default_value(tIntLenMax.def)
				->notifier(boost::bind(&CommandLineParsing::validate_tIntLenMax,this,_1))
			, std::string("interaction site : maximal window size to be considered for"
					" interaction (and thus accessible) within the target"
					" (arg in range ["+toString(tIntLenMax.min)+","+toString(tIntLenMax.max)+"];"
					" 0 defaults to the full sequence length)."
					" If --tAccW is provided, the smaller window size of both is used."
					).c_str())
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

		("seedBP"
			, value<int>(&(seedBP.val))
				->default_value(seedBP.def)
				->notifier(boost::bind(&CommandLineParsing::validate_seedBP,this,_1))
			, std::string("number of inter-molecular base pairs within the seed region (arg in range ["+toString(seedBP.min)+","+toString(seedBP.max)+"])").c_str())
		("seedMaxUP"
			, value<int>(&(seedMaxUP.val))
				->default_value(seedMaxUP.def)
				->notifier(boost::bind(&CommandLineParsing::validate_seedMaxUP,this,_1))
			, std::string("maximal overall number (query+target) of unpaired bases within the seed region (arg in range ["+toString(seedMaxUP.min)+","+toString(seedMaxUP.max)+"])").c_str())
		;
	opts_cmdline_short.add(opts_seed);
	opts_seed.add_options()
		("seedMaxUPq"
			, value<int>(&(seedMaxUPq.val))
				->default_value(seedMaxUPq.def)
				->notifier(boost::bind(&CommandLineParsing::validate_seedMaxUPq,this,_1))
			, std::string("maximal number of unpaired bases within the query's seed region (arg in range ["+toString(seedMaxUPq.min)+","+toString(seedMaxUPq.max)+"]); if -1 the value of seedMaxUP is used.").c_str())
		("seedMaxUPt"
			, value<int>(&(seedMaxUPt.val))
				->default_value(seedMaxUPt.def)
				->notifier(boost::bind(&CommandLineParsing::validate_seedMaxUPt,this,_1))
			, std::string("maximal number of unpaired bases within the target's seed region (arg in range ["+toString(seedMaxUPt.min)+","+toString(seedMaxUPt.max)+"]); if -1 the value of seedMaxUP is used.").c_str())
		("seedMaxE"
			, value<E_type>(&(seedMaxE.val))
				->default_value(seedMaxE.def)
				->notifier(boost::bind(&CommandLineParsing::validate_seedMaxE,this,_1))
			, std::string("maximal energy a seed region may have (arg in range ["+toString(seedMaxE.min)+","+toString(seedMaxE.max)+"]).").c_str())
		("seedMinPu"
			, value<E_type>(&(seedMinPu.val))
				->default_value(seedMinPu.def)
				->notifier(boost::bind(&CommandLineParsing::validate_seedMinPu,this,_1))
			, std::string("minimal unpaired probability (per sequence) a seed region may have (arg in range ["+toString(seedMinPu.min)+","+toString(seedMinPu.max)+"]).").c_str())
		("seedRangeq"
			, value<std::string>(&(seedRangeq))
				->notifier(boost::bind(&CommandLineParsing::validate_seedRangeq,this,_1))
			, std::string("interval(s) in the query to search for seeds in format 'from1-to1,from2-to2,...' (Note, only for single query)").c_str())
		("seedRanget"
			, value<std::string>(&(seedRanget))
				->notifier(boost::bind(&CommandLineParsing::validate_seedRanget,this,_1))
			, std::string("interval(s) in the target to search for seeds in format 'from1-to1,from2-to2,...' (Note, only for single target)").c_str())
		;

	////  INTERACTION/ENERGY OPTIONS  ////////////////////////

	opts_inter.add_options()
		("mode,m"
			, value<int>(&(predMode.val))
				->default_value(predMode.def)
				->notifier(boost::bind(&CommandLineParsing::validate_predMode,this,_1))
			, std::string("prediction mode : "
					+toString(PredictionMode::HEURISTIC)+"= MFE-heuristic (n^2:S|T), "
					+toString(PredictionMode::SPACEEFFICIENT)+"= MFE (n^2:S,n^4:T), "
					+toString(PredictionMode::FULL)+"= MFE (n^4:S|T), "
					+toString(PredictionMode::MAXPROB)+"= MaxProb (n^4:S|T)").c_str())
		;
	opts_cmdline_short.add(opts_inter);
	opts_inter.add_options()
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


	////  OUTPUT OPTIONS  ////////////////////////////////////

	opts_output.add_options()
		("out"
			, value<std::string>(&(out))
				->default_value(out)
				->notifier(boost::bind(&CommandLineParsing::validate_out,this,_1))
			, std::string("output : provide a file name for output or 'STDOUT/STDERR' to write to the according stream").c_str())
		("outMode"
			, value<int>(&(outMode.val))
				->default_value(outMode.def)
				->notifier(boost::bind(&CommandLineParsing::validate_outMode,this,_1))
			, std::string("output mode : "
					+toString(OutputMode::DETAILED)+"= detailed, "
					+toString(OutputMode::V1_DETAILED)+"= v1-detailed"
					+toString(OutputMode::CSV)+"= CSV"
					).c_str())
	    ("outNumber,n"
			, value<int>(&(outNumber.val))
				->default_value(outNumber.def)
				->notifier(boost::bind(&CommandLineParsing::validate_outNumber,this,_1))
			, std::string("maximal overall number (query+target) of unpaired bases within the seed region (arg in range ["+toString(outNumber.min)+","+toString(outNumber.max)+"])").c_str())
	    ("outOverlap"
			, value<int>(&(outOverlap.val))
				->default_value(outOverlap.def)
				->notifier(boost::bind(&CommandLineParsing::validate_outOverlap,this,_1))
			, std::string("suboptimal output : interactions can overlap ("
					+toString(OutputConstraint::OVERLAP_NONE)+") in none of the sequences, ("
					+toString(OutputConstraint::OVERLAP_SEQ1)+") in the target only, ("
					+toString(OutputConstraint::OVERLAP_SEQ2)+") in the query only, ("
					+toString(OutputConstraint::OVERLAP_BOTH)+") in both sequences").c_str())
		;
	opts_cmdline_short.add(opts_output);
	opts_output.add_options()
	    ("outMaxE"
			, value<double>(&(outMaxE.val))
				->default_value(outMaxE.def)
				->notifier(boost::bind(&CommandLineParsing::validate_outMaxE,this,_1))
			, std::string("only interactions with E <= maxE are reported").c_str())
	    ("outDeltaE"
			, value<double>(&(outDeltaE.val))
				->default_value(outDeltaE.def)
				->notifier(boost::bind(&CommandLineParsing::validate_outDeltaE,this,_1))
			, std::string("suboptimal output : only interactions with E <= (minE+deltaE) are reported").c_str())
		("outCsvCols"
			, value<std::string>(&(outCsvCols))
				->default_value(outCsvCols)
				->notifier(boost::bind(&CommandLineParsing::validate_outCsvCols,this,_1))
			, std::string("output : comma separated list of CSV column IDs to print if outMode=CSV."
					" An empty argument prints all possible columns from the following available ID list: "
					+ boost::replace_all_copy(OutputHandlerCsv::list2string(OutputHandlerCsv::string2list("")), ",", ", ")
					).c_str())
		("outAccFileq"
			, value<std::string>(&(outAccFileq))
				->notifier(boost::bind(&CommandLineParsing::validate_outAccFileq,this,_1))
			, std::string("output : writes the query's ED values to the given file/stream"
					" in a format similar to RNAplfold unpaired probability output."
					" Use STDOUT/STDERR to write to the respective output stream."
					).c_str())
		("outAccFilet"
			, value<std::string>(&(outAccFilet))
				->notifier(boost::bind(&CommandLineParsing::validate_outAccFilet,this,_1))
			, std::string("output : writes the target's ED values to the given file/stream"
					" in a format similar to RNAplfold unpaired probability output."
					" Use STDOUT/STDERR to write to the respective output stream."
					).c_str())
		("outPuFileq"
			, value<std::string>(&(outPuFileq))
				->notifier(boost::bind(&CommandLineParsing::validate_outPuFileq,this,_1))
			, std::string("output : writes the query's unpaired probabilities used for ED values to the given file/stream"
					" in RNAplfold unpaired probability output format."
					" Use STDOUT/STDERR to write to the respective output stream."
					).c_str())
		("outPuFilet"
			, value<std::string>(&(outPuFilet))
				->notifier(boost::bind(&CommandLineParsing::validate_outPuFilet,this,_1))
			, std::string("output : writes the target's unpaired probabilities used for ED values to the given file/stream"
					" in RNAplfold unpaired probability output format."
					" Use STDOUT/STDERR to write to the respective output stream."
					).c_str())
	    ("verbose,v", "verbose output") // handled via easylogging++
//	    (logFile_argument.c_str(), "name of log file to be used for output")
	    ;

	////  GENERAL OPTIONS  ////////////////////////////////////

	opts_general.add_options()
	    ("threads"
			, value<int>(&(threads.val))
				->default_value(threads.def)
				->notifier(boost::bind(&CommandLineParsing::validate_threads,this,_1))
			, std::string("maximal number of threads to be used for parallel computation of query-target-combinations."
					" Note, the number of threads multiplies the required memory needed for computation!"
					" (arg in range ["+toString(threads.min)+","+toString(threads.max)+"])").c_str())
	    ("version", "print version")
	    ("help,h", "show the help page for basic parameters")
	    ("fullhelp", "show the extended help page for all available parameters")
	    ;
	opts_cmdline_short.add(opts_general);

	////  GENERAL OPTIONS  ////////////////////////////////////

	opts_cmdline_all.add(opts_query).add(opts_target).add(opts_seed).add(opts_inter).add(opts_output).add(opts_general);


}

////////////////////////////////////////////////////////////////////////////

CommandLineParsing::~CommandLineParsing() {

	CLEANUP(seedConstraint);

	if (outStream != &std::cout && outStream != &std::cerr) {
		std::fstream *outFileStream = dynamic_cast<std::fstream*>(outStream);
		assert(outFileStream != NULL);
		// flush and close file stream
		outFileStream->flush();
		outFileStream->close();
		// delete file handler
		CLEANUP(outFileStream);
	}
	// reset output stream
	outStream = & std::cout;

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
		int parseStyle =
				  command_line_style::style_t::allow_long
				| command_line_style::style_t::long_allow_adjacent
				| command_line_style::style_t::long_allow_next
				| command_line_style::style_t::allow_short
				| command_line_style::style_t::allow_dash_for_short
				| command_line_style::style_t::short_allow_next
				| command_line_style::style_t::case_insensitive
				;
		store( parse_command_line(argc, argv, opts_cmdline_all, parseStyle), vm);
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
				<<"\nThe following basic program arguments are supported:\n"
				<< opts_cmdline_short
				<< "\n"
				<< "Run --fullhelp for the extended list of parameters\n"
				<< "\n";
			parsingCode = ReturnCode::STOP_ALL_FINE;
			return parsingCode;
		}
		if (vm.count("fullhelp")) {
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

			// open output stream
			{
				// get output selection upper case
				std::string outUpperCase = boost::to_upper_copy<std::string>(out,std::locale());
				// check if standard stream
				if (boost::iequals(out,"STDOUT")) {
					outStream = & std::cout;
				} else
				if (boost::iequals(out,"STDERR")) {
					outStream = & std::cerr;
				} else {
					// open file stream
					std::fstream * outFileStream = new std::fstream();
					outFileStream->open( out.c_str(), std::ios_base::out );
					if (!outFileStream->is_open()) {
						delete outFileStream;
						LOG(ERROR) <<"could not open output file --out='"<<out << "' for writing";
						updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
					} else {
						// set output stream
						outStream = outFileStream;
					}
				}


			}

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
					NOTIMPLEMENTED("--qAccConstr only supported for single sequence input");
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
					NOTIMPLEMENTED("--tAccConstr only supported for single sequence input");
				}
			} else {
				// generate empty constraint
				tAccConstr = std::string(target.at(0).size(),'.');
			}

			// check sanity of accessibility setup
			switch(qAcc.val) {
			case 'C' : {
				if (!qAccFile.empty()) LOG(INFO) <<"qAcc = "<<qAcc.val<<" : ignoring --qAccFile";
				break;
			}
			case 'E' : // drop to next handling
			case 'P' : {
				if (qAccFile.empty()) LOG(INFO) <<"qAcc = "<<qAcc.val<<" but no --qAccFile given";
				if (!qAccConstr.empty()) LOG(INFO) <<"qAcc = "<<qAcc.val<<" : accessibility constraints (--qAccConstr) possibly not used in computation of loaded ED values";
				if (getQuerySequences().size()>1) throw std::runtime_error("qAcc = "+toString(qAcc.val)+" only supported for single query sequence input");
			}	// drop to next handling
			case 'N' : {
				if (qAccL.val != qAccL.def) LOG(INFO) <<"qAcc = "<<qAcc.val<<" : ignoring --qAccL";
				if (qAccW.val != qAccW.def) LOG(INFO) <<"qAcc = "<<qAcc.val<<" : ignoring --qAccW";
				if (qAcc.val != 'N' && !qAccFile.empty()) LOG(INFO) <<"qAcc = "<<qAcc.val<<" : ignoring --qAccFile";
				break;
			}
			} // switch
			switch(tAcc.val) {
			case 'C' : {
				if (!tAccFile.empty()) LOG(INFO) <<"tAcc = "<<tAcc.val<<" : ignoring --tAccFile";
				break;
			}
			case 'E' : // drop to next handling
			case 'P' : {
				if (tAccFile.empty()) LOG(INFO) <<"tAcc = "<<tAcc.val<<" but no --tAccFile given";
				if (!tAccConstr.empty()) LOG(INFO) <<"tAcc = "<<tAcc.val<<" : accessibility constraints (--tAccConstr) possibly not used in computation of loaded ED values";
				if (getTargetSequences().size()>1) throw std::runtime_error("tAcc = "+toString(tAcc.val)+" only supported for single target sequence input");
			}	// drop to next handling
			case 'N' : {
				if (tAccL.val != tAccL.def) LOG(INFO) <<"tAcc = "<<tAcc.val<<" : ignoring --tAccL";
				if (tAccW.val != tAccW.def) LOG(INFO) <<"tAcc = "<<tAcc.val<<" : ignoring --tAccW";
				if (tAcc.val=='N' && !tAccFile.empty()) LOG(INFO) <<"tAcc = "<<tAcc.val<<" : ignoring --tAccFile";
				break;
			}
			} // switch

			// check energy setup
			if (vm.count("energyVRNA") > 0 && energy.val != 'F') {
				throw error("--energyVRNA provided but no (F)ull energy computation requested (--energy = "+toString(energy.val)+")");
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
				if (seedMinPu.val != seedMinPu.def) LOG(INFO) <<"no seed constraint wanted, but seedMinPu provided (will be ignored)";
				if (!seedRangeq.empty()) LOG(INFO) <<"no seed constraint wanted, but seedRangeq provided (will be ignored)";
				if (!seedRanget.empty()) LOG(INFO) <<"no seed constraint wanted, but seedRanget provided (will be ignored)";
			} else {
				// check query search ranges
				if (!seedRangeq.empty()) {
					if (query.size()!=1) {
						LOG(ERROR) <<"seedRangeq given but not only one query sequence provided";
						updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
					} else {
						validate_indexRangeList("seedRangeq",seedRangeq, 1, query.begin()->size());
					}
				}
				// check target search ranges
				if (!seedRanget.empty()) {
					if (target.size()!=1) {
						LOG(ERROR) <<"seedRanget given but not only one target sequence provided";
						updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
					} else {
						validate_indexRangeList("seedRanget",seedRanget, 1, target.begin()->size());
					}
				}
			}

			// check qAcc upper bound
			if (qAccL.val > qAccW.val && qAccW.val != 0) {
				LOG(ERROR) <<"qAccL = " <<qAccL.val <<" : has to be <= qAccW (=" <<qAccW.val<<")";
				updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
			}

			// check qAcc upper bound
			if (tAccL.val > tAccW.val && tAccW.val != 0) {
				LOG(ERROR) <<"tAccL = " <<tAccL.val <<" : has to be <= tAccW (=" <<tAccW.val<<")";
				updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
			}

			// check CSV stuff
			if (outCsvCols != outCsvCols_default && outMode.val != OutputMode::CSV) {
				LOG(ERROR) <<"outCsvCols set but outMode != "<<OutputMode::CSV;
				updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
			}

			// check ED/Pu output sanity
			if (!outAccFileq.empty() && getQuerySequences().size()>1) throw std::runtime_error("--outAccFileq only supported for single query sequence input");
			if (!outAccFilet.empty() && getTargetSequences().size()>1) throw std::runtime_error("--outAccFilet only supported for single target sequence input");
			if (!outPuFileq.empty() && getQuerySequences().size()>1) throw std::runtime_error("--outPuFileq only supported for single query sequence input");
			if (!outPuFilet.empty() && getTargetSequences().size()>1) throw std::runtime_error("--outPuFilet only supported for single target sequence input");

			// check if multi-threading
			if (threads.val > 1 && getTargetSequences().size() > 1) {
				// warn if 4D space prediction enabled
				if (predMode.val > PredictionMode::SPACEEFFICIENT) {
					LOG(WARNING) <<"Multi-threading enabled in high-mem-prediction mode : ensure you have enough memory available!";
				}
				if (outMode.val == OutputMode::V1_NORMAL || outMode.val == OutputMode::V1_DETAILED) {
					throw std::runtime_error("Multi-threading not supported for IntaRNA v1 output");
				}
			}

			// trigger initial output handler output
			initOutputHandler();

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
	if ( boost::iequals(value,"STDIN") ) {
		setStdinUsed();
	} else {
		// check if it is a sequence
		if ( RnaSequence::isValidSequenceIUPAC(value) ) {

			// do nothing, all fine

		} else { // it is NOT a sequence! can be a file name!

			// check if a file of this name exists and is readable
			if ( ! validateFile( value ) ) {
				LOG(ERROR) <<""<<name<<" '"<<value <<"' : is neither STDIN, a file name, or a sequence!";
				updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
			}
		}

	}

}


////////////////////////////////////////////////////////////////////////////

void
CommandLineParsing::
validate_indexRangeList(const std::string & argName, const std::string & value
		, const size_t indexMin , const size_t indexMax)
{
	assert(indexMin <= indexMax);

	// check if empty
	if (value.empty()) {
		return;
	}
	// check if matched by regex
	if ( ! boost::regex_match(value,IndexRangeList::regex, boost::match_perl) ) {
		LOG(ERROR) <<argName<<" : '"<<value<<"' does not match the format 'from1-to1,from2-to2,..'";
		updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
	} else {
		// create range list for further checking
		IndexRangeList r(value);
		for (IndexRangeList::const_iterator i=r.begin(); i!=r.end(); i++) {
			// ensure range is ascending
			if (!i->isAscending()) {
				LOG(ERROR)  <<argName<<" : subrange " <<*i <<" is not ascending";
				updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
				return;
			}
			// check if boundaries in range (given they are ascending)
			if (i->to < indexMin || i->to > indexMax) {
				LOG(ERROR)  <<argName<<" : subrange " <<*i <<" is out of bounds [1,"<<indexMax<<"]";
				updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
				return;
			}
		}
	} // matches regex
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

void CommandLineParsing::validate_outCsvCols(const std::string & value) {
	try {
		// try to parse
		OutputHandlerCsv::string2list( value );
	} catch (const std::runtime_error & e ) {
		LOG(ERROR) <<"--outCsvCols : " <<e.what();
		updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
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
	if (boost::regex_match( value, IndexRangeList::regex, boost::match_perl )) {
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
	if (boost::regex_match( value, IndexRangeList::regex, boost::match_perl )) {
		// ensure single sequence input
		if(sequences.size() != 1) {
			LOG(ERROR) <<argName <<" : string range list encoding provided but more than one sequence present.";
			updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
			return;
		}
		// validate range encodings
		validate_indexRangeList(argName, value, 1, sequences.begin()->size());
		// ensure range list size sufficient
		rangeList.resize(1);
		// fill range list from string but shift by -1
		rangeList[0] = IndexRangeList( value ).shift(-1, sequences.begin()->size()-1);
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

	// create temporary constraint object (will be copied)
	AccessibilityConstraint accConstraint(qAccConstr);
	bool computeES = false;
#if IN_DEBUG_MODE
	// TODO replace based on predictor selection
	computeES = qAccW.val==0;
#endif
	// construct selected accessibility object
	switch(qAcc.val) {

	case 'N' : // no accessibility
		return new AccessibilityDisabled( seq
										, qIntLenMax.val
										, &accConstraint );

	case 'E' : // drop to next handling
	case 'P' : { // VRNA RNAplfold unpaired probability file output
		std::istream * accStream = NULL;
		std::ifstream * accFileStream = NULL;
		if ( boost::iequals(qAccFile,"STDIN") ) {
			accStream = &(std::cin);
		} else {
			// file support
			accFileStream = new std::ifstream(qAccFile);
			try {
				if(!accFileStream->good()){
					accFileStream->close();
					CLEANUP(accFileStream);
					throw std::runtime_error("accessibility parsing of --qAccFile : could not open file '"+qAccFile+"'");
				}
			} catch (std::exception & ex) {
				accFileStream->close();
				CLEANUP(accFileStream);
				throw std::runtime_error("accessibility parsing of --qAccFile : error while opening '"+qAccFile+"' : "+ex.what());
			}
			// set file stream as input stream
			accStream = accFileStream;
		}
		Accessibility * acc = new AccessibilityFromStream( seq
										, qIntLenMax.val
										, &accConstraint
										, *accStream
										, (qAcc.val == 'P' ? AccessibilityFromStream::Pu_RNAplfold_Text : AccessibilityFromStream::ED_RNAplfold_Text)
										, vrnaHandler.getRT() );
		// cleanup
		if ( accFileStream != NULL ) {
			accFileStream->close();
			CLEANUP( accFileStream );
		}
		return acc;
	}

	case 'C' : // compute VRNA-based accessibilities
		return new AccessibilityVrna( seq
									, std::min( qIntLenMax.val == 0 ? seq.size() : qIntLenMax.val
												, qAccW.val == 0 ? seq.size() : qAccW.val )
									, &accConstraint
									, vrnaHandler
									, qAccW.val
									, qAccL.val
									, computeES);
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
	// create temporary constraint object (will be copied)
	AccessibilityConstraint accConstraint(tAccConstr);
	bool computeES = false;
#if IN_DEBUG_MODE
	// TODO replace based on predictor selection
	computeES = tAccW.val==0;
#endif
	const RnaSequence& seq = getTargetSequences().at(sequenceNumber);
	switch(tAcc.val) {

	case 'N' : // no accessibility
		return new AccessibilityDisabled( seq
										, tIntLenMax.val
										, &accConstraint );

	case 'E' : // drop to next handling
	case 'P' : { // VRNA RNAplfold unpaired probability file output
		std::istream * accStream = NULL;
		std::ifstream * accFileStream = NULL;
		// select stream to read from
		if ( boost::iequals(tAccFile,"STDIN") ) {
			accStream = &(std::cin);
		} else {
			// file support
			accFileStream = new std::ifstream(tAccFile);
			try {
				if(!accFileStream->good()){
					accFileStream->close();
					CLEANUP(accFileStream);
					throw std::runtime_error("accessibility parsing of --tAccFile : could not open file '"+tAccFile+"'");
				}
			} catch (std::exception & ex) {
				accFileStream->close();
				CLEANUP(accFileStream);
				throw std::runtime_error("accessibility parsing of --tAccFile : error while opening '"+tAccFile+"' : "+ex.what());
			}
			// set file stream as input stream
			accStream = accFileStream;
		}
		// read data
		Accessibility * acc = new AccessibilityFromStream( seq
										, tIntLenMax.val
										, &accConstraint
										, *accStream
										, ( tAcc.val == 'P' ? AccessibilityFromStream::Pu_RNAplfold_Text : AccessibilityFromStream::ED_RNAplfold_Text )
										, vrnaHandler.getRT() );
		// cleanup
		if ( accFileStream != NULL ) {
			accFileStream->close();
			CLEANUP( accFileStream );
		}
		return acc;
	}
	case 'C' : // compute VRNA-based accessibilities
		return new AccessibilityVrna( seq
									, std::min( tIntLenMax.val == 0 ? seq.size() : tIntLenMax.val
											, tAccW.val == 0 ? seq.size() : tAccW.val )
									, &accConstraint
									, vrnaHandler
									, tAccW.val
									, tAccL.val
									, computeES);
	default :
		NOTIMPLEMENTED("CommandLineParsing::getTargetAccessibility : tAcc = '"+toString(tAcc.val)+"' is not supported");
	}
	return NULL;
}

////////////////////////////////////////////////////////////////////////////

InteractionEnergy*
CommandLineParsing::
getEnergyHandler( const Accessibility& accTarget, const ReverseAccessibility& accQuery ) const
{
	checkIfParsed();
	switch( energy.val ) {
	case 'B' : return new InteractionEnergyBasePair( accTarget, accQuery, tIntLoopMax.val, qIntLoopMax.val );
	case 'F' : return new InteractionEnergyVrna( accTarget, accQuery, vrnaHandler, tIntLoopMax.val, qIntLoopMax.val );
	default :
		NOTIMPLEMENTED("CommandLineParsing::getEnergyHandler : energy = '"+toString(energy.val)+"' is not supported");
	}
	return NULL;

}

////////////////////////////////////////////////////////////////////////////

OutputConstraint
CommandLineParsing::
getOutputConstraint()  const
{
	checkIfParsed();
	return OutputConstraint(
			  outNumber.val
			, static_cast<OutputConstraint::ReportOverlap>(outOverlap.val)
			, static_cast<E_type>(outMaxE.val)
			, static_cast<E_type>(outDeltaE.val)
			);
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
	if (boost::iequals(paramArg,"STDIN")) {
		parseSequencesFasta(paramName, std::cin, sequences);
	} else
	if (RnaSequence::isValidSequenceIUPAC(paramArg)) {
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
	while( std::getline( input, line ) ) {
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
		if (! RnaSequence::isValidSequenceIUPAC(sequences.at(i).asString())) {
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
		switch ( (PredictionMode)predMode.val ) {
		case HEURISTIC :  return new PredictorMfe2dHeuristic( energy, output );
		case SPACEEFFICIENT :  return new PredictorMfe2d( energy, output );
		case FULL :  return new PredictorMfe4d( energy, output );
		case MAXPROB :  return new PredictorMaxProb( energy, output );
		}
	} else {
		switch ( predMode.val ) {
		case HEURISTIC :  return new PredictorMfe2dHeuristicSeed( energy, output, getSeedConstraint( energy ) );
		case SPACEEFFICIENT :  return new PredictorMfe2dSeed( energy, output, getSeedConstraint( energy ) );
		case FULL :  return new PredictorMfe4dSeed( energy, output, getSeedConstraint( energy ) );
		case MAXPROB :  NOTIMPLEMENTED("mode "+toString(predMode.val)+" not implemented for seed constraint (try --noSeed)"); return NULL;
		}
	}
}

////////////////////////////////////////////////////////////////////////////

std::ostream &
CommandLineParsing::
getOutputStream() const
{
	return *outStream;
}

////////////////////////////////////////////////////////////////////////////

void
CommandLineParsing::
initOutputHandler()
{
	// check if output mode == IntaRNA1-detailed
	switch ((OutputMode)outMode.val) {
	case V1_DETAILED :
		getOutputStream()
		<<"-------------------------" <<"\n"
		<<"INPUT" <<"\n"
		<<"-------------------------" <<"\n"
		<<"number of base pairs in seed                                  : "<<seedBP.val <<"\n"
		<<"max. number of unpaired bases in the seed region of seq. 1    : "<<(seedMaxUPt.val<0 ? seedMaxUP.val : seedMaxUPt.val) <<"\n"
		<<"max. number of unpaired bases in the seed region of seq. 2    : "<<(seedMaxUPq.val<0 ? seedMaxUP.val : seedMaxUPq.val) <<"\n"
		<<"max. number of unpaired bases in the seed region of both seq's: "<<seedMaxUP.val <<"\n"
		<<"RNAup used                                                    : "<<((qAccConstr.empty() && (qAccW.val ==0))?"true":"false") <<"\n"
		<<"RNAplfold used                                                : "<<(((! tAccConstr.empty()) || (tAccW.val !=0))?"true":"false") <<"\n"
		<<"sliding window size                                           : "<<tAccW.val<<"\n" //(tAccW.val!=0?tAccW.val:energy.size1()) <<"\n"
		<<"max. length of unpaired region                                : "<<tAccW.val<<"\n" //(tAccW.val!=0?tAccW.val:energy.size1()) <<"\n"
		<<"max. distance of two paired bases                             : "<<tAccL.val<<"\n" //(tAccL.val!=0?tAccL.val:energy.size1()) <<"\n"
		<<"weight for ED values of target RNA in energy                  : 1" <<"\n"
		<<"weight for ED values of binding RNA in energy                 : 1" <<"\n"
		<<"temperature                                                   : "<<temperature.val <<" Celsius" <<"\n"
		<<"max. number of subopt. results                                : "<<(getOutputConstraint().reportMax-1) <<"\n"
		<<"Heuristic for hybridization end used                          : "<<((PredictionMode)predMode.val==HEURISTIC?"true":"false") <<"\n"
		<<"\n"
		<<"-------------------------" <<"\n"
		<<"OUTPUT" <<"\n"
		<<"-------------------------" <<"\n"
		; break;
	case CSV :
		getOutputStream()
		<<OutputHandlerCsv::getHeader( OutputHandlerCsv::string2list( outCsvCols ) )
		; break;
	}

}

////////////////////////////////////////////////////////////////////////////

OutputHandler*
CommandLineParsing::
getOutputHandler( const InteractionEnergy & energy ) const
{
	switch ((OutputMode)outMode.val) {
	case DETAILED :
		return new OutputHandlerText( getOutputStream(), energy );
	case CSV :
		return new OutputHandlerCsv( getOutputStream(), energy, OutputHandlerCsv::string2list( outCsvCols ));
	case V1_NORMAL :
		return new OutputHandlerIntaRNA1( getOutputStream(), energy, false );
	case V1_DETAILED :
		return new OutputHandlerIntaRNA1( getOutputStream(), energy, true );
	}
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
getSeedConstraint( const InteractionEnergy & energy ) const
{
	if (seedConstraint == NULL) {
		// setup according to user data
		seedConstraint = new SeedConstraint(
							  seedBP.val
							, seedMaxUP.val
							, seedMaxUPt.val<0 ? seedMaxUP.val : seedMaxUPt.val
							, seedMaxUPq.val<0 ? seedMaxUP.val : seedMaxUPq.val
							, seedMaxE.val
							, (seedMinPu.val>0 ? std::min<E_type>(Accessibility::ED_UPPER_BOUND, - energy.getRT() * std::log( seedMinPu.val )) : Accessibility::ED_UPPER_BOUND) // transform unpaired prob to ED value
							// shift ranges to start counting with 0
							, IndexRangeList( seedRanget ).shift(-1,energy.size1()-1)
							, IndexRangeList( seedRangeq ).shift(-1,energy.size2()-1).reverse(energy.size2())
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

void
CommandLineParsing::
writeAccessibility( const Accessibility& acc, const std::string fileOrStream, const bool writeED ) const
{
	if (fileOrStream.empty())
		return;

	// setup output stream
	std::ostream * out = NULL;
	std::fstream * outFile = NULL;
	if ( boost::iequals(fileOrStream,"STDOUT")) {
		out = &std::cout;
	} else
	if ( boost::iequals(fileOrStream,"STDERR")) {
		out = &std::cerr;
	} else {
		// open file
		outFile = new std::fstream();
		outFile->open( fileOrStream.c_str(), std::ios_base::out );
		if (!outFile->is_open()) {
			CLEANUP(outFile);
			throw std::runtime_error("could not open output file '"+fileOrStream +"' for "+(writeED?"accessibility":"unpaired probability")+" output");
		} else {
			out = outFile;
		}
	}

	// write data to stream
	if (writeED) {
		acc.writeRNAplfold_ED_text( *out );
	} else {
		acc.writeRNAplfold_Pu_text( *out, vrnaHandler.getRT() );
	}

	// clean up
	if (outFile != NULL) { outFile->close(); }
	CLEANUP(outFile);
}

////////////////////////////////////////////////////////////////////////////



