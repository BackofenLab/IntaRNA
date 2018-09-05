
#include "CommandLineParsing.h"

#include "IntaRNA/general.h"

#include <cmath>
#include <stdexcept>
#include <fstream>

#if INTARNA_MULITHREADING
	#include <omp.h>
#endif

#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/assign/list_of.hpp>


extern "C" {
	#include <ViennaRNA/vrna_config.h>
}

#define _STRINGIFY(s) #s
#define STRINGIFY(s) _STRINGIFY(s)


#include "IntaRNA/AccessibilityConstraint.h"

#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/AccessibilityFromStream.h"
#include "IntaRNA/AccessibilityVrna.h"
#include "IntaRNA/AccessibilityBasePair.h"

#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/InteractionEnergyVrna.h"

#include "IntaRNA/PredictorMfe2dHeuristic.h"
#include "IntaRNA/PredictorMfe2d.h"
#include "IntaRNA/PredictorMfe4d.h"
#include "IntaRNA/PredictorMaxProb.h"

#include "IntaRNA/PredictorMfe2dHeuristicSeed.h"
#include "IntaRNA/PredictorMfe2dSeed.h"
#include "IntaRNA/PredictorMfe4dSeed.h"

#include "IntaRNA/PredictionTracker.h"
#include "IntaRNA/PredictionTrackerHub.h"
#include "IntaRNA/PredictionTrackerPairMinE.h"
#include "IntaRNA/PredictionTrackerProfileMinE.h"
#include "IntaRNA/PredictionTrackerSpotProb.h"
#include "IntaRNA/PredictionTrackerProfileSpotProb.h"

#include "IntaRNA/SeedHandlerMfe.h"

#include "IntaRNA/OutputHandlerCsv.h"
#include "IntaRNA/OutputHandlerIntaRNA1.h"
#include "IntaRNA/OutputHandlerText.h"


using namespace IntaRNA;


////////////////////////////////////////////////////////////////////////////

const std::string CommandLineParsing::outCsvCols_default = "id1,start1,end1,id2,start2,end2,subseqDP,hybridDP,E";

////////////////////////////////////////////////////////////////////////////

CommandLineParsing::CommandLineParsing()
	:
	stdinUsed(false),
	opts_query("Query"),
	opts_target("Target"),
	opts_seed("Seed"),
	opts_shape("SHAPE"),
	opts_inter("Interaction"),
	opts_general("General"),
	opts_output("Output"),
	opts_cmdline_all(),
	opts_cmdline_short(),

	parsingCode(NOT_PARSED_YET),

	queryArg(""),
	query(),
	qSet(),
	qSetString(""),
	qAcc("NCPE",'C'),
	qAccW( 0, 99999, 150),
	qAccL( 0, 99999, 100),
	qAccConstr(""),
	qAccFile(""),
	qIntLenMax( 0, 99999, 0),
	qIntLoopMax( 0, 30, 16),
	qRegionString(""),
	qRegion(),
	qRegionLenMax( 0, 99999, 0),
	qShape(""),
	qShapeMethod("Zb0.89"),
	qShapeConversion("Os1.6i-2.29"),

	targetArg(""),
	target(),
	tSet(),
	tSetString(""),
	tAcc("NCPE",'C'),
	tAccW( 0, 99999, 150),
	tAccL( 0, 99999, 100),
	tAccConstr(""),
	tAccFile(""),
	tIntLenMax( 0, 99999, 0),
	tIntLoopMax( 0, 30, 16),
	tRegionString(""),
	tRegion(),
	tRegionLenMax( 0, 99999, 0),
	tShape(""),
	tShapeMethod("Zb0.89"),
	tShapeConversion("Os1.6i-2.29"),

	noSeedRequired(false),
	seedTQ(""),
	seedBP(2,20,7),
	seedMaxUP(0,20,0),
	seedQMaxUP(-1,20,-1),
	seedTMaxUP(-1,20,-1),
	seedMaxE(-999,+999,0),
	seedMinPu(0,1,0),
	seedQRange(""),
	seedTRange(""),
	seedConstraint(NULL),

	temperature(0,100,37),

	pred( "SP", 'S'),
	predMode( "HME", 'H'),
#if INTARNA_MULITHREADING
	threads( 0, omp_get_max_threads(), 1),
#endif
	windowWidth(0,99999,0),
	windowOverlap(0,99999,0),

	energy("BV",'V'),
	energyFile(""),

	out(),
	outPrefix2streamName(),
	outStream(&(std::cout)),
	outMode( "NDC1O", 'N' ),
	outNumber( 0, 1000, 1),
	outOverlap( "NTQB", 'Q' ),
	outDeltaE( 0.0, 100.0, 100.0),
	outMaxE( -999.0, +999.0, 0.0),
	outMinPu( 0.0, 1.0, 0.0),
	outCsvCols(outCsvCols_default),
	outPerRegion(false),
	outSpotProbSpots(""),

	logFileName(""),

	vrnaHandler()

{
	using namespace boost::program_options;

	////  REMAINING INITIALIZATIONS  /////////////////////////////////

	// register all allowed prefix codes for the --out argument
	outPrefix2streamName[OutPrefixCode::OP_EMPTY] = "STDOUT";
	for (int c=1; c<OutPrefixCode::OP_UNKNOWN; c++) {
		outPrefix2streamName[(OutPrefixCode)c] = "";
	}


	////  QUERY SEQUENCE OPTIONS  ////////////////////////////////////

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
			, std::string("accessibility computation :"
					"\n 'N' no accessibility contributions"
					"\n 'C' computation of accessibilities"
					"\n 'P' unpaired probabilities in RNAplfold format from --qAccFile"
					"\n 'E' ED values in RNAplfold Pu-like format from --qAccFile"
					).c_str())
		("qAccW"
			, value<int>(&(qAccW.val))
				->default_value(qAccW.def)
				->notifier(boost::bind(&CommandLineParsing::validate_qAccW,this,_1))
			, std::string("accessibility computation : sliding window size for query accessibility computation"
					" (arg in range ["+toString(qAccW.min)+","+toString(qAccW.max)+"];"
					" 0 will use to the full sequence length)."
					" Note, this also restricts the maximal interaction length (see --qIntLenMax)."
					).c_str())
		("qAccL"
			, value<int>(&(qAccL.val))
				->default_value(qAccL.def)
				->notifier(boost::bind(&CommandLineParsing::validate_qAccL,this,_1))
			, std::string("accessibility computation : maximal loop length (base pair span) for query accessibility computation"
					" (arg in range ["+toString(qAccL.min)+","+toString(qAccL.max)+"]; 0 will use to sliding window size 'qAccW')").c_str())
		;
	opts_cmdline_short.add(opts_query);
	opts_query.add_options()
		("qSet"
			, value<std::string>(&(qSetString))
				->notifier(boost::bind(&CommandLineParsing::validate_qSet,this,_1))
			, std::string("query subset : List of sequence indices to consider for prediction in the format 'from1-to1,from2-to2,..' assuming indexing starts with 1").c_str())
		("qAccConstr"
			, value<std::string>(&(qAccConstr))
				->notifier(boost::bind(&CommandLineParsing::validate_qAccConstr,this,_1))
			, std::string("accessibility computation : structure constraint :"
					"\n EITHER a string of query sequence length encoding for each position:"
					" '.' no constraint,"
					" '"+toString(AccessibilityConstraint::dotBracket_accessible)+"' unpaired,"
					" '"+toString(AccessibilityConstraint::dotBracket_paired)+"' paired (intramolecularly), or"
					" '"+toString(AccessibilityConstraint::dotBracket_blocked)+"' blocked."
					" Note, blocked positions are excluded from interaction prediction and constrained to be unpaired!"
					"\n OR an index range based encoding that is prefixed by the according constraint letter and a colon,"
					" e.g. 'b:3-4,33-40,p:1-2,12-20'"
					).c_str())
		("qAccFile"
			, value<std::string>(&(qAccFile))
			, std::string("accessibility computation : the file/stream to be parsed, if --qAcc is to be read from file. Used 'STDIN' if to read from standard input stream.").c_str())
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
			, std::string("interaction site : maximal number of unpaired bases between neighbored interacting bases to be considered in interactions within the query"
					" (arg in range ["+toString(qIntLoopMax.min)+","+toString(qIntLoopMax.max)+"]; 0 enforces stackings only)").c_str())
		("qRegion"
			, value<std::string>(&(qRegionString))
				->notifier(boost::bind(&CommandLineParsing::validate_qRegion,this,_1))
			, std::string("interaction site : query regions to be considered for"
					" interaction prediction. Either given as BED file (for"
					" multi-sequence FASTA input) or in the format"
					" 'from1-to1,from2-to2,..' assuming indexing starts with 1."
					" Consider '--qRegionLenMax' for automatic region setup for"
					" long sequences."
					).c_str())
		("qRegionLenMax"
			, value<int>(&(qRegionLenMax.val))
				->default_value(qRegionLenMax.def)
				->notifier(boost::bind(&CommandLineParsing::validate_qRegionLenMax,this,_1))
			, std::string("interaction site : maximal length of highly accessible regions"
					" to be automatically identified. To this end, most inaccessible regions"
					" of length '--seedBP' are iteratively removed from the available indices"
					" until only regions below the given maximal length or remaining."
					" (arg in range ["+toString(qRegionLenMax.min)+","+toString(qRegionLenMax.max)+"];"
					" 0 defaults to no automatic range detection)"
					).c_str())
		;

	////  TARGET SEQUENCE OPTIONS  ////////////////////////////////////

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
			, std::string("accessibility computation :"
					"\n 'N' no accessibility contributions"
					"\n 'C' computation of accessibilities"
					"\n 'P' unpaired probabilities in RNAplfold format from --tAccFile"
					"\n 'E' ED values in RNAplfold Pu-like format from --tAccFile"
					).c_str())
		("tAccW"
			, value<int>(&(tAccW.val))
				->default_value(tAccW.def)
				->notifier(boost::bind(&CommandLineParsing::validate_tAccW,this,_1))
			, std::string("accessibility computation : sliding window size for query accessibility computation"
					" (arg in range ["+toString(tAccW.min)+","+toString(tAccW.max)+"];"
					" 0 will use the full sequence length)"
					" Note, this also restricts the maximal interaction length (see --tIntLenMax)."
					).c_str())
		("tAccL"
			, value<int>(&(tAccL.val))
				->default_value(tAccL.def)
				->notifier(boost::bind(&CommandLineParsing::validate_tAccL,this,_1))
			, std::string("accessibility computation : maximal loop size (base pair span) for query accessibility computation"
					" (arg in range ["+toString(tAccL.min)+","+toString(tAccL.max)+"]; 0 will use the sliding window size 'tAccW')").c_str())
		;
	opts_cmdline_short.add(opts_target);
	opts_target.add_options()
		("tSet"
			, value<std::string>(&(tSetString))
				->notifier(boost::bind(&CommandLineParsing::validate_tSet,this,_1))
			, std::string("target subset : List of sequence indices to consider for prediction in the format 'from1-to1,from2-to2,..' assuming indexing starts with 1").c_str())
		("tAccConstr"
			, value<std::string>(&(tAccConstr))
				->notifier(boost::bind(&CommandLineParsing::validate_tAccConstr,this,_1))
			, std::string("accessibility computation : structure constraint :"
					"\n EITHER a string of target sequence length encoding for each position:"
					" '.' no constraint,"
					" '"+toString(AccessibilityConstraint::dotBracket_accessible)+"' unpaired,"
					" '"+toString(AccessibilityConstraint::dotBracket_paired)+"' paired (intramolecularly), or"
					" '"+toString(AccessibilityConstraint::dotBracket_blocked)+"' blocked."
					" Note, blocked positions are excluded from interaction prediction and constrained to be unpaired!"
					"\n OR an index range based encoding that is prefixed by the according constraint letter and a colon,"
					" e.g. 'b:3-4,33-40,p:1-2,12-20'"
					).c_str())
		("tAccFile"
			, value<std::string>(&(tAccFile))
			, std::string("accessibility computation : the file/stream to be parsed, if --tAcc is to be read from file. Used 'STDIN' if to read from standard input stream.").c_str())
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
			, std::string("interaction site : maximal number of unpaired bases between neighbored interacting bases to be considered in interactions within the target"
					" (arg in range ["+toString(tIntLoopMax.min)+","+toString(tIntLoopMax.max)+"]; 0 enforces stackings only)").c_str())
		("tRegion"
			, value<std::string>(&(tRegionString))
				->notifier(boost::bind(&CommandLineParsing::validate_tRegion,this,_1))
			, std::string("interaction site : target regions to be considered for"
					" interaction prediction. Either given as BED file (for"
					" multi-sequence FASTA input) or in the format"
					" 'from1-to1,from2-to2,..' assuming indexing starts with 1."
					" Consider '--tRegionLenMax' for automatic region setup for"
					" long sequences."
					).c_str())
		("tRegionLenMax"
			, value<int>(&(tRegionLenMax.val))
				->default_value(tRegionLenMax.def)
				->notifier(boost::bind(&CommandLineParsing::validate_tRegionLenMax,this,_1))
			, std::string("interaction site : maximal length of highly accessible regions"
					" to be automatically identified. To this end, most inaccessible regions"
					" of length '--seedBP' are iteratively removed from the available indices"
					" until only regions below the given maximal length or remaining."
					" (arg in range ["+toString(tRegionLenMax.min)+","+toString(tRegionLenMax.max)+"];"
					" 0 defaults to no automatic range detection)"
					).c_str())
		;

	////  SEED OPTIONS  ////////////////////////////////////


	opts_seed.add_options()
	    ("noSeed", "if present, no seed is enforced within the predicted interactions")

		("seedTQ"
			, value<std::string>(&(seedTQ))
				->notifier(boost::bind(&CommandLineParsing::validate_seedTQ,this,_1))
			, std::string("comma separated list of explicit seed base pair encoding(s) in the format startTbpsT&startQbpsQ, e.g. '3|||.|&7||.||', where startT/Q are the indices of the 5' seed ends in target/query sequence and 'bps' the dot-bar base pair encodings. This disables all other seed constraints and seed identification.").c_str())
		("seedBP"
			, value<int>(&(seedBP.val))
				->default_value(seedBP.def)
				->notifier(boost::bind(&CommandLineParsing::validate_seedBP,this,_1))
			, std::string("number of inter-molecular base pairs within the seed region"
					" (arg in range ["+toString(seedBP.min)+","+toString(seedBP.max)+"])").c_str())
		("seedMaxUP"
			, value<int>(&(seedMaxUP.val))
				->default_value(seedMaxUP.def)
				->notifier(boost::bind(&CommandLineParsing::validate_seedMaxUP,this,_1))
			, std::string("maximal overall number (query+target) of unpaired bases within the seed region"
					" (arg in range ["+toString(seedMaxUP.min)+","+toString(seedMaxUP.max)+"])").c_str())
		;
	opts_cmdline_short.add(opts_seed);
	opts_seed.add_options()
		("seedQMaxUP"
			, value<int>(&(seedQMaxUP.val))
				->default_value(seedQMaxUP.def)
				->notifier(boost::bind(&CommandLineParsing::validate_seedQMaxUP,this,_1))
			, std::string("maximal number of unpaired bases within the query's seed region"
					" (arg in range ["+toString(seedQMaxUP.min)+","+toString(seedQMaxUP.max)+"]); if -1 the value of seedMaxUP is used.").c_str())
		("seedTMaxUP"
			, value<int>(&(seedTMaxUP.val))
				->default_value(seedTMaxUP.def)
				->notifier(boost::bind(&CommandLineParsing::validate_seedTMaxUP,this,_1))
			, std::string("maximal number of unpaired bases within the target's seed region"
					" (arg in range ["+toString(seedTMaxUP.min)+","+toString(seedTMaxUP.max)+"]); if -1 the value of seedMaxUP is used.").c_str())
		("seedMaxE"
			, value<E_type>(&(seedMaxE.val))
				->default_value(seedMaxE.def)
				->notifier(boost::bind(&CommandLineParsing::validate_seedMaxE,this,_1))
			, std::string("maximal energy a seed region may have"
					" (arg in range ["+toString(seedMaxE.min)+","+toString(seedMaxE.max)+"]).").c_str())
		("seedMinPu"
			, value<E_type>(&(seedMinPu.val))
				->default_value(seedMinPu.def)
				->notifier(boost::bind(&CommandLineParsing::validate_seedMinPu,this,_1))
			, std::string("minimal unpaired probability (per sequence) a seed region may have"
					" (arg in range ["+toString(seedMinPu.min)+","+toString(seedMinPu.max)+"]).").c_str())
		("seedQRange"
			, value<std::string>(&(seedQRange))
				->notifier(boost::bind(&CommandLineParsing::validate_seedQRange,this,_1))
			, std::string("interval(s) in the query to search for seeds in format 'from1-to1,from2-to2,...' (Note, only for single query)").c_str())
		("seedTRange"
			, value<std::string>(&(seedTRange))
				->notifier(boost::bind(&CommandLineParsing::validate_seedTRange,this,_1))
			, std::string("interval(s) in the target to search for seeds in format 'from1-to1,from2-to2,...' (Note, only for single target)").c_str())
		;

	////  SHAPE OPTIONS  ////////////////////////

	opts_shape.add_options()
		("qShape"
			, value<std::string>(&qShape)
				->notifier(boost::bind(&CommandLineParsing::validate_qShape,this,_1))
			, "file name from where to read the query sequence's SHAPE reactivity data to guide its accessibility computation")
		("tShape"
			, value<std::string>(&tShape)
				->notifier(boost::bind(&CommandLineParsing::validate_tShape,this,_1))
			, "SHAPE: file name from where to read the target sequence's SHAPE reactivity data to guide its accessibility computation")
		("qShapeMethod"
				, value<std::string>(&qShapeMethod)
					->notifier(boost::bind(&CommandLineParsing::validate_qShapeMethod,this,_1))
				, std::string("SHAPE: method how to integrate SHAPE reactivity data into query accessibility computation via pseudo energies:"
					"\n"
					" 'D': Convert by using a linear equation according to Deigan et al. (2009). The calculated pseudo energies will "
					"be applied for every nucleotide involved in a stacked pair. "
					"The slope 'm' and the intercept 'b' can be set using e.g. 'Dm1.8b-0.6' (=defaults for 'D')."
					"\n"
					" 'Z':  Convert according to Zarringhalam et al. (2012) via pairing probabilities by using linear mapping. "
					"Aberration from the observed pairing probabilities will be penalized, which can be adjusted by the factor beta e.g. 'Zb0.89' (=default for 'Z')."
					"\n"
					" 'W': Apply a given vector of perturbation energies to unpaired nucleotides according to Washietl et al. (2012). "
					"Perturbation vectors can be calculated by using RNApvmin."
					).c_str())
		("tShapeMethod"
			, value<std::string>(&tShapeMethod)
				->notifier(boost::bind(&CommandLineParsing::validate_tShapeMethod,this,_1))
			, std::string("SHAPE: method how to integrate SHAPE reactivity data into target accessibility computation via pseudo energies.\n"
					"[for encodings see --qShapeMethod]"
					).c_str())
		("qShapeConversion"
				, value<std::string>(&qShapeConversion)
					->notifier(boost::bind(&CommandLineParsing::validate_qShapeConversion,this,_1))
				, std::string("SHAPE: method how to convert SHAPE reactivities to pairing probabilities for query accessibility computation. "
					"This parameter is useful when dealing with the SHAPE incorporation according to Zarringhalam et al. (2012). "
					"The following methods can be used to convert SHAPE reactivities into the probability for a certain nucleotide to be unpaired:"
					"\n"
					" 'M': Linear mapping according to Zarringhalam et al. (2012)"
					"\n"
					" 'C': Use a cutoff-approach to divide into paired and unpaired nucleotides, e.g. 'C0.25' (= default for 'C')"
					"\n"
					" 'S': Skip the normalizing step since the input data already represents probabilities for being unpaired rather than raw reactivity values"
					"\n"
					" 'L': Linear model to convert reactivity into a probability for being unpaired, e.g. 'Ls0.68i0.2' for slope of 0.68 and intercept of 0.2 (=default for 'L')"
					"\n"
					" 'O': Linear model to convert the log reactivity into a probability for being unpaired, e.g. 'Os1.6i-2.29' to use a slope of 1.6 and an intercept of -2.29 (=default for 'O')"
					).c_str())
		("tShapeConversion"
				, value<std::string>(&tShapeConversion)
					->notifier(boost::bind(&CommandLineParsing::validate_tShapeConversion,this,_1))
				, std::string("SHAPE: method how to convert SHAPE reactivities to pairing probabilities for target accessibility computation.\n"
					"[for encodings see --qShapeConversion]"
					).c_str())
			;

	////  INTERACTION/ENERGY OPTIONS  ////////////////////////

	opts_inter.add_options()
		("mode,m"
			, value<char>(&(predMode.val))
				->default_value(predMode.def)
				->notifier(boost::bind(&CommandLineParsing::validate_predMode,this,_1))
			, std::string("prediction mode : "
					"\n 'H' = heuristic (fast and low memory), "
					"\n 'M' = exact and low memory, "
					"\n 'E' = exact (high memory)"
					).c_str())
		;
	opts_cmdline_short.add(opts_inter);
	opts_inter.add_options()
		("pred"
			, value<char>(&(pred.val))
				->default_value(pred.def)
				->notifier(boost::bind(&CommandLineParsing::validate_pred,this,_1))
			, std::string("prediction target : "
					"\n 'S' = single-site minimum-free-energy interaction (interior loops only), "
					"\n 'P' = single-site maximum-probability interaction (interior loops only)"
					).c_str())
		("energy,e"
			, value<char>(&(energy.val))
				->default_value(energy.def)
				->notifier(boost::bind(&CommandLineParsing::validate_energy,this,_1))
			, std::string("energy computation :"
					"\n 'B' base pair energy -1 (Nussinov-like, i.e. independently of context), or"
					"\n 'V' VRNA-based computation (Nearest-Neighbor model, see also --energVRNA)").c_str())
		("energyVRNA"
			, value<std::string>(&energyFile)
				->notifier(boost::bind(&CommandLineParsing::validate_energyFile,this,_1))
			, std::string("energy parameter file of VRNA package to be used. If not provided, the default parameter set of the linked Vienna RNA package is used.").c_str())
		("temperature"
			, value<T_type>(&(temperature.val))
				->default_value(temperature.def)
				->notifier(boost::bind(&CommandLineParsing::validate_temperature,this,_1))
			, std::string("temperature in Celsius to setup the VRNA energy parameters"
					" (arg in range ["+toString(temperature.min)+","+toString(temperature.max)+"])").c_str())
		("windowWidth"
			, value<int>(&(windowWidth.val))
				->default_value(windowWidth.def)
				->notifier(boost::bind(&CommandLineParsing::validate_windowWidth,this,_1))
			, std::string("Window-based computation: width of the window to be used; 0 disables window-based computation"
					" (arg in range ["+toString(windowWidth.min)+","+toString(windowWidth.max)+"])").c_str())
		("windowOverlap"
			, value<int>(&(windowOverlap.val))
				->default_value(windowOverlap.def)
				->notifier(boost::bind(&CommandLineParsing::validate_windowOverlap,this,_1))
			, std::string("Window-based computation: overlap of the window to be used."
					" Has to be smaller than --windowWidth and greater or equal than the maximal interaction length (see --q|tIntMaxLen)."
					" (arg in range ["+toString(windowOverlap.min)+","+toString(windowOverlap.max)+"])").c_str())
		;


	////  OUTPUT OPTIONS  ////////////////////////////////////

	opts_output.add_options()
		("out"
			, value< std::vector<std::string> >(&(out))
				->composing()
				->default_value(boost::assign::list_of(outPrefix2streamName.at(OutPrefixCode::OP_EMPTY)), outPrefix2streamName.at(OutPrefixCode::OP_EMPTY).c_str())
				->notifier(boost::bind(&CommandLineParsing::validate_out,this,_1))
			, std::string("output (multi-arg) : provide a file name for output (will be overwritten)"
					" or 'STDOUT/STDERR' to write to the according stream (according to --outMode)."
					"\nUse one of the following PREFIXES (colon-separated) to generate"
					" ADDITIONAL output:"
					"\n 'qMinE:' (query) for each position the minimal energy of any interaction covering the position (CSV format)"
					"\n 'qSpotProb:' (query) for each position the probability that is is covered by an interaction covering (CSV format)"
					"\n 'qAcc:' (query) ED accessibility values ('qPu'-like format)."
					"\n 'qPu:' (query) unpaired probabilities values (RNAplfold format)."
					"\n 'tMinE:' (target) for each position the minimal energy of any interaction covering the position (CSV format)"
					"\n 'tSpotProb:' (target) for each position the probability that is is covered by an interaction covering (CSV format)"
					"\n 'tAcc:' (target) ED accessibility values ('tPu'-like format)."
					"\n 'tPu:' (target) unpaired probabilities values (RNAplfold format)."
					"\n 'pMinE:' (query+target) for each index pair the minimal energy of any interaction covering the pair (CSV format)"
					"\n 'spotProb:' (query+target) tracks for a given set of interaction spots their probability to be covered by an interaction. Spots are encoded by comma-separated 'idx1&idx2' pairs. For each spot a probability is provided in concert with the probability that none of the spots (encoded by '0&0') is covered (CSV format). The spot encoding is followed colon-separated by the output stream/file name, eg. '--out=\"spotProb:3&76,59&2:STDERR\"'. NOTE: value has to be quoted due to '&' symbol!"
					"\nFor each, provide a file name or STDOUT/STDERR to write to the respective output stream."
					).c_str())
		("outMode"
			, value<char>(&(outMode.val))
				->default_value(outMode.def)
				->notifier(boost::bind(&CommandLineParsing::validate_outMode,this,_1))
			, std::string("output mode :"
					"\n 'N' normal output (ASCII char + energy),"
					"\n 'D' detailed output (ASCII char + energy/position details),"
					"\n 'C' CSV output (see --outCsvCols),"
					"\n '1' backward compatible IntaRNA v1.* normal output,"
					"\n 'O' backward compatible IntaRNA v1.* detailed output (former -o)"
					).c_str())
	    ("outNumber,n"
			, value<int>(&(outNumber.val))
				->default_value(outNumber.def)
				->notifier(boost::bind(&CommandLineParsing::validate_outNumber,this,_1))
			, std::string("number of (sub)optimal interactions to report"
					" (arg in range ["+toString(outNumber.min)+","+toString(outNumber.max)+"])").c_str())
	    ("outOverlap"
			, value<char>(&(outOverlap.val))
				->default_value(outOverlap.def)
				->notifier(boost::bind(&CommandLineParsing::validate_outOverlap,this,_1))
			, std::string("suboptimal output : interactions can overlap "
					"\n 'N' in none of the sequences, "
					"\n 'T' in the target only, "
					"\n 'Q' in the query only, "
					"\n 'B' in both sequences").c_str())
		;
	opts_cmdline_short.add(opts_output);
	opts_output.add_options()
	    ("outMaxE"
			, value<double>(&(outMaxE.val))
				->default_value(outMaxE.def)
				->notifier(boost::bind(&CommandLineParsing::validate_outMaxE,this,_1))
			, std::string("only interactions with E <= maxE are reported"
					" (arg in range ["+toString(outMaxE.min)+","+toString(outMaxE.max)+"])").c_str())
	    ("outMinPu"
			, value<double>(&(outMinPu.val))
				->default_value(outMinPu.def)
				->notifier(boost::bind(&CommandLineParsing::validate_outMinPu,this,_1))
			, std::string("only interactions where all individual positions of both interacting sites have an unpaired probability >= minPu are reported"
					" (arg in range ["+toString(outMinPu.min)+","+toString(outMinPu.max)+"])").c_str())
	    ("outDeltaE"
			, value<double>(&(outDeltaE.val))
				->default_value(outDeltaE.def)
				->notifier(boost::bind(&CommandLineParsing::validate_outDeltaE,this,_1))
			, std::string("suboptimal output : only interactions with E <= (minE+deltaE) are reported"
					" (arg in range ["+toString(outDeltaE.min)+","+toString(outDeltaE.max)+"])").c_str())
		("outCsvCols"
			, value<std::string>(&(outCsvCols))
				->default_value(outCsvCols,"see text")
				->notifier(boost::bind(&CommandLineParsing::validate_outCsvCols,this,_1))
			, std::string("output : comma separated list of CSV column IDs to print if outMode=CSV."
					" An empty argument (using '') prints all possible columns from the following available ID list: "
					+ boost::replace_all_copy(OutputHandlerCsv::list2string(OutputHandlerCsv::string2list("")), ",", ", ")+"."
					+ "\nDefault = '"+outCsvCols+"'."
					).c_str())
	    ("outPerRegion", "output : if given, best interactions are reported independently"
	    		" for all region combinations; otherwise only the best for each query-target combination")
	    ("verbose,v", "verbose output") // handled via easylogging++
	    ("default-log-file", value<std::string>(&(logFileName)), "name of file to be used for log output (INFO, WARNING, VERBOSE, DEBUG)")
	    ;

	////  GENERAL OPTIONS  ////////////////////////////////////

	opts_general.add_options()
#if INTARNA_MULITHREADING
	    ("threads"
			, value<int>(&(threads.val))
				->default_value(threads.def)
				->notifier(boost::bind(&CommandLineParsing::validate_threads,this,_1))
			, std::string("maximal number of threads to be used for parallel computation of query-target combinations."
					" A value of 0 requests all available CPUs."
					" Note, the number of threads multiplies the required memory used for computation!"
					" (arg in range ["+toString(threads.min)+","+toString(threads.max)+"])").c_str())
#endif
	    ("version", "print version")
	    ("help,h", "show the help page for basic parameters")
	    ("fullhelp", "show the extended help page for all available parameters")
	    ;
	opts_cmdline_short.add(opts_general);

	////  GENERAL OPTIONS  ////////////////////////////////////

	opts_cmdline_all.add(opts_query).add(opts_target).add(opts_seed).add(opts_shape).add(opts_inter).add(opts_output).add(opts_general);


}

////////////////////////////////////////////////////////////////////////////

CommandLineParsing::~CommandLineParsing() {

	 INTARNA_CLEANUP(seedConstraint);

	// reset output stream
	deleteOutputStream( outStream );
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
			std::cout
					<<INTARNA_PACKAGE_STRING
					<< "\n"
					<<" using Vienna RNA package "
						<<STRINGIFY(VRNA_VERSION)
					<<" and boost "
						<<(BOOST_VERSION / 100000)
						<<'.'<<(BOOST_VERSION / 100 % 1000)
						<<'.'<<(BOOST_VERSION%100)
					<< "\n";
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
				// open according stream
				outStream = newOutputStream( outPrefix2streamName.at(OutPrefixCode::OP_EMPTY) );
				// check success
				if (outStream == NULL) {
					throw error("could not open output file --out='"+toString(outPrefix2streamName.at(OutPrefixCode::OP_EMPTY))+ "' for writing");
				}
			}
			outPerRegion = vm.count("outPerRegion") > 0;

			// parse the sequences
			parseSequences("query",queryArg,query,qSet);
			parseSequences("target",targetArg,target,tSet);

			// valide accessibility input from file (requires parsed sequences)
			validate_qAccFile( qAccFile );
			validate_tAccFile( tAccFile );

			// check seed setup
			noSeedRequired = vm.count("noSeed") > 0;
			if (noSeedRequired) {
				// input sanity check : maybe seed constraints defined -> warn
				if (!seedTQ.empty()) LOG(INFO) <<"no seed constraint wanted, but explicit seedTQ provided (will be ignored)";
				if (seedBP.val != seedBP.def) LOG(INFO) <<"no seed constraint wanted, but seedBP provided (will be ignored)";
				if (seedMaxUP.val != seedMaxUP.def) LOG(INFO) <<"no seed constraint wanted, but seedMaxUP provided (will be ignored)";
				if (seedQMaxUP.val != seedQMaxUP.def) LOG(INFO) <<"no seed constraint wanted, but seedQMaxUP provided (will be ignored)";
				if (seedTMaxUP.val != seedTMaxUP.def) LOG(INFO) <<"no seed constraint wanted, but seedTMaxUP provided (will be ignored)";
				if (seedMaxE.val != seedMaxE.def) LOG(INFO) <<"no seed constraint wanted, but seedMaxE provided (will be ignored)";
				if (seedMinPu.val != seedMinPu.def) LOG(INFO) <<"no seed constraint wanted, but seedMinPu provided (will be ignored)";
				if (!seedQRange.empty()) LOG(INFO) <<"no seed constraint wanted, but seedQRange provided (will be ignored)";
				if (!seedTRange.empty()) LOG(INFO) <<"no seed constraint wanted, but seedTRange provided (will be ignored)";
			} else {
				// check query search ranges
				if (!seedQRange.empty()) {
					if (query.size()!=1) {
						throw error("seedQRange given but not only one query sequence provided");
					} else {
						validate_indexRangeList("seedQRange",seedQRange, 1, query.begin()->size());
					}
				}
				// check target search ranges
				if (!seedTRange.empty()) {
					if (target.size()!=1) {
						throw error("seedTRange given but not only one target sequence provided");
					} else {
						validate_indexRangeList("seedTRange",seedTRange, 1, target.begin()->size());
					}
				}

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

				// check for explicit seed constraints
				if (!seedTQ.empty()) {
					if (target.size()>1 || query.size() > 1) {
						throw error("explicit seed definition only for single query/target available");
					}
					// input sanity check : maybe seed constraints defined -> warn
					if (seedBP.val != seedBP.def) LOG(INFO) <<"explicit seeds defined, but seedBP provided (will be ignored)";
					if (seedMaxUP.val != seedMaxUP.def) LOG(INFO) <<"explicit seeds defined, but seedMaxUP provided (will be ignored)";
					if (seedQMaxUP.val != seedQMaxUP.def) LOG(INFO) <<"explicit seeds defined, but seedQMaxUP provided (will be ignored)";
					if (seedTMaxUP.val != seedTMaxUP.def) LOG(INFO) <<"explicit seeds defined, but seedTMaxUP provided (will be ignored)";
					if (seedMaxE.val != seedMaxE.def) LOG(INFO) <<"explicit seeds defined, but seedMaxE provided (will be ignored)";
					if (seedMinPu.val != seedMinPu.def) LOG(INFO) <<"explicit seeds defined, but seedMinPu provided (will be ignored)";
					if (!seedQRange.empty()) LOG(INFO) <<"explicit seeds defined, but seedQRange provided (will be ignored)";
					if (!seedTRange.empty()) LOG(INFO) <<"explicit seeds defined, but seedTRange provided (will be ignored)";
				}
			}

			///////////////  PARSE AND PREPARE PREDICTION RANGES  //////////////

			// check regions to be used for interaction prediction
			if ( qRegionLenMax.val > 0 && qAcc.val == 'N' ) {
				throw error("automatic query accessible region identification requested (--qRegionLenMax) but accessibility computation disabled (--qAcc=N)");
			}
			if ( tRegionLenMax.val > 0 && tAcc.val == 'N' ) {
				throw error("automatic query accessible region identification requested (--tRegionLenMax) but accessibility computation disabled (--tAcc=N)");
			}
			if (qRegionLenMax.val > 0 && !qRegionString.empty()) {
				throw error("automatic query accessible region identification requested (--qRegionLenMax) but manual regions provided via (--qRegion)");
			}
			if (tRegionLenMax.val > 0 && !tRegionString.empty()) {
				throw error("automatic target accessible region identification requested (--tRegionLenMax) but manual regions provided via (--tRegion)");
			}
			// parse region string if available
			parseRegion( "qRegion", qRegionString, query, qRegion );
			parseRegion( "tRegion", tRegionString, target, tRegion );


			//////////////// WINDOW-BASED COMPUTATION ///////////////////

			// check if window-based computation enabled
			if (windowWidth.val > 0) {
				// minimal window width
				if (windowWidth.val < 10) {
					throw error("window-based computation: --windowWidth should be at least 10 but is "+toString(windowWidth.val));
				}
				// ensure interaction length is restricted
				size_t maxIntLength = 0;
				if (qAcc.val != 'N') { // check if accessibility computation enabled
					if (qIntLenMax.val == 0 && qAccW.val == 0) {
						throw error("window-based computation: maximal query interaction length has to be restricted either via --qAccW or --qIntLenMax");
					} else { // update maximum
						maxIntLength = std::max(qIntLenMax.val, qAccW.val);
					}
				} else { // max interaction length of query
					if (qIntLenMax.val == 0) {
						throw error("window-based computation: maximal query interaction length has to be restricted via --qIntLenMax");
					} else {
						maxIntLength = qIntLenMax.val;
					}
				}
				if (tAcc.val != 'N') {
					if (tIntLenMax.val == 0 && tAccW.val == 0) {
						throw error("window-based computation: maximal target interaction length has to be restricted either via --tAccW or --tIntLenMax");
					} else { // update maximum
						maxIntLength = maxIntLength < std::max(tAccW.val,tIntLenMax.val) ? std::max(tAccW.val,tIntLenMax.val) : maxIntLength;
					}
				} else {
					if (tIntLenMax.val == 0) {
						throw error("window-based computation: maximal target interaction length has to be restricted via --tIntLenMax");
					} else { // update maximum
						maxIntLength = maxIntLength < tIntLenMax.val ? tIntLenMax.val : maxIntLength;
					}
				}
				if (windowOverlap.val < maxIntLength) {
					throw error("window-based computation: --windowOverlap ("+toString(windowOverlap.val)+") has to be at least as large as the maximum of --q|tAccW or --q|tIntLenMax ("+toString(maxIntLength)+")");
				}
				if (windowWidth.val <= windowOverlap.val) {
					throw error("window-based computation: --windowWidth ("+toString(windowWidth.val)+") has to exceed --windowOverlap ("+toString(windowOverlap.val)+")");
				}
			}

			//////////////// ACCESSIBILITY CONSTRAINTS ///////////////////

			// check qAccConstr - query sequence compatibility
			if (vm.count("qAccConstr") > 0) {
				// only for single sequence input supported
				if (!validateSequenceNumber("qAccConstr",query,1,1)) {
					// TODO report error
					INTARNA_NOT_IMPLEMENTED("--qAccConstr only supported for single sequence input");
				}
			} else {
				// generate empty constraint
				qAccConstr = std::string(query.at(0).size(),'.');
			}
			// check tAccConstr - target sequence compatibility
			if (vm.count("tAccConstr") > 0) {
				// only for single sequence input supported
				if (!validateSequenceNumber("tAccConstr",target,1,1)) {
					// TODO report error
					INTARNA_NOT_IMPLEMENTED("--tAccConstr only supported for single sequence input");
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
				if (vm.count("qAccConstr")>0) LOG(INFO) <<"qAcc = "<<qAcc.val<<" : accessibility constraints (--qAccConstr) possibly not used in computation of loaded ED values";
				break;
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
				if (vm.count("tAccConstr")>0) LOG(INFO) <<"tAcc = "<<tAcc.val<<" : accessibility constraints (--tAccConstr) possibly not used in computation of loaded ED values";
				break;
			}	// drop to next handling
			case 'N' : {
				if (tAccL.val != tAccL.def) LOG(INFO) <<"tAcc = "<<tAcc.val<<" : ignoring --tAccL";
				if (tAccW.val != tAccW.def) LOG(INFO) <<"tAcc = "<<tAcc.val<<" : ignoring --tAccW";
				if (tAcc.val=='N' && !tAccFile.empty()) LOG(INFO) <<"tAcc = "<<tAcc.val<<" : ignoring --tAccFile";
				break;
			}
			} // switch

			// check energy setup
			if (vm.count("energyVRNA") > 0 && energy.val != 'V') {
				throw error("--energyVRNA provided but no VRNA energy computation (V) requested (--energy = "+toString(energy.val)+")");
			}

			// check qAcc upper bound
			if (qAccL.val > qAccW.val && qAccW.val != 0) {
				throw error("qAccL = " +toString(qAccL.val) + " : has to be <= qAccW (=" +toString(qAccW.val) + ")");
			}

			// check qAcc upper bound
			if (tAccL.val > tAccW.val && tAccW.val != 0) {
				throw error("tAccL = " +toString(tAccL.val)+" : has to be <= tAccW (=" +toString(tAccW.val)+")");
			}

			// check CSV stuff
			if (outCsvCols != outCsvCols_default && outMode.val != 'C') {
				throw error("outCsvCols set but outMode != C ("+toString(outMode.val)+")");
			}

			// check output sanity
			{	// check for duplicates
				bool noDuplicate = true;
				for (auto c1=outPrefix2streamName.begin(); noDuplicate && c1!=outPrefix2streamName.end(); c1++) {
					// skip empty entries
					if (c1->second.empty()) { continue; }
					for (auto c2=c1; noDuplicate && (++c2)!=outPrefix2streamName.end();) {
						if ( ! c2->second.empty() && boost::iequals( c1->second, c2->second ) ) {
							noDuplicate = false;
							throw error("--out argument shows multiple times '"+toString(c1->second)+"' as target file/stream.");
						}
					}
				}
			}

#if INTARNA_MULITHREADING
			// check if multi-threading
			if (threads.val != 1 && getTargetSequences().size() > 1) {
				// warn if >= 4D space prediction enabled
				if (pred.val != 'S' || predMode.val == 'E') {
					LOG(WARNING) <<"Multi-threading enabled in high-mem-prediction mode : ensure you have enough memory available!";
				}
				if (outMode.val == '1' || outMode.val == 'O') {
					throw error("Multi-threading not supported for IntaRNA v1 output");
				}
			}
#endif

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
	if ( energy.val == 'V') {
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
	// check if valid encoding
	if (!boost::regex_match( value, AccessibilityConstraint::regex, boost::match_perl )) {
		LOG(ERROR) <<"constraint "<<name<<" = '"<<value <<"' is not correctly encoded!";
		updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
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
		try {
			// test if parsable as ascending, non-overlapping range list
			IndexRangeList tmpList(value);
		} catch (std::exception & ex) {
			return false;
		}
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
		INTARNA_NOT_IMPLEMENTED("BED file input for index range list not implemented");
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

const size_t
CommandLineParsing::
getWindowWidth() const
{
	// check if disabled (== 0)
	return windowWidth.val == 0 ? std::numeric_limits<size_t>::max() : windowWidth.val;
}

////////////////////////////////////////////////////////////////////////////

const size_t
CommandLineParsing::
getWindowOverlap() const
{
	return windowOverlap.val;
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
	AccessibilityConstraint accConstraint(seq.size(),0,"","","");
	try {
		// try parsing
		accConstraint = AccessibilityConstraint(seq.size(), qAccConstr, qAccL.val, qShape, qShapeMethod, qShapeConversion);
	} catch (std::exception & ex) {
		throw std::runtime_error(toString("query accessibility constraint : ")+ex.what());
	}
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
			// file support : add sequence-specific prefix (for multi-sequence input)
			accFileStream = new std::ifstream( getFullFilename(qAccFile, NULL, &(seq)) );
			try {
				if(!accFileStream->good()){
					accFileStream->close();
					 INTARNA_CLEANUP(accFileStream);
					throw std::runtime_error("accessibility parsing of --qAccFile : could not open file '"+qAccFile+"'");
				}
			} catch (std::exception & ex) {
				accFileStream->close();
				 INTARNA_CLEANUP(accFileStream);
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
			 INTARNA_CLEANUP( accFileStream );
		}
		return acc;
	}

	case 'C' : // compute VRNA-based accessibilities
		switch( energy.val ) {

		case 'B' : // base-pair based accessibility
    		return new AccessibilityBasePair(
    						seq
    						, std::min( qIntLenMax.val == 0 ? seq.size() : qIntLenMax.val
    									, qAccW.val == 0 ? seq.size() : qAccW.val )
    						, &accConstraint
    						);

		case 'V' : // VRNA-based accessibilities
			return new AccessibilityVrna(
							seq
							, std::min( qIntLenMax.val == 0 ? seq.size() : qIntLenMax.val
										, qAccW.val == 0 ? seq.size() : qAccW.val )
							, &accConstraint
							, vrnaHandler
							, qAccW.val
							);
		default :
			INTARNA_NOT_IMPLEMENTED("query accessibility computation not implemented for energy = '"+toString(energy.val)+"'. Disable via --qAcc=N.");
		} break;
	default :
		INTARNA_NOT_IMPLEMENTED("CommandLineParsing::getQueryAccessibility : qAcc = '"+toString(qAcc.val)+"' is not supported");
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
	// create temporary constraint object (will be copied)
	AccessibilityConstraint accConstraint(seq.size(), 0, "","","");
	try {
		accConstraint = AccessibilityConstraint(seq.size(), tAccConstr, tAccL.val, tShape, tShapeMethod, tShapeConversion);
	} catch (std::exception & ex) {
		throw std::runtime_error(toString("target accessibility constraint : ")+ex.what());
	}

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
			// file support : add sequence-specific prefix (for multi-sequence input)
			accFileStream = new std::ifstream( getFullFilename(tAccFile, &(seq), NULL) );
			try {
				if(!accFileStream->good()){
					accFileStream->close();
					 INTARNA_CLEANUP(accFileStream);
					throw std::runtime_error("accessibility parsing of --tAccFile : could not open file '"+tAccFile+"'");
				}
			} catch (std::exception & ex) {
				accFileStream->close();
				 INTARNA_CLEANUP(accFileStream);
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
			 INTARNA_CLEANUP( accFileStream );
		}
		return acc;
	}
	case 'C' : // compute accessibilities
		switch( energy.val ) {

		case 'B' : // base-pair based accessibility
			return new AccessibilityBasePair(
								seq
								, std::min( qIntLenMax.val == 0 ? seq.size() : qIntLenMax.val
										, qAccW.val == 0 ? seq.size() : qAccW.val )
								, &accConstraint
								);

		case 'V' : // VRNA-based accessibilities
			return new AccessibilityVrna(
								seq
								, std::min( tIntLenMax.val == 0 ? seq.size() : tIntLenMax.val
										, tAccW.val == 0 ? seq.size() : tAccW.val )
								, &accConstraint
								, vrnaHandler
								, tAccW.val
								);
		default :
			INTARNA_NOT_IMPLEMENTED("target accessibility computation not implemented for energy = '"+toString(energy.val)+"'. Disable via --tAcc=N.");
		} break;
	default :
		INTARNA_NOT_IMPLEMENTED("CommandLineParsing::getTargetAccessibility : tAcc = '"+toString(tAcc.val)+"' is not supported");
	}
	return NULL;
}

////////////////////////////////////////////////////////////////////////////

InteractionEnergy*
CommandLineParsing::
getEnergyHandler( const Accessibility& accTarget, const ReverseAccessibility& accQuery ) const
{
	checkIfParsed();

	// check whether to compute ES values (for multi-site predictions)
	const bool initES = std::string("M").find(pred.val) != std::string::npos;

	switch( energy.val ) {
	case 'B' : return new InteractionEnergyBasePair( accTarget, accQuery, tIntLoopMax.val, qIntLoopMax.val, initES );
	case 'V' : return new InteractionEnergyVrna( accTarget, accQuery, vrnaHandler, tIntLoopMax.val, qIntLoopMax.val, initES );
	default :
		INTARNA_NOT_IMPLEMENTED("CommandLineParsing::getEnergyHandler : energy = '"+toString(energy.val)+"' is not supported");
	}
	return NULL;

}

////////////////////////////////////////////////////////////////////////////

OutputConstraint
CommandLineParsing::
getOutputConstraint()  const
{
	checkIfParsed();
	OutputConstraint::ReportOverlap overlap = OutputConstraint::ReportOverlap::OVERLAP_BOTH;
	switch(outOverlap.val) {
	case 'N' : overlap = OutputConstraint::ReportOverlap::OVERLAP_NONE; break;
	case 'T' : overlap = OutputConstraint::ReportOverlap::OVERLAP_SEQ1; break;
	case 'Q' : overlap = OutputConstraint::ReportOverlap::OVERLAP_SEQ2; break;
	case 'B' : overlap = OutputConstraint::ReportOverlap::OVERLAP_BOTH; break;
	default : throw std::runtime_error("CommandLineParsing::getOutputConstraint() : unsupported outOverlap value "+toString(outOverlap.val));
	}
	return OutputConstraint(
			  outNumber.val
			, overlap
			, static_cast<E_type>(outMaxE.val)
			, static_cast<E_type>(outDeltaE.val)
			);
}

////////////////////////////////////////////////////////////////////////////

void
CommandLineParsing::
parseSequences(const std::string & paramName,
					const std::string& paramArg,
					RnaSequenceVec& sequences,
					const IndexRangeList & seqSubset )
{

	// clear sequence container
	sequences.clear();

	// read FASTA from STDIN stream
	if (boost::iequals(paramArg,"STDIN")) {
		parseSequencesFasta(paramName, std::cin, sequences, seqSubset);
	} else
	if (RnaSequence::isValidSequenceIUPAC(paramArg)) {
		// check if sequence is to be stored
		if (seqSubset.empty() || seqSubset.covers(1)) {
			// direct sequence input
			sequences.push_back(RnaSequence(paramName,paramArg));
		} else {
			// would result in no sequence -> error
			LOG(ERROR) <<"Parsing of "<<paramName<<" : only single sequence given but not listed in sequence subset '"<<seqSubset<<"'";
			updateParsingCode( ReturnCode::STOP_PARSING_ERROR );
		}
		// check if sequence index range is within number of sequences
		if (!seqSubset.empty()
				&& seqSubset.rbegin()->to < IndexRange::LAST_INDEX
				&& seqSubset.rbegin()->to > 1)
		{
			// provide user warning of maybe wrongly defined sequence subset
			LOG(WARNING) <<"Sequence subset definition "<<seqSubset<<" exceeds sequence number "<<1<<" for parameter "<<paramName;
		}

	} else
	{
		// open file handle
		std::ifstream infile(paramArg);
		try {
			if(!infile.good()){
				LOG(ERROR) <<"FASTA parsing of "<<paramName<<" : could not open FASTA file  '"<<paramArg<<"'";
				updateParsingCode( ReturnCode::STOP_PARSING_ERROR );
			} else {
				parseSequencesFasta(paramName, infile, sequences, seqSubset);
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
	valid = valid && validateSequenceNumber(paramName, sequences, 1, 999999);
	// validate alphabet
	valid = valid && validateSequenceAlphabet(paramName, sequences);

}


////////////////////////////////////////////////////////////////////////////

void
CommandLineParsing::
parseSequencesFasta( const std::string & paramName,
					std::istream& input,
					RnaSequenceVec& sequences,
					const IndexRangeList & seqSubset)
{
	// temporary variables
	std::string line, name, sequence;
	int trimStart = 0;
	size_t seqNumber = 0;

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
					// update sequence counter
					seqNumber++;
					// check if sequence is to be stored
					if (seqSubset.empty() || seqSubset.covers(seqNumber)) {
						// store sequence
						sequences.push_back( RnaSequence( name, sequence ) );
					}
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
		// update sequence counter
		seqNumber++;
		// check if sequence is to be stored
		if (seqSubset.empty() || seqSubset.covers(seqNumber)) {
			// store sequence
			sequences.push_back( RnaSequence( name, sequence ) );
		}
	}

	// check if sequence index range is within number of sequences
	if (!seqSubset.empty()
			&& seqSubset.rbegin()->to < IndexRange::LAST_INDEX
			&& seqSubset.rbegin()->to > seqNumber)
	{
		// provide user warning of maybe wrongly defined sequence subset
		LOG(WARNING) <<"Sequence subset definition "<<seqSubset<<" exceeds sequence number "<<seqNumber<<" for parameter "<<paramName;
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
	// set up hub for prediction tracking (if needed)
	PredictionTrackerHub * predTracker = new PredictionTrackerHub();

	// check if minE-profile is to be generated
	if (!outPrefix2streamName.at(OutPrefixCode::OP_tMinE).empty() || !outPrefix2streamName.at(OutPrefixCode::OP_qMinE).empty()) {
		predTracker->addPredictionTracker(
				new PredictionTrackerProfileMinE( energy
						// add sequence-specific prefix for output file
						, getFullFilename( outPrefix2streamName.at(OutPrefixCode::OP_tMinE)
								, &(energy.getAccessibility1().getSequence())
								, &(energy.getAccessibility2().getAccessibilityOrigin().getSequence()))
						// add sequence-specific prefix for output file
						, getFullFilename( outPrefix2streamName.at(OutPrefixCode::OP_qMinE)
								, &(energy.getAccessibility1().getSequence())
								, &(energy.getAccessibility2().getAccessibilityOrigin().getSequence()))
						, "NA") );
	}

	// check if spotProb-profile is to be generated
	if (!outPrefix2streamName.at(OutPrefixCode::OP_tSpotProb).empty() || !outPrefix2streamName.at(OutPrefixCode::OP_qSpotProb).empty()) {
		predTracker->addPredictionTracker(
				new PredictionTrackerProfileSpotProb( energy
						// add sequence-specific prefix for output file
						, getFullFilename( outPrefix2streamName.at(OutPrefixCode::OP_tSpotProb)
								, &(energy.getAccessibility1().getSequence())
								, &(energy.getAccessibility2().getAccessibilityOrigin().getSequence()))
						// add sequence-specific prefix for output file
						, getFullFilename( outPrefix2streamName.at(OutPrefixCode::OP_qSpotProb)
								, &(energy.getAccessibility1().getSequence())
								, &(energy.getAccessibility2().getAccessibilityOrigin().getSequence()))
						, "0") );
	}

	// check if minE-pairs are to be generated
	if (!outPrefix2streamName.at(OutPrefixCode::OP_pMinE).empty()) {
		predTracker->addPredictionTracker(
				new PredictionTrackerPairMinE( energy
						// add sequence-specific prefix for output file
						, getFullFilename( outPrefix2streamName.at(OutPrefixCode::OP_pMinE)
								, &(energy.getAccessibility1().getSequence())
								, &(energy.getAccessibility2().getAccessibilityOrigin().getSequence()))
						, "NA") );
	}

	// check if spotProbs are to be tracked
	if (!outPrefix2streamName.at(OutPrefixCode::OP_spotProb).empty()) {
		predTracker->addPredictionTracker(
				new PredictionTrackerSpotProb( energy
								// get encoding
								, outSpotProbSpots
								, outPrefix2streamName.at(OutPrefixCode::OP_spotProb) )
							);
	}

	// check if any tracker registered
	if (predTracker->empty()) {
		// cleanup to avoid overhead
		 INTARNA_CLEANUP(predTracker);
		predTracker == NULL;
	}

	if (noSeedRequired) {
		// predictors without seed constraint
		switch( pred.val ) {
		// single-site mfe interactions (contain only interior loops)
		case 'S' : {
			switch ( predMode.val ) {
			case 'H' :  return new PredictorMfe2dHeuristic( energy, output, predTracker );
			case 'M' :  return new PredictorMfe2d( energy, output, predTracker );
			case 'E' :  return new PredictorMfe4d( energy, output, predTracker );
			default :  INTARNA_NOT_IMPLEMENTED("mode "+toString(predMode.val)+" not implemented for prediction target "+toString(pred.val));
			}
		} break;
		// single-site max-prob interactions (contain only interior loops)
		case 'P' : {
			switch ( predMode.val ) {
			case 'E' :  return new PredictorMaxProb( energy, output, predTracker );
			default :  INTARNA_NOT_IMPLEMENTED("mode "+toString(predMode.val)+" not implemented for prediction target "+toString(pred.val)+" : try --mode=E");
			}
		} break;
		// multi-site mfe interactions (contain interior and multi-loops loops)
		case 'M' : {
			switch ( predMode.val ) {
			default :  INTARNA_NOT_IMPLEMENTED("mode "+toString(predMode.val)+" not implemented for prediction target "+toString(pred.val));
			}
		} break;
		default : INTARNA_NOT_IMPLEMENTED("mode "+toString(predMode.val)+" not implemented");
		}
	} else {
		// seed-constrained predictors
		switch( pred.val ) {
		// single-site mfe interactions (contain only interior loops)
		case 'S' : {
			switch ( predMode.val ) {
			case 'H' :  return new PredictorMfe2dHeuristicSeed( energy, output, predTracker, getSeedHandler( energy ) );
			case 'M' :  return new PredictorMfe2dSeed( energy, output, predTracker, getSeedHandler( energy ) );
			case 'E' :  return new PredictorMfe4dSeed( energy, output, predTracker, getSeedHandler( energy ) );
			}
		} break;
		// single-site max-prob interactions (contain only interior loops)
		case 'P' : {
			switch ( predMode.val ) {
			case 'E' :  INTARNA_NOT_IMPLEMENTED("mode "+toString(predMode.val)+" not implemented for seed constraint (try --noSeed)"); return NULL;
			default :  INTARNA_NOT_IMPLEMENTED("mode "+toString(predMode.val)+" not implemented for prediction target "+toString(pred.val));
			}
		} break;
		// multi-site mfe interactions (contain interior and multi-loops loops)
		case 'M' : {
			switch ( predMode.val ) {
			default :  INTARNA_NOT_IMPLEMENTED("mode "+toString(predMode.val)+" not implemented for prediction target "+toString(pred.val));
			}
		} break;
		default : INTARNA_NOT_IMPLEMENTED("mode "+toString(predMode.val)+" not implemented");
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
	switch (outMode.val) {
	case 'O' :
		getOutputStream()
		<<"-------------------------" <<"\n"
		<<"INPUT" <<"\n"
		<<"-------------------------" <<"\n"
		<<"number of base pairs in seed                                  : "<<seedBP.val <<"\n"
		<<"max. number of unpaired bases in the seed region of seq. 1    : "<<(seedTMaxUP.val<0 ? seedMaxUP.val : seedTMaxUP.val) <<"\n"
		<<"max. number of unpaired bases in the seed region of seq. 2    : "<<(seedQMaxUP.val<0 ? seedMaxUP.val : seedQMaxUP.val) <<"\n"
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
		<<"Heuristic for hybridization end used                          : "<<(predMode.val=='H'?"true":"false") <<"\n"
		<<"\n"
		<<"-------------------------" <<"\n"
		<<"OUTPUT" <<"\n"
		<<"-------------------------" <<"\n"
		; break;
	case 'C' :
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
	switch (outMode.val) {
	case 'N' :
		return new OutputHandlerText( getOutputStream(), energy, 10, false );
	case 'D' :
		return new OutputHandlerText( getOutputStream(), energy, 10, true );
	case 'C' :
		return new OutputHandlerCsv( getOutputStream(), energy, OutputHandlerCsv::string2list( outCsvCols ));
	case '1' :
		return new OutputHandlerIntaRNA1( getOutputStream(), energy, false );
	case 'O' :
		return new OutputHandlerIntaRNA1( getOutputStream(), energy, true );
	default :
		INTARNA_NOT_IMPLEMENTED("Output mode "+toString(outMode.val)+" not implemented yet");
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
							, seedTMaxUP.val<0 ? seedMaxUP.val : seedTMaxUP.val
							, seedQMaxUP.val<0 ? seedMaxUP.val : seedQMaxUP.val
							, seedMaxE.val
							, (seedMinPu.val>0 ? std::min<E_type>(Accessibility::ED_UPPER_BOUND, energy.getE( seedMinPu.val )) : Accessibility::ED_UPPER_BOUND) // transform unpaired prob to ED value
							// shift ranges to start counting with 0
							, IndexRangeList( seedTRange ).shift(-1,energy.size1()-1)
							, IndexRangeList( seedQRange ).shift(-1,energy.size2()-1).reverse(energy.size2())
							, seedTQ
						);
	}
	return *seedConstraint;
}

////////////////////////////////////////////////////////////////////////////

SeedHandler *
CommandLineParsing::
getSeedHandler( const InteractionEnergy & energy ) const
{
	// get seed constraint
	const SeedConstraint & seedConstr = getSeedConstraint( energy );

	if (!seedTQ.empty()) {
		// create new seed handler for explicit seed definitions
		return new SeedHandlerExplicit( energy, seedConstr );
	} else {
		// create new seed handler using mfe computation
		return new SeedHandlerMfe( energy, seedConstr );
	}
}

////////////////////////////////////////////////////////////////////////////

const IndexRangeList&
CommandLineParsing::
getQueryRanges( const InteractionEnergy & energy, const size_t sequenceNumber ) const
{
	checkIfParsed();
#if INTARNA_IN_DEBUG_MODE
	if (sequenceNumber>=qRegion.size())
		throw std::runtime_error("CommandLineParsing::getQueryRanges("+toString(sequenceNumber)+") out of bounds");
	if (qRegion.at(sequenceNumber).empty())
		throw std::runtime_error("CommandLineParsing::getQueryRanges("+toString(sequenceNumber)+") is empty");
#endif

	// check if ranges are to be computed
	if (qRegionLenMax.val > 0) {
		// check if computation is needed
		if (qRegion.at(sequenceNumber).begin()->to - qRegion.at(sequenceNumber).begin()->from +1 > qRegionLenMax.val) {
			// compute highly accessible regions using ED-window-size = seedBP and minRangeLength = seedBP
			qRegion.at(sequenceNumber) = getQueryAccessibility( sequenceNumber )->decomposeByMaxED( qRegionLenMax.val, seedBP.val, seedBP.val);
			// inform user
			VLOG(1) <<"detected accessible regions for query '"<<getQuerySequences().at(sequenceNumber).getId()<<"' : "<<qRegion.at(sequenceNumber);
		}
	}

	// decompose ranges based in minimal unpaired probability value per position
	// since all ranges covering a position will have a lower unpaired probability
	getQueryAccessibility( sequenceNumber )->decomposeByMinPu( qRegion.at(sequenceNumber), outMinPu.val, energy.getRT() );

	return qRegion.at(sequenceNumber);
}

////////////////////////////////////////////////////////////////////////////

const IndexRangeList&
CommandLineParsing::
getTargetRanges( const InteractionEnergy & energy, const size_t sequenceNumber ) const
{
	checkIfParsed();
#if INTARNA_IN_DEBUG_MODE
	if (sequenceNumber>=tRegion.size())
		throw std::runtime_error("CommandLineParsing::getTargetRanges("+toString(sequenceNumber)+") out of bounds");
	if (tRegion.at(sequenceNumber).empty())
		throw std::runtime_error("CommandLineParsing::getTargetRanges("+toString(sequenceNumber)+") is empty");
#endif

	// check if to be computed
	if (tRegionLenMax.val > 0) {
		// check if computation is needed
		if (tRegion.at(sequenceNumber).begin()->to - tRegion.at(sequenceNumber).begin()->from +1 > tRegionLenMax.val) {
			// compute highly accessible regions using ED-window-size = seedBP and minRangeLength = seedBP
			tRegion.at(sequenceNumber) = getTargetAccessibility( sequenceNumber )->decomposeByMaxED( tRegionLenMax.val, seedBP.val, seedBP.val);
			// inform user
			VLOG(1) <<"detected accessible regions for target '"<<getTargetSequences().at(sequenceNumber).getId()<<"' : "<<tRegion.at(sequenceNumber);
		}
	}

	// decompose ranges based in minimal unpaired probability value per position
	// since all ranges covering a position will have a lower unpaired probability
	getTargetAccessibility( sequenceNumber )->decomposeByMinPu( tRegion.at(sequenceNumber), outMinPu.val, energy.getRT() );


	return tRegion.at(sequenceNumber);
}

////////////////////////////////////////////////////////////////////////////

void
CommandLineParsing::
writeAccessibility( const Accessibility& acc, const std::string & fileOrStream, const bool writeED ) const
{
	if (fileOrStream.empty())
		return;

	// setup output stream
	std::ostream * out = newOutputStream( fileOrStream );
	if (out == NULL) {
		throw std::runtime_error("could not open output file '"+fileOrStream +"' for "+(writeED?"accessibility":"unpaired probability")+" output");
	}

	// write data to stream
	if (writeED) {
		acc.writeRNAplfold_ED_text( *out );
	} else {
		acc.writeRNAplfold_Pu_text( *out, vrnaHandler.getRT() );
	}

	// clean up
	deleteOutputStream( out );
}

////////////////////////////////////////////////////////////////////////////



