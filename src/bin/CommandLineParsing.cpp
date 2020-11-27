
#include "CommandLineParsing.h"

#include "IntaRNA/general.h"

#include <cmath>
#include <stdexcept>
#include <fstream>
#include <cstdio>

#if INTARNA_MULITHREADING
	#include <omp.h>
#endif

#include <boost/bind/bind.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/assign/list_of.hpp>


extern "C" {
	#include <ViennaRNA/vrna_config.h>
}


#include "IntaRNA/AccessibilityConstraint.h"

#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/AccessibilityFromStream.h"
#include "IntaRNA/AccessibilityVrna.h"
#include "IntaRNA/AccessibilityBasePair.h"

#include "IntaRNA/HelixHandler.h"

#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/InteractionEnergyVrna.h"

#include "IntaRNA/PredictorMfe2dHeuristic.h"
#include "IntaRNA/PredictorMfe2dHelixBlockHeuristic.h"
#include "IntaRNA/PredictorMfe2d.h"

#include "IntaRNA/PredictorMfeSeedOnly.h"
#include "IntaRNA/PredictorMfe2dHeuristicSeed.h"
#include "IntaRNA/PredictorMfe2dHelixBlockHeuristicSeed.h"
#include "IntaRNA/PredictorMfe2dSeed.h"
#include "IntaRNA/PredictorMfe2dSeedExtension.h"
#include "IntaRNA/PredictorMfe2dSeedExtensionRIblast.h"
#include "IntaRNA/PredictorMfe2dHeuristicSeedExtension.h"

#include "IntaRNA/PredictorMfeEnsSeedOnly.h"
#include "IntaRNA/PredictorMfeEns2d.h"
#include "IntaRNA/PredictorMfeEns2dHeuristic.h"
#include "IntaRNA/PredictorMfeEns2dSeedExtension.h"
#include "IntaRNA/PredictorMfeEns2dHeuristicSeedExtension.h"

#include "IntaRNA/PredictionTracker.h"
#include "IntaRNA/PredictionTrackerHub.h"
#include "IntaRNA/PredictionTrackerPairMinE.h"
#include "IntaRNA/PredictionTrackerProfileMinE.h"
#include "IntaRNA/PredictionTrackerSpotProb.h"
#include "IntaRNA/PredictionTrackerSpotProbAll.h"
#include "IntaRNA/PredictionTrackerProfileSpotProb.h"

#include "IntaRNA/SeedHandlerMfe.h"
#include "IntaRNA/SeedHandlerNoBulge.h"

#include "IntaRNA/OutputStreamHandlerSortedCsv.h"

#include "IntaRNA/OutputHandlerCsv.h"
#include "IntaRNA/OutputHandlerEnsemble.h"
#include "IntaRNA/OutputHandlerText.h"


using namespace IntaRNA;


////////////////////////////////////////////////////////////////////////////

const std::string CommandLineParsing::outCsvCols_default = "id1,start1,end1,id2,start2,end2,subseqDP,hybridDP,E";

////////////////////////////////////////////////////////////////////////////

const std::string CommandLineParsing::outCsvLstSep = ":";

////////////////////////////////////////////////////////////////////////////

CommandLineParsing::CommandLineParsing( const Personality personality  )
	:
	personality( personality ),
	personalityParamValue(""),
	stdinUsed(false),
	opts_query("Query"),
	opts_target("Target"),
	opts_helix("Helix (only if --model=B)"),
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
	qIdxPos0("qIdxPos0",-9999999,9999999,1),
	qSet(),
	qSetString(""),
	qAcc("qAcc", "NCPE", 'C'),
	qAccW("qAccW", 0, 99999, 150),
	qAccL("qAccL", 0, 99999, 100),
	qAccConstr(""),
	qAccFile(""),
	qIntLenMax("qIntLenMax", 0, 99999, 0),
	qIntLoopMax("qIntLoopMax", 0, 30, 16),
	qRegionString(""),
	qRegion(),
	qRegionLenMax("qRegionLenMax", 0, 99999, 0),
	qShape(""),
	qShapeMethod("Zb0.89"),
	qShapeConversion("Os1.6i-2.29"),

	targetArg(""),
	target(),
	tIdxPos0("tIdxPos0",-9999999,9999999,1),
	tSet(),
	tSetString(""),
	tAcc("tAcc","NCPE", 'C'),
	tAccW("tAccW", 0, 99999, 150),
	tAccL("tAccL", 0, 99999, 100),
	tAccConstr(""),
	tAccFile(""),
	tIntLenMax("tIntLenMax", 0, 99999, 0),
	tIntLoopMax("tIntLoopMax", 0, 30, 10),
	tRegionString(""),
	tRegion(),
	tRegionLenMax("tRegionLenMax", 0, 99999, 0),
	tShape(""),
	tShapeMethod("Zb0.89"),
	tShapeConversion("Os1.6i-2.29"),

	// meta constraints
	acc("acc","NC", 'C'),
	accW("accW", 0, 99999, 150),
	accL("accL", 0, 99999, 100),
	intLenMax("intLenMax", 0, 99999, 0),
	intLoopMax("intLoopMax", 0, 30, 10),

	// Helix Constraints
	helixMinBP("helixMinBP",2,4,2),
	helixMaxBP("helixMaxBP",2,20,10),
	helixMaxIL("helixMaxIL",0,2,0),
	helixMinPu("helixMinPu",0,1,0),
	helixMaxE("helixMaxE",-999,+999,0),
	helixFullE(false),
	helixConstraint(NULL),

	noSeedRequired(false),
	seedTQ(""),
	seedBP("seedBP",2,20,7),
	seedMaxUP("seedMaxUP",0,20,0),
	seedQMaxUP("seedQMaxUP",-1,20,-1),
	seedTMaxUP("seedTMaxUP",-1,20,-1),
	seedMaxE("seedMaxE",-999,+999,0),
	seedMinPu("seedMinPu",0,1,0),
	seedMaxEhybrid("seedMaxEhybrid",-999,+999,999),
	seedNoGU(false),
	seedNoGUend(false),
	seedQRange(""),
	seedTRange(""),
	seedConstraint(NULL),

	temperature("temperature",0,100,37),

	model("model", "SPBX", 'X'),
	mode("mode", "HMSR", 'H'),  // R for RIblast heuristic only
#if INTARNA_MULITHREADING
	threads("threads", 0, omp_get_max_threads(), 1),
#endif
	windowWidth("windowWidth",0,99999,0),
	windowOverlap("windowOverlap",0,99999,150),

	energy("energy","BV",'V'),
	energyFile(""),
	energyAdd("energyAdd",-999,+999,0),
	energyNoDangles(false),

	accNoLP(false),
	accNoGUend(false),

	out(),
	outPrefix2streamName(),
	outMode("outMode", "NDCE", 'N' ),
	outNumber("outNumber", 0, 1000, 1),
	outOverlap("outOverlap", "NTQB", 'B' ),
	outDeltaE("outDeltaE", 0.0, 100.0, 100.0),
	outMaxE("outMaxE", -999.0, +999.0, 0.0),
	outMinPu("outMinPu", 0.0, 1.0, 0.0),
	outBestSeedOnly(false),
	outNoLP(false),
	outNoGUend(false),
	outSep(";"),
	outCsvCols(outCsvCols_default),
	outPerRegion(false),
	outSpotProbSpots(""),
	outNeedsZall(false),
	outNeedsBPs(true),

	logFileName(""),
	configFileName(""),

	vrnaHandler(),
	outStreamHandler(NULL)

{
	// report the used personality
	VLOG(1) <<"You called me "<<getPersonalityName(personality)<<" ..";

	switch (personality) {
	case IntaRNA :
	case IntaRNA3 :
		// no changes
		break;
	case IntaRNA1 :
		resetParamDefault<>(model, 'S');
		resetParamDefault<>(mode, 'H');
		resetParamDefault<>(outOverlap, 'Q');
		resetParamDefault<>(qAccW, 0);
		resetParamDefault<>(qAccL, 0);
		resetParamDefault<>(qIntLenMax, 0);
		resetParamDefault<>(qIntLoopMax, 16);
		resetParamDefault<>(tAccW, 0);
		resetParamDefault<>(tAccL, 0);
		resetParamDefault<>(tIntLoopMax, 16);
		resetParamDefault<>(accW, 0);
		resetParamDefault<>(accL, 0);
		resetParamDefault<>(intLenMax, 0);
		resetParamDefault<>(intLoopMax, 16);
#if INTARNA_MULITHREADING
		resetParamDefault<>(threads, 1);
#endif
		resetParamDefault<>(outBestSeedOnly, true, "outBestSeedOnly");
		resetParamDefault<>(seedBP, 6);
		resetParamDefault<>(seedMaxUP, 0);
		resetParamDefault<>(seedQMaxUP, -1);
		resetParamDefault<>(seedTMaxUP, -1);
		resetParamDefault<>(seedMaxE, 999);
		resetParamDefault<>(seedMinPu, 0);
		resetParamDefault<>(seedMaxEhybrid, 999);
		resetParamDefault<>(seedNoGU, false, "seedNoGU");
		resetParamDefault<>(seedNoGUend, false, "seedNoGUend");
		resetParamDefault<>(outNoLP, false, "outNoLP");
		resetParamDefault<>(outNoGUend, false, "outNoGUend");
		resetParamDefault<>(accNoLP, false, "accNoLP");
		resetParamDefault<>(accNoGUend, false, "accNoGUend");
		resetParamDefault<>(energyFile, std::string(VrnaHandler::Turner99), "energyVRNA");
		break;
	case IntaRNA2 :
		// IntaRNA v2 parameters
		resetParamDefault<>(model, 'S');
		resetParamDefault<>(mode, 'H');
		resetParamDefault<>(outOverlap, 'Q');
		resetParamDefault<>(qAccW, 150);
		resetParamDefault<>(qAccL, 100);
		resetParamDefault<>(qIntLenMax, 0);
		resetParamDefault<>(qIntLoopMax, 16);
		resetParamDefault<>(tAccW, 150);
		resetParamDefault<>(tAccL, 100);
		resetParamDefault<>(tIntLoopMax, 16);
		resetParamDefault<>(accW, 150);
		resetParamDefault<>(accL, 100);
		resetParamDefault<>(intLenMax, 0);
		resetParamDefault<>(intLoopMax, 16);
		resetParamDefault<>(energyFile, std::string(VrnaHandler::Turner04), "energyVRNA");
#if INTARNA_MULITHREADING
		resetParamDefault<>(threads, 1);
#endif
		resetParamDefault<>(outBestSeedOnly, true, "outBestSeedOnly");
		resetParamDefault<>(seedBP, 7);
		resetParamDefault<>(seedMaxUP, 0);
		resetParamDefault<>(seedQMaxUP, -1);
		resetParamDefault<>(seedTMaxUP, -1);
		resetParamDefault<>(seedMaxE, 0);
		resetParamDefault<>(seedMinPu, 0);
		resetParamDefault<>(seedMaxEhybrid, 999);
		resetParamDefault<>(seedNoGU, false, "seedNoGU");
		resetParamDefault<>(seedNoGUend, false, "seedNoGUend");
		resetParamDefault<>(outNoLP, false, "outNoLP");
		resetParamDefault<>(outNoGUend, false, "outNoGUend");
		resetParamDefault<>(accNoLP, false, "accNoLP");
		resetParamDefault<>(accNoGUend, false, "accNoGUend");
		break;
	case IntaRNAens :
		// ensemble-based predictions
		resetParamDefault<>(model, 'P');
		break;
	case IntaRNAhelix :
		// helix-block-based predictions
		resetParamDefault<>(model, 'B');
		break;
	case IntaRNAduplex :
		// RNAhybrid/RNAduplex-like optimizing hybrid only
		resetParamDefault<>(qAcc, 'N');
		resetParamDefault<>(tAcc, 'N');
		break;
	case IntaRNAexact :
		// RNAup-like exact predictions
		resetParamDefault<>(model, 'X');
		resetParamDefault<>(mode, 'M');
		resetParamDefault<>(outOverlap, 'B');
		resetParamDefault<>(qAccW, 0);
		resetParamDefault<>(qAccL, 0);
		resetParamDefault<>(qIntLenMax, 60);
		resetParamDefault<>(tAccW, 0);
		resetParamDefault<>(tAccL, 0);
		resetParamDefault<>(tIntLenMax, 60);
		resetParamDefault<>(accW, 0);
		resetParamDefault<>(accL, 0);
		resetParamDefault<>(intLenMax, 60);
		break;
	case IntaRNAseed :
		// seed-only prediction
		resetParamDefault<>(mode, 'S');
		break;
	case IntaRNAsTar :
		resetParamDefault<>(outOverlap, 'Q');
		// optimized parameters for sRNA-target prediction
		resetParamDefault<>(seedNoGU, true, "seedNoGU");
		resetParamDefault<>(seedMinPu, 0.001);
		resetParamDefault<>(tIntLenMax, 60);
		resetParamDefault<>(tIntLoopMax, 8);
		resetParamDefault<>(qIntLenMax, 60);
		resetParamDefault<>(qIntLoopMax, 8);
		resetParamDefault<>(intLenMax, 60);
		resetParamDefault<>(intLoopMax, 8);
		resetParamDefault<>(outMinPu, 0.001);
		// new features added with 3.1.0
		resetParamDefault<>(outNoLP, true, "outNoLP");
		resetParamDefault<>(outNoGUend, true, "outNoGUend");
		// ensure CSV output
		resetParamDefault<>(outMode, 'C');
		// avoid interaction traceback to speedup
		resetParamDefault<>(outCsvCols, std::string("id1,id2,start1,end1,start2,end2,E"), "outCsvCols");
		break;
	default : // no changes
		break;
	}

	using namespace boost::program_options;
	using namespace boost::placeholders;

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
				->notifier(boost::bind(&CommandLineParsing::validate_sequenceArgument,this,"query",_1))
			, "either an RNA sequence or the stream/file name from where to read"
				" the query sequences (should be the shorter sequences to increase efficiency);"
				" use 'STDIN' to read from standard input stream; sequences have to use"
				" IUPAC nucleotide encoding; output alias is [seq2]")
		;
	opts_cmdline_short.add(opts_query);
	opts_query.add_options()
		("qId"
			, value<std::string>(&qId)
				->default_value("")
				->notifier(boost::bind(&CommandLineParsing::validate_id,this,"qId",_1))
			, std::string("id (FASTA-prefix) to be used for query sequence naming.").c_str())
		(qIdxPos0.name.c_str()
			, value<long>(&(qIdxPos0.val))
				->default_value(qIdxPos0.def)
				->notifier(boost::bind(&CommandLineParsing::validate_numberArgument<long>,this,qIdxPos0,_1))
			, std::string("index of first (5') sequence position of all query sequences"
					" (arg in range ["+toString(qIdxPos0.min)+","+toString(qIdxPos0.max)+"];").c_str())
		("qSet"
			, value<std::string>(&(qSetString))
				->notifier(boost::bind(&CommandLineParsing::validate_qSet,this,_1))
			, std::string("query subset : List of sequence indices to consider for prediction in the format 'from1-to1,from2-to2,..' assuming indexing starts with 1").c_str())
		(qAcc.name.c_str()
			, value<char>(&(qAcc.val))
				->default_value(qAcc.def)
				->notifier(boost::bind(&CommandLineParsing::validate_charArgument,this,qAcc,_1))
			, std::string("accessibility computation :"
					"\n 'N' no accessibility contributions"
					"\n 'C' computation of accessibilities"
					"\n 'P' unpaired probabilities in RNAplfold format from --qAccFile"
					"\n 'E' ED values in RNAplfold Pu-like format from --qAccFile"
					).c_str())
		(qAccW.name.c_str()
			, value<int>(&(qAccW.val))
				->default_value(qAccW.def)
				->notifier(boost::bind(&CommandLineParsing::validate_numberArgumentExcludeRange<int>,this,qAccW,_1,1,2))
			, std::string("accessibility computation : sliding window size for query accessibility computation"
					" (arg in range ["+toString(qAccW.min)+","+toString(qAccW.max)+"];"
					" 0 will use to the full sequence length)."
					" Note, this also restricts the maximal interaction length (see --qIntLenMax)."
					).c_str())
		(qAccL.name.c_str()
			, value<int>(&(qAccL.val))
				->default_value(qAccL.def)
				->notifier(boost::bind(&CommandLineParsing::validate_numberArgumentExcludeRange<int>,this,qAccL,_1,1,2))
			, std::string("accessibility computation : maximal loop length (base pair span) for query accessibility computation"
					" (arg in range ["+toString(qAccL.min)+","+toString(qAccL.max)+"]; 0 will use to sliding window size 'qAccW')").c_str())
		("qAccConstr"
			, value<std::string>(&(qAccConstr))
				->notifier(boost::bind(&CommandLineParsing::validate_structureConstraintArgument,this,"qAccConstr",_1))
			, std::string("accessibility computation : structure constraint :"
					"\n EITHER a string of query sequence length encoding for each position:"
					" '.' no constraint,"
					" '"+toString(AccessibilityConstraint::dotBracket_accessible)+"' unpaired,"
					" '"+toString(AccessibilityConstraint::dotBracket_paired)+"' paired (intra-molecularly), or"
					" '"+toString(AccessibilityConstraint::dotBracket_blocked)+"' blocked."
					" Note, blocked positions are excluded from interaction prediction and constrained to be unpaired!"
					"\n OR an index range based encoding that is prefixed by the according constraint letter and a colon,"
					" e.g. 'b:3-4,33-40,p:1-2,12-20'.\n"
					"You might also want to check --energyAdd to correct the computed energies."
					).c_str())
		("qAccFile"
			, value<std::string>(&(qAccFile))
			, std::string("accessibility computation : the file/stream to be parsed, if --qAcc is to be read from file. Used 'STDIN' if to read from standard input stream.").c_str())
		(qIntLenMax.name.c_str()
			, value<int>(&(qIntLenMax.val))
				->default_value(qIntLenMax.def)
				->notifier(boost::bind(&CommandLineParsing::validate_numberArgument<int>,this,qIntLenMax,_1))
			, std::string("interaction site : maximal window size to be considered"
					" for interaction (and thus accessible) within the query"
					" (arg in range ["+toString(qIntLenMax.min)+","+toString(qIntLenMax.max)+"];"
					" 0 defaults to the full sequence length)"
					" If --qAccW is provided, the smaller window size of both is used."
					).c_str())
		(qIntLoopMax.name.c_str()
			, value<int>(&(qIntLoopMax.val))
				->default_value(qIntLoopMax.def)
				->notifier(boost::bind(&CommandLineParsing::validate_numberArgument<int>,this,qIntLoopMax,_1))
			, std::string("interaction site : maximal number of unpaired bases between neighbored interacting bases to be considered in interactions within the query"
					" (arg in range ["+toString(qIntLoopMax.min)+","+toString(qIntLoopMax.max)+"]; 0 enforces stackings only)").c_str())
		("qRegion"
			, value<std::string>(&(qRegionString))
				->notifier(boost::bind(&CommandLineParsing::validateRegion,this,"qRegion",_1))
			, std::string("interaction site : query regions to be considered for"
					" interaction prediction. Format ="
					" 'from1-to1,from2-to2,..' where indexing starts with 'qIdxPos0'."
					" Consider '--qRegionLenMax' for automatic region setup for"
					" long sequences."
					).c_str())
		(qRegionLenMax.name.c_str()
			, value<int>(&(qRegionLenMax.val))
				->default_value(qRegionLenMax.def)
				->notifier(boost::bind(&CommandLineParsing::validate_numberArgument<int>,this,qRegionLenMax,_1))
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
				->notifier(boost::bind(&CommandLineParsing::validate_sequenceArgument,this,"target",_1))
				, "either an RNA sequence or the stream/file name from where to read the target sequences (should be the longer sequences to increase efficiency); use 'STDIN' to read from standard input stream; sequences have to use IUPAC nucleotide encoding; output alias is [seq1]")
		;
	opts_cmdline_short.add(opts_target);
	opts_target.add_options()
		("tId"
			, value<std::string>(&tId)
				->default_value("")
				->notifier(boost::bind(&CommandLineParsing::validate_id,this,"tId",_1))
			, std::string("id (FASTA-prefix) to be used for target sequence naming.").c_str())
		(tIdxPos0.name.c_str()
			, value<long>(&(tIdxPos0.val))
				->default_value(tIdxPos0.def)
				->notifier(boost::bind(&CommandLineParsing::validate_numberArgument<long>,this,tIdxPos0,_1))
			, std::string("index of first (5') sequence position of all target sequences"
					" (arg in range ["+toString(tIdxPos0.min)+","+toString(tIdxPos0.max)+"];").c_str())
		("tSet"
			, value<std::string>(&(tSetString))
				->notifier(boost::bind(&CommandLineParsing::validate_tSet,this,_1))
			, std::string("target subset : List of sequence indices to consider for prediction in the format 'from1-to1,from2-to2,..' assuming indexing starts with 1").c_str())
		(tAcc.name.c_str()
			, value<char>(&(tAcc.val))
				->default_value(tAcc.def)
				->notifier(boost::bind(&CommandLineParsing::validate_charArgument,this,tAcc,_1))
			, std::string("accessibility computation :"
					"\n 'N' no accessibility contributions"
					"\n 'C' computation of accessibilities"
					"\n 'P' unpaired probabilities in RNAplfold format from --tAccFile"
					"\n 'E' ED values in RNAplfold Pu-like format from --tAccFile"
					).c_str())
		(tAccW.name.c_str()
			, value<int>(&(tAccW.val))
				->default_value(tAccW.def)
				->notifier(boost::bind(&CommandLineParsing::validate_numberArgumentExcludeRange<int>,this,tAccW,_1,1,2))
			, std::string("accessibility computation : sliding window size for target accessibility computation"
					" (arg in range ["+toString(tAccW.min)+","+toString(tAccW.max)+"];"
					" 0 will use the full sequence length)"
					" Note, this also restricts the maximal interaction length (see --tIntLenMax)."
					).c_str())
		(tAccL.name.c_str()
			, value<int>(&(tAccL.val))
				->default_value(tAccL.def)
				->notifier(boost::bind(&CommandLineParsing::validate_numberArgumentExcludeRange<int>,this,tAccL,_1,1,2))
			, std::string("accessibility computation : maximal loop size (base pair span) for target accessibility computation"
					" (arg in range ["+toString(tAccL.min)+","+toString(tAccL.max)+"]; 0 will use the sliding window size 'tAccW')").c_str())
		("tAccConstr"
			, value<std::string>(&(tAccConstr))
				->notifier(boost::bind(&CommandLineParsing::validate_structureConstraintArgument,this,"tAccConstr",_1))
			, std::string("accessibility computation : structure constraint :"
					"\n EITHER a string of target sequence length encoding for each position:"
					" '.' no constraint,"
					" '"+toString(AccessibilityConstraint::dotBracket_accessible)+"' unpaired,"
					" '"+toString(AccessibilityConstraint::dotBracket_paired)+"' paired (intra-molecularly), or"
					" '"+toString(AccessibilityConstraint::dotBracket_blocked)+"' blocked."
					" Note, blocked positions are excluded from interaction prediction and constrained to be unpaired!"
					"\n OR an index range based encoding that is prefixed by the according constraint letter and a colon,"
					" e.g. 'b:3-4,33-40,p:1-2,12-20'.\n"
					"You might also want to check --energyAdd to correct the computed energies."
					).c_str())
		("tAccFile"
			, value<std::string>(&(tAccFile))
			, std::string("accessibility computation : the file/stream to be parsed, if --tAcc is to be read from file. Used 'STDIN' if to read from standard input stream.").c_str())
		(tIntLenMax.name.c_str()
			, value<int>(&(tIntLenMax.val))
				->default_value(tIntLenMax.def)
				->notifier(boost::bind(&CommandLineParsing::validate_numberArgument<int>,this,tIntLenMax,_1))
			, std::string("interaction site : maximal window size to be considered for"
					" interaction (and thus accessible) within the target"
					" (arg in range ["+toString(tIntLenMax.min)+","+toString(tIntLenMax.max)+"];"
					" 0 defaults to the full sequence length)."
					" If --tAccW is provided, the smaller window size of both is used."
					).c_str())
		(tIntLoopMax.name.c_str()
			, value<int>(&(tIntLoopMax.val))
				->default_value(tIntLoopMax.def)
				->notifier(boost::bind(&CommandLineParsing::validate_numberArgument<int>,this,tIntLoopMax,_1))
			, std::string("interaction site : maximal number of unpaired bases between neighbored interacting bases to be considered in interactions within the target"
					" (arg in range ["+toString(tIntLoopMax.min)+","+toString(tIntLoopMax.max)+"]; 0 enforces stackings only)").c_str())
		("tRegion"
			, value<std::string>(&(tRegionString))
				->notifier(boost::bind(&CommandLineParsing::validateRegion,this,"tRegion",_1))
			, std::string("interaction site : target regions to be considered for"
					" interaction prediction. Format ="
					" 'from1-to1,from2-to2,..' where indexing starts with 'tIdxPos0'."
					" Consider '--tRegionLenMax' for automatic region setup for"
					" long sequences."
					).c_str())
		(tRegionLenMax.name.c_str()
			, value<int>(&(tRegionLenMax.val))
				->default_value(tRegionLenMax.def)
				->notifier(boost::bind(&CommandLineParsing::validate_numberArgument<int>,this,tRegionLenMax,_1))
			, std::string("interaction site : maximal length of highly accessible regions"
					" to be automatically identified. To this end, most inaccessible regions"
					" of length '--seedBP' are iteratively removed from the available indices"
					" until only regions below the given maximal length or remaining."
					" (arg in range ["+toString(tRegionLenMax.min)+","+toString(tRegionLenMax.max)+"];"
					" 0 defaults to no automatic range detection)"
					).c_str())
		;


	////  HELIX OPTIONS  ///////////////////////////////////

	opts_helix.add_options()
		(helixMinBP.name.c_str()
			, value<int>(&(helixMinBP.val))
				->default_value(helixMinBP.def)
				->notifier(boost::bind(&CommandLineParsing::validate_numberArgument<int>,this,helixMinBP,_1))
			, std::string("minimal number of base pairs inside a helix"
						  " (arg in range ["+toString(helixMinBP.min)+","+toString(helixMinBP.max)+"])").c_str())
		(helixMaxBP.name.c_str()
			, value<int>(&(helixMaxBP.val))
				->default_value(helixMaxBP.def)
				->notifier(boost::bind(&CommandLineParsing::validate_numberArgument<int>,this,helixMaxBP,_1))
			, std::string("maximal number of base pairs inside a helix"
						  " (arg in range ["+toString(helixMaxBP.min)+","+toString(helixMaxBP.max)+"])").c_str())
		(helixMaxIL.name.c_str()
			, value<int>(&(helixMaxIL.val))
				->default_value(helixMaxIL.def)
				->notifier(boost::bind(&CommandLineParsing::validate_numberArgument<int>,this,helixMaxIL,_1))
			, std::string("maximal size for each internal loop size in a helix"
						  " (arg in range ["+toString(helixMaxIL.min)+","+toString(helixMaxIL.max)+"]).").c_str())
		(helixMinPu.name.c_str()
			, value<Z_type>(&(helixMinPu.val))
				->default_value(helixMinPu.def)
				->notifier(boost::bind(&CommandLineParsing::validate_numberArgument<Z_type>,this,helixMinPu,_1))
			, std::string("minimal unpaired probability (per sequence) of considered helices"
						  " (arg in range ["+toString(helixMinPu.min)+","+toString(helixMinPu.max)+"]).").c_str())
		(helixMaxE.name.c_str()
			, value<E_kcal_type>(&(helixMaxE.val))
				->default_value(helixMaxE.def)
				->notifier(boost::bind(&CommandLineParsing::validate_numberArgument<E_kcal_type>,this,helixMaxE,_1))
			, std::string("maximal energy (excluding) a helix may have"
						  " (arg in range ["+toString(helixMaxE.min)+","+toString(helixMaxE.max)+"]).").c_str())
		("helixFullE"
			, value<bool>(&helixFullE)
				->default_value(helixFullE)
				->implicit_value(true)
			, std::string("if given (or true), the overall energy of a helix (including E_init, ED, dangling ends, ..)"
						" will be used for helixMaxE checks; otherwise only loop-terms are considered.").c_str())
		;
	opts_cmdline_short.add(opts_helix);

	////  SEED OPTIONS  ////////////////////////////////////


	opts_seed.add_options()
	    ("noSeed"
	    	, value<bool>(&noSeedRequired)
				->default_value(noSeedRequired)
				->implicit_value(true)
			, "if given (or true), no seed is enforced within the predicted interactions")

		("seedTQ"
			, value<std::string>(&(seedTQ))
				->notifier(boost::bind(&CommandLineParsing::validate_seedTQ,this,_1))
			, std::string("comma separated list of explicit seed base pair encoding(s) in the format"
					" startTbpsT&startQbpsQ, e.g. '3|||.|&7||.||', where 'startT/Q' are the"
					" indices of the 5' seed ends in target/query sequence and 'bpsT/Q'"
					" the respective dot-bar base pair encodings. This disables all other seed"
					" constraints and seed identification.").c_str())
		(seedBP.name.c_str()
			, value<int>(&(seedBP.val))
				->default_value(seedBP.def)
				->notifier(boost::bind(&CommandLineParsing::validate_numberArgument<int>,this,seedBP,_1))
			, std::string("number of inter-molecular base pairs within the seed region"
					" (arg in range ["+toString(seedBP.min)+","+toString(seedBP.max)+"])").c_str())
		;
	opts_cmdline_short.add(opts_seed);
	opts_seed.add_options()
		(seedMaxUP.name.c_str()
			, value<int>(&(seedMaxUP.val))
				->default_value(seedMaxUP.def)
				->notifier(boost::bind(&CommandLineParsing::validate_numberArgument<int>,this,seedMaxUP,_1))
			, std::string("maximal overall number (query+target) of unpaired bases within the seed region"
					" (arg in range ["+toString(seedMaxUP.min)+","+toString(seedMaxUP.max)+"])").c_str())
		(seedQMaxUP.name.c_str()
			, value<int>(&(seedQMaxUP.val))
				->default_value(seedQMaxUP.def)
				->notifier(boost::bind(&CommandLineParsing::validate_numberArgument<int>,this,seedQMaxUP,_1))
			, std::string("maximal number of unpaired bases within the query's seed region"
					" (arg in range ["+toString(seedQMaxUP.min)+","+toString(seedQMaxUP.max)+"]); if -1 the value of seedMaxUP is used.").c_str())
		(seedTMaxUP.name.c_str()
			, value<int>(&(seedTMaxUP.val))
				->default_value(seedTMaxUP.def)
				->notifier(boost::bind(&CommandLineParsing::validate_numberArgument<int>,this,seedTMaxUP,_1))
			, std::string("maximal number of unpaired bases within the target's seed region"
					" (arg in range ["+toString(seedTMaxUP.min)+","+toString(seedTMaxUP.max)+"]); if -1 the value of seedMaxUP is used.").c_str())
		(seedMaxE.name.c_str()
			, value<E_kcal_type>(&(seedMaxE.val))
				->default_value(seedMaxE.def)
				->notifier(boost::bind(&CommandLineParsing::validate_numberArgument<E_kcal_type>,this,seedMaxE,_1))
			, std::string("maximal energy a seed region may have"
					" (arg in range ["+toString(seedMaxE.min)+","+toString(seedMaxE.max)+"]).").c_str())
		(seedMaxEhybrid.name.c_str()
			, value<E_kcal_type>(&(seedMaxEhybrid.val))
				->default_value(seedMaxEhybrid.def)
				->notifier(boost::bind(&CommandLineParsing::validate_numberArgument<E_kcal_type>,this,seedMaxEhybrid,_1))
			, std::string("maximal hybridization energy (including E_init) a seed region may have"
					" (arg in range ["+toString(seedMaxEhybrid.min)+","+toString(seedMaxEhybrid.max)+"]).").c_str())
		(seedMinPu.name.c_str()
			, value<Z_type>(&(seedMinPu.val))
				->default_value(seedMinPu.def)
				->notifier(boost::bind(&CommandLineParsing::validate_numberArgument<Z_type>,this,seedMinPu,_1))
			, std::string("minimal unpaired probability (per sequence) a seed region may have"
					" (arg in range ["+toString(seedMinPu.min)+","+toString(seedMinPu.max)+"]).").c_str())
		("seedQRange"
			, value<std::string>(&(seedQRange))
				->notifier(boost::bind(&CommandLineParsing::validate_seedRange,this,"seedQRange",_1))
			, std::string("interval(s) in the query to search for seeds in format 'from1-to1,from2-to2,...' (Note, only for single query)").c_str())
		("seedTRange"
			, value<std::string>(&(seedTRange))
				->notifier(boost::bind(&CommandLineParsing::validate_seedRange,this,"seedTRange",_1))
			, std::string("interval(s) in the target to search for seeds in format 'from1-to1,from2-to2,...' (Note, only for single target)").c_str())
	    ("seedNoGU", value<bool>(&seedNoGU)
						->default_value(seedNoGU)
						->implicit_value(true)
	    		, "if given (or true), no GU base pairs are allowed within seeds")
	    ("seedNoGUend", value<bool>(&seedNoGUend)
						->default_value(seedNoGUend)
						->implicit_value(true)
	    		, "if given (or true), no GU base pairs are allowed at seed ends")
		;

	////  SHAPE OPTIONS  ////////////////////////

	opts_shape.add_options()
		("qShape"
			, value<std::string>(&qShape)
				->notifier(boost::bind(&CommandLineParsing::validate_shape,this,"qShape",_1))
			, "file name from where to read the query sequence's SHAPE reactivity data to guide its accessibility computation")
		("tShape"
			, value<std::string>(&tShape)
				->notifier(boost::bind(&CommandLineParsing::validate_shape,this,"tShape",_1))
			, "SHAPE: file name from where to read the target sequence's SHAPE reactivity data to guide its accessibility computation")
		("qShapeMethod"
				, value<std::string>(&qShapeMethod)
					->notifier(boost::bind(&CommandLineParsing::validate_shapeMethod,this,"qShapeMethod",_1))
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
				->notifier(boost::bind(&CommandLineParsing::validate_shapeMethod,this,"tShapeMethod",_1))
			, std::string("SHAPE: method how to integrate SHAPE reactivity data into target accessibility computation via pseudo energies.\n"
					"[for encodings see --qShapeMethod]"
					).c_str())
		("qShapeConversion"
				, value<std::string>(&qShapeConversion)
					->notifier(boost::bind(&CommandLineParsing::validate_shapeConversion,this,"qShapeConversion",_1))
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
					->notifier(boost::bind(&CommandLineParsing::validate_shapeConversion,this,"tShapeConversion",_1))
				, std::string("SHAPE: method how to convert SHAPE reactivities to pairing probabilities for target accessibility computation.\n"
					"[for encodings see --qShapeConversion]"
					).c_str())
			;

	////  INTERACTION/ENERGY OPTIONS  ////////////////////////

	opts_inter.add_options()
		((mode.name+",m").c_str()
			, value<char>(&(mode.val))
				->default_value(mode.def)
				->notifier(boost::bind(&CommandLineParsing::validate_charArgument,this,mode,_1))
			, std::string("prediction mode : "
					"\n 'H' = heuristic (fast and low memory), "
					"\n 'M' = exact (slow), "
					"\n 'S' = seed-only"
					).c_str())
		(model.name.c_str()
			, value<char>(&(model.val))
				->default_value(model.def)
				->notifier(boost::bind(&CommandLineParsing::validate_charArgument,this,model,_1))
			, std::string("interaction model : "
					"\n 'S' = single-site, minimum-free-energy interaction (interior loops only), "
					"\n 'X' = single-site, minimum-free-energy interaction via seed-extension (interior loops only), "
					"\n 'B' = single-site, helix-block-based, minimum-free-energy interaction (blocks of stable helices and interior loops only), "
					"\n 'P' = single-site interaction with minimal free ensemble energy per site (interior loops only)"
					).c_str())
		(acc.name.c_str()
			, value<char>(&(acc.val))
				->default_value(acc.def)
				->notifier(boost::bind(&CommandLineParsing::validate_charArgument,this,acc,_1))
			, std::string("accessibility computation :"
					"\n 'N' no accessibility contributions"
					"\n 'C' computation of accessibilities (see --accW and --accL)"
					).c_str())
		(intLenMax.name.c_str()
			, value<int>(&(intLenMax.val))
				->default_value(intLenMax.def)
				->notifier(boost::bind(&CommandLineParsing::validate_numberArgument<int>,this,intLenMax,_1))
			, std::string("interaction site : maximal window size to be considered for"
					" interaction"
					" (arg in range ["+toString(intLenMax.min)+","+toString(intLenMax.max)+"];"
					" 0 refers to the full sequence length)."
					" If --accW is provided, the smaller window size of both is used."
					).c_str())
		(intLoopMax.name.c_str()
			, value<int>(&(intLoopMax.val))
				->default_value(intLoopMax.def)
				->notifier(boost::bind(&CommandLineParsing::validate_numberArgument<int>,this,intLoopMax,_1))
			, std::string("interaction site : maximal number of unpaired bases between neighbored interacting bases to be considered in interactions"
					" (arg in range ["+toString(intLoopMax.min)+","+toString(intLoopMax.max)+"]; 0 enforces stackings only)").c_str())
		;
	opts_cmdline_short.add(opts_inter);
	opts_inter.add_options()
		(accW.name.c_str()
			, value<int>(&(accW.val))
				->default_value(accW.def)
				->notifier(boost::bind(&CommandLineParsing::validate_numberArgumentExcludeRange<int>,this,accW,_1,1,2))
			, std::string("accessibility computation : sliding window size for accessibility computation"
					" (arg in range ["+toString(accW.min)+","+toString(accW.max)+"];"
					" 0 will use the full sequence length)"
					" Note, this also restricts the maximal interaction length (see --intLenMax)."
					).c_str())
		(accL.name.c_str()
			, value<int>(&(accL.val))
				->default_value(accL.def)
				->notifier(boost::bind(&CommandLineParsing::validate_numberArgumentExcludeRange<int>,this,accL,_1,1,2))
			, std::string("accessibility computation : maximal loop size (base pair span) for accessibility computation"
					" (arg in range ["+toString(tAccL.min)+","+toString(accL.max)+"]; 0 will use the sliding window size 'accW')").c_str())
		((energy.name+",e").c_str()
			, value<char>(&(energy.val))
				->default_value(energy.def)
				->notifier(boost::bind(&CommandLineParsing::validate_charArgument,this,energy,_1))
			, std::string("energy computation :"
					"\n 'B' base pair energy -1 (Nussinov-like, i.e. independently of context), or"
					"\n 'V' VRNA-based computation (Nearest-Neighbor model, see also --energVRNA)").c_str())
		("energyVRNA"
			, value<std::string>(&energyFile)
				->default_value("Turner04")
				->notifier(boost::bind(&CommandLineParsing::validate_energyFile,this,_1))
			, std::string("energy parameter file of VRNA package to be used."
					" Directly provided models are \n"
					" - Turner99\n"
					" - Turner04\n"
					" - Andronescu07\n"
					" If not provided, the 'Turner04' parameter set of the linked Vienna RNA package is used.").c_str())
		("energyNoDangles"
				, value<bool>(&(energyNoDangles))
					->default_value(energyNoDangles)
					->implicit_value(true)
				, std::string("if given (or true), no dangling end contributions are considered within the overall interaction energy").c_str())
		("accNoLP"
				, value<bool>(&(accNoLP))
					->default_value(accNoLP)
					->implicit_value(true)
				, std::string("if given (or true), no lonely base pairs are considered for accessibility computation").c_str())
		("accNoGUend"
				, value<bool>(&(accNoGUend))
					->default_value(accNoGUend)
					->implicit_value(true)
				, std::string("if given (or true), no GU-helix-ends are considered for accessibility computation").c_str())
		(energyAdd.name.c_str()
			, value<E_kcal_type>(&(energyAdd.val))
				->default_value(energyAdd.def)
				->notifier(boost::bind(&CommandLineParsing::validate_numberArgument<E_kcal_type>,this,energyAdd,_1))
			, std::string("energy computation :"
					" if provided, this term is added to compute the overall energy of an interaction."
					" This is useful to incorporate the energy shift of applied accessibility constraints.").c_str())
		(temperature.name.c_str()
			, value<Z_type>(&(temperature.val))
				->default_value(temperature.def)
				->notifier(boost::bind(&CommandLineParsing::validate_numberArgument<Z_type>,this,temperature,_1))
			, std::string("temperature in Celsius to setup the VRNA energy parameters"
					" (arg in range ["+toString(temperature.min)+","+toString(temperature.max)+"])").c_str())
		(windowWidth.name.c_str()
			, value<int>(&(windowWidth.val))
				->default_value(windowWidth.def)
				->notifier(boost::bind(&CommandLineParsing::validate_numberArgument<int>,this,windowWidth,_1))
			, std::string("Window-based computation: width of the window to be used; 0 disables window-based computation"
					" (arg in range ["+toString(windowWidth.min)+","+toString(windowWidth.max)+"])").c_str())
		(windowOverlap.name.c_str()
			, value<int>(&(windowOverlap.val))
				->default_value(windowOverlap.def)
				->notifier(boost::bind(&CommandLineParsing::validate_numberArgument<int>,this,windowOverlap,_1))
			, std::string("Window-based computation: overlap of the window to be used."
					" Has to be smaller than --windowWidth and greater or equal than the maximal interaction length (see --q|tIntLenMax)."
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
					"\n 'pMinE:' (target+query) for each index pair the minimal energy of any interaction covering the pair (CSV format)"
					"\n 'spotProb:' (target+query) tracks for a given set of interaction spots their probability to be covered by an interaction. If no spots are provided, probabilities for all index combinations are computed. Spots are encoded by comma-separated 'idxT&idxQ' pairs (target-query). For each spot a probability is provided in concert with the probability that none of the spots (encoded by '0&0') is covered (CSV format). The spot encoding is followed colon-separated by the output stream/file name, eg. '--out=\"spotProb:3&76,59&2:STDERR\"'. NOTE: value has to be quoted due to '&' symbol!"
					"\nFor each, provide a file name or STDOUT/STDERR to write to the respective output stream."
					).c_str())
		("outSep"
			, value< std::string >(&(outSep))
				->composing()
				->default_value(outSep.c_str())
				->notifier(boost::bind(&CommandLineParsing::validate_outSep,this,_1))
			, std::string("column separator to be used in tabular CSV output").c_str())
		(outMode.name.c_str()
			, value<char>(&(outMode.val))
				->default_value(outMode.def)
				->notifier(boost::bind(&CommandLineParsing::validate_charArgument,this,outMode,_1))
			, std::string("output mode :"
					"\n 'N' normal output (ASCII char + energy),"
					"\n 'D' detailed output (ASCII char + energy/position details),"
					"\n 'C' CSV output (see --outCsvCols),"
					"\n 'E' ensemble information"
					).c_str())
	    ((outNumber.name+",n").c_str()
			, value<int>(&(outNumber.val))
				->default_value(outNumber.def)
				->notifier(boost::bind(&CommandLineParsing::validate_numberArgument<int>,this,outNumber,_1))
			, std::string("number of (sub)optimal interactions to report"
					" (arg in range ["+toString(outNumber.min)+","+toString(outNumber.max)+"])").c_str())
	    (outOverlap.name.c_str()
			, value<char>(&(outOverlap.val))
				->default_value(outOverlap.def)
				->notifier(boost::bind(&CommandLineParsing::validate_charArgument,this,outOverlap,_1))
			, std::string("suboptimal output : interactions can overlap "
					"\n 'N' in none of the sequences, "
					"\n 'T' in the target only, "
					"\n 'Q' in the query only, "
					"\n 'B' in both sequences").c_str())
		;
	opts_cmdline_short.add(opts_output);
	opts_output.add_options()
	    (outMaxE.name.c_str()
			, value<E_kcal_type>(&(outMaxE.val))
				->default_value(outMaxE.def)
				->notifier(boost::bind(&CommandLineParsing::validate_numberArgument<E_kcal_type>,this,outMaxE,_1))
			, std::string("only interactions with E < maxE are reported"
					" (arg in range ["+toString(outMaxE.min)+","+toString(outMaxE.max)+"])").c_str())
	    (outMinPu.name.c_str()
			, value<Z_type>(&(outMinPu.val))
				->default_value(outMinPu.def)
				->notifier(boost::bind(&CommandLineParsing::validate_numberArgument<Z_type>,this,outMinPu,_1))
			, std::string("only interactions where all individual positions of both interacting sites have an unpaired probability >= minPu are reported"
					" (arg in range ["+toString(outMinPu.min)+","+toString(outMinPu.max)+"])").c_str())
	    (outDeltaE.name.c_str()
			, value<E_kcal_type>(&(outDeltaE.val))
				->default_value(outDeltaE.def)
				->notifier(boost::bind(&CommandLineParsing::validate_numberArgument<E_kcal_type>,this,outDeltaE,_1))
			, std::string("suboptimal output : only interactions with E <= (minE+deltaE) are reported"
					" (arg in range ["+toString(outDeltaE.min)+","+toString(outDeltaE.max)+"])").c_str())
	    ("outBestSeedOnly"
			, value<bool>(&(outBestSeedOnly))
				->default_value(outBestSeedOnly)
				->implicit_value(true)
			, std::string("if given (or true), only the energetically best putative seed is reported").c_str())
	    ("outNoLP"
			, value<bool>(&(outNoLP))
				->default_value(outNoLP)
				->implicit_value(true)
			, std::string("if given (or true), no lonely (non-stacked) inter-molecular base pairs are allowed in predictions").c_str())
	    ("outNoGUend"
			, value<bool>(&(outNoGUend))
				->default_value(outNoGUend)
				->implicit_value(true)
			, std::string("if given (or true), no GU inter-molecular base pairs are allowed at interaction ends and interior loops (helix ends)").c_str())
		("outCsvCols"
			, value<std::string>(&(outCsvCols))
				->default_value(outCsvCols,"see text")
				->notifier(boost::bind(&CommandLineParsing::validate_outCsvCols,this,_1))
			, std::string("output : comma separated list of CSV column IDs to print if outMode=C."
					" Using '*' or an empty argument ('') prints all possible columns from the following available ID list: "
					+ OutputHandlerCsv::list2string(OutputHandlerCsv::string2list(""),", ")+"."
					+ "\nDefault = '"+outCsvCols+"'."
					).c_str())
		("outCsvSort"
			, value<std::string>(&(outCsvSort))
				->notifier(boost::bind(&CommandLineParsing::validate_outCsvSort,this,_1))
			, std::string("output : column ID from [outCsvCols] to be used for CSV row sorting if outMode=C.").c_str())
	    ("outPerRegion"
	    		, value<bool>(&outPerRegion)
						->default_value(outPerRegion)
						->implicit_value(true)
	    		, "output : if given (or true), best interactions are reported independently"
	    		" for all region combinations; otherwise only the best for each query-target combination")
	    ("verbose,v", "verbose output") // handled via easylogging++
	    ("default-log-file", value<std::string>(&(logFileName)), "file to be used for log output (INFO, WARNING, VERBOSE, DEBUG)")
	    ;

	////  GENERAL OPTIONS  ////////////////////////////////////

	opts_general.add_options()
#if INTARNA_MULITHREADING
	    (threads.name.c_str()
			, value<int>(&(threads.val))
				->default_value(threads.def)
				->notifier(boost::bind(&CommandLineParsing::validate_numberArgument<int>,this,threads,_1))
			, std::string("maximal number of threads to be used for parallel computation of query-target combinations."
					" A value of 0 requests all available CPUs."
					" Note, the number of threads multiplies the required memory used for computation!"
					" (arg in range ["+toString(threads.min)+","+toString(threads.max)+"])").c_str())
#endif
	    ("version", "print version")
	    ("personality", value<std::string>(&(personalityParamValue))
			, "IntaRNA personality to be used, which defines default values, available program arguments and tool behavior")
	    ("parameterFile", value<std::string>(&(configFileName))
			, "file from where to read additional command line arguments")
	    ("help,h", "show the help page for basic parameters")
	    ("fullhelp", "show the extended help page for all available parameters")
	    ;
	opts_cmdline_short.add(opts_general);

	////  GENERAL OPTIONS  ////////////////////////////////////

	opts_cmdline_all.add(opts_query).add(opts_target).add(opts_seed).add(opts_shape).add(opts_inter).add(opts_helix).add(opts_output).add(opts_general);


}

////////////////////////////////////////////////////////////////////////////

CommandLineParsing::~CommandLineParsing() {

	INTARNA_CLEANUP(helixConstraint);
	INTARNA_CLEANUP(seedConstraint);
	INTARNA_CLEANUP(outStreamHandler);

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
		// parse and store
		store( parse_command_line(argc, argv, opts_cmdline_all, parseStyle), vm);
		// parsing fine so far
		parsingCode = ReturnCode::KEEP_GOING;

		// check if we have to parse parameters from file
		if (vm.count("parameterFile") > 0 ) {
			const std::string paramFileName = vm.at("parameterFile").as<std::string>();
			validate_configFileName(paramFileName);
			if (parsingCode == ReturnCode::KEEP_GOING) {
				// open file handle
				std::ifstream paramfile(paramFileName);
				try {
					if(!paramfile.good()){
						LOG(ERROR) <<"Parsing of 'parameterFile' : could not open file  '"<<paramFileName<<"'";
						updateParsingCode( ReturnCode::STOP_PARSING_ERROR );
					} else {
						store( parse_config_file<char>(paramfile, opts_cmdline_all), vm);
					}
				} catch (std::exception & ex) {
					LOG(ERROR) <<"error while parsing of 'parameterFile="<<paramFileName<<"' : "<<ex.what();
					updateParsingCode( ReturnCode::STOP_PARSING_ERROR );
				}
				// close stream
				paramfile.close();
			}
		}
	} catch (error& e) {
		LOG(ERROR) <<e.what() << " : run with '--help' for allowed arguments";
		updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
	}

	// if parsing was successful, check for help request
	if (parsingCode == ReturnCode::KEEP_GOING) {
		if ( vm.count("help") > 0 ) {
			std::cout
				<<"\nIntaRNA predicts RNA-RNA interactions.\n"
				<<"\nThe following basic program arguments are supported:\n"
				<< opts_cmdline_short
				<< "\n"
				<< "Run --fullhelp for the extended list of parameters\n"
				<< "\n";
			if (personality != IntaRNA) {
				std::cout
					<<"\n############################################################################\n"
					<<"\nNote, you are using the '"
					<<getPersonalityName(personality)
					<<"' personality, which alters some default values!\n"
					<<"\n############################################################################\n"
					<< "\n";
			}
			std::cout
				<<"You can find more detailed documentation of IntaRNA on \n"
				<<"\n  https://backofenlab.github.io/IntaRNA/\n"
				<<std::endl
				;
			parsingCode = ReturnCode::STOP_ALL_FINE;
			return parsingCode;
		}
		if ( vm.count("fullhelp") > 0 ) {
			std::cout
				<<"\nIntaRNA predicts RNA-RNA interactions.\n"
				<<"\nThe following program arguments are supported:\n"
				<< opts_cmdline_all
				<< "\n";
			if (personality != IntaRNA) {
				std::cout
					<<"\n############################################################################\n"
					<<"\nNote, you are using personality '"
					<<getPersonalityName(personality)
					<<"', which alters some default values!\n(see -v output)\n"
					<<"\n############################################################################\n"
					<< "\n";
			}
			parsingCode = ReturnCode::STOP_ALL_FINE;
			return parsingCode;
		}
		if ( vm.count("version") > 0 ) {
			std::cout
					<<INTARNA_PACKAGE_STRING
					<< "\n"
					<<" using Vienna RNA package "
						<<(VRNA_VERSION)
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


			// parsing escape literals
			outSep = unescaped_string<std::string::const_iterator>::getUnescaped( outSep );

			// open output stream
			{
				// open according stream
				std::ostream* outStream = newOutputStream( outPrefix2streamName.at(OutPrefixCode::OP_EMPTY) );
				// check success
				if (outStream == NULL) {
					throw error("could not open output --out='"+toString(outPrefix2streamName.at(OutPrefixCode::OP_EMPTY))+ "' for writing");
				}
				// create final output handler
				outStreamHandler = new OutputStreamHandler(outStream);
			}

			// parse the sequences
			parseSequences("query",qId,queryArg,query,qSet,qIdxPos0.val);
			parseSequences("target",tId,targetArg,target,tSet,tIdxPos0.val);

			// validate accessibility input from file (requires parsed sequences)
			validate_qAccFile( qAccFile );
			validate_tAccFile( tAccFile );

			//////////////// INTERACTION CONSTRAINTS ///////////////////

			// interaction length: meta vs. special
			if (intLenMax.isSet()) {
				if (qIntLenMax.isSet() && intLenMax.val != qIntLenMax.val) { throw error("--intLenMax and --qIntLenMax are set to different values"); }
				if (tIntLenMax.isSet() && intLenMax.val != tIntLenMax.val) { throw error("--intLenMax and --tIntLenMax are set to different values"); }
				qIntLenMax.val = intLenMax.val;
				tIntLenMax.val = intLenMax.val;
			}

			// interior loop length: meta vs. special
			if (intLoopMax.isSet()) {
				if (qIntLoopMax.isSet() && intLoopMax.val != qIntLoopMax.val) { throw error("--intLoopMax and --qIntLoopMax are set to different values"); }
				if (tIntLoopMax.isSet() && intLoopMax.val != tIntLoopMax.val) { throw error("--intLoopMax and --tIntLoopMax are set to different values"); }
				qIntLoopMax.val = intLoopMax.val;
				tIntLoopMax.val = intLoopMax.val;
			}

			////////////////   SEED  ////////////////////////

			// check seed setup
			if (noSeedRequired) {
				// reset model if needed
				if (model.val == 'X') {
					LOG(INFO) <<"Since no seed constraint needed: resetting model from 'X' to 'S'";
					model.val='S';
				}
				if (mode.val == 'S') throw error("mode=S not applicable for non-seed predictions!");
				// input sanity check : maybe seed constraints defined -> warn
				if (!seedTQ.empty()) LOG(INFO) <<"no seed constraint wanted, but explicit seedTQ provided (will be ignored)";
				if (seedBP.isSet()) LOG(INFO) <<"no seed constraint wanted, but seedBP provided (will be ignored)";
				if (seedMaxUP.isSet()) LOG(INFO) <<"no seed constraint wanted, but seedMaxUP provided (will be ignored)";
				if (seedQMaxUP.isSet()) LOG(INFO) <<"no seed constraint wanted, but seedQMaxUP provided (will be ignored)";
				if (seedTMaxUP.isSet()) LOG(INFO) <<"no seed constraint wanted, but seedTMaxUP provided (will be ignored)";
				if (seedMaxE.isSet()) LOG(INFO) <<"no seed constraint wanted, but seedMaxE provided (will be ignored)";
				if (seedMaxEhybrid.isSet()) LOG(INFO) <<"no seed constraint wanted, but seedMaxEhybrid provided (will be ignored)";
				if (seedMinPu.isSet()) LOG(INFO) <<"no seed constraint wanted, but seedMinPu provided (will be ignored)";
				if (!seedQRange.empty()) LOG(INFO) <<"no seed constraint wanted, but seedQRange provided (will be ignored)";
				if (!seedTRange.empty()) LOG(INFO) <<"no seed constraint wanted, but seedTRange provided (will be ignored)";
				if (seedNoGU) LOG(INFO) <<"no seed constraint wanted, but seedNoGU provided (will be ignored)";
				if (seedNoGUend) LOG(INFO) <<"no seed constraint wanted, but seedNoGUend provided (will be ignored)";
			} else {
				// check query search ranges
				if (!seedQRange.empty()) {
					for (auto q : query) {
						validate_indexRangeList("seedQRange",seedQRange, q);
					}
				}
				// check target search ranges
				if (!seedTRange.empty()) {
					for (auto t : target) {
						validate_indexRangeList("seedTRange",seedTRange, t);
					}
				}


				// check for explicit seed constraints
				if (seedTQ.empty()) {
					// check for minimal sequence length (>=seedBP)
					for( auto q : query) {
						if (q.size() < seedBP.val) {
							throw error("length of query "+q.getId()+" is below minimal number of seed base pairs (seedBP="+toString(seedBP.val)+")");
						}
					}
					for( auto t : target) {
						if (t.size() < seedBP.val) {
							throw error("length of target sequence "+t.getId()+" is below minimal number of seed base pairs (seedBP="+toString(seedBP.val)+")");
						}
					}
				} else {
					if (target.size()>1 || query.size() > 1) {
						throw error("explicit seed definition only for single query/target available");
					}
					// input sanity check : maybe seed constraints defined -> warn
					if (seedBP.isSet()) LOG(INFO) <<"explicit seeds defined, but seedBP provided (will be ignored)";
					if (seedMaxUP.isSet()) LOG(INFO) <<"explicit seeds defined, but seedMaxUP provided (will be ignored)";
					if (seedQMaxUP.isSet()) LOG(INFO) <<"explicit seeds defined, but seedQMaxUP provided (will be ignored)";
					if (seedTMaxUP.isSet()) LOG(INFO) <<"explicit seeds defined, but seedTMaxUP provided (will be ignored)";
					if (seedMaxE.isSet()) LOG(INFO) <<"explicit seeds defined, but seedMaxE provided (will be ignored)";
					if (seedMaxEhybrid.isSet()) LOG(INFO) <<"explicit seeds defined, but seedMaxEhybrid provided (will be ignored)";
					if (seedMinPu.isSet()) LOG(INFO) <<"explicit seeds defined, but seedMinPu provided (will be ignored)";
					if (!seedQRange.empty()) LOG(INFO) <<"explicit seeds defined, but seedQRange provided (will be ignored)";
					if (!seedTRange.empty()) LOG(INFO) <<"explicit seeds defined, but seedTRange provided (will be ignored)";
					if (seedNoGU) LOG(INFO) <<"explicit seeds defined, but seedNoGU provided (will be ignored)";
					if (seedNoGUend) LOG(INFO) <<"explicit seeds defined, but seedNoGUend provided (will be ignored)";
				}
				// compare maximal interaction length with minimal seed length
				if (qIntLenMax.val > 0 && qIntLenMax.val < seedBP.val) { throw error("maximal query interaction length < seedBP"); }
				if (tIntLenMax.val > 0 && tIntLenMax.val < seedBP.val) { throw error("maximal target interaction length < seedBP"); }
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

			//////////////// ACCESSIBILITY CONSTRAINTS ///////////////////

			// accessibility setup: meta vs. special
			if (acc.isSet()) {
				if (qAcc.isSet() && acc.val != qAcc.val) { throw error("--acc and --qAcc are set to different values"); }
				if (tAcc.isSet() && acc.val != tAcc.val) { throw error("--acc and --tAcc are set to different values"); }
				qAcc.val = acc.val;
				tAcc.val = acc.val;
			}
			if (accW.isSet()) {
				if (qAccW.isSet() && accW.val != qAccW.val) { throw error("--accW and --qAccW are set to different values"); }
				if (tAccW.isSet() && accW.val != tAccW.val) { throw error("--accW and --tAccW are set to different values"); }
				qAccW.val = accW.val;
				tAccW.val = accW.val;
			}
			if (accL.isSet()) {
				if (qAccL.isSet() && accL.val != qAccL.val) { throw error("--accL and --qAccL are set to different values"); }
				if (tAccL.isSet() && accL.val != tAccL.val) { throw error("--accL and --tAccL are set to different values"); }
				qAccL.val = accL.val;
				tAccL.val = accL.val;
			}

			// check qAccConstr - query sequence compatibility
			if (vm.count("qAccConstr") > 0) {
				// only for single sequence input supported
				if (!validateSequenceNumber("qAccConstr",query,1,1)) {
					// report error
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
					// report error
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
				if (energy.val != 'V') {
					if (accNoLP) LOG(INFO) <<"ignoring --accNoLP, not supported for energy=="<<energy.val;
					if (accNoGUend) LOG(INFO) <<"ignoring --accNoGUend, not supported for energy=="<<energy.val;
				}
				break;
			}
			case 'E' : // drop to next handling
			case 'P' : {
				if (qAccFile.empty()) LOG(INFO) <<"qAcc = "<<qAcc.val<<" but no --qAccFile given";
				if (vm.count("qAccConstr")>0) LOG(INFO) <<"qAcc = "<<qAcc.val<<" : accessibility constraints (--qAccConstr) possibly not used in computation of loaded ED values";
				break;
			}	// drop to next handling
			case 'N' : {
				if (qAccL.isSet()) LOG(INFO) <<"qAcc = "<<qAcc.val<<" : ignoring --qAccL";
				if (qAccW.isSet()) LOG(INFO) <<"qAcc = "<<qAcc.val<<" : ignoring --qAccW";
				if (qAcc.val != 'N' && !qAccFile.empty()) LOG(INFO) <<"qAcc = "<<qAcc.val<<" : ignoring --qAccFile";
				break;
			}
			} // switch
			switch(tAcc.val) {
			case 'C' : {
				if (!tAccFile.empty()) LOG(INFO) <<"tAcc = "<<tAcc.val<<" : ignoring --tAccFile";
				if (energy.val != 'V') {
					if (accNoLP) LOG(INFO) <<"ignoring --accNoLP, not supported for energy=="<<energy.val;
					if (accNoGUend) LOG(INFO) <<"ignoring --accNoGUend, not supported for energy=="<<energy.val;
				}
				break;
			}
			case 'E' : // drop to next handling
			case 'P' : {
				if (tAccFile.empty()) LOG(INFO) <<"tAcc = "<<tAcc.val<<" but no --tAccFile given";
				if (vm.count("tAccConstr")>0) LOG(INFO) <<"tAcc = "<<tAcc.val<<" : accessibility constraints (--tAccConstr) possibly not used in computation of loaded ED values";
				break;
			}	// drop to next handling
			case 'N' : {
				if (tAccL.isSet()) LOG(INFO) <<"tAcc = "<<tAcc.val<<" : ignoring --tAccL";
				if (tAccW.isSet()) LOG(INFO) <<"tAcc = "<<tAcc.val<<" : ignoring --tAccW";
				if (tAcc.val=='N' && !tAccFile.empty()) LOG(INFO) <<"tAcc = "<<tAcc.val<<" : ignoring --tAccFile";
				break;
			}
			} // switch


			// check qAcc upper bound
			if (qAccL.val > qAccW.val && qAccW.val != 0) {
				throw error("qAccL = " +toString(qAccL.val) + " : has to be <= qAccW (=" +toString(qAccW.val) + ")");
			}

			// check qAcc upper bound
			if (tAccL.val > tAccW.val && tAccW.val != 0) {
				throw error("tAccL = " +toString(tAccL.val)+" : has to be <= tAccW (=" +toString(tAccW.val)+")");
			}


			/////////////////  OUTPUT MODE  //////////////

			// check ensemble sanity
			switch(outMode.val) {
			case 'E' :
				// check single sequence input
				if (target.size() > 1 || query.size() > 1) {
					throw error("outmode=E allows only single sequence input for query and target");
				}
				// ensure separator is not used within sequence name
				if (qId.find(outSep) != std::string::npos) { throw error("--qId contains --outSep separator..."); }
				if (tId.find(outSep) != std::string::npos) { throw error("--tId contains --outSep separator..."); }
				break;
			case 'C' :
				// ensure separator is not used within sequence name
				if (qId.find(outSep) != std::string::npos) { throw error("--qId contains --outSep separator..."); }
				if (tId.find(outSep) != std::string::npos) { throw error("--tId contains --outSep separator..."); }
				break;
			default:
				break;
			}

			// check CSV stuff
			if (outCsvCols != outCsvCols_default && outMode.val != 'C') {
				throw error("outCsvCols set but outMode != C ("+toString(outMode.val)+")");
			}
			if (!outCsvSort.empty()) {
				if (outMode.val != 'C') {
					throw error("outCsvSort set but outMode != C ("+toString(outMode.val)+")");
				}
				// check if column to sort is within output list
				std::string curCsvCols = (outCsvCols.empty() || outCsvCols=="*" ? OutputHandlerCsv::list2string(OutputHandlerCsv::string2list(""),",") : outCsvCols);
				if ( ("," + curCsvCols + ",").find(","+outCsvSort+",") == std::string::npos ) {
					throw error("outCsvSort column ID '"+outCsvSort+"' is not within outCsvCols list");
				}
				// get index of column to sort
				OutputHandlerCsv::ColTypeList csvCols = OutputHandlerCsv::string2list( curCsvCols );
				OutputHandlerCsv::ColType outCsvColType = *(OutputHandlerCsv::string2list( outCsvSort ).begin());
				auto it = std::find(csvCols.begin(), csvCols.end(), outCsvColType);
				size_t outCsvSortIdx = std::distance(csvCols.begin(), it);
				// check whether to do lex or numeric ordering
				bool sortLexOrder = (std::find( OutputHandlerCsv::colTypeNumericSort.begin(), OutputHandlerCsv::colTypeNumericSort.end(), outCsvColType) == OutputHandlerCsv::colTypeNumericSort.end());
				// setup sorted CSV output
				OutputStreamHandler * tmpOSH = outStreamHandler;
				outStreamHandler = new OutputStreamHandlerSortedCsv( tmpOSH, outCsvSortIdx, sortLexOrder, outSep, true, outCsvLstSep );
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
					throw error("window-based computation: --windowOverlap ("+toString(windowOverlap.val)+") has to be at least as large as the maximal interaction length, i.e. max of --q|tAccW and --q|tIntLenMax ("+toString(maxIntLength)+")");
				}
				if (windowWidth.val <= windowOverlap.val) {
					throw error("window-based computation: --windowWidth ("+toString(windowWidth.val)+") has to exceed --windowOverlap ("+toString(windowOverlap.val)+")");
				}
				if (outNumber.val > 1 && outOverlap.val != 'B') {
					throw error("window-based computation: non-overlapping subopt output (-n > 1) only supported for --outOverlap=B");
				}
			}


			//////////////// MODEL ////////////////////////

			// final checks if parameter set compatible with model
			switch (model.val) {
			// ensemble based predictions
			case 'P' : {
				// ensure Zall is computed
				outNeedsZall = true;
				// no window decomposition of regions (overlapping regions break overall partition function computation)
				if (windowWidth.val != 0)  throw error("windowWidth not supported for --model=P");
				break;}
			case 'B' : {
				// check compatibility of seed constraint with helix setup
				if (!noSeedRequired) {
					// check if helixMaxBP >= seedBP
					if (helixMaxBP.val < ( seedTQ.empty() ? seedBP.val : SeedHandlerExplicit::getSeedMaxBP(seedTQ) ) ) {
						throw error("maximum number of seed base pairs ("
								+toString(( seedTQ.empty() ? seedBP.val : SeedHandlerExplicit::getSeedMaxBP(seedTQ) ))
								+") exceeds the maximally allowed number of helix base pairs ("+toString(helixMaxBP.val)+")");
					}
				}
				
				// check for minimal sequence length
				for(size_t i=0; i<query.size(); i++) {
					if (query.at(i).size() < helixMinBP.val) {
						throw error("length of query sequence "+query.at(i).getId()+" is below minimal number of helix base pairs (helixMinBP="+toString(helixMinBP.val)+")");
					}
				}
				for(size_t i=0; i<target.size(); i++) {
					if (target.at(i).size() < helixMinBP.val) {
						throw error("length of target sequence "+target.at(i).getId()+" is below minimal number of helix base pairs (helixMinBP="+toString(helixMinBP.val)+")");
					}
				}
				// Ensure that min is smaller than max.
				if (helixMinBP.val > helixMaxBP.val) {
					throw error("the minimum number of base pairs (" +toString(helixMinBP.val)+") is higher than the maximum number of base pairs (" +toString(helixMaxBP.val)+")");
				}
				break;}
			default:
				break;
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
	if ( energy.val == 'V') {
		vrnaHandler = VrnaHandler( temperature.val, (energyFile.size() > 0 ? energyFile : "Turner04"), accNoGUend, accNoLP );
	}


	// return validate_* dependent parsing code
	return parsingCode;
}

////////////////////////////////////////////////////////////////////////////

void
CommandLineParsing::
validate_charArgument(const CommandLineParsing::CharParameter& param, const char & value)
{
	// alphabet check
	if ( ! param.isInAlphabet(value) ) {
		LOG(ERROR) <<""<<param.name<<" = " <<value <<" : has to be one of '" <<param.alphabet <<"'";
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
		, const RnaSequence & seq )
{
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
		try {
			IndexRangeList r(value, false, &seq);
			for (IndexRangeList::const_iterator i=r.begin(); i!=r.end(); i++) {
				// ensure range is ascending
				if (!i->isAscending()) {
					LOG(ERROR)  <<argName<<" : subrange " <<*i <<" is not ascending";
					updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
					return;
				}
				// check if boundaries in range (given they are ascending)
				if (i->to >= seq.size() ) {
					LOG(ERROR)  <<argName<<" : subrange " <<*i <<" is out of bounds [,"<<(seq.size()-1)<<"] of sequence "<<seq.getId();
					updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
					return;
				}
			}
		} catch (std::runtime_error & e) {
			LOG(ERROR)  <<argName<<" : '" <<value <<"' can not be parsed for sequence index range ["<<seq.getInOutIndex(0)<<","<<seq.getInOutIndex(seq.size()-1)<<"] of sequence "<<seq.getId();
			updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
			return;
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

void CommandLineParsing::validate_outCsvSort(const std::string & value) {
	OutputHandlerCsv::ColTypeList lst;
	try {
		// try to parse
		lst = OutputHandlerCsv::string2list( value );
	} catch (const std::runtime_error & e ) {
		LOG(ERROR) <<"--outCsvSort : " <<e.what();
		updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
	}
	if (lst.size()>1) {
		LOG(ERROR) <<"--outCsvSort : more than one column IDs given";
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
		LOG(ERROR) <<"'"<<filename<<"' : is no file!";
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
			// NOTE: non-overlapping is needed for ensemble-based prediction to ensure correctness of overall partition function!!!
			IndexRangeList tmpList(value, false);
		} catch (std::exception & ex) {
			return false;
		}
		return true;
	} else
	// no valid input
	{
		LOG(ERROR) <<"the argument for "<<argName<<" is no valid range string encoding";
		updateParsingCode(ReturnCode::STOP_PARSING_ERROR);
	}
	return false;
}

////////////////////////////////////////////////////////////////////////////

void
CommandLineParsing::
parseRegion( const std::string & argName, const std::string & value, const RnaSequenceVec & sequences, IndexRangeListVec & rangeList )
{
	// ensure range list size sufficient
	rangeList.resize( sequences.size() );
	// check if nothing given
	if (value.empty()) {
		for( size_t i=0; i<sequences.size(); i++) {
			// clear if any datathere
			rangeList[i].clear();
			// set sequence-specific full-length range
			rangeList[i].push_back( IndexRange(0,sequences.at(i).size()-1) );
		}
		return;
	} else
	// check direct range input
	if (boost::regex_match( value, IndexRangeList::regex, boost::match_perl )) {
		for( size_t i=0; i<sequences.size(); i++) {
			// validate range encodings for current sequence
			validate_indexRangeList(argName, value, sequences.at(i));
			// fill range list from string using index correction from first sequence
			rangeList[i] = IndexRangeList( value, false, &(sequences.at(i)) );
		}
		return;
	} else {
		throw boost::program_options::error(argName+" is not a comma-separated list of index ranges.");
	}
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
		accConstraint = AccessibilityConstraint(seq, qAccConstr, qAccL.val, qShape, qShapeMethod, qShapeConversion);
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
		std::istream * accStream = newInputStream( getFullFilename(qAccFile, NULL, &(seq)) );
		if (accStream == NULL) {
			throw std::runtime_error("accessibility parsing of --qAccFile : could not open file '"+qAccFile+"'");
		}
		Accessibility * acc = new AccessibilityFromStream( seq
										, qIntLenMax.val
										, &accConstraint
										, *accStream
										, (qAcc.val == 'P' ? AccessibilityFromStream::Pu_RNAplfold_Text : AccessibilityFromStream::ED_RNAplfold_Text)
										, vrnaHandler.getRT() );
		// cleanup
		deleteInputStream( accStream );
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
		accConstraint = AccessibilityConstraint(seq, tAccConstr, tAccL.val, tShape, tShapeMethod, tShapeConversion);
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
		std::istream * accStream = newInputStream( getFullFilename(tAccFile, &(seq), NULL) );
		if (accStream == NULL) {
			throw std::runtime_error("accessibility parsing of --tAccFile : could not open file '"+tAccFile+"'");
		}
		// read data
		Accessibility * acc = new AccessibilityFromStream( seq
										, tIntLenMax.val
										, &accConstraint
										, *accStream
										, ( tAcc.val == 'P' ? AccessibilityFromStream::Pu_RNAplfold_Text : AccessibilityFromStream::ED_RNAplfold_Text )
										, vrnaHandler.getRT() );
		// cleanup
		deleteInputStream( accStream );
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
	const bool initES = std::string("M").find(model.val) != std::string::npos;

	switch( energy.val ) {
	case 'B' : return new InteractionEnergyBasePair( accTarget, accQuery
						, tIntLoopMax.val, qIntLoopMax.val
						, initES, Z_type(1.0), Ekcal_2_E(-1), 3
						, Ekcal_2_E(energyAdd.val), !energyNoDangles, !outNoGUend );
	case 'V' : return new InteractionEnergyVrna( accTarget, accQuery, vrnaHandler, tIntLoopMax.val, qIntLoopMax.val, initES, Ekcal_2_E(energyAdd.val), !energyNoDangles, !outNoGUend );
	default :
		INTARNA_NOT_IMPLEMENTED("CommandLineParsing::getEnergyHandler : energy = '"+toString(energy.val)+"' is not supported");
	}
	return NULL;

}

////////////////////////////////////////////////////////////////////////////

OutputConstraint
CommandLineParsing::
getOutputConstraint( const InteractionEnergy & energy )  const
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
			, Ekcal_2_E(outMaxE.val)
			, Ekcal_2_E(outDeltaE.val)
			, outBestSeedOnly
			, outNoLP
			, outNoGUend
			, outNeedsZall
			, outNeedsBPs
			, (outMinPu.val>0 ? std::min<E_type>(Accessibility::ED_UPPER_BOUND, energy.getE( outMinPu.val )) : Accessibility::ED_UPPER_BOUND)
			);
}

////////////////////////////////////////////////////////////////////////////

void
CommandLineParsing::
parseSequences(const std::string & paramName,
				const std::string & idPrefix,
					const std::string& paramArg,
					RnaSequenceVec& sequences,
					const IndexRangeList & seqSubset,
					const long idxPos0 )
{

	// clear sequence container
	sequences.clear();

	// read FASTA from STDIN stream
	if (boost::iequals(paramArg,"STDIN")) {
		parseSequencesFasta(paramName,idPrefix, std::cin, sequences, seqSubset,idxPos0);
	} else
	if (RnaSequence::isValidSequenceIUPAC(paramArg)) {
		// check if sequence is to be stored
		if (seqSubset.empty() || seqSubset.covers(1)) {
			// direct sequence input
			sequences.push_back(RnaSequence(idPrefix.empty()?paramName:idPrefix,paramArg,idxPos0));
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
		std::istream * infile = newInputStream( paramArg );
		try {
			if (infile == NULL) {
				LOG(ERROR) <<"FASTA parsing of "<<paramName<<" : could not open FASTA file  '"<<paramArg<<"'";
				updateParsingCode( ReturnCode::STOP_PARSING_ERROR );
			} else {
				parseSequencesFasta(paramName,idPrefix, *infile, sequences, seqSubset,idxPos0);
			}
		} catch (std::exception & ex) {
			LOG(ERROR) <<"error while FASTA parsing of "<<paramName<<" : "<<ex.what();
			updateParsingCode( ReturnCode::STOP_PARSING_ERROR );
		}
		// close stream
		deleteInputStream( infile );
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
					const std::string& idPrefix,
					std::istream& input,
					RnaSequenceVec& sequences,
					const IndexRangeList & seqSubset,
					const long idxPos0 )
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
						sequences.push_back( RnaSequence( idPrefix+name, sequence, idxPos0 ) );
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
				// name = prefix up to first whitespace
				name = line.substr( 0, line.find_first_of(" \t\n\r"));
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
			sequences.push_back( RnaSequence( idPrefix+name, sequence, idxPos0 ) );
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
			LOG(ERROR) <<"sequence " <<sequences.at(i).getId()<<" for parameter "<<paramName<<" is not valid!";
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

Z_type
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
						, "NA", outSep ) );
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
						, "0", outSep ) );
	}

	// check if minE-pairs are to be generated
	if (!outPrefix2streamName.at(OutPrefixCode::OP_pMinE).empty()) {
		predTracker->addPredictionTracker(
				new PredictionTrackerPairMinE( energy
						// add sequence-specific prefix for output file
						, getFullFilename( outPrefix2streamName.at(OutPrefixCode::OP_pMinE)
								, &(energy.getAccessibility1().getSequence())
								, &(energy.getAccessibility2().getAccessibilityOrigin().getSequence()))
						, "NA", outSep ) );
	}

	// check if specific spotProbs are to be tracked
	if (!outPrefix2streamName.at(OutPrefixCode::OP_spotProb).empty()) {
		// track only specific spots
		predTracker->addPredictionTracker(
				new PredictionTrackerSpotProb( energy
								// get encoding
								, outSpotProbSpots
								, outPrefix2streamName.at(OutPrefixCode::OP_spotProb)
								, outSep )
							);
	}

	// check if all spotProbs are to be tracked
	if (!outPrefix2streamName.at(OutPrefixCode::OP_spotProbAll).empty()) {
		// track all spots
		predTracker->addPredictionTracker(
				new PredictionTrackerSpotProbAll( energy
						// add sequence-specific prefix for output file
						, getFullFilename( outPrefix2streamName.at(OutPrefixCode::OP_spotProbAll)
								, &(energy.getAccessibility1().getSequence())
								, &(energy.getAccessibility2().getAccessibilityOrigin().getSequence()))
						, "0", outSep ) );
	}

	// check if any tracker registered
	if (predTracker->empty()) {
		// cleanup to avoid overhead
		 INTARNA_CLEANUP(predTracker);
		predTracker == NULL;
	}

	if (noSeedRequired) {
		// predictors without seed constraint
		switch( model.val ) {
		case 'B':  {
			switch ( mode.val ) {
			case 'H' :	return new PredictorMfe2dHelixBlockHeuristic( energy, output, predTracker, getHelixConstraint(energy));
			default :  INTARNA_NOT_IMPLEMENTED("mode "+toString(mode.val)+" not available for model "+toString(model.val));
			}
		} break;
		// single-site mfe interactions (contain only interior loops)
		case 'S' : {
			switch ( mode.val ) {
			case 'H' :  return new PredictorMfe2dHeuristic( energy, output, predTracker );
			case 'M' :  return new PredictorMfe2d( energy, output, predTracker );
			default :  INTARNA_NOT_IMPLEMENTED("mode "+toString(mode.val)+" not available for model "+toString(model.val));
			}
		} break;
		// single-site mfe ensemble interactions (contain only interior loops)
		case 'P' : {
			switch ( mode.val ) {
			case 'H' :  return new PredictorMfeEns2dHeuristic( energy, output, predTracker );
			case 'M' :  return new PredictorMfeEns2d( energy, output, predTracker );
			default :  INTARNA_NOT_IMPLEMENTED("mode "+toString(mode.val)+" not available for model "+toString(model.val));
			}
		} break;
		// multi-site mfe interactions (contain interior and multi-loops loops)
		case 'M' : {
			switch ( mode.val ) {
			default :  INTARNA_NOT_IMPLEMENTED("mode "+toString(mode.val)+" not available for model "+toString(model.val));
			}
		} break;
		default : INTARNA_NOT_IMPLEMENTED("no model "+toString(model.val)+" available for mode "+toString(mode.val));
		}
	} else {
		// seed-constrained predictors
		switch( model.val ) {
		case 'B' : {
			switch  ( mode.val ) {
				case 'H' : return new PredictorMfe2dHelixBlockHeuristicSeed(energy, output, predTracker, getHelixConstraint(energy), getSeedHandler(energy));
				default :  INTARNA_NOT_IMPLEMENTED("mode "+toString(mode.val)+" not available for model "+toString(model.val));
			}
		} break;
		// single-site mfe interactions (contain only interior loops)
		case 'S' : {
			switch ( mode.val ) {
			case 'H' :  return new PredictorMfe2dHeuristicSeed( energy, output, predTracker, getSeedHandler( energy ) );
			case 'M' :  return new PredictorMfe2dSeed( energy, output, predTracker, getSeedHandler( energy ) );
			case 'S' :  return new PredictorMfeSeedOnly( energy, output, predTracker, getSeedHandler( energy ) );
			default :  INTARNA_NOT_IMPLEMENTED("mode "+toString(mode.val)+" not available for model "+toString(model.val));
			}
		} break;
		// single-site min energy interactions via seed extension (contain only interior loops)
		case 'X' : {
			switch ( mode.val ) {
			case 'H' :  return new PredictorMfe2dHeuristicSeedExtension( energy, output, predTracker, getSeedHandler( energy ) );
			case 'M' :  return new PredictorMfe2dSeedExtension( energy, output, predTracker, getSeedHandler( energy ) );
			case 'R' :  return new PredictorMfe2dSeedExtensionRIblast( energy, output, predTracker, getSeedHandler( energy ) );
			case 'S' :  return new PredictorMfeSeedOnly( energy, output, predTracker, getSeedHandler( energy ) );
			default  :  INTARNA_NOT_IMPLEMENTED("mode "+toString(mode.val)+" not implemented"); return NULL;
			}
		} break;
		// single-site max-prob interactions (contain only interior loops)
		case 'P' : {
			switch ( mode.val ) {
			case 'H' :  return new PredictorMfeEns2dHeuristicSeedExtension( energy, output, predTracker, getSeedHandler( energy ) );
			case 'M' :  return new PredictorMfeEns2dSeedExtension( energy, output, predTracker, getSeedHandler( energy ) );
			case 'S' :  return new PredictorMfeEnsSeedOnly( energy, output, predTracker, getSeedHandler(energy) );
			default :  INTARNA_NOT_IMPLEMENTED("mode "+toString(mode.val)+" not available for model "+toString(model.val));
			}
		} break;
		// multi-site mfe interactions (contain interior and multi-loops loops)
		case 'M' : {
			switch ( mode.val ) {
			default :  INTARNA_NOT_IMPLEMENTED("mode "+toString(mode.val)+" not available for model "+toString(model.val));
			}
		} break;
		default : INTARNA_NOT_IMPLEMENTED("no model "+toString(model.val)+" available for mode "+toString(mode.val));
		}
	}
}


////////////////////////////////////////////////////////////////////////////

void
CommandLineParsing::
initOutputHandler()
{
	// check if initial output needed
	switch (outMode.val) {
	case 'C' :
		outStreamHandler->getOutStream()
		<<OutputHandlerCsv::getHeader( OutputHandlerCsv::string2list( outCsvCols ), outSep )
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
		return new OutputHandlerText( getOutputConstraint(energy), outStreamHandler->getOutStream(), energy, 10, false );
	case 'D' :
		return new OutputHandlerText( getOutputConstraint(energy), outStreamHandler->getOutStream(), energy, 10, true );
	case 'E' :
		// ensure that Zall is computed
		outNeedsZall = true;
		// no interaction details needed
		outNeedsBPs = false;
		return new OutputHandlerEnsemble( getOutputConstraint(energy), outStreamHandler->getOutStream(), energy );
	case 'C' :
		// ensure that Zall is computed if needed
		outNeedsZall = outNeedsZall || OutputHandlerCsv::needsZall(OutputHandlerCsv::string2list( outCsvCols ));
		// check whether interaction details are needed
		outNeedsBPs = OutputHandlerCsv::needBPs(OutputHandlerCsv::string2list( outCsvCols ));;
		// create output handler
		return new OutputHandlerCsv( getOutputConstraint(energy), outStreamHandler->getOutStream(), energy, OutputHandlerCsv::string2list( outCsvCols ), outSep, false, outCsvLstSep );
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

const HelixConstraint &
CommandLineParsing::
getHelixConstraint(const InteractionEnergy &energy) const
{
	if (helixConstraint == NULL) {
		// setup according to user data
		helixConstraint = new HelixConstraint(
				  helixMinBP.val
				, helixMaxBP.val
			    , helixMaxIL.val
			    , ( (helixMinPu.val > Z_type(0)) ? std::min<E_type>(Accessibility::ED_UPPER_BOUND, energy.getE(helixMinPu.val)) : Accessibility::ED_UPPER_BOUND )
			    , Ekcal_2_E(helixMaxE.val)
				, helixFullE
		);
	}
	return *helixConstraint;
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
							, Ekcal_2_E(seedMaxE.val)
							, (seedMinPu.val>0 ? std::min<E_type>(Accessibility::ED_UPPER_BOUND, energy.getE( seedMinPu.val )) : Accessibility::ED_UPPER_BOUND) // transform unpaired prob to ED value
							, Ekcal_2_E(seedMaxEhybrid.val)
							// shift ranges to start counting with 0
							, IndexRangeList( seedTRange, false, &(energy.getAccessibility1().getSequence()) )
							, IndexRangeList( seedQRange, false, &(energy.getAccessibility2().getAccessibilityOrigin().getSequence()) ).reverse(energy.size2())
							, seedTQ
							, seedNoGU
							, seedNoGUend
							, outNoLP
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
		// check if we have to allow for bulges in seed
		if (seedConstr.getMaxUnpaired1()+seedConstr.getMaxUnpaired2()+seedConstr.getMaxUnpairedOverall() > 0) {
			// create new seed handler using mfe computation
			return new SeedHandlerMfe( energy, seedConstr );
		} else {
			// create new bulge-free seed handler
			return new SeedHandlerNoBulge( energy, seedConstr );
		}
	}
}

////////////////////////////////////////////////////////////////////////////

const IndexRangeList&
CommandLineParsing::
getQueryRanges( const InteractionEnergy & energy, const size_t sequenceNumber, const Accessibility & acc ) const
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
			qRegion[sequenceNumber] = acc.decomposeByMaxED( qRegionLenMax.val, seedBP.val, seedBP.val);
			// inform user
			VLOG(1) <<"detected accessible regions for query '"<<getQuerySequences().at(sequenceNumber).getId()<<"' : "<<qRegion.at(sequenceNumber);
		}
	}

	if (outMinPu.val > Z_type(0) && !Z_equal(outMinPu.val, Z_type(0))) {
		// decompose ranges based in minimal unpaired probability value per position
		// since all ranges covering a position will have a lower unpaired probability
		acc.decomposeByMaxED( qRegion[sequenceNumber], energy.getE( outMinPu.val ), (noSeedRequired ? RnaSequence::lastPos : seedBP.val ) );
	}

	return qRegion.at(sequenceNumber);
}

////////////////////////////////////////////////////////////////////////////

const IndexRangeList&
CommandLineParsing::
getTargetRanges( const InteractionEnergy & energy, const size_t sequenceNumber, const Accessibility & acc ) const
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
			tRegion[sequenceNumber] = acc.decomposeByMaxED( tRegionLenMax.val, seedBP.val, seedBP.val);
			// inform user
			VLOG(1) <<"detected accessible regions for target '"<<getTargetSequences().at(sequenceNumber).getId()<<"' : "<<tRegion.at(sequenceNumber);
		}
	}

	if (outMinPu.val > Z_type(0) && !Z_equal(outMinPu.val, Z_type(0))) {
		// decompose ranges based in minimal unpaired probability value per position
		// since all ranges covering a position will have a lower unpaired probability
		acc.decomposeByMaxED( tRegion[sequenceNumber], energy.getE( outMinPu.val ), (noSeedRequired ? RnaSequence::lastPos : seedBP.val ) );
	}

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

CommandLineParsing::Personality
CommandLineParsing::
getPersonality( int argc, char ** argv )
{
	// default : check via call name
	std::string value(argv[0]);
	// strip path if present
	size_t cutPos = value.find_last_of(R"(/\)");
	if (cutPos != std::string::npos) {
		value = value.substr(cutPos+1);
	}

	// check if respective parameter provided
	// --> overwrites default from call name
	const boost::regex paramNameRegex("^--personality.*", boost::regex::icase);
	const boost::regex paramRegex(R"(^--personality=\S+$)", boost::regex::icase);
	bool setViaParameter = false;
	for (int i=1; i<argc; i++) {
		std::string argi(argv[i]);
		// found parameter
		if (boost::regex_match(argi,paramNameRegex, boost::match_perl)) {
			if (boost::regex_match(argi,paramRegex, boost::match_perl)) {
				// cut off argument
				value = argi.substr( argi.find('=')+1 );
				setViaParameter = true;
			} else {
				// handle missing value
				throw std::runtime_error("--personality provided without argument (no spaces allowed!)");
			}
			break;
		}

	}

	// parse personality
	if (boost::regex_match(value,boost::regex("IntaRNAexact"), boost::match_perl)) {
		return Personality::IntaRNAexact;
	}
	if (boost::regex_match(value,boost::regex("IntaRNAduplex"), boost::match_perl)) {
		return Personality::IntaRNAduplex;
	}
	if (boost::regex_match(value,boost::regex("IntaRNAhelix"), boost::match_perl)) {
		return Personality::IntaRNAhelix;
	}
	if (boost::regex_match(value,boost::regex("IntaRNAseed"), boost::match_perl)) {
		return Personality::IntaRNAseed;
	}
	if (boost::regex_match(value,boost::regex("IntaRNAens"), boost::match_perl)) {
		return Personality::IntaRNAens;
	}
	if (boost::regex_match(value,boost::regex("IntaRNA3"), boost::match_perl)) {
		return Personality::IntaRNA3;
	}
	if (boost::regex_match(value,boost::regex("IntaRNA2"), boost::match_perl)) {
		return Personality::IntaRNA2;
	}
	if (boost::regex_match(value,boost::regex("IntaRNA1"), boost::match_perl)) {
		return Personality::IntaRNA1;
	}
	if (boost::regex_match(value,boost::regex("IntaRNAsTar"), boost::match_perl)) {
		return Personality::IntaRNAsTar;
	}

	if (setViaParameter) {
		// handle wrong value
		throw std::runtime_error("Personality '"+value+"' does not share my head ...");
	}

	// default behaviour
	return Personality::IntaRNA;
}


////////////////////////////////////////////////////////////////////////////



