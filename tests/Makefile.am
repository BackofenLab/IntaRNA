
# sub directories to check for Makefiles
SUBDIRS = data .

# ensure parallel test build and run
AUTOMAKE_OPTIONS = parallel-tests

# setup environment needed for the tests
AM_TESTS_ENVIRONMENT = export INTARNABINPATH=$(abs_top_builddir);

# list of tests to run
TESTS = runApiTests$(EXEEXT) $(dist_check_SCRIPTS)
# how to handle the tests
TEST_EXTENSIONS = $(EXEEXT) .sh
SH_LOG_COMPILER = $(SHELL)

# the script needed for tests
dist_check_SCRIPTS = runIntaRNA.sh

# the program to build
check_PROGRAMS = runApiTests

# test sources
runApiTests_SOURCES =	\
					catch.hpp \
					testEasyLoggingSetup.icc \
					AccessibilityConstraint_test.cpp \
					AccessibilityFromStream_test.cpp \
					AccessibilityBasePair_test.cpp \
					AccessibilityVrna_test.cpp \
					HelixConstraint_test.cpp \
					HelixHandlerNoBulgeMax_test.cpp \
					HelixHandlerNoBulgeMaxSeed_test.cpp \
					HelixHandlerNoBulgeMaxIdxOffset_test.cpp \
					HelixHandlerNoBulgeMaxSeedIdxOffset_test.cpp \
					HelixHandlerUnpaired_test.cpp \
					HelixHandlerUnpairedIdxOffset_test.cpp \
					HelixHandlerUnpairedSeed_test.cpp \
					HelixHandlerUnpairedSeedIdxOffset_test.cpp \
					IndexRange_test.cpp  \
					IndexRangeList_test.cpp  \
					Interaction_test.cpp  \
					InteractionEnergyBasePair_test.cpp  \
					InteractionRange_test.cpp  \
					PredictionTrackerBasePairProb_test.cpp \
					PredictionTrackerProfileMinE_test.cpp \
					PredictionTrackerSpotProb_test.cpp \
					PredictorMfe2dHelixBlockHeuristic_test.cpp \
					PredictorMfe2dHelixBlockHeuristicSeed_test.cpp \
					PredictorMfeEns2d_test.cpp \
					NussinovHandler_test.cpp \
					RnaSequence_test.cpp \
					OutputStreamHandlerSortedCsv_test.cpp \
					OutputHandlerInteractionList_test.cpp \
					SeedHandlerExplicit_test.cpp \
					SeedHandlerNoBulge_test.cpp \
					SeedHandlerMfe_test.cpp \
					SeedHandlerIdxOffset_test.cpp \
					runApiTests.cpp

# add IntaRNA lib for linking
LIBS= 	-L$(top_builddir)/src/IntaRNA -lIntaRNA \
		-L$(top_builddir)/src -leasylogging \
		@LIBS@

# add source include for compilation
runApiTests_CXXFLAGS =	-I$(top_builddir)/src \
						@AM_CXXFLAGS@ @CXXFLAGS@ \
						-DELPP_NO_LOG_TO_FILE
