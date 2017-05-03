
# ensure serial test build and run
AUTOMAKE_OPTIONS = serial-tests

# list of tests to run
TESTS = runTests


# the program to build
check_PROGRAMS = runTests 

# test sources
runTests_SOURCES =	catch.hpp \
					AccessibilityConstraint_test.cpp \
					AccessibilityFromStream_test.cpp \
					AccessibilityBasePair_test.cpp \
					IndexRange_test.cpp  \
					IndexRangeList_test.cpp  \
					Interaction_test.cpp  \
					InteractionEnergyBasePair_test.cpp  \
					InteractionRange_test.cpp  \
					PredictionTrackerProfileMinE_test.cpp \
					NussinovHandler_test.cpp \
					RnaSequence_test.cpp \
					OutputHandlerRangeOnly_test.cpp \
					runTests.cpp

# add IntaRNA lib for linking
LIBS= -L$(top_builddir)/src/IntaRNA -lIntaRNA @LIBS@

# add source include for compilation
runTests_CXXFLAGS =	-I$(top_builddir)/src @AM_CXXFLAGS@ @CXXFLAGS@
