
### std-options
# check that installed programs support the --help and --version options
# basic check that the program's run-time dependencies are satisfied
### subdir-objects
# Place generated object files (.o) into the same directory as their source
# files, in order to avoid collisions when non-recursive make is used.
AUTOMAKE_OPTIONS = std-options subdir-objects

AM_DEFAULT_SOURCE_EXT = .cpp

###############################################################################
# THE INTARNA LIBRARY
###############################################################################

# the library names to build (note we are building static libs only)
lib_LIBRARIES = libIntaRNA.a

# where to install the headers on the system
libIntaRNA_adir = $(includedir)/IntaRNA

# the list of header files that belong to the library (to be installed later)
libIntaRNA_a_HEADERS =    \
					intarna_config.h \
					general.h \
					Accessibility.h \
					AccessibilityConstraint.h \
					AccessibilityDisabled.h \
					AccessibilityFromStream.h \
					AccessibilityVrna.h \
					AccessibilityBasePair.h \
					IndexRange.h \
					IndexRangeList.h \
					Interaction.h \
					InteractionEnergy.h \
					InteractionEnergyBasePair.h \
					InteractionEnergyIdxOffset.h \
					InteractionEnergyVrna.h \
					InteractionRange.h \
					NussinovHandler.h \
					OutputConstraint.h \
					OutputHandler.h \
					OutputHandlerCsv.h \
					OutputHandlerHub.h \
					OutputHandlerIntaRNA1.h \
					OutputHandlerInteractionList.h \
					OutputHandlerRangeOnly.h \
					OutputHandlerText.h \
					PredictionTracker.h \
					PredictionTrackerHub.h \
					PredictionTrackerPairMinE.h \
					PredictionTrackerProfileMinE.h \
					PredictionTrackerSpotProb.h \
					PredictionTrackerProfileSpotProb.h \
					Predictor.h \
					PredictorMaxProb.h \
					PredictorMfe.h \
					PredictorMfe2d.h \
					PredictorMfe2dSeed.h \
					PredictorMfe2dHeuristic.h \
					PredictorMfe2dHeuristicSeed.h \
					PredictorMfe4d.h \
					PredictorMfe4dSeed.h \
					ReverseAccessibility.h \
					RnaSequence.h \
					SeedConstraint.h \
					SeedHandler.h \
					SeedHandlerExplicit.h \
					SeedHandlerIdxOffset.h \
					SeedHandlerMfe.h \
					VrnaHandler.h

# the sources to add to the library and to add to the source distribution
libIntaRNA_a_SOURCES =    \
                    $(libmylib_a_HEADERS) \
					general.cpp \
					Accessibility.cpp \
					AccessibilityConstraint.cpp \
					AccessibilityFromStream.cpp \
					AccessibilityVrna.cpp \
					AccessibilityBasePair.cpp \
					IndexRange.cpp \
					IndexRangeList.cpp \
					Interaction.cpp \
					InteractionEnergy.cpp \
					InteractionEnergyBasePair.cpp \
					InteractionEnergyIdxOffset.cpp \
					InteractionEnergyVrna.cpp \
					InteractionRange.cpp \
					NussinovHandler.cpp \
					OutputConstraint.cpp \
					OutputHandler.cpp \
					OutputHandlerCsv.cpp \
					OutputHandlerIntaRNA1.cpp \
					OutputHandlerInteractionList.cpp \
					OutputHandlerRangeOnly.cpp \
					OutputHandlerText.cpp \
					PredictionTrackerPairMinE.cpp \
					PredictionTrackerProfileMinE.cpp \
					PredictionTrackerSpotProb.cpp \
					PredictionTrackerProfileSpotProb.cpp \
					PredictorMaxProb.cpp \
					PredictorMfe.cpp \
					PredictorMfe2d.cpp \
					PredictorMfe2dSeed.cpp \
					PredictorMfe2dHeuristic.cpp \
					PredictorMfe2dHeuristicSeed.cpp \
					PredictorMfe4d.cpp \
					PredictorMfe4dSeed.cpp \
					ReverseAccessibility.cpp \
					RnaSequence.cpp \
					SeedConstraint.cpp \
					SeedHandlerExplicit.cpp \
					SeedHandlerMfe.cpp \
					VrnaHandler.cpp


# generated intaRNA sources 
nodist_libIntaRNA_a_SOURCES = 	\
					../config.h \
					intarna_config.h
						
################################################################
# pkg-config
################################################################

if enable_pkg_config
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = $(top_srcdir)/IntaRNA.pc
endif

###############################################################################

