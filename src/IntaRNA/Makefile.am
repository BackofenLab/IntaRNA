
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
					Accessibility.h \
					AccessibilityConstraint.h \
					AccessibilityDisabled.h \
					AccessibilityFromStream.h \
					AccessibilityVrna.h \
					AccessibilityBasePair.h \
					HelixConstraint.h \
					HelixHandler.h \
					HelixHandlerIdxOffset.h \
					HelixHandlerNoBulgeMax.h \
					HelixHandlerUnpaired.h \
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
					OutputStreamHandler.h \
					OutputStreamHandlerSortedCsv.h \
					OutputHandler.h \
					OutputHandlerCsv.h \
					OutputHandlerEnsemble.h \
					OutputHandlerHub.h \
					OutputHandlerInteractionList.h \
					OutputHandlerText.h \
					PredictionTracker.h \
					PredictionTrackerHub.h \
					PredictionTrackerPairMinE.h \
					PredictionTrackerProfileMinE.h \
					PredictionTrackerSpotProb.h \
					PredictionTrackerSpotProbAll.h \
					PredictionTrackerProfileSpotProb.h \
					Predictor.h \
					PredictorMfe.h \
					PredictorMfeSeedOnly.h \
					PredictorMfe2d.h \
					PredictorMfe2dSeed.h \
					PredictorMfe2dSeedExtension.h \
					PredictorMfe2dSeedExtensionRIblast.h \
					PredictorMfe2dHeuristic.h \
					PredictorMfe2dHeuristicSeed.h \
					PredictorMfe2dHelixBlockHeuristic.h \
					PredictorMfe2dHelixBlockHeuristicSeed.h \
					PredictorMfe2dHeuristicSeedExtension.h \
					PredictorMfeEns.h \
					PredictorMfeEnsSeedOnly.h \
					PredictorMfeEns2d.h \
					PredictorMfeEns2dHeuristic.h \
					PredictorMfeEns2dSeedExtension.h \
					PredictorMfeEns2dHeuristicSeedExtension.h \
					ReverseAccessibility.h \
					RnaSequence.h \
					SeedConstraint.h \
					SeedHandler.h \
					SeedHandlerExplicit.h \
					SeedHandlerIdxOffset.h \
					SeedHandlerMfe.h \
					SeedHandlerNoBulge.h \
					VrnaHandler.h \
					general.h


# the sources to add to the library and to add to the source distribution
libIntaRNA_a_SOURCES =    \
                    $(libmylib_a_HEADERS) \
					Accessibility.cpp \
					AccessibilityConstraint.cpp \
					AccessibilityFromStream.cpp \
					AccessibilityVrna.cpp \
					AccessibilityBasePair.cpp \
					HelixHandler.cpp \
					HelixHandlerNoBulgeMax.cpp \
				    HelixHandlerUnpaired.cpp \
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
					OutputStreamHandlerSortedCsv.cpp \
					OutputHandler.cpp \
					OutputHandlerEnsemble.cpp \
					OutputHandlerCsv.cpp \
					OutputHandlerInteractionList.cpp \
					OutputHandlerText.cpp \
					PredictionTrackerPairMinE.cpp \
					PredictionTrackerProfileMinE.cpp \
					PredictionTrackerSpotProb.cpp \
					PredictionTrackerSpotProbAll.cpp \
					PredictionTrackerProfileSpotProb.cpp \
					PredictorMfe.cpp \
					PredictorMfeSeedOnly.cpp \
					PredictorMfe2d.cpp \
					PredictorMfe2dSeed.cpp \
					PredictorMfe2dSeedExtension.cpp \
					PredictorMfe2dSeedExtensionRIblast.cpp \
					PredictorMfe2dHeuristic.cpp \
					PredictorMfe2dHeuristicSeed.cpp \
					PredictorMfe2dHelixBlockHeuristic.cpp \
					PredictorMfe2dHelixBlockHeuristicSeed.cpp \
					PredictorMfe2dHeuristicSeedExtension.cpp \
					PredictorMfeEns.cpp \
					PredictorMfeEns2d.cpp \
					PredictorMfeEns2dHeuristic.cpp \
					PredictorMfeEnsSeedOnly.cpp \
					PredictorMfeEns2dSeedExtension.cpp \
					PredictorMfeEns2dHeuristicSeedExtension.cpp \
					ReverseAccessibility.cpp \
					RnaSequence.cpp \
					SeedHandler.cpp \
					SeedHandlerExplicit.cpp \
					SeedHandlerMfe.cpp \
					SeedHandlerNoBulge.cpp \
					VrnaHandler.cpp \
					general.cpp


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
