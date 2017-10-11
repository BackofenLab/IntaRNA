
### std-options
# check that installed programs support the --help and --version options
# basic check that the program's run-time dependencies are satisfied
### subdir-objects
# Place generated object files (.o) into the same directory as their source
# files, in order to avoid collisions when non-recursive make is used.
AUTOMAKE_OPTIONS = std-options subdir-objects

AM_DEFAULT_SOURCE_EXT = .cpp

###############################################################################

# the programs to be build
bin_PROGRAMS = IntaRNA

###############################################################################
# THE INTARNA BINARY
###############################################################################

# generated intaRNA sources 
nodist_IntaRNA_SOURCES = 	\
					../config.h \
					IntaRNA/intarna_config.h

# intaRNA sources
IntaRNA_SOURCES =	IntaRNA.cpp \
					CommandLineParsing.h \
					CommandLineParsing.cpp
					
IntaRNA_LDADD = 	$(top_builddir)/src/IntaRNA/libIntaRNA.a \
					$(top_builddir)/src/libeasylogging.a

###############################################################################

