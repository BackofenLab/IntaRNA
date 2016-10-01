
# the program to build
bin_PROGRAMS = runTests 

# test sources
runTests_SOURCES =	catch.hpp \
					$(wildcard *_test.cpp) \
					runTests.cpp

# add source objects for linking
runTests_LDFLAGS = `ls ../src/*.o | tr " " "\n" | grep -v "intaRNA" | tr "\n" " "` @AM_LDFLAGS@ @LDFLAGS@

# add source include for compilation
runTests_CXXFLAGS =	-I$(top_srcdir)/src @AM_CXXFLAGS@ @CXXFLAGS@

# disable general build
all:

# global test directive
tests-local: all-am
	@./runTests
