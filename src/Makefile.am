
# sub directories to check for Makefiles
SUBDIRS = . IntaRNA bin

################################################################
# extra files to be added to distribution
################################################################

EXTRA_DIST = easylogging++.LICENCE.txt

###############################################################################
# THE INTARNA LIBRARY
###############################################################################

# the library names to build (note we are building static libs only)
lib_LIBRARIES = libeasylogging.a

# where to install the headers on the system
libeasylogging_adir = $(includedir)

# the list of header files that belong to the library (to be installed later)
libeasylogging_a_HEADERS = easylogging++.h

# the sources to add to the library and to add to the source distribution
libeasylogging_a_SOURCES = easylogging++.cc

###############################################################################

