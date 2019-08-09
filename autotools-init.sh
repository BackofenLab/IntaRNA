#!/usr/bin/env bash
# 
# Run this before configure
#

set -e -o pipefail

#test -d config || mkdir config
# Produce aclocal.m4, so autoconf gets the automake macros it needs
echo "Creating aclocal.m4..."
aclocal || exit $?;

# Produce config.h raw file
autoheader || exit $?;

# Produce all the `Makefile.in's, verbosely, and create neat missing things
# like `libtool', `install-sh', etc.
automake --add-missing --gnu --copy || exit $?;

# If there's a config.cache file, we may need to delete it.  
# If we have an existing configure script, save a copy for comparison.
if [ -f config.cache ] && [ -f configure ]; then
  cp configure configure.$$.tmp
fi

# Produce ./configure
echo "Creating configure..."
autoconf || exit $?;

echo ""
echo "You can run ./configure [--prefix=$HOME] now."
echo ""
