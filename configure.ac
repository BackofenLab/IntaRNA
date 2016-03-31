
AC_PREREQ(2.59)

AC_INIT([IntaRNA], [1.2.5], [intarna@informatik.uni-freiburg.de])

AC_CONFIG_SRCDIR(src/energy.h)
AC_CONFIG_AUX_DIR([.])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

AM_INIT_AUTOMAKE

AC_LANG_CPLUSPLUS

###############################################################################
# BEGIN PATH-SUPPORT CHECK
###############################################################################

# Vienna RNA package library path support, if not installed in usual directories
AC_ARG_WITH(RNA,
	AC_HELP_STRING(
	    [--with-RNA=PREFIX],
	    [alternative prefix path to Vienna RNA library]
	  ),
	  RNAPATHSET=1,
	  RNAPATHSET=0
)
if test  $RNAPATHSET = 1 ; then
		CXXFLAGS="-I$with_RNA/include $CXXFLAGS"
		LDFLAGS="-L$with_RNA/lib $LDFLAGS"
fi

# add the Vienna RNA lib to the libraries used
LIBS="-lRNA $LIBS"

###############################################################################
# END PATH-SUPPORT CHECK
###############################################################################

# add all libraries to the LDFLAGS for linking
LDFLAGS="$LDFLAGS $LIBS"

###############################################################################
# BEGIN MISC CHECKS
###############################################################################

# Checks for header files.
AC_HEADER_STDC

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T

# Check if everything should be compiled in DEBUG mode
AC_MSG_CHECKING(whether to build in debug mode with assertions)
debuger=no
AC_ARG_ENABLE(debug,
	AC_HELP_STRING(
	    [--enable-debug],
	    [enable debug mode and assertions (default=disabled)]
	),
  	debuger="$enableval"
)
AC_MSG_RESULT($debuger)
if test x"$debuger" = x"yes"; then
	AC_DEFINE(DEBUG)
	CXXFLAGS="$CXXFLAGS -g -Wall -Wno-uninitialized -O0"
	LDFLAGS="$LDFLAGS -g -Wall -Wno-uninitialized -O0"
else
	AC_DEFINE(NDEBUG)
	CXXFLAGS="-O2 $CXXFLAGS"
fi

###############################################################################
# END MISC CHECKS
###############################################################################

###############################################################################
# BEGIN VIENNA CHECK
###############################################################################

# check for Vienna RNA headers
AC_MSG_CHECKING([for the Vienna RNA C header])
AC_COMPILE_IFELSE(
     [ #include <ViennaRNA/fold_vars.h> ],
     [ 
		AC_MSG_RESULT(yes)
		AC_CHECK_LIB(RNA, [energy_of_struct], [RNANOTFOUND=0], [ RNANOTFOUND=1])
     ],
     [
       AC_MSG_RESULT(no)
       RNANOTFOUND=1;
     ]
)

###############################################################################
# END VIENNA CHECK
###############################################################################


###############################################################################
# BEGIN DO NOT COMPILE CHECKS
###############################################################################

# flag for final abort if one of the dependency was not found but all were checked
DEPENDENCYNOTFOUND=0;

# error output if ViennaRNA not found
if test "$RNANOTFOUND" = "1"; then
	AC_MSG_NOTICE()
	AC_MSG_NOTICE(the Vienna RNA C library is required for building the program.)
	AC_MSG_NOTICE(it can be obtained from http://www.tbi.univie.ac.at/~ivo/RNA/.)
	AC_MSG_NOTICE()
 	if test "$RNAPATHSET" = "1"; then
		AC_MSG_NOTICE(cant find Vienna RNA library in given path '$with_RNA'.)
 	else
		AC_MSG_NOTICE(if Vienna RNA is installed elsewhere use --with-RNA=PREFIX.)
 	fi
 	DEPENDENCYNOTFOUND=1;
fi

# error ABORT if on of the libraries was not found
if test "$DEPENDENCYNOTFOUND" = "1"; then
	AC_MSG_NOTICE()
	AC_MSG_ERROR(library not found! Try './configure --help'.)
fi

###############################################################################
# END DO NOT COMPILE CHECKS
###############################################################################



AC_CONFIG_FILES([ Makefile ])
AC_CONFIG_FILES([ src/Makefile ])

AC_OUTPUT
