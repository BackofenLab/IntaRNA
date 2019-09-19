#!/usr/bin/env bash


if [[ -z "${GENOUT}" ]]; then
	GENERATE_OUTPUT=false
else
	GENERATE_OUTPUT=true
fi

# subdirectory of test data
datadir=data

# list of all test names (parameter file in data subfolder) 
IntaRNA_tests=`cd $datadir; for f in *.parameter; do printf "${f%.*} "; done`

DIFFERENCES=false

## ----------------------------------------
## call test and compare output
##
## @param $1 the name of the test; shall not contain special symbols
## and white space, since this is used to derive file names
## @param $2 target directory, which is wiped after the test
## @param $3 options for diff
##
function calltest {
    name="$1"
    diffopts="${2:-"--normal"}"

    resultfile=$datadir/$name.testout
    reference_resultsfile=$datadir/$name.testresult
    
    difftmp=testDiff.tmp

    shift 2
    
    testCall="$INTARNABINPATH/src/bin/IntaRNA --parameterFile=$datadir/$name.parameter --default-log-file=/dev/null"

    echo " IntaRNA TEST $name"
	#	echo CALL $testCall
    #echo
    
    $testCall 2>&1 | grep -v INFO > $resultfile

	if [ -e "$reference_resultsfile" ] ; then
	    if ! diff "$reference_resultsfile" "$resultfile" "${diffopts}" > $difftmp; then
	        #	        echo "==================== OK"
	    #else
	        DIFFERENCES=true
		    echo "============================================================"
		    cat $difftmp
	        echo "==================== DIFFERENT"
	    fi
	    rm -f $difftmp
	else
	    echo "============================================================"
	    echo "WARNING: file '$reference_resultsfile' does not exist!"
	    echo "==================== NO_REFERENCE"
		DIFFERENCES=true
	
        if $GENERATE_OUTPUT ; then
            echo "Write new reference '$reference_resultsfile'."
    		\cp $resultfile "$reference_resultsfile"
        fi
	fi

}



# run all tests
for t in $IntaRNA_tests; do
	calltest $t
done

# check if differences in diffs
if $DIFFERENCES ; then
    exit -1
fi
