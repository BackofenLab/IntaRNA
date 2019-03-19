#!/usr/bin/env bash


if [[ -z "${GENOUT}" ]]; then
	GENERATE_OUTPUT=false
else
	GENERATE_OUTPUT=true
fi

# list of all test names (parameter and output file in data subfolder) 
IntaRNA_tests="\
				noLP-exact \
				noLP-exact-seed \
				"

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
    
    datadir=data

    resultfile=$datadir/$name.testout
    reference_resultsfile=$datadir/$name.testresult

    shift 2
    
    testCall="$INTARNABINPATH/src/bin/IntaRNA --parameterFile=$datadir/$name.parameter"

    echo "============================================================"
    echo TEST $name
	#	echo CALL $testCall
    #echo
    
    $testCall 2>&1 > $resultfile

	if [ -e "$reference_resultsfile" ] ; then
	    if diff "$resultfile" "$reference_resultsfile" "${diffopts}"; then
	        echo "==================== OK"
	    else
	        DIFFERENCES=true
	        echo "==================== DIFFERENT"
	    fi
	else
	    echo "WARNING: file '$reference_resultsfile' does not exist!"
	    echo "==================== NO_REFERENCE"
	            DIFFERENCES=true
	        fi
	
	        if $GENERATE_OUTPUT ; then
	            echo "Write new reference '$reference_resultsfile'."
	    \cp $resultfile "$reference_resultsfile"
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