#!/usr/bin/env bash

$CXX --version

which $CXX

echo "#### activating build-IntaRNA"

#conda activate build-IntaRNA
source ${CONDA_PATH}/bin/activate build-IntaRNA

echo "#### build-IntaRNA activated"

$CXX --version

which $CXX

source autotools-init.sh 

./configure --prefix=$HOME/IntaRNA --with-vrna=${CONDA_PATH}

make -j 2 && make tests -j 2 && make install