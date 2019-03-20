
conda activate build-IntaRNA

$CXX --version

which $CXX

autotools-init.sh 

./configure --prefix=$HOME/IntaRNA --with-vrna=${CONDA_PATH}

make -j 2 && make tests -j 2 && make install