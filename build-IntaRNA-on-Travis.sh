
$CXX --version

which $CXX

bash autotools-init.sh 

./configure --prefix=$HOME/IntaRNA --with-vrna=${CONDA_PATH}

make -j 2 && make tests -j 2 && make install