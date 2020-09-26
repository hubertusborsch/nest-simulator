#!/bin/bash/
#git_hash=$(git rev-parse --short=7 HEAD)
. /home/bouhadjar/miniconda3/etc/profile.d/conda.sh
conda activate nest_sim
mkdir -p $PWD/build
cd $PWD/build
 
# build --> '4' 
# hash --> begin CMAKE $PWD/../git_hash
# source --> last CMAKE ../

#cmake -DCMAKE_INSTALL_PREFIX=$PWD/ -Dwith-python=3 -Dwith-mpi=ON -DCMAKE_CXX_COMPILER=mpic++ -DCMAKE_C_COMPILER=mpicc ../../nest-simulator/ && make -j8 && make install
#cmake -DCMAKE_INSTALL_PREFIX=$PWD/ -Dwith-mpi=ON -DCMAKE_CXX_COMPILER=mpic++ -DCMAKE_C_COMPILER=mpicc ../../nest-simulator/ && make -j8 && make install && make installcheck
cmake -DCMAKE_INSTALL_PREFIX=$PWD/ -Dwith-gsl=$CONDA_PREFIX -DREADLINE_ROOT_DIR=$CONDA_PREFIX -DLTDL_ROOT_DIR=$CONDA_PREFIX -Dwith-mpi=ON ../ && make -j8 && make install

