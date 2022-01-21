#!/bin/bash/
. /home/hubertus/miniconda3/etc/profile.d/conda.sh
conda activate spiking-htm
#git_hash=$(git rev-parse --short=7 HEAD)
git_branch=$(git symbolic-ref --short HEAD)
mkdir -p $PWD/build/$git_branch
cd $PWD/build/$git_branch
cmake -DCMAKE_INSTALL_PREFIX=$PWD/ -Dwith-gsl=$CONDA_PREFIX -DREADLINE_ROOT_DIR=$CONDA_PREFIX -DLTDL_ROOT_DIR=$CONDA_PREFIX -Dwith-mpi=ON ../../ && make -j8 && make install
