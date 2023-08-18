#!/bin/bash

cd build || exit
module load compiler/gnu/10.2
module load devel/cmake
module load mpi/openmpi/4.1
cmake -DCMAKE_BUILD_TYPE=Release ..
build
module purge
