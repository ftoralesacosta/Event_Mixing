#!/bin/bash
rm parallel_pair_gale_shapley_hdf5
make parallel_pair_gale_shapley_hdf5
export HDF5_USE_FILE_LOCKING=FALSE
echo "./parallel_pair_gale_shapley_hdf5 15o_246001.root 15o_246001.hdf5 $1 $2 0"
./parallel_pair_gale_shapley_hdf5 15o_246001.root 15o_246001.hdf5 $1 $2 0
