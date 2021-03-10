#!/bin/bash

export HDF5_USE_FILE_LOCKING=FALSE
#./parallel_pair_gale_shapley_hdf5 ../InputData/13f_new_skimClusterMinE12.root ../InputData/13f_MB_0GeV.hdf5 0 300 0
./parallel_pair_gale_shapley_hdf5 ../InputData/13f_new_skimClusterMinE12.root ../InputData/13f_MB_4GeV.hdf5 0 300 4

