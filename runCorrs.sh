#!/bin/bash

if [[ $# -eq 0 ]] ; then
    echo 'please give [{13d, 13e, 17q,...}.root [ "" or Parallel]] as input'
    exit 0
fi

date

module load root
module load hdf5

n_mix=300

for p in {0,4}

do
    name=${1}
    #name=$(basename ../${1}..root)
    echo $name

    if [ "17q" = "$name" ]; then #account for 17q-17p mixing
	echo "./Parallel_Mix_Correlations ../InputData/${name}_${p}GeVTrack_paired.root ../InputData/17p_MB_${p}GeV.hdf5 0 $n_mix $p"
	./Parallel_Mix_Correlations ../InputData/${name}_${p}GeVTrack_paired.root ../InputData/17p_MB_${p}GeV.hdf5 0 $n_mix $p
    else
	echo "./Parallel_Mix_Correlations ../InputData/${name}_${p}GeVTrack_paired.root ../InputData/${name}_MB_${p}GeV.hdf5 0 $n_mix $p"
	./Parallel_Mix_Correlations ../InputData/${name}_${p}GeVTrack_paired.root ../InputData/${name}_MB_${p}GeV.hdf5 0 $n_mix $p
    fi
    echo "$nmix Mixed Events DATASET =  $name"
    echo
done
