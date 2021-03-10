#!/bin/bash

date

module load root
export HDF5_USE_FILE_LOCKING=FALSE
n_mix=300

for dataset in {13fnew_skimClusterMinE12,13f,13d,13e,17q} #13f first since troublesome (only uses half block size)
#for dataset in {13f_new_skimClusterMinE12,13}
#dataset="13f_new_skimClusterMinE12"
#for dataset in {13f_new_skimClusterMinE12,17q}
do
    for i in {1}
			
    do
	
	for p in {0,4}
		 
	do
	    name=${dataset}
	    echo $name
	    
	    if [ "17q" = "$name" ]; then #account for 17q-17p mixing
		echo "./Parallel_Mix_Correlations_pWeight ../InputData/${name}_${p}GeVTrack_paired.root ../InputData/17p_MB_${p}GeV.hdf5 0 $n_mix $p pp"
		./Parallel_Mix_Correlations_pWeight ../InputData/${name}_${p}GeVTrack_paired.root ../InputData/17p_MB_${p}GeV.hdf5 0 $n_mix $p pp
		
	    elif [ "13fnew_skimClusterMinE12" = "$name" ];
	    then
		echo "./Parallel_Mix_Correlations_pWeight ../InputData/${name}_${p}GeVTrack_paired.root ../InputData/13f_MB_${p}GeV.hdf5 0 $n_mix $p pPb"
                ./Parallel_Mix_Correlations_pWeight ../InputData/${name}_${p}GeVTrack_paired.root ../InputData/13f_MB_${p}GeV.hdf5 0 $n_mix $p pPb
	    else
		echo "./Parallel_Mix_Correlations_pWeight ../InputData/${name}_${p}GeVTrack_paired.root ../InputData/${name}_MB_${p}GeV.hdf5 0 $n_mix $p pPb"
		./Parallel_Mix_Correlations_pWeight ../InputData/${name}_${p}GeVTrack_paired.root ../InputData/${name}_MB_${p}GeV.hdf5 0 $n_mix $p pPb
	    fi
	    echo "$nmix Mixed Events DATASET =  $name"
	    echo
	done
	
    done
done
mail -s "Mixing Done" ftoralesacosta@lbl.gov < ./
