#!/bin/bash

module load root
make Divide_Same_Mix

for data in {"13d","13e","13f","13fnew","17q"}
#for data in {"13d","13e","13f","17q_wSDD"}
#for data in {"13d"}
do
    echo "${1}${data}_SE_L0_Correlation.root"
    ./Divide_Same_Mix ${1}${data}_SE_L0_Correlation.root
done

cd $1
hadd -f pPb_SE_L0_Correlation_GMB_Ratio.root 13*_SE_L0_Correlation_GMB_Ratio.root
cp 17q_SE_L0_Correlation_GMB_Ratio.root pp_SE_L0_Correlation_GMB_Ratio.root

echo "DIVIDE CORRELATIONS AND CREATED pPb and pp Files"
