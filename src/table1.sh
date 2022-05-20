#!/usr/bin/env bash

DATA_DIR=$1
TE=$2
WE=$3
OUTPUT_FILE=$4

rm -f tmp0_tauQ ${OUTPUT_FILE}
for j in 2 4 6 8; do
	for i in $(grep "^${j}" ${DATA_DIR}/DATA_FILES | cut -d, -f5)
    do
		python src/Topo_binning.py ${TE} ${WE} $i >> tmp0_tauQ
		grep "${i}" DATA_FILES | cut -d, -f4 | paste -d" " tmp0_tauQ - > tmp1_tauQ
		awk '{print "$"$1"$ & $"$2"$ & $"$3"$ & $"$4"$ \\\\"}' tmp1_tauQ >> OUTPUT_FILE
		rm tmp0_tauQ tmp1_tauQ
	done
done
