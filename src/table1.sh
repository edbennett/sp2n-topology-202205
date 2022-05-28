#!/usr/bin/env bash

DATA_DIR=$1
PROC_DIR=$2
TE=$3
WE=$4
OUTPUT_FILE=$5

rm -f tmp0_tauQ ${OUTPUT_FILE}

cat > ${OUTPUT_FILE} <<EOF
\\begin{tabular}{|c|c|c|c|c|c|}
\\hline
\\hline
\$N_c\$ & \$L/a\$ & \$\\beta\$ & \$N_{\\mathrm{tot}}\$ & \$N_{\\mathrm{sw}}\$ & \$\\tau_Q\$ \\\\
EOF
for j in 2 4 6 8; do
	echo '\hline' >> ${OUTPUT_FILE}
	for i in $(grep "^${j}" ${DATA_DIR}/DATA_FILES | cut -d, -f5)
    do
		python src/Topo_binning.py ${TE} ${WE} ${PROC_DIR}/$i >> tmp0_tauQ
		grep "${i}\$" ${DATA_DIR}/DATA_FILES | cut -d, -f4 | paste -d" " tmp0_tauQ <(grep -c 0\\.000000 ${PROC_DIR}/${i}) - > tmp1_tauQ
		awk '{print "$"$1"$ & $"$2"$ & $"$3"$ & $"$5"$ & $"$6"$ & $"$4"$ \\\\"}' tmp1_tauQ >> ${OUTPUT_FILE}
		rm tmp0_tauQ tmp1_tauQ
	done
done

cat >> ${OUTPUT_FILE} <<EOF
\\hline
\\end{tabular}
EOF
