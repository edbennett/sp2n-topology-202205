DDIR=/home/davide/documenti/physics/montecarlo/data_WF
TE=0.5
WE=0.5

#TABLE 1
rm tmp0_tauQ
for j in 2 4 6 8; do
	for i in `grep "^${j}" ${DDIR}/DATA_FILES | cut -d, -f5`; do 
		python scripts/Topo_binning.py ${TE} ${WE} $i >> tmp0_tauQ
		grep "${i}" DATA_FILES | cut -d, -f4 | paste -d" " tmp0_tauQ - > tmp1_tauQ
		awk '{print "$"$1"$ & $"$2"$ & $"$3"$ & $"$4"$ \\\\"}' tmp1_tauQ
		rm tmp0_tauQ tmp1_tauQ
	done
done

python produce_bs_sample.py WF_2_20_2.55 


#FIGURE 1
python scripts/vis_WF_W_E.py ${TE} ${WE} ${DDIR}/WF_6_{16_15.9,16_16.1}
#python scripts/vis_WF_W_E.py ${TE} ${WE} ${DDIR}/WF_2_{20_2.55,24_2.60,32_2.65,32_2.70}
#python scripts/vis_WF_W_E.py ${TE} ${WE} ${DDIR}/WF_4_{20_{7.7,7.72,7.76,7.78,7.80,7.85},24_8.2}
#python scripts/vis_WF_W_E.py ${TE} ${WE} ${DDIR}/WF_6_{18_15.75,16_15.9,16_16.1,20_16.3}
#python scripts/vis_WF_W_E.py ${TE} ${WE} ${DDIR}/WF_8_16_{26.5,26.7,27.0,27.2}

# FIGURE 2, 12 PLUS TABLES 2, 3, 7, 8, 9
python scripts/vis_WF_Scale.py ${TE} ${WE} ${DDIR}/WF_2_{20_2.55,24_2.60,32_2.65,32_2.70} > Scale_2.dat
python scripts/vis_WF_Scale.py ${TE} ${WE} ${DDIR}/WF_4_20_{7.7,7.72,7.76,7.78,7.80,7.85} > Scale_4.dat
python scripts/vis_WF_Scale.py ${TE} ${WE} ${DDIR}/WF_6_{18_15.75,16_15.9,16_16.1,20_16.3} > Scale_6.dat
python scripts/vis_WF_Scale.py ${TE} ${WE} ${DDIR}/WF_8_16_{26.5,26.7,27.0,27.2} > Scale_8.dat

for i in 2 4 6 8; do
	awk '{print "$",$3,"$&$",$4,"$&$",$5,"$&$",$6,"$ \\\\"}' Scale_${i}.dat > Scale_${i}_t0.dat
	sed -i -e '$s/\\\\//' Scale_${i}_t0.dat
	awk '{print "$",$3,"$&$",$7,"$&$",$8,"$&$",$9,"$ \\\\"}' Scale_${i}.dat > Scale_${i}_w0.dat
	sed -i -e '$s/\\\\//' Scale_${i}_w0.dat
done

# FIGURE 3, 4
python scripts/vis_WF_N_scaling.py ${TE} ${WE} ${DDIR}/WF_{{2_20_2.55,2_32_2.70},{4_20_7.7,4_24_8.2},{6_18_15.75,6_16.3},{8_16_26.5,8_16_27.2}}

# FIGURE 5
python scripts/Topo_t.py ${TE} ${WE} ${DDIR}/WF_8_16_26.7

# FIGURE 6-9
python scripts/TC_hist.py ${DDIR}/WF_2_{20_2.55,24_2.60}
python scripts/TC_hist.py ${DDIR}/WF_4_{20_7.7,24_8.2}
python scripts/TC_hist.py ${DDIR}/WF_6_{18_15.75,20_16.3}
python scripts/TC_hist.py ${DDIR}/WF_8_16_{26.5,27.2}

# FIGURE 10
python scripts/Topo_tauQ.py tauQ_vs_t0_0.5_0.5.dat

# TABLE 4, 5, 6, FIG 11
python scripts/Topology.py ${TE} ${WE} ${DDIR}/WF_2_{20_2.55,24_2.60,32_2.65,32_2.70} \
	${DDIR}/WF_4_{20_7.7,20_7.72,20_7.76,20_7.78,20_7.80,20_7.85,24_8.2} \
	${DDIR}/WF_6_{18_15.75,16_15.9,16_16.1,20_16.3} \
	${DDIR}/WF_8_16_{26.5,26.7,27.0,27.2} > chi_vs_t0_${TE}_w0_${WE}_scaled.dat
python scripts/Topology_contlim.py chi_vs_t0_${TE}_w0_${WE}_scaled.dat
