# Values of the reference scales used for the full analysis
TE = 0.225
WE = 0.225

# Values of the reference scales used for the demonstationstrative figures
TE_DEMO = 0.5
WE_DEMO = 0.5

# Set this to where the raw output files are
DATA_DIR = data

# These data are included in this code package
QUOTED_DIR = quoted_data

# These will be created automatically
PICKLE_DIR = pkl_flows_bs
TABLES_DIR = tables
PLOT_DIR = plots
PROC_DIR = processed_data
SHORT_PLOT_DIR = short_plots
SHORT_TABLES_DIR = short_tables

# Required to enable some of the more intricate array lookups
.SECONDEXPANSION :

# Keep all intermediary files
.SECONDARY :

${PICKLE_DIR} ${TABLES_DIR} ${PLOT_DIR} ${PROC_DIR} ${SHORT_PLOT_DIR} ${SHORT_TABLES_DIR} :
	mkdir -p $@

${PROC_DIR}/WF_% : ${DATA_DIR}/out_% | ${PROC_DIR}
	grep "\[WILSONFLOW\]" $^ | cut -d"=" -f2 | awk '{if($$1==0.0) i++; print i,$$0;}' > $@

${PROC_DIR}/WF_%_plaq : ${DATA_DIR}/out_% | ${PROC_DIR}
	grep "Plaquette" $^ | cut -d"=" -f2 > $@

# Bootstrap samples
$(foreach SUFFIX, t_E w_E t_symE w_symE, ${PICKLE_DIR}/pkl_bs_%_${SUFFIX}) &: ${PROC_DIR}/WF_% | ${PICKLE_DIR}
	python src/produce_bs_sample.py $^

# Table 1
TAB1_OUTPUT = ${TABLES_DIR}/table_tauQ_t0_${TE}_w0_${WE}.tex
TAB1_DEPS = $(foreach SUFFIX, 2_20_2.55 2_24_2.60 2_32_2.65 2_32_2.70 4_20_7.7 4_20_7.72 4_20_7.76 4_20_7.78 4_20_7.80 4_20_7.85 4_24_8.2 6_18_15.75 6_16_15.9 6_16_16.1 6_20_16.3 8_16_26.5 8_16_26.7 8_16_27.0 8_16_27.2, ${PROC_DIR}/WF_${SUFFIX})
${TABLES_DIR}/table_tauQ_t0_%.tex ${PROC_DIR}/tauQ_vs_t0_%.dat &: ${TAB1_DEPS} | ${TABLES_DIR}
	PROC_DIR=${PROC_DIR} bash src/table1.sh ${DATA_DIR} ${PROC_DIR} $(subst _w0_, ,$*) ${TAB1_OUTPUT}

# Figure 1
NON_SCALED_SUFFIXES_6 = 16_15.9 16_16.1
FIG1_OUTPUT = ${PLOT_DIR}/non_scaled_flows_6.pdf

${PLOT_DIR}/non_scaled_flows_%.pdf : $$(foreach SUFFIX, $${NON_SCALED_SUFFIXES_$$*}, ${PROC_DIR}/WF_$$*_$${SUFFIX}) | $$(foreach SUFFIX, $${NON_SCALED_SUFFIXES_$$*}, $$(foreach FLOW, t_E w_E t_symE w_symE, ${PICKLE_DIR}/pkl_bs_$$*_$${SUFFIX}_$${FLOW})) ${PLOT_DIR}
	PLOT_DIR=${PLOT_DIR} python src/vis_WF_W_E.py ${TE_DEMO} ${WE_DEMO} $^ --absolute_scales --outfile $@

# Figure 2, 12-17; tables 2, 3, 7-10
FIG2_OUTPUT = $(foreach OBS, w t, ${PLOT_DIR}/Scale_6_${OBS}0.pdf)
FIG12_OUTPUT = ${PLOT_DIR}/Scale_2_t0.pdf
FIG13_OUTPUT = ${PLOT_DIR}/Scale_2_w0.pdf
FIG14_OUTPUT = ${PLOT_DIR}/Scale_4_t0.pdf
FIG15_OUTPUT = ${PLOT_DIR}/Scale_4_w0.pdf
FIG16_OUTPUT = ${PLOT_DIR}/Scale_8_t0.pdf
FIG17_OUTPUT = ${PLOT_DIR}/Scale_8_w0.pdf
TAB2_OUTPUT = ${TABLES_DIR}/Scale_6_t0.dat
TAB3_OUTPUT = ${TABLES_DIR}/Scale_6_w0.dat
TAB7_OUTPUT = ${TABLES_DIR}/Scale_2_t0.dat
TAB8_OUTPUT = ${TABLES_DIR}/Scale_2_w0.dat
TAB9_OUTPUT = ${TABLES_DIR}/Scale_4_t0.dat
TAB10_OUTPUT = ${TABLES_DIR}/Scale_4_w0.dat
TAB11_OUTPUT = ${TABLES_DIR}/Scale_8_t0.dat
TAB12_OUTPUT = ${TABLES_DIR}/Scale_8_w0.dat

FIG2_NC2_SUFFIXES = 20_2.55 24_2.60 32_2.65 32_2.70
FIG2_NC4_SUFFIXES = $(foreach BETA, 7.7 7.72 7.76 7.78 7.80 7.85, 20_${BETA})
FIG2_NC6_SUFFIXES = 18_15.75 16_15.9 16_16.1 20_16.3
FIG2_NC8_SUFFIXES = $(foreach BETA, 26.5 26.7 27.0 27.2, 16_${BETA})

# Thanks to Mad Scientist on StackOverflow for explaining this:
# https://stackoverflow.com/questions/72324675/can-gnu-make-use-pattern-matching-to-look-up-variables/72324829
FIG2_RAW_DEPS = $(foreach SUFFIX,${FIG2_NC$1_SUFFIXES},${PROC_DIR}/WF_$1_${SUFFIX})
FIG2_PICKLE_DEPS = $(foreach SUFFIX,${FIG2_NC$1_SUFFIXES},$(foreach VAR,t_E w_E t_symE w_symE,${PICKLE_DIR}/pkl_bs_$1_${SUFFIX}_${VAR}))

${PROC_DIR}/Scale_%.dat ${PLOT_DIR}/Scale_%_t0.pdf ${PLOT_DIR}/Scale_%_w0.pdf &: $$(call FIG2_RAW_DEPS,$$*) $$(call FIG2_PICKLE_DEPS,$$*) | ${PROC_DIR} ${PLOT_DIR}
	python src/vis_WF_Scale.py $(call FIG2_RAW_DEPS,$*) --t0_plot_filename ${PLOT_DIR}/Scale_$*_t0.pdf --w0_plot_filename ${PLOT_DIR}/Scale_$*_w0.pdf > ${PROC_DIR}/Scale_$*.dat

${TABLES_DIR}/Scale_%_t0.dat ${TABLES_DIR}/Scale_%_w0.dat &: ${PROC_DIR}/Scale_%.dat | ${TABLES_DIR}
	bash src/tabulate_scale.sh $< ${TABLES_DIR}/Scale_$*_t0.dat ${TABLES_DIR}/Scale_$*_w0.dat

# Figure 3, 4
FIG3_OUTPUT = ${PLOT_DIR}/flows_${TE_DEMO}_${WE_DEMO}_scaled.pdf
FIG4_OUTPUT = ${PLOT_DIR}/flows_${TE}_${WE}_scaled.pdf
FIG3_4_SUFFIXES = 2_20_2.55 2_32_2.70 4_20_7.7 4_24_8.2 6_18_15.75 6_20_16.3 8_16_26.5 8_16_27.2
FIG3_4_ARGS = $(foreach SUFFIX,${FIG3_4_SUFFIXES},${PROC_DIR}/WF_${SUFFIX})
FIG3_4_REQS = $(foreach SUFFIX,${FIG3_4_SUFFIXES},${PROC_DIR}/WF_${SUFFIX}_plaq ${PICKLE_DIR}/pkl_bs_${SUFFIX}_t_symE ${PICKLE_DIR}/pkl_bs_${SUFFIX}_w_symE)

${PLOT_DIR}/flows_%_scaled.pdf : ${FIG3_4_ARGS} ${FIG3_4_REQS} | ${PLOT_DIR}
	PLOT_DIR=${PLOT_DIR} python src/vis_WF_N_scaling.py $$(echo $* | sed "s/_/ /") ${FIG3_4_ARGS}

# Figure 5
FIG5_OUTPUT = ${PLOT_DIR}/TCvst_8_16_26.7.pdf

${PLOT_DIR}/TCvst_%.pdf : ${PROC_DIR}/WF_% | ${PLOT_DIR}
	PLOT_DIR=${PLOT_DIR} python src/Topo_t.py ${TE} ${WE} $^

# Figure 6, 7, 8, 9
FIG6789_NC2_SUFFIXES = 32_2.65 32_2.70
FIG6789_NC4_SUFFIXES = 20_7.7 24_8.2
FIG6789_NC6_SUFFIXES = 18_15.75 20_16.3
FIG6789_NC8_SUFFIXES = 16_26.5 16_27.2
FIG6789_DEPS = $(foreach SUFFIX, ${FIG6789_NC$1_SUFFIXES}, ${PROC_DIR}/WF_$1_${SUFFIX})

FIG6789_OUTPUT = $(foreach NUM, 2 4 6 8, ${PLOT_DIR}/TC_hist_${NUM}.pdf)

${PLOT_DIR}/TC_hist_%.pdf : $$(call FIG6789_DEPS,$$*) | ${PLOT_DIR}
	PLOT_DIR=${PLOT_DIR} python src/TC_hist.py $^

# Figure 10
FIG10_OUTPUT = ${PLOT_DIR}/tauQ_vs_t0_${TE}.pdf
${PLOT_DIR}/tauQ_vs_t0_%.pdf : ${PROC_DIR}/tauQ_vs_t0_$$*_w0_$$*.dat | ${PLOT_DIR}
	PLOT_DIR=${PLOT_DIR} python src/Topo_tauQ.py $^ $@

# Table 4, 5, 6; Figure 11
FIG11_NC2_SUFFIXES = 20_2.55 24_2.60 32_2.65 32_2.70
FIG11_NC4_SUFFIXES = $(foreach BETA, 7.7 7.72 7.76 7.78 7.80 7.85, 20_${BETA}) 24_8.2
FIG11_NC6_SUFFIXES = 18_15.75 16_15.9 16_16.1 20_16.3
FIG11_NC8_SUFFIXES = $(foreach BETA, 26.5 26.7 27.0 27.2, 16_${BETA})


TAB45_OUTPUT = ${TABLES_DIR}/table_chi_t0_${TE}_w0_${WE}.tex
${TABLES_DIR}/table_chi_t0_%.tex ${PROC_DIR}/chi_vs_t0_%_scaled.dat &: $(foreach NC, 2 4 6 8, $(foreach SUFFIX,${FIG11_NC${NC}_SUFFIXES},${PROC_DIR}/WF_${NC}_${SUFFIX})) | $(foreach SUFFIX,${FIG11_NC${NC}_SUFFIXES}, $(foreach FLOW, t_E w_E t_symE w_symE, ${PICKLE_DIR}/pkl_bs_${NC}_${SUFFIX}_${FLOW})) ${TABLES_DIR}
	TABLES_DIR=${TABLES_DIR} python src/Topology.py $(subst _w0_, ,$*) $^ > ${PROC_DIR}/chi_vs_t0_$*_scaled.dat

FIG11_OUTPUT = $(foreach SUFFIX,.pdf _w0.pdf,$(foreach SCALE,${TE_DEMO} ${TE},${PLOT_DIR}/SPN_Topology_contlim_${SCALE}_${SCALE}_scaled${SUFFIX}))
${PLOT_DIR}/SPN_Topology_contlim_%_scaled.pdf ${PLOT_DIR}/SPN_Topology_contlim_%_scaled_w0.pdf ${PROC_DIR}/clim_SP_%.dat : ${PROC_DIR}/chi_vs_t0_$$(subst _,_w0_,$$*)_scaled.dat | ${PLOT_DIR}
	PLOT_DIR=${PLOT_DIR} python src/Topology_contlim.py $^ ${PROC_DIR}/clim_SP_$*.dat

SHORT_FIGS = ${SHORT_PLOT_DIR}/ScaledChi.pdf ${SHORT_PLOT_DIR}/NONScaledChi.pdf
SHORT_TABLE = ${SHORT_TABLES_DIR}/summary_table.tex
${SHORT_FIGS} ${SHORT_TABLE} &: $(foreach AUTHORS,AA_MT DD_EV_HP CB_MD CB_CB_MD BL_MT,quoted_data/clim_SU_${AUTHORS}.dat) ${PROC_DIR}/clim_SP_${TE}_${WE}.dat | ${SHORT_PLOT_DIR} ${SHORT_TABLES_DIR}
	PLOT_DIR=${SHORT_PLOT_DIR} TABLES_DIR=${SHORT_TABLES_DIR} PROC_DIR=${PROC_DIR} QUOTED_DIR=${QUOTED_DIR} python src/fits_largeN_universality.py

datapackage.h5 : ${TAB1_DEPS} $(foreach PREFIX, ${TAB1_DEPS}, ${PREFIX}_plaq)
	python src/package_data.py $^ --hdf5_filename=$@

long-figures : $(foreach FIGNUM, 1 2 3 4 5 6789 10 11 12 13 14 15 16 17,${FIG${FIGNUM}_OUTPUT}) ${SHORT_FIGS}
long-tables : $(foreach TABNUM, 1 2 3 45 6 7 8 9 10 11 12, ${TAB${TABNUM}_OUTPUT})
short-figures : ${SHORT_FIGS}
short-tables : ${SHORT_TABLE}
all-figures : long-figures short-figures
all-tables : long-tables short-tables

clean : clean-processed clean-plots clean-tables

clean-plots :
	rm -f ${PLOT_DIR}/* ${SHORT_PLOT_DIR}/*

clean-tables :
	rm -f ${TABLES_DIR}/* ${SHORT_TABLES_DIR}/*

clean-processed :
	rm -f ${PROC_DIR}/*

clean-pickles :
	rm -f ${PICKLE_DIR}/*

clean-all : clean clean-pickles

all : all-figures all-tables

datapackage : datapackage.h5

.DEFAULT_GOAL := all
