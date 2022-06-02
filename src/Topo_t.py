import os
import sys
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import MaxNLocator

import lib_topology as es
import pickle as pkl

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Computer Modern Roman"]
plt.rcParams["text.usetex"] = True
plt.rcParams["lines.linewidth"] = 0.5
plt.rcParams["figure.figsize"] = [3.4, 3.4]
plt.rcParams["errorbar.capsize"] = 2
plt.rcParams["lines.markersize"] = 2
plt.matplotlib.rc("font", size=11)

outdir = os.environ.get("PLOT_DIR", ".")

TE = float(sys.argv[1])
WE = float(sys.argv[2])

for fname in sys.argv[3:]:
    N, L, beta, TCdata = es.topo_load_raw_data(fname)

    Nconf = int(np.max(TCdata["nconf"]))
    TE_scaled = TE * es.Casimir_SP(N)
    WE_scaled = WE * es.Casimir_SP(N)

    fn_bs = "pkl_flows_bs/pkl_bs_" + N + "_" + L + "_" + beta + "_"
    infile = open(fn_bs + "t_E", "rb")
    bs_flow_E = pkl.load(infile)
    infile.close()
    infile = open(fn_bs + "w_E", "rb")
    w0_flow_E = pkl.load(infile)
    infile.close()
    infile = open(fn_bs + "t_symE", "rb")
    bs_flow_symE = pkl.load(infile)
    infile.close()
    infile = open(fn_bs + "w_symE", "rb")
    w0_flow_symE = pkl.load(infile)
    infile.close()

    t0_tmp_symE = es.find_t0(bs_flow_symE, TE_scaled)
    TC = es.find_TC(TCdata, t0_tmp_symE[0])
    a_min = es.topo_find_alpha(TC)
    TCdata_t = np.zeros(len(TC), dtype=[("nconf", "i8"), ("TC", "f8")])
    TCdata_t["nconf"] = range(1, len(TC) + 1)
    TCdata_t["TC"] = np.rint(a_min * TC)

    w0_tmp_symE = es.find_w0(w0_flow_symE, WE_scaled)
    TC_w = es.find_TC(TCdata, w0_tmp_symE[0] ** 2)
    a_min = es.topo_find_alpha(TC_w)
    TCdata_w = np.zeros(len(TC_w), dtype=[("nconf", "i8"), ("TC", "f8")])
    TCdata_w["nconf"] = range(1, len(TC_w) + 1)
    TCdata_w["TC"] = np.rint(a_min * TC_w)

    ax = plt.figure().gca()
    plt.xlabel(r"$t/a^2$")
    plt.ylabel(r"$Q_L$")
    plt.ylim(-10, 10)
    color0 = next(plt.gca()._get_lines.prop_cycler)["color"]
    lab = r"$N_c=" + str(int(N)) + "$, $\\beta=" + str(beta) + "$"
    print(lab)
    plt.plot([np.nan], [np.nan], label=lab, color=color0)
    for i in np.arange(1, 101):
        dplot = np.compress(TCdata["nconf"] == i, TCdata)
        plt.plot(dplot["t"], dplot["TC"], alpha=0.3, color=color0)
    xmax_red = np.max(np.unique(dplot["t"]))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    TC_red = np.unique(TCdata_t["TC"][0:100])
    x_TC = np.repeat(t0_tmp_symE[0], len(TC_red))
    plt.plot(x_TC, TC_red, "r+")
    for y in TC_red:
        plt.hlines(
            y,
            xmin=t0_tmp_symE[0],
            xmax=xmax_red,
            alpha=0.4,
            linestyle="dotted",
            color="r",
        )
    plt.tight_layout()
    plt.legend(loc=2, frameon=False)
    fname_out = outdir + "/TCvst_" + str(N) + "_" + str(L) + "_" + str(beta) + ".pdf"
    plt.savefig(fname_out)
    plt.close(plt.gcf())
