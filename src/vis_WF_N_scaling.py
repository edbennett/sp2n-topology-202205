import argparse
import numpy as np
import matplotlib.pyplot as plt

from uncertainties import ufloat, unumpy

import lib_topology as es
import pickle as pkl

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Computer Modern Roman"]
plt.rcParams["text.usetex"] = True
plt.rcParams["lines.linewidth"] = 1
plt.rcParams["figure.figsize"] = [3.4, 6.8]
plt.rcParams["errorbar.capsize"] = 2
plt.rcParams["lines.markersize"] = 2
plt.matplotlib.rc("font", size=11)

parser = argparse.ArgumentParser()
parser.add_argument("TE", type=int)
parser.add_argument("WE", type=int)
parser.add_argument("faddrs", nargs="+")
parser.add_argument("--plot_dir", default=".")
parser.add_argument("--pickle_dir", default="pkl_flows_bs")
parser.add_argument("--num_bs", type=int, default=es.DEFAULT_NUM_BS)
args = parser.parse_args()

fig, ax = plt.subplots(2, 1)
for faddr in args.faddrs:

    N, L, beta, rawdata = es.topo_load_raw_data(faddr)
    TE_scaled = args.TE * es.Casimir_SP(N)
    WE_scaled = args.WE * es.Casimir_SP(N)

    plaq_fname = faddr + "_plaq"
    plaq_data = np.genfromtxt(plaq_fname)

    # Seed an RNG so that the bootstrap is reproducible
    plaq_rng = es.get_rng(faddr)
    plaq_bs = es.bstrap(plaq_data, len(plaq_data), rng=plaq_rng)
    plaq_u = ufloat(plaq_bs[0], plaq_bs[1])
    print(N, L, beta, "{:.2uS}".format(plaq_u))
    lthooft = es.d_g_SPN(int(N)) * plaq_u / float(beta)
    print(N, L, beta, "{:.2uS}".format(lthooft))

    print("loading flow files")
    fn_bs = args.pickle_dir + "/pkl_bs_" + N + "_" + L + "_" + beta + "_"
    infile = open(fn_bs + "t_symE", "rb")
    bs_flow_symE = pkl.load(infile)
    infile.close()
    infile = open(fn_bs + "w_symE", "rb")
    w0_flow_symE = pkl.load(infile)
    infile.close()

    # Seed this RNG compatibly with `vis_WF_Scale.py`
    tw_rng = es.get_rng(f"{TE_scaled}_{faddr}_sym")
    t0_tmp_symE = es.find_t0(bs_flow_symE, TE_scaled, rng=tw_rng, num_bs=args.num_bs)
    w0_tmp_symE = es.find_w0(w0_flow_symE, WE_scaled, rng=tw_rng, num_bs=args.num_bs)

    ax[0].set_xlabel(r"$t/t_0$")
    ax[0].set_ylabel(r"$\mathcal{E}(t)/C_2(F)$")
    ax[0].axhline(y=TE, linestyle="dotted", alpha=0.7)
    color1 = next(plt.gca()._get_lines.prop_cycler)["color"]
    lab1 = (
        r"$N_c="
        + str(int(N))
        + "$, $\\tilde{\lambda}="
        + str(np.around(lthooft.n, 3))
        + "$, cl."
    )
    plot_flow = es.avg_flows(bs_flow_symE)
    ax[0].fill_between(
        plot_flow["t"] / t0_tmp_symE[0],
        y1=(plot_flow["flow"] + plot_flow["err"]) / (es.Casimir_SP(N)),
        y2=(plot_flow["flow"] - plot_flow["err"]) / (es.Casimir_SP(N)),
        label=lab1,
        color=color1,
        alpha=0.5,
    )
    ax[0].legend(loc=2, prop={"size": 7}, frameon=False)

    ax[1].set_xlabel(r"$t/w_0^2$")
    ax[1].set_ylabel(r"$\mathcal{W}(t)/C_2(F)$")
    ax[1].axhline(y=WE, linestyle="dotted", alpha=0.7)
    plot_flow = es.avg_flows(w0_flow_symE)
    ax[1].fill_between(
        plot_flow["t"] / w0_tmp_symE[0] ** 2,
        y1=(plot_flow["flow"] + plot_flow["err"]) / (es.Casimir_SP(N)),
        y2=(plot_flow["flow"] - plot_flow["err"]) / (es.Casimir_SP(N)),
        label=lab1,
        color=color1,
        alpha=0.5,
    )
    ax[1].legend(loc=2, prop={"size": 7}, frameon=False)
plt.tight_layout()
plt.savefig(args.plot_dir + "/flows_" + str(TE) + "_" + str(WE) + "_scaled.pdf")
