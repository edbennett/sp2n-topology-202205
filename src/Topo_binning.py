import csv
import argparse
import numpy as np
import matplotlib.pyplot as plt
import re
import matplotlib.ticker as ticker
from matplotlib import gridspec
from matplotlib.ticker import AutoLocator
from uncertainties import ufloat
from scipy.interpolate import interp1d, UnivariateSpline, CubicSpline
import itertools
import lib_topology as es
import pickle as pkl

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Computer Modern Roman"]
plt.rcParams["text.usetex"] = True
plt.rcParams["lines.linewidth"] = 1
plt.rcParams["figure.figsize"] = [1.8 * 3.375, 1.8 * 3.375]
plt.rcParams["errorbar.capsize"] = 2
plt.rcParams["lines.markersize"] = 2
plt.matplotlib.rc("font", size=11)


def binned_std(din, bin_size, rng=np.random.default_rng(), num_bs=es.DEFAULT_NUM_BS):
    bin1 = []
    for i in range(int(len(din) / bin_size)):
        bin1.append(np.average(din[i * bin_size : (i + 1) * bin_size]))

    resampled1 = []
    for j in range(num_bs):
        sam = rng.integers(0, len(bin1), size=len(bin1))
        avg1 = np.average([bin1[i] for i in sam])
        resampled1.append(avg1)

    a = np.std(resampled1)

    return a


parser = argparse.ArgumentParser()
parser.add_argument("--data_dir", default=".")
parser.add_argument("--proc_dir", default=".")
parser.add_argument("--pickle_dir", default="pkl_flows_bs")
parser.add_argument("TE", type=float)
parser.add_argument("WE", type=float)
parser.add_argument("fnames", nargs="+")
parser.add_argument("--num_bs", type=int, default=es.DEFAULT_NUM_BS)
args = parser.parse_args()


def get_n_sw(N, L, beta):
    with open(args.data_dir + "/DATA_FILES", "r") as f:
        reader = csv.reader(f)
        for row in reader:
            if row[0] == N and row[1] == L and row[2] == beta:
                return int(row[3])
    raise ValueError(f"{N}, {L}, {beta} not found in DATA_FILES")


bin_range = np.arange(1, 100, 1)

outf = args.proc_dir + "/tauQ_vs_t0_" + str(args.TE) + "_w0_" + str(args.WE) + ".dat"
f = open(outf, "a")
for fname in args.fnames:
    N, L, beta, TCdata = es.topo_load_raw_data(fname)
    TE_scaled = args.TE
    WE_scaled = args.WE

    fn_bs = args.pickle_dir + "/pkl_bs_" + N + "_" + L + "_" + beta + "_"
    infile = open(fn_bs + "t_symE", "rb")
    bs_flow_symE = pkl.load(infile)
    infile.close()

    # Seed RNG compatibly with other calls
    rng = es.get_rng(f"{TE_scaled}_{fname}_sym")

    t0_tmp_symE = es.find_t0(bs_flow_symE, TE_scaled, rng=rng, num_bs=args.num_bs)

    TC = es.find_TC(TCdata, t0_tmp_symE[0])
    autC, autC_err = es.autocorr(TC)
    tmp = ufloat(autC, autC_err)
    n_sw = get_n_sw(N, L, beta)
    print(N, L, beta, n_sw, t0_tmp_symE[0], t0_tmp_symE[1], autC, autC_err, file=f)
    print(N, L, beta, n_sw, "{:.2uS}".format(tmp))

    plt.figure()
    plt.ylim(0.0, 5.0)
    plt.xlabel(r"b")
    plt.ylabel(r"$\frac{\sigma_b}{\sigma}$")
    name_postfix = "_"

    rng = es.get_rng(f"{fname}_{TE}_bin")
    std_1 = binned_std(TC, 1, num_bs=args.num_bs)
    binned_res = [(binned_std(TC, i, num_bs=args.num_bs) / std_1) ** 2 for i in bin_range]
    lab = r"$\beta=" + str(beta) + "$"
    plt.plot(bin_range, binned_res, label=lab)
    name_postfix = name_postfix + beta
    plt.legend(frameon=False)
    fname_out = "binning_aut_" + str(N) + name_postfix + ".pdf"
    plt.savefig(args.proc_dir + "/" + fname_out)
    plt.close(plt.gcf())
