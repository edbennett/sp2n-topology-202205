from argparse import ArgumentParser, FileType
import csv
import numpy as np
import matplotlib.pyplot as plt
from uncertainties import ufloat, unumpy
import lib_topology as es
import pickle as pkl

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Computer Modern Roman"]
plt.rcParams["text.usetex"] = True
plt.rcParams["lines.linewidth"] = 1
plt.rcParams["figure.figsize"] = [3.4, 3.4]
plt.rcParams["errorbar.capsize"] = 1
plt.rcParams["lines.markersize"] = 2
plt.matplotlib.rc("font", size=11)

parser = ArgumentParser()
parser.add_argument("files_to_process", nargs="+")
parser.add_argument("--t0_plot_filename", default="t0_scale.pdf")
parser.add_argument("--w0_plot_filename", default="w0_scale.pdf")
parser.add_argument("--csv_file", type=FileType("w"), default=None)
parser.add_argument("--pickle_dir", default="pkl_flows_bs")
parser.add_argument("--num_bs", type=int, default=es.DEFAULT_NUM_BS)
args = parser.parse_args()

if args.csv_file:
    csv_writer = csv.writer(args.csv_file)

fdtype = np.dtype(
    [
        ("N", "i"),
        ("L", "i"),
        ("beta", "d"),
        ("obs", "U"),
        ("E0", "d"),
        ("t", "d"),
        ("st", "d"),
        ("W0", "d"),
        ("w2", "d"),
        ("sw2", "d"),
    ]
)
final_E = np.empty(0, dtype=fdtype)
final_symE = np.empty(0, dtype=fdtype)


e0_list = [0.3, 0.4, 0.5, 0.6]
for TE in e0_list:
    for faddr in args.files_to_process:
        N, L, beta, rawdata = es.topo_load_raw_data(faddr)
        TE_scaled = TE
        WE_scaled = TE

        fn_bs = args.pickle_dir + "/pkl_bs_" + N + "_" + L + "_" + beta + "_"
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

        # Seed an RNG seeded by our current parameters so that
        # the results reproduce.
        tw_plaq_rng = es.get_rng(f"{TE}_{faddr}_plaq")

        try:
            t0_tmp_E = es.find_t0(
                bs_flow_E, TE_scaled, rng=tw_plaq_rng, num_bs=args.num_bs
            )
            w0_tmp_E = es.find_w0(
                w0_flow_E, WE_scaled, rng=tw_plaq_rng, num_bs=args.num_bs
            )
            t0_tmp_E_u = ufloat(t0_tmp_E[0], t0_tmp_E[1])
            w0_tmp_E_u = ufloat(w0_tmp_E[0], w0_tmp_E[1])
            tmp_fin_E = np.array(
                [
                    (
                        N,
                        L,
                        beta,
                        "t2E",
                        TE_scaled,
                        t0_tmp_E[0],
                        t0_tmp_E[1],
                        WE_scaled,
                        w0_tmp_E[0],
                        w0_tmp_E[1],
                    )
                ],
                dtype=fdtype,
            )
            final_E = np.append(final_E, tmp_fin_E)
        except:
            t0_tmp_E = np.nan
            w0_tmp_E = np.nan

        # Seed a fresh RNG so that we get results compatible with
        # vis_WF_N_scaling.py
        tw_sym_rng = es.get_rng(f"{TE}_{faddr}_sym")

        try:
            t0_tmp_symE = es.find_t0(
                bs_flow_symE, TE_scaled, rng=tw_sym_rng, num_bs=args.num_bs
            )
            w0_tmp_symE = es.find_w0(
                w0_flow_symE, WE_scaled, rng=tw_sym_rng, num_bs=args.num_bs
            )
            t0_tmp_symE_u = ufloat(t0_tmp_symE[0], t0_tmp_symE[1])
            w0_tmp_symE_u = ufloat(w0_tmp_symE[0], w0_tmp_symE[1])

            tmp_fin_symE = np.array(
                [
                    (
                        N,
                        L,
                        beta,
                        "t2symE",
                        TE_scaled,
                        t0_tmp_symE[0],
                        t0_tmp_symE[1],
                        WE_scaled,
                        w0_tmp_symE[0],
                        w0_tmp_symE[1],
                    )
                ],
                dtype=fdtype,
            )
            final_symE = np.append(final_symE, tmp_fin_symE)
        except:
            t0_tmp_symE = np.nan
            w0_tmp_symE = np.nan

        print(N, L, beta, TE_scaled, end=" ")
        to_print = []
        csv_row = [N, L, beta, TE_scaled]
        try:
            for v in [t0_tmp_E, t0_tmp_symE]:
                pr = ufloat(np.sqrt(v[0]), v[1] / (2.0 * np.sqrt(v[0])))
                to_print.append("{:.2uS}".format(pr))
                csv_row.extend([pr.n, pr.s])
        except:
            to_print = ["-", "-"]
            csv_row.extend([None, None, None, None])
        finally:
            print(" ".join(to_print), end=" ")

        print(WE_scaled, end=" ")
        to_print = []
        try:
            for v in [w0_tmp_E, w0_tmp_symE]:
                pr = ufloat(np.sqrt(v[0]), v[1] / (2.0 * np.sqrt(v[0])))
                to_print.append("{:.2uS}".format(pr))
                csv_row.extend([pr.n, pr.s])
        except:
            to_print = ["-", "-"]
            csv_row.extend([None, None, None, None])
        finally:
            print(" ".join(to_print), end=" ")
        print("")
        if args.csv_file:
            csv_writer.writerow(csv_row)


plt.figure()
plt.xlabel(r"$\beta$")
plt.ylabel(r"$\sqrt{t_0}/a$")
color0 = next(plt.gca()._get_lines.prop_cycler)["color"]
for e0 in e0_list:
    arr = np.compress(final_E["E0"] == e0, final_E)
    xdata = arr["beta"]
    ydata = np.sqrt(arr["t"])
    ydata_err = arr["st"] / (2.0 * np.sqrt(arr["t"]))
    lab = r"$\mathcal{E}_0=" + str(e0) + "$"

    color1 = next(plt.gca()._get_lines.prop_cycler)["color"]
    plt.errorbar(
        xdata, ydata, yerr=ydata_err, linestyle="None", marker="^", color=color1
    )
    plt.fill_between([np.nan], [np.nan], [np.nan], label=lab, color=color1)
    arr = np.compress(final_symE["E0"] == e0, final_symE)
    xdata = arr["beta"]
    ydata = np.sqrt(arr["t"])
    ydata_err = arr["st"] / (2.0 * np.sqrt(arr["t"]))
    lab = r"$\mathcal{E}_0=" + str(e0) + "$"

    plt.errorbar(
        xdata, ydata, yerr=ydata_err, linestyle="None", marker="o", color=color1
    )
plt.errorbar(
    [np.nan], [np.nan], yerr=[np.nan], color="black", fmt="^", label=r"plaquette"
)
plt.errorbar([np.nan], [np.nan], yerr=[np.nan], color="black", fmt="o", label=r"clover")

plt.tight_layout()
plt.legend(loc="upper left", prop={"size": 6}, frameon=False)
plt.savefig(args.t0_plot_filename)
plt.close(plt.gcf())

plt.figure()
plt.xlabel(r"$\beta$")
plt.ylabel(r"$w_0/a$")
color0 = next(plt.gca()._get_lines.prop_cycler)["color"]
for e0 in e0_list:
    arr = np.compress(final_E["W0"] == e0, final_E)
    xdata = arr["beta"]
    ydata = np.sqrt(arr["w2"])
    ydata_err = arr["sw2"] / (2.0 * arr["w2"])
    lab = r"$\mathcal{W}_0=" + str(e0) + "$"
    color1 = next(plt.gca()._get_lines.prop_cycler)["color"]
    plt.errorbar(
        xdata, ydata, yerr=ydata_err, linestyle="None", marker="^", color=color1
    )
    arr = np.compress(final_symE["W0"] == e0, final_symE)
    xdata = arr["beta"]
    ydata = np.sqrt(arr["w2"])
    ydata_err = arr["sw2"] / (2.0 * arr["w2"])
    lab = r"$\mathcal{W}_0=" + str(e0) + "$"
    plt.errorbar(
        xdata, ydata, yerr=ydata_err, linestyle="None", marker="o", color=color1
    )
    plt.fill_between([np.nan], [np.nan], [np.nan], label=lab, color=color1)
plt.errorbar(
    [np.nan], [np.nan], yerr=[np.nan], color="black", fmt="^", label=r"plaquette"
)
plt.errorbar([np.nan], [np.nan], yerr=[np.nan], color="black", fmt="o", label=r"clover")
plt.tight_layout()
plt.legend(loc="upper left", prop={"size": 6}, frameon=False)
plt.savefig(args.w0_plot_filename)
