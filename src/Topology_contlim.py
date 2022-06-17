import argparse
import csv
import numpy as np
import matplotlib.pyplot as plt
import re
import matplotlib.ticker as ticker
from matplotlib import gridspec
from matplotlib.ticker import AutoLocator
import itertools
import lib_topology as es
from scipy.optimize import curve_fit
from uncertainties import ufloat, umath

from tables_base import table_start, table_end

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Computer Modern Roman"]
plt.rcParams["text.usetex"] = True
plt.rcParams["lines.linewidth"] = 1
plt.rcParams["figure.figsize"] = [3.74, 3.54]
plt.rcParams["errorbar.capsize"] = 2
plt.rcParams["lines.markersize"] = 2
plt.matplotlib.rc("font", size=8)

parser = argparse.ArgumentParser()
parser.add_argument("summary_data_filename")
parser.add_argument("--t0_plot_filename", default=None)
parser.add_argument("--w0_plot_filename", default=None)
parser.add_argument("--beta_table_file", type=argparse.FileType("w"), default="-")
parser.add_argument("--beta_table_data_file", type=argparse.FileType("w"), default=None)
parser.add_argument("--cont_t0_table_file", type=argparse.FileType("w"), default="-")
parser.add_argument("--cont_w0_table_file", type=argparse.FileType("w"), default="-")
parser.add_argument("--clim_data_file", type=argparse.FileType("w"), default="-")
parser.add_argument("--clim_csv_file", type=argparse.FileType("w"), default=None)
args = parser.parse_args()

clim_summary_data = {}
if args.clim_csv_file:
    csv_writer = csv.writer(args.clim_csv_file, lineterminator="\n")

marks = itertools.cycle(("v", "s", "^"))


def f(x, a, b):
    return a + b * x


def f2(x, a, b):
    return a + b * x * x


def g(x, a, b, c):
    return a + b * x + c * x**2


patt = re.compile(r"chi_vs_t0_([0-9]+.[0-9]+)_w0_([0-9]+.[0-9]+)_([a-z]+).dat")
(
    t0,
    w0,
    labf,
) = patt.search(args.summary_data_filename).groups()

chi_SPN = np.genfromtxt(
    args.summary_data_filename,
    usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15),
    names=[
        "N",
        "L",
        "nconf",
        "beta",
        "sqrtS",
        "sqrtS_err",
        "TE",
        "tscale",
        "tscale_err",
        "chi_t",
        "schi_t",
        "WE",
        "wscale",
        "wscale_err",
        "chi_w",
        "schi_w",
    ],
)

plt.figure()
plt.xlabel(r"$a^2/t_0$")
plt.ylabel(r"$\chi_L t_0^2$")

dtclim = np.dtype(
    [
        ("N", "i"),
        ("chi", "d"),
        ("chi_err", "d"),
        ("sl", "d"),
        ("sl_err", "d"),
        ("chi2", "d"),
    ]
)
clim = np.empty(0, dtype=dtclim)

t_exponent = -int(
    np.floor(
        np.log10(np.mean(chi_SPN["chi_t"] * chi_SPN["tscale"] ** 2 / chi_SPN["L"] ** 4))
    )
)
w_exponent = -int(
    np.floor(
        np.log10(np.mean(chi_SPN["chi_w"] * chi_SPN["wscale"] ** 4 / chi_SPN["L"] ** 4))
    )
)

print(table_start.format(columnspec="|cccccc|"), file=args.beta_table_file)
print(
    (
        r"$N_c$ & $\beta$ & $\sigma t_0$ & $\chi_L t_0^2 \cdot 10^{t_exponent}$ &"
        + r"$\sigma w_0^2$ & $\chi_L w_0^4 \cdot 10^{w_exponent}$ \\"
    ).format(t_exponent=t_exponent, w_exponent=w_exponent),
    file=args.beta_table_file,
)

for i in np.unique(chi_SPN["N"]):
    pl_data = np.compress(chi_SPN["N"] == i, chi_SPN)
    xdata_t0 = 1.0 / pl_data["tscale"]
    ydata_t0 = pl_data["chi_t"] * pl_data["tscale"] ** 2 / pl_data["L"] ** 4

    ydata_t0_err = (
        np.sqrt(
            (2.0 * pl_data["chi_t"] * pl_data["tscale"] * pl_data["tscale_err"]) ** 2
            + (pl_data["schi_t"] * pl_data["tscale"] ** 2) ** 2
        )
        / pl_data["L"] ** 4
    )
    xdata_w0 = 1.0 / pl_data["wscale"] ** 2
    ydata_w0 = pl_data["chi_w"] * pl_data["wscale"] ** 4 / pl_data["L"] ** 4
    ydata_w0_err = (
        np.sqrt(
            (4.0 * pl_data["chi_w"] * pl_data["wscale"] ** 3 * pl_data["wscale_err"])
            ** 2
            + (pl_data["schi_w"] * pl_data["wscale"] ** 4) ** 2
        )
        / pl_data["L"] ** 4
    )

    print(r"\hline", file=args.beta_table_file)
    for l in range(len(xdata_t0)):
        val_t0 = 10**t_exponent * ufloat(ydata_t0[l], ydata_t0_err[l])
        uval_sqrtS_t0 = ufloat(
            pl_data["sqrtS"][l] ** 2.0 * pl_data["tscale"][l],
            np.sqrt(
                (
                    2.0
                    * pl_data["sqrtS"][l]
                    * pl_data["sqrtS_err"][l]
                    * pl_data["tscale"][l]
                )
                ** 2
                + (pl_data["sqrtS"][l] ** 2 * pl_data["tscale_err"][l]) ** 2
            ),
        )
        val_w0 = 10**w_exponent * ufloat(ydata_w0[l], ydata_w0_err[l])
        uval_sqrtS_w0 = ufloat(
            pl_data["sqrtS"][l] ** 2.0 * pl_data["wscale"][l] ** 2,
            np.sqrt(
                (
                    2.0
                    * pl_data["sqrtS"][l]
                    * pl_data["sqrtS_err"][l]
                    * pl_data["wscale"][l] ** 2
                )
                ** 2
                + (
                    pl_data["sqrtS"][l] ** 2
                    * 2.0
                    * pl_data["wscale"][l]
                    * pl_data["wscale_err"][l]
                )
                ** 2
            ),
        )
        print(
            "$",
            int(pl_data["N"][l]),
            "$ & $",
            pl_data["beta"][l],
            "$ & $",
            "{:.2uS}".format(uval_sqrtS_t0),
            "$ & $",
            "{:.2uS}".format(val_t0),
            "$ & $",
            "{:.2uS}".format(uval_sqrtS_w0),
            "$ & $",
            "{:.2uS}".format(val_w0),
            "$ \\\\",
            file=args.beta_table_file,
        )
        if args.beta_table_data_file:
            print(
                int(pl_data["N"][l]),
                pl_data["beta"][l],
                int(pl_data["nconf"][l]),
                uval_sqrtS_t0.n,
                uval_sqrtS_t0.s,
                val_t0.n,
                val_t0.s,
                uval_sqrtS_w0.n,
                uval_sqrtS_w0.s,
                val_w0.n,
                val_w0.s,
                file=args.beta_table_data_file,
            )
print(table_end, file=args.beta_table_file)
args.beta_table_file.close()

# Continuum extrapolations using the t scale
for i in np.unique(chi_SPN["N"]):
    pl_data = np.compress(chi_SPN["N"] == i, chi_SPN)
    marker1 = next(marks)

    xdata = 1.0 / pl_data["tscale"]
    ydata = pl_data["chi_t"] * pl_data["tscale"] ** 2 / pl_data["L"] ** 4

    ydata_err = (
        np.sqrt(
            (2.0 * pl_data["chi_t"] * pl_data["tscale"] * pl_data["tscale_err"]) ** 2
            + (pl_data["schi_t"] * pl_data["tscale"] ** 2) ** 2
        )
        / pl_data["L"] ** 4
    )

    popt, pcov = curve_fit(f, xdata, ydata, sigma=ydata_err, absolute_sigma=True)
    chi2 = np.sum((ydata - f(xdata, *popt)) ** 2 / ydata_err**2) / (
        len(ydata) - len(popt)
    )
    perr = np.sqrt(np.diag(pcov))
    lab = (
        r"$N_c= "
        + str(int(i))
        + "$, $\\tilde{\mathcal{X}}^2="
        + str(np.around(chi2, 2))
        + "$"
    )

    color1 = next(plt.gca()._get_lines.prop_cycler)["color"]
    xr = np.arange(0.0, np.max(xdata) * 1.15, 0.01)
    clim_tmp = np.array([(i, popt[0], perr[0], popt[1], perr[1], chi2)], dtype=dtclim)
    clim = np.append(clim, clim_tmp)
    plt.plot(xr, f(xr, *popt), color=color1, linestyle="dashed")
    plt.errorbar(
        0.0, popt[0], yerr=perr[0], linestyle="None", marker=marker1, color=color1
    )
    plt.errorbar(
        xdata, ydata, yerr=ydata_err, linestyle="None", marker=marker1, color=color1
    )
    plt.errorbar(
        [np.nan], [np.nan], yerr=[np.nan], color=color1, label=lab, linestyle="dashed"
    )

plt.legend(bbox_to_anchor=(0.0, 1.0), loc="lower left", ncol=2, frameon=False)
plt.tight_layout()
if args.t0_plot_filename:
    plt.savefig(args.t0_plot_filename)
    plt.close(plt.gcf())
else:
    plt.show()

for i in clim["N"]:
    table_fits = np.compress(clim["N"] == i, clim)
    val = ufloat(table_fits["chi"], table_fits["chi_err"])
    c2 = np.around(table_fits["chi2"], 2)[0]
    print("& ${:.2uS}$ & ${:.2f}$".format(val, c2), file=args.cont_t0_table_file)
    clim_summary_data[i] = [val.n, val.s, table_fits["chi2"][0]]


plt.figure()
plt.xlabel(r"$a^2/w_0^2$")
plt.ylabel(r"$\chi_L w_0^4$")


clim = np.empty(0, dtype=dtclim)
for i in np.unique(chi_SPN["N"]):
    pl_data = np.compress(chi_SPN["N"] == i, chi_SPN)
    marker1 = next(marks)

    xdata = 1.0 / pl_data["wscale"] ** 2
    ydata = pl_data["chi_w"] * pl_data["wscale"] ** 4 / pl_data["L"] ** 4

    ydata_err = (
        np.sqrt(
            (4.0 * pl_data["chi_w"] * pl_data["wscale"] ** 3 * pl_data["wscale_err"])
            ** 2
            + (pl_data["schi_w"] * pl_data["wscale"] ** 4) ** 2
        )
        / pl_data["L"] ** 4
    )

    popt, pcov = curve_fit(f, xdata, ydata, sigma=ydata_err, absolute_sigma=True)
    chi2 = np.sum((ydata - f(xdata, *popt)) ** 2 / ydata_err**2) / (
        len(ydata) - len(popt)
    )
    perr = np.sqrt(np.diag(pcov))
    lab = (
        r"$N_c= "
        + str(int(i))
        + "$, $\\tilde{\mathcal{X}}^2="
        + str(np.around(chi2, 2))
        + "$"
    )

    color1 = next(plt.gca()._get_lines.prop_cycler)["color"]
    xr = np.arange(0.0, np.max(xdata) * 1.15, 0.01)
    clim_tmp = np.array([(i, popt[0], perr[0], popt[1], perr[1], chi2)], dtype=dtclim)
    clim = np.append(clim, clim_tmp)
    plt.plot(xr, f(xr, *popt), color=color1, linestyle="dashed")
    plt.errorbar(
        0.0, popt[0], yerr=perr[0], linestyle="None", marker=marker1, color=color1
    )
    plt.errorbar(
        xdata, ydata, yerr=ydata_err, linestyle="None", marker=marker1, color=color1
    )
    plt.errorbar(
        [np.nan], [np.nan], yerr=[np.nan], color=color1, label=lab, linestyle="dashed"
    )

for i in clim["N"]:
    table_fits = np.compress(clim["N"] == i, clim)
    print("N = ", i)
    val = ufloat(table_fits["chi"], table_fits["chi_err"])
    print("     chi*scale (a=0 ) = ", "{:.2uS}".format(val))
    val = ufloat(table_fits["sl"], table_fits["sl_err"])
    print("     coefficient = ", "{:.2uS}".format(val))
    print("     chi^2 = ", np.around(table_fits["chi2"], 2))

for i in clim["N"]:
    table_fits = np.compress(clim["N"] == i, clim)
    val = ufloat(table_fits["chi"], table_fits["chi_err"])
    c2 = np.around(table_fits["chi2"], 2)[0]
    print("& ${:.2uS}$ & ${:.2f}$".format(val, c2), file=args.cont_w0_table_file)
    clim_summary_data[i].extend([val.n, val.s, table_fits["chi2"][0]])

plt.legend(bbox_to_anchor=(0.0, 1.0), loc="lower left", ncol=2, frameon=False)
plt.tight_layout()
if args.w0_plot_filename:
    plt.savefig(args.w0_plot_filename)
    plt.close(plt.gcf())
else:
    plt.show()

plt.figure()
plt.xlabel(r"$\sigma a^2$")
plt.ylabel(r"$\chi_L /\sigma^2$")

clim = np.empty(0, dtype=dtclim)
for i in np.unique(chi_SPN["N"]):
    pl_data = np.compress(chi_SPN["N"] == i, chi_SPN)
    marker1 = next(marks)

    xdata = pl_data["sqrtS"] ** 2
    ydata = pl_data["chi_t"] / pl_data["sqrtS"] ** 4 / pl_data["L"] ** 4

    ydata_err = (
        np.sqrt(
            (4.0 * pl_data["chi_t"] / pl_data["sqrtS"] ** 5 * pl_data["sqrtS_err"]) ** 2
            + (pl_data["schi_t"] / pl_data["sqrtS"] ** 4) ** 2
        )
        / pl_data["L"] ** 4
    )

    popt, pcov = curve_fit(f, xdata, ydata, sigma=ydata_err, absolute_sigma=True)
    chi2 = np.sum((ydata - f(xdata, *popt)) ** 2 / ydata_err**2) / (
        len(ydata) - len(popt)
    )
    perr = np.sqrt(np.diag(pcov))
    lab = (
        r"$N_c= "
        + str(int(i))
        + "$, $\\tilde{\mathcal{X}}^2="
        + str(np.around(chi2, 2))
        + "$"
    )

    color1 = next(plt.gca()._get_lines.prop_cycler)["color"]
    xr = np.arange(0.0, np.max(xdata) * 1.15, 0.01)
    clim_tmp = np.array([(i, popt[0], perr[0], popt[1], perr[1], chi2)], dtype=dtclim)
    clim = np.append(clim, clim_tmp)
    plt.plot(xr, f(xr, *popt), color=color1, linestyle="dashed")
    plt.errorbar(
        0.0, popt[0], yerr=perr[0], linestyle="None", marker=marker1, color=color1
    )
    plt.errorbar(
        xdata, ydata, yerr=ydata_err, linestyle="None", marker=marker1, color=color1
    )
    plt.errorbar(
        [np.nan], [np.nan], yerr=[np.nan], color=color1, label=lab, linestyle="dashed"
    )

plt.legend(bbox_to_anchor=(0.0, 1.0), loc="lower left", ncol=2, frameon=False)

for i in clim["N"]:
    table_fits = np.compress(clim["N"] == i, clim)
    print("N = ", i)
    val = ufloat(table_fits["chi"], table_fits["chi_err"])
    print("     chi / sig^2 (a=0 ) = ", "{:.2uS}".format(val))
    val = ufloat(table_fits["sl"], table_fits["sl_err"])
    print("     coefficient = ", "{:.2uS}".format(val))
    print("     chi^2 = ", np.around(table_fits["chi2"], 2))

for i in clim["N"]:
    table_fits = np.compress(clim["N"] == i, clim)
    val = ufloat(table_fits["chi"], table_fits["chi_err"])
    c2 = np.around(table_fits["chi2"], 2)[0]
    print("${}$ & ${:.2uS}$ & ${:.2f}$ \\\\".format(i, val, c2))
    print(
        i // 2, table_fits["chi"][0], table_fits["chi_err"][0], file=args.clim_data_file
    )

if args.clim_csv_file:
    csv_writer.writerows([v for _, v in sorted(clim_summary_data.items())])
