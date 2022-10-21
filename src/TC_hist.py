import os
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.ticker as ticker
from matplotlib.ticker import AutoLocator
from collections import Counter
from uncertainties import ufloat, unumpy
import numpy as np
import re
import sys
import lib_topology as es

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Computer Modern Roman"]
plt.rcParams["text.usetex"] = True
plt.rcParams["lines.linewidth"] = 0.3
plt.rcParams["figure.figsize"] = [3.54, 5.08]
plt.rcParams["errorbar.capsize"] = 2
plt.rcParams["lines.markersize"] = 2
plt.matplotlib.rc("font", size=9)

outdir = os.environ.get("PLOT_DIR", ".")

patt = re.compile(r"WF_([0-9])_([0-9]+)_([0-9]+.[0-9]+)")


def f(x, A, m, s):
    return A * np.exp(-((x - m) ** 2) / (2.0 * s**2))


@ticker.FuncFormatter
def major_formatter(x, pos):
    return "%d" % x


ngraf = len(sys.argv[1:])
i = 0
fig = plt.figure()
gs = gridspec.GridSpec(ngraf, 2, width_ratios=[2.3, 0.6])
gs.update(wspace=0.05, hspace=0.1)

for fname in sys.argv[1:]:
    rawdata = np.genfromtxt(
        fname,
        usecols=(0, 1, 2, 3, 4, 5, 6),
        names=["nconf", "t", "E", "t2E", "tsym", "t2symE", "TC"],
    )
    (
        N,
        L,
        beta,
    ) = patt.search(fname).groups()
    print(N, L, beta)

    tmax = np.max(rawdata["t"])
    TC = np.compress(rawdata["t"] == tmax, rawdata)
    TC = np.compress(TC["nconf"] < 4000, TC)
    maxTC = np.max(np.abs(TC["TC"])) + 5

    tc, stc = es.bs_avg_binned(TC["TC"])

    print(tmax)

    axs0 = plt.subplot(gs[i, 0])
    axs0.yaxis.set_major_formatter(major_formatter)
    axs0.set_ylim(-maxTC, maxTC)
    axs0.set_ylabel(r"$Q_L(t_0)$")
    lab = r"$N_c=" + str(int(N)) + "$, $\\beta=" + str(beta) + "$, $L=" + str(L) + "a$"
    axs0.step(TC["nconf"], TC["TC"], label=lab)
    axs0.legend(frameon=False, fontsize=9, loc="upper right")
    axs0.set_xticks([])
    axs1 = plt.subplot(gs[i, 1])
    axs1.set_ylim(-maxTC, maxTC)
    axs1.set_xticks([])
    axs1.set_yticklabels([])
    # axs1.hist(TC["TC"], histtype="step", orientation="horizontal")

    Q_bins = Counter(np.rint(TC["TC"]))
    range_min = min(min(Q_bins), -max(Q_bins)) - 1
    # Add one to be inclusive
    range_max = -range_min + 1
    Q_range = np.arange(range_min, range_max)

    Q_counts = [Q_bins[Q] for Q in Q_range]

    axs1.step(Q_counts, Q_range - 0.5)

    xdata = Q_range
    ydata = Q_counts
    ydata_err = np.sqrt(Q_counts)
    popt, pcov = curve_fit(f, xdata, ydata)
    smooth_Q_range = np.linspace(range_min - 0.5, range_max - 0.5, 1000)
    axs1.plot(f(smooth_Q_range, *popt), smooth_Q_range)
    textstr = r"$\langle Q_L \rangle = {:0.2uS}$".format(ufloat(tc, stc))
    axs0.text(
        0.05,
        0.10,
        textstr,
        transform=axs0.transAxes,
        fontsize=10,
        verticalalignment="top",
    )

    i = i + 1

axs0.xaxis.set_major_locator(AutoLocator())
axs0.set_xlabel(r"Trajectory")
axs1.set_xlabel(r"Count")

gs.tight_layout(fig)
outname = outdir + "/TC_hist_" + N + ".pdf"
plt.savefig(outname)
