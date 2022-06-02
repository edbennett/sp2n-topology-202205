import numpy as np
import re
from scipy.interpolate import interp1d


def autocorr(indata):

    L = len(indata)
    if L == 0:
        print("Empty data set")
        autocorr = 0.0
        err = 0.0
        return autocorr, err
    avr = 0.0
    f2 = 0.0
    C0 = 0.0
    for j in range(L):
        avr = avr + indata[j]
        f2 = f2 + indata[j] ** 2
    f = 2.0 * avr
    avr = avr / (1.0 * L)
    C0 = f2 / (1.0 * L) - avr ** 2
    tint = 0.5
    val = False
    for M in range(1, L):
        avr = 0.0
        f = 0.0
        f2 = 0.0
        for i in range(L - M):
            f2 = f2 + indata[i] * indata[i + M]
            f = f + indata[i] + indata[i + M]
            avr = avr + indata[i]
        avr = avr / (1.0 * (L - M))
        Ct = (f2 / (1.0 * (L - M))) + avr * (avr - f / (1.0 * (L - M)))
        rho = Ct / C0
        tint = tint + rho
        if 4.0 * tint < M:
            val = True
            break
    if val:
        autocorr = tint
        err = np.sqrt(2.0 * (2.0 * M + 1) / (1.0 * L)) * tint
    else:
        autocorr = 0.0
        err = 0.0
    return autocorr, err


def topo_find_alpha(din):
    a_min = 0.4
    delta = np.average((a_min * din - np.rint(a_min * din)) ** 2)
    for ia in np.arange(a_min, 2.0, 0.000025):
        delta_tmp = np.average((ia * din - np.rint(ia * din)) ** 2)
        if delta_tmp < delta:
            a_min = ia
            delta = delta_tmp
    return a_min


def topo_load_raw_data(fname):
    patt = re.compile(r"WF_([0-9])_([0-9]+)_([0-9]+.[0-9]+)")
    N, L, beta, = patt.search(fname).groups()
    out_data = np.genfromtxt(
        fname,
        usecols=(0, 1, 2, 3, 4, 5, 6),
        dtype=[
            ("nconf", "int"),
            ("t", "f8"),
            ("E", "f8"),
            ("t2E", "f8"),
            ("Esym", "f8"),
            ("t2symE", "f8"),
            ("TC", "f8"),
        ],
    )
    return N, L, beta, out_data


def bs_avg_err_TC(din, bin_size=20):
    bin1 = []
    bin2 = []
    for i in range(int(len(din) / bin_size)):
        bin1.append(np.average(din[i * bin_size : (i + 1) * bin_size]))
        bin2.append(np.average(din[i * bin_size : (i + 1) * bin_size] ** 2))

    resampled1 = []
    resampled2 = []
    for j in range(100):
        sam = np.random.randint(0, len(bin1), size=len(bin1))
        avg = np.average([bin1[i] for i in sam])
        avg2 = np.average([bin2[i] for i in sam])
        resampled1.append(avg)
        resampled2.append(avg2 - avg ** 2)
    return np.average(resampled2), np.std(resampled2)


def bstrap(a, nbstrap):
    stat1 = np.zeros(nbstrap)
    stat2 = np.zeros(nbstrap)
    lena = len(a)
    for i in range(nbstrap):
        atmp = np.random.choice(a, lena, replace=True)
        stat1[i] = np.average(atmp)
        stat2[i] = np.std(atmp, ddof=1)
    s1 = np.average(stat1)
    e1 = np.std(stat1)
    s2 = np.average(stat2)
    e2 = np.std(stat2)
    return s1, e1, s2, e2


def find_t0(indata, TE):
    t0 = []
    for i in np.unique(indata["bs"]):
        data = np.compress(indata["bs"] == i, indata)
        idx = np.abs(data["flow"] - TE).argmin()
        traj = data[idx - 5 : idx + 5]
        f = interp1d(traj["flow"], traj["t"])
        t0 = np.append(t0, f(TE))
    a, b, c, d = bstrap(t0, 100)
    return a, c


def find_w0(indata, TE):
    t0 = []
    for i in np.unique(indata["bs"]):
        data = np.compress(indata["bs"] == i, indata)
        idx = np.abs(data["flow"] - TE).argmin()
        traj = data[idx - 5 : idx + 5]
        f = interp1d(traj["flow"], traj["t"])
        t0 = np.append(t0, np.sqrt(f(TE)))
    a, b, c, d = bstrap(t0, 100)
    return a, c


rng = np.random.default_rng()


def flows(rawdata, N_bs, obs, Ntherm=500, rng=rng):
    Nconf = np.max(np.unique(rawdata["nconf"]))
    rawdata["t"] = np.around(rawdata["t"], 8)
    ts = np.unique(rawdata["t"])
    delta_t = ts[1] - ts[0]
    bs_flow = np.empty(N_bs * len(ts), dtype=[("bs", "i"), ("t", "f8"), ("flow", "f8")])
    w0_flow = np.empty(
        N_bs * (len(ts) - 1), dtype=[("bs", "i"), ("t", "f8"), ("flow", "f8")]
    )
    tmp_bs_flow = np.empty(
        Nconf * len(ts),
        dtype=[
            ("nconf", "int"),
            ("t", "f8"),
            ("E", "f8"),
            ("t2E", "f8"),
            ("Esym", "f8"),
            ("t2symE", "f8"),
            ("TC", "f8"),
        ],
    )
    per_conf_flows = {}
    for s in set(rawdata["nconf"]):
        current_conf_locs = np.where(rawdata["nconf"] == s)[0]
        conf_lbound = current_conf_locs[0]
        conf_ubound = current_conf_locs[-1] + 1
        per_conf_flows[s] = rawdata[conf_lbound:conf_ubound]
    for i_bs in range(N_bs):
        sam = rng.integers(Ntherm, Nconf + 1, size=Nconf)
        i_c = 0
        for s in sam:
            tmp_bs_flow[i_c * len(ts) : (i_c + 1) * len(ts)] = per_conf_flows[s]
            i_c = i_c + 1
        for i_t in range(len(ts)):
            t0_tmp = np.around(ts[i_t], 8)
            f0tmp = np.compress(tmp_bs_flow["t"] == t0_tmp, tmp_bs_flow)
            f0_avg = np.average(f0tmp[obs])
            bs_fl_tmp = np.array(
                [(i_bs, t0_tmp, f0_avg)],
                dtype=[("bs", "i"), ("t", "f8"), ("flow", "f8")],
            )
            bs_flow[i_bs * len(ts) + i_t] = bs_fl_tmp
            if i_t < len(ts) - 1:
                t1_tmp = np.around(ts[i_t + 1], 8)
                f1tmp = np.compress(tmp_bs_flow["t"] == t1_tmp, tmp_bs_flow)
                f1_avg = np.average(f1tmp[obs])
                delta_t = np.around(t1_tmp - t0_tmp, 8)
                w0_fl_tmp = np.array(
                    [(i_bs, t0_tmp, t0_tmp * (f1_avg - f0_avg) / delta_t)],
                    dtype=[("bs", "i"), ("t", "f8"), ("flow", "f8")],
                )
                w0_flow[i_bs * (len(ts) - 1) + i_t] = w0_fl_tmp
    return bs_flow, w0_flow


def avg_flows(in_flow):
    ts = np.unique(in_flow["t"])
    Lts = len(ts)
    out_flow = np.empty(Lts, dtype=[("t", "f8"), ("flow", "f8"), ("err", "f8")])
    for i_t in range(Lts):
        out_flow["t"][i_t] = ts[i_t]
        toavg = np.compress(in_flow["t"] == ts[i_t], in_flow["flow"])
        out_flow["flow"][i_t] = np.average(toavg)
        out_flow["err"][i_t] = np.std(toavg)
    return out_flow


def Casimir_SP(N):
    return (float(N) + 1.0) / 4.0


def Casimir_SUN(n):
    return (float(n) ** 2 - 1.0) / (2.0 * float(n))


def find_TC(rawdata, t):
    confs = np.unique(rawdata["nconf"])
    Nconf = len(confs)
    TC = np.empty(Nconf, dtype="d")
    for i in confs:
        data = np.compress(rawdata["nconf"] == i, rawdata)
        idx = np.abs(data["t"] - t).argmin()
        traj = data[idx - 3 : idx + 3]
        f = interp1d(data["t"], data["TC"])
        TC[i - 1] = f(t)
    return TC


def d_g_SPN(n):
    return n * (2.0 * n + 1.0)


def d_g_SUN(n):
    return n ** 2 - 1.0
