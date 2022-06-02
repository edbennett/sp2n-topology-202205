import argparse
import numpy as np
from uncertainties import ufloat
import lib_topology as es
import pickle as pkl

parser = argparse.ArgumentParser()
parser.add_argument("TE", type=float)
parser.add_argument("WE", type=float)
parser.add_argument("input_datafiles", nargs="+")
parser.add_argument("--output_data", type=argparse.FileType("w"), default="-")
parser.add_argument("--output_tex", type=argparse.FileType("w"), default="-")
parser.add_argument("--sqrt_sigma_filename", default="./sqrts_vs_beta.dat")
args = parser.parse_args()

sqrts_data = np.genfromtxt(
    args.sqrt_sigma_filename,
    usecols=(0, 1, 3, 4, 5),
    dtype=[("N", "i"), ("L", "i"), ("beta", "d"), ("sqrts", "d"), ("sqrts_err", "d")],
)

fdtype = np.dtype(
    [
        ("N", "i"),
        ("L", "i"),
        ("Nconf", "i8"),
        ("beta", "d"),
        ("TE", "f8"),
        ("scale", "d"),
        ("scale_err", "d"),
        ("chiTC", "d"),
        ("chiTC_err", "d"),
    ]
)
final = np.empty(0, dtype=fdtype)

for fname in args.input_datafiles:
    iN, iL, iB, rawdata = es.topo_load_raw_data(fname)

    tmp_sqrts = np.compress(sqrts_data["beta"] == float(iB), sqrts_data)
    sqrtS = tmp_sqrts["sqrts"][0]
    sqrtS_err = tmp_sqrts["sqrts_err"][0]

    TE_scaled = args.TE * es.Casimir_SP(iN)
    WE_scaled = args.WE * es.Casimir_SP(iN)

    fn_bs = "pkl_flows_bs/pkl_bs_" + iN + "_" + iL + "_" + iB + "_"
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
    TC = es.find_TC(rawdata, t0_tmp_symE[0])
    a_min = es.topo_find_alpha(TC)
    TCdata = np.zeros(len(TC), dtype=[("nconf", "i8"), ("TC", "f8")])
    TCdata["nconf"] = range(1, len(TC) + 1)
    TCdata["TC"] = np.rint(a_min * TC)

    w0_tmp_symE = es.find_w0(w0_flow_symE, WE_scaled)
    TC_w = es.find_TC(rawdata, w0_tmp_symE[0] ** 2)
    a_min = es.topo_find_alpha(TC_w)
    TCdata_w = np.zeros(len(TC_w), dtype=[("nconf", "i8"), ("TC", "f8")])
    TCdata_w["nconf"] = range(1, len(TC_w) + 1)
    TCdata_w["TC"] = np.rint(a_min * TC_w)

    s_TC_avg, s_TC_err = es.bs_avg_err_TC(TCdata["TC"], 20)
    s_TC_w_avg, s_TC_w_err = es.bs_avg_err_TC(TCdata_w["TC"], 20)

    Nconfs = len(TCdata["nconf"])

    print(
        iN,
        iL,
        Nconfs,
        iB,
        sqrtS,
        sqrtS_err,
        TE_scaled,
        t0_tmp_symE[0],
        t0_tmp_symE[1],
        s_TC_avg,
        s_TC_err,
        WE_scaled,
        w0_tmp_symE[0],
        w0_tmp_symE[1],
        s_TC_w_avg,
        s_TC_w_err,
        file=args.output_data,
    )

    print(
        "$",
        iN,
        "$  &  $",
        iL,
        "$  &  $",
        Nconfs,
        "$  &  $",
        iB,
        "$  &  $",
        sqrtS,
        "$  &  $",
        sqrtS_err,
        "$  &  $",
        TE_scaled,
        "$  &  $",
        "{:.2uS}".format(ufloat(t0_tmp_symE[0], t0_tmp_symE[1])),
        "$  &  $",
        "{:.2uS}".format(ufloat(s_TC_avg, s_TC_err)),
        "$  &  $",
        WE_scaled,
        "$  &  $",
        "{:.2uS}".format(ufloat(w0_tmp_symE[0], w0_tmp_symE[1])),
        "$  &  $",
        "{:.2uS}".format(ufloat(s_TC_w_avg, s_TC_w_err)),
        "$ \\\\",
        file=args.output_tex,
    )
