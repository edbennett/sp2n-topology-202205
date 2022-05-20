import sys
import numpy as np
import matplotlib.pyplot as plt
import re
import matplotlib.ticker as ticker
#from matplotlib import gridspec
#from matplotlib.ticker import AutoLocator
from uncertainties import ufloat, unumpy
from scipy.interpolate import interp1d,UnivariateSpline, CubicSpline
import itertools
import lib_topology as es
import pickle as pkl

#plt.rcParams['font.family'] = 'serif'
#plt.rcParams['font.serif'] = ['Computer Modern Roman']
#plt.rcParams['text.usetex'] = True
#plt.rcParams['lines.linewidth'] = 1
#plt.rcParams['figure.figsize']=[1.8*3.375,1.8*3.375]
##colours = plt.rcParams['axes.color_cycle']
#plt.rcParams['errorbar.capsize'] = 2
#plt.rcParams['lines.markersize'] = 2
#plt.matplotlib.rc('font', size=11)
##plt.style.use('classic')
#
#marks = itertools.cycle( ("v" , "s" , "^" ))
#
#@ticker.FuncFormatter
#def major_formatter(x, pos):
#    return "%d" % x
#data_nconf = np.genfromtxt('NCONF_vs_setup.txt', usecols=(0,1,2,3), dtype=[("N",'i'), ("L",'i'),("beta", 'd'), ("Nconf",'i8')])

sqrts_data = np.genfromtxt('sqrts_vs_beta.dat', usecols=(0,1,3,4,5), dtype=[("N",'i'), ("L",'i'), ("beta",'d'),('sqrts','d'),('sqrts_err','d')]) 

fdtype = np.dtype([("N",'i'), ("L",'i'), ("Nconf",'i8'),("beta", 'd'),("TE",'f8'),("scale",'d'),("scale_err",'d'),("chiTC",'d'),("chiTC_err",'d')])
final = np.empty(0, dtype=fdtype)

TE=float(sys.argv[1])
WE=float(sys.argv[2])
for fname in sys.argv[3:]:
    iN, iL, iB, rawdata = es.topo_load_raw_data(fname)

    tmp_sqrts = np.compress( sqrts_data['beta']==float(iB), sqrts_data)
    sqrtS = tmp_sqrts['sqrts'][0]
    sqrtS_err = tmp_sqrts['sqrts_err'][0]
    
    TE_scaled = TE*es.Casimir_SP(iN)
    WE_scaled = WE*es.Casimir_SP(iN)

    fn_bs='pkl_flows_bs/pkl_bs_'+iN+'_'+iL+'_'+iB+'_'
    infile = open(fn_bs+'t_E','rb')
    bs_flow_E = pkl.load(infile)
    infile.close()
    infile = open(fn_bs+'w_E','rb')
    w0_flow_E = pkl.load(infile)
    infile.close()
    infile = open(fn_bs+'t_symE','rb')
    bs_flow_symE = pkl.load(infile)
    infile.close()
    infile = open(fn_bs+'w_symE','rb')
    w0_flow_symE = pkl.load(infile)
    infile.close()


    t0_tmp_symE = es.find_t0(bs_flow_symE, TE_scaled)
    TC = es.find_TC(rawdata, t0_tmp_symE[0])
    a_min = es.topo_find_alpha(TC)
    TCdata = np.zeros(len(TC) , dtype=[("nconf",'i8'),("TC",'f8')])
    TCdata['nconf'] = range(1,len(TC)+1)
    TCdata['TC'] = np.rint(a_min*TC)

    w0_tmp_symE = es.find_w0(w0_flow_symE, WE_scaled)
    TC_w = es.find_TC(rawdata, w0_tmp_symE[0]**2)
    a_min = es.topo_find_alpha(TC_w)
    TCdata_w = np.zeros(len(TC_w) , dtype=[("nconf",'i8'),("TC",'f8')])
    TCdata_w['nconf'] = range(1,len(TC_w)+1)
    TCdata_w['TC'] = np.rint(a_min*TC_w)

    s_TC_avg, s_TC_err=es.bs_avg_err_TC(TCdata['TC'],20)
    s_TC_w_avg, s_TC_w_err=es.bs_avg_err_TC(TCdata_w['TC'],20)
    
    Nconfs = len(TCdata['nconf'])

    print(iN,iL,Nconfs, iB, sqrtS, sqrtS_err, 
                        TE_scaled,
                         t0_tmp_symE[0], t0_tmp_symE[1],
                         s_TC_avg, s_TC_err,
                        WE_scaled,
                         w0_tmp_symE[0], w0_tmp_symE[1],
                         s_TC_w_avg, s_TC_w_err)
    print("$",iN,"$  &  $",iL,"$  &  $",Nconfs,"$  &  $", iB,"$  &  $", 
                        sqrtS,"$  &  $",sqrtS_err, "$  &  $",
                        TE_scaled,"$  &  $",
                        '{:.2uS}'.format(ufloat(t0_tmp_symE[0], t0_tmp_symE[1])),"$  &  $",
                         '{:.2uS}'.format(ufloat(s_TC_avg, s_TC_err)),"$  &  $",
                        WE_scaled,"$  &  $",
                         '{:.2uS}'.format(ufloat(w0_tmp_symE[0], w0_tmp_symE[1])),"$  &  $",
                         '{:.2uS}'.format(ufloat(s_TC_w_avg, s_TC_w_err)),"$ \\\\", 
                        file=open("table_chi_t0_w0.tex", "a"))
