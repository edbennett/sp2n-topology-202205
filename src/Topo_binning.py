import csv
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import re
import matplotlib.ticker as ticker
from matplotlib import gridspec
from matplotlib.ticker import AutoLocator
from uncertainties import ufloat
from scipy.interpolate import interp1d,UnivariateSpline, CubicSpline
import itertools
import lib_topology as es
import pickle as pkl

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern Roman']
plt.rcParams['text.usetex'] = True
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['figure.figsize']=[1.8*3.375,1.8*3.375]
#colours = plt.rcParams['axes.color_cycle']
plt.rcParams['errorbar.capsize'] = 2
plt.rcParams['lines.markersize'] = 2
plt.matplotlib.rc('font', size=11)
#plt.style.use('classic')

#def find_TC(rawdata, t):
#    confs = np.unique(rawdata['nconf'])
#    Nconf = len(confs)
#    TC = np.empty(Nconf, dtype='d')
#    for i in confs:
#        data = np.compress(rawdata['nconf'] == i, rawdata)
#        idx = np.abs(data['t']-t).argmin()
#        traj = data[idx-3:idx+3]
#        f=interp1d(data['t'], data['TC'])
#        TC[i-1]=f(t)
#        #print(i, TC[i-1])
#        #TC = np.append(TC, f(t))
#    return TC

def binned_std(din, bin_size):
    bin1 = []
    for i in range(int(len(din)/bin_size)):
        bin1.append(np.average(din[ i*bin_size : (i+1)*bin_size ]))

    resampled1 = []
    for j in range(200):
    	sam = np.random.randint(0,len(bin1),size=len(bin1))
    	avg1 = np.average( [ bin1[i] for i in sam ] ) 
    	resampled1.append( avg1 )
    
    #print(N, L, beta)
    a=np.std(resampled1)
    
    return a

#def Casimir_SP(N):
#    return (float(N)+1.)/4.

#for i in bin_range:
#    s_TC_avg, s_TC_err=es.bs_avg_err_TC(TC,i )
#    s_TC_avg0, s_TC_err0=es.bs_avg_err_TC(TC,1)
#    print(i, s_TC_avg/s_TC_avg0, s_TC_err/s_TC_err0)
#
#print(binned_res)

DATA_DIR = os.environ.get('DATA_DIR', '.')
PROC_DIR = os.environ.get('PROC_DIR', '.')

def get_n_sw(N, L, beta):
    with open(DATA_DIR + '/DATA_FILES', 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if row[0] == N and row[1] == L and row[2] == beta:
                return int(row[3])
    raise ValueError(f'{N}, {L}, {beta} not found in DATA_FILES')


bin_range=np.arange(1,100,1)

TE=float(sys.argv[1])
WE=float(sys.argv[2])
outf = PROC_DIR + '/tauQ_vs_t0_'+str(TE)+'_w0_'+str(WE)+'.dat'
f=open(outf,'a')
for fname in sys.argv[3:]:
    N,L,beta,TCdata =es.topo_load_raw_data(fname)
    TE_scaled = TE
    WE_scaled = WE

    #print("loading flow files")
    fn_bs='pkl_flows_bs/pkl_bs_'+N+'_'+L+'_'+beta+'_'
    infile = open(fn_bs+'t_symE','rb')
    bs_flow_symE = pkl.load(infile)
    infile.close()
    #infile = open(fn_bs+'w_symE','rb')
    #w0_flow_symE = pkl.load(infile)
    #infile.close()
    
    t0_tmp_symE = es.find_t0(bs_flow_symE, TE_scaled)
    #w0_tmp_symE = es.find_w0(w0_flow_symE, WE_scaled)


    TC = es.find_TC(TCdata, t0_tmp_symE[0])
    autC, autC_err = es.autocorr(TC)
    tmp=ufloat(autC, autC_err)
    n_sw = get_n_sw(N, L, beta)
    print(N,L,beta, n_sw, t0_tmp_symE[0], t0_tmp_symE[1], autC,autC_err, file=f)
    print(N,L,beta, n_sw, '{:.2uS}'.format(tmp))

    plt.figure()
    plt.ylim(0.0,5.0)
    plt.xlabel(r'b')
    plt.ylabel(r'$\frac{\sigma_b}{\sigma}$')
    name_postfix='_'

    #autC, autC_err = es.autocorr(TC)
    #tmp=ufloat(autC, autC_err)
    #print(N,L,beta, len(TC), '{:.2uS}'.format(tmp))
    #print(N,L,beta, len(TC), autC, autC_err)
    std_1=binned_std(TC,1)
    binned_res=[ (binned_std(TC,i)/std_1)**2 for i in bin_range]
    lab=r'$\beta='+str(beta)+'$'
    plt.plot( bin_range, binned_res, label=lab)
    name_postfix=name_postfix+beta
    plt.legend(frameon=False)
    fname_out='binning_aut_'+str(N)+name_postfix+'.pdf'
    plt.savefig(PROC_DIR + '/' + fname_out)
    #plt.show()



