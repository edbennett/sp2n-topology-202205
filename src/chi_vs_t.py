import argparse
import numpy as np
import matplotlib.pyplot as plt
import re
import matplotlib.ticker as ticker
from matplotlib import gridspec
from matplotlib.ticker import AutoLocator, MaxNLocator
from scipy.interpolate import interp1d,UnivariateSpline, CubicSpline
import itertools
import lib_topology as es
import pickle as pkl

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern Roman']
plt.rcParams['text.usetex'] = True
plt.rcParams['lines.linewidth'] = .5
plt.rcParams['figure.figsize']=[3.4,3.4]
#colours = plt.rcParams['axes.color_cycle']
plt.rcParams['errorbar.capsize'] = 2
plt.rcParams['lines.markersize'] = 2
plt.matplotlib.rc('font', size=11)
#plt.style.use('classic')

parser = argparse.ArgumentParser()
parser.add_argument("TE", type=float)
parser.add_argument("WE", type=float)
parser.add_argument("fnames", nargs="+")
parser.add_argument("--plot_dir", default=".")
parser.add_argument("--pickle_dir", default="pkl_flows_bs")
parser.add_argument("--num_bs", type=int, default=es.DEFAULT_NUM_BS)
args = parser.parse_args()

#data_nconf = np.genfromtxt('NCONF_vs_setup.txt', usecols=(0,1,2,3), 
#        dtype=[("N",'i'), ("L",'i'),("beta", 'd'), ("Nconf",'i8')])

def find_TC(rawdata, t):
    confs = np.unique(rawdata['nconf'])
    Nconf = len(confs)
    TC = np.empty(Nconf, dtype='d')
    for i in confs:
        data = np.compress(rawdata['nconf'] == i, rawdata)
        idx = np.abs(data['t']-t).argmin()
        traj = data[idx-3:idx+3]
        f=interp1d(data['t'], data['TC'])
        TC[i-1]=f(t)
        #print(i, TC[i-1])
        #TC = np.append(TC, f(t))
    return TC

def Casimir_SP(N):
    return (float(N)+1.)/4.

for fname in args.fnames:

    # parse the parameters of the ensemble and
    # the T charge of the flowed configs
    N,L,beta,TCdata =es.topo_load_raw_data(fname)
    
    Nconf= int( np.max(TCdata['nconf']))
    #NCONF = np.compress( data_nconf['beta']==float(beta), data_nconf['Nconf'])
    #rawdata = np.compress( TCdata['nconf'] < NCONF+1, TCdata)
    #TE_scaled = TE*Casimir_SP(iN)
    TE_scaled = args.TE * es.Casimir_SP(N)
    WE_scaled = args.WE * es.Casimir_SP(N)

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


    # the following for cycle computes the top susceptibility
    # at every value of the flow-time. Note that there is
    # alpha rounding (see line 83)
    chivst = np.empty(0,dtype=[('t','d'),('amin','d'),('chi','d'),('chierr','d')])
    trange = np.unique(TCdata['t'])
    for it in trange:
        TC_tmp = np.compress(TCdata['t'] == it, TCdata)
        TC = TC_tmp['TC']
        a_min = es.topo_find_alpha(TC)
        TCf = np.rint(a_min*TC)

        chi, err_chi = es.bs_avg_err_TC(TCf)
        chi_tmp = np.array([(it,a_min, chi, err_chi)], dtype=[('t','d'),('amin','d'),('chi','d'),('chierr','d')])
        chivst = np.append(chivst, chi_tmp)
        


    plt.figure()
    plt.xlabel(r'$t/a^2$')
    plt.ylabel(r'$\chi_L a^4$')
    lab=r'$N_c='+str(int(N))+'$, $\\beta='+str(beta)+'$'
    #print(lab)
    
    #find the t_0 corresponding to the TE WE parameters
    t0_tmp_symE = es.find_t0(bs_flow_symE, TE_scaled)
    TC = find_TC(TCdata, t0_tmp_symE[0])
    # plot a vertical line corresponding to the value of t_0
    plt.axvline(x=t0_tmp_symE[0], linestyle='dotted',
                color='r')
    plt.plot([np.nan],[np.nan],label=r'$c_e='+str(args.TE)+'$', color='r')

    color0 = next(plt.gca()._get_lines.prop_cycler)['color']
    # compute and plot the susceptibility corresponding to 
    # that value of t_0
    a_min = es.topo_find_alpha(TC)
    TCdata = np.zeros(len(TC) , dtype=[("nconf",'i8'),("TC",'f8')])
    TCdata['nconf'] = range(1,len(TC)+1)
    TCdata['TC'] = np.rint(a_min*TC)
    s_TC_avg, s_TC_err=es.bs_avg_err_TC(TCdata['TC'],args.num_bs)
    plt.errorbar(t0_tmp_symE[0], s_TC_avg, yerr=s_TC_err, marker='o',
                 color='r')
    
    #Plot the chi as a function of flow-time t
    plt.plot([np.nan],[np.nan],label=lab, color='r')
    plt.plot(chivst['t'], chivst['chi'], color='r')
    plt.fill_between(chivst['t'], chivst['chi']-chivst['chierr'], chivst['chi']+chivst['chierr'], 
            alpha=0.3, color='r')


    plt.legend(loc=1, frameon=False)
    
    fname_out = (
        args.plot_dir + "/chivst_"+str(N)+"_"+str(L)+"_"+str(beta)+".pdf"
    )
    plt.tight_layout()
    plt.savefig(fname_out)
