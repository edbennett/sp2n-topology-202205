import os
import sys
import numpy as np
import matplotlib.pyplot as plt
#import re
#import matplotlib.ticker as ticker
#from scipy.interpolate import interp1d,UnivariateSpline, CubicSpline,PchipInterpolator
#from matplotlib import gridspec
#from matplotlib.ticker import AutoLocator
from uncertainties import ufloat, unumpy
#import itertools
import lib_topology as es
import pickle as pkl

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern Roman']
plt.rcParams['text.usetex'] = True
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['figure.figsize']=[3.4 ,6.8]
#colours = plt.rcParams['axes.color_cycle']
plt.rcParams['errorbar.capsize'] = 2
plt.rcParams['lines.markersize'] = 2
plt.matplotlib.rc('font', size=11)
#plt.style.use('classic')

outdir = os.environ.get('PLOT_DIR', '.')

#marks = itertools.cycle( ("v" , "s" , "^" ))


#def find_t0(indata, TE):
#    t0 = []
#    #try:
#    for i in np.unique( indata['bs'] ):
#        data = np.compress( indata['bs'] == i , indata)
#        idx = np.abs( data['flow']-TE ).argmin()
#        traj = data[idx-5:idx+5]
#        f = interp1d( traj['flow'], traj['t'])
#        t0 = np.append(t0,f(TE))
#    a,b,c,d =es.bstrap( t0, 200)
#    #except:
#    #    a=np.nan
#    #    c=np.nan
#    return a,c
#
#def find_w0(indata, TE):
#    t0 = []
#    #try:
#    for i in np.unique( indata['bs'] ):
#        data = np.compress( indata['bs'] == i , indata)
#        idx = np.abs( data['flow']-TE ).argmin()
#        traj = data[idx-5:idx+5]
#        f = interp1d( traj['flow'], traj['t'])
#        t0 = np.append(t0,np.sqrt(f(TE)))
#    a,b,c,d =es.bstrap( t0, 200)
#    #except:
#    #    a=np.nan
#    #    c=np.nan
#    return a,c
#
#def flows(rawdata, N_bs, obs):
#    Nconf = np.max(np.unique(rawdata['nconf']))
#    ts = np.unique(rawdata['t'])
#    delta_t = ts[1]-ts[0]
#    bs_flow = np.empty(0, dtype=[('bs','i'),('t','f8'),('flow','f8')])
#    w0_flow = np.empty(0, dtype=[('bs','i'),('t','f8'),('flow','f8')])
#    for i_bs in range(N_bs):
#        tmp_bs_flow = np.empty(0, dtype=[('nconf','int'),('t','f8'),('E','f8'),('t2E','f8'),('tsym','f8'),('t2symE','f8'),('TC','f8')])
#        sam = np.random.randint(0,Nconf+1,size=100)
#        for s in sam:
#            tmp = np.compress(rawdata['nconf'] == s, rawdata)
#            tmp_bs_flow = np.append(tmp_bs_flow, tmp)
#        for i_t in range(len(ts)-1):
#            f0tmp = np.compress( tmp_bs_flow['t'] == ts[i_t], tmp_bs_flow)
#            t0_tmp = ts[i_t]
#            f0_avg = np.average( f0tmp[obs] )
#            f1tmp = np.compress( tmp_bs_flow['t'] == ts[i_t+1], tmp_bs_flow)
#            t1_tmp = ts[i_t+1]
#            f1_avg = np.average( f1tmp[obs] )
#            delta_t = t1_tmp - t0_tmp
#            bs_fl_tmp = np.array([(i_bs, t0_tmp, f0_avg)], dtype=[('bs','i'),('t','f8'),('flow','f8')])
#            bs_flow = np.append( bs_flow, bs_fl_tmp)
#            w0_fl_tmp = np.array([(i_bs, t0_tmp, t0_tmp*(f1_avg-f0_avg)/delta_t)], dtype=[('bs','i'),('t','f8'),('flow','f8')])
#            w0_flow = np.append( w0_flow, w0_fl_tmp)
#    return bs_flow, w0_flow
#
#def avg_flows(in_flow):
#    ts = np.unique(in_flow['t'])
#    Lts = len(ts)
#    out_flow = np.empty(Lts, dtype=[('t','f8'),('flow','f8'),('err','f8')])
#    for i_t in range(Lts):
#        out_flow['t'][i_t] = ts[i_t]
#        toavg = np.compress(in_flow['t']==ts[i_t], in_flow['flow'])
#        out_flow['flow'][i_t] = np.average(toavg)
#        out_flow['err'][i_t] = np.std(toavg)
#    return out_flow



TE=float(sys.argv[1])
WE=float(sys.argv[2])


fig, ax = plt.subplots(2,1)
for faddr in sys.argv[3:]:
    

    N,L,beta,rawdata=es.topo_load_raw_data(faddr)
    TE_scaled = TE*es.Casimir_SP(N)
    WE_scaled = WE*es.Casimir_SP(N)
    
    plaq_fname=faddr+'_plaq'
    plaq_data = np.genfromtxt(plaq_fname)
    plaq_bs=es.bstrap(plaq_data,len(plaq_data))
    plaq_u = ufloat(plaq_bs[0], plaq_bs[1])
    print(N,L,beta,'{:.2uS}'.format(plaq_u))
    lthooft = es.d_g_SPN(int(N))*plaq_u/float(beta)
    print(N,L,beta,'{:.2uS}'.format(lthooft))
    
    print("loading flow files")
    fn_bs='pkl_flows_bs/pkl_bs_'+N+'_'+L+'_'+beta+'_'
    infile = open(fn_bs+'t_symE','rb')
    bs_flow_symE = pkl.load(infile)
    infile.close()
    infile = open(fn_bs+'w_symE','rb')
    w0_flow_symE = pkl.load(infile)
    infile.close()
    
    t0_tmp_symE = es.find_t0(bs_flow_symE, TE_scaled)
    w0_tmp_symE = es.find_w0(w0_flow_symE, WE_scaled)


    ax[0].set_xlabel(r'$t/t_0$')
    ax[0].set_ylabel(r'$\mathcal{E}(t)/C_2(F)$')
    ax[0].axhline(y=TE, linestyle='dotted', alpha=0.7)
    color1 = next(plt.gca()._get_lines.prop_cycler)['color']
    lab1=r'$N_c='+str(int(N))+'$, $\\tilde{\lambda}='+str(np.around(lthooft.n,3))+'$, cl.'
    plot_flow = es.avg_flows(bs_flow_symE)
    ax[0].fill_between(plot_flow['t']/t0_tmp_symE[0], 
                     y1=(plot_flow['flow']+plot_flow['err'])/(es.Casimir_SP(N)), 
                     y2=(plot_flow['flow']-plot_flow['err'])/(es.Casimir_SP(N)), 
                     label=lab1, color=color1, alpha=0.5)
    ax[0].legend(loc=2, prop={'size':7}, frameon=False)
    
    ax[1].set_xlabel(r'$t/w_0^2$')
    ax[1].set_ylabel(r'$\mathcal{W}(t)/C_2(F)$')
    ax[1].axhline(y=WE, linestyle='dotted', alpha=0.7)
    plot_flow = es.avg_flows(w0_flow_symE)
    ax[1].fill_between(plot_flow['t']/w0_tmp_symE[0]**2, 
                     y1=(plot_flow['flow']+plot_flow['err'])/(es.Casimir_SP(N)), 
                     y2=(plot_flow['flow']-plot_flow['err'])/(es.Casimir_SP(N)), 
                     label=lab1, color=color1, alpha=0.5)
    ax[1].legend(loc=2, prop={'size':7}, frameon=False )
plt.tight_layout()
#plt.savefig('non_scaled_flows_'+str(int(N))+'_'+str(beta)+'.pdf')
plt.savefig(outdir + '/flows_'+str(TE)+'_'+str(WE)+'_scaled.pdf')
#plt.show()
