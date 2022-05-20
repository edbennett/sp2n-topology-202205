import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import lib_topology as es
import pickle as pkl

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern Roman']
plt.rcParams['text.usetex'] = True
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['figure.figsize']=[4.4 ,6.8]
#colours = plt.rcParams['axes.color_cycle']
plt.rcParams['errorbar.capsize'] = 2
plt.rcParams['lines.markersize'] = 2
plt.matplotlib.rc('font', size=11)
#plt.style.use('classic')

outdir = os.environ.get('PLOT_DIR', '.')

TE=float(sys.argv[1])
WE=float(sys.argv[2])

fig, ax = plt.subplots(2,1, sharex=True)
for faddr in sys.argv[3:]:
    N,L,beta,rawdata=es.topo_load_raw_data(faddr)
    TE_scaled = TE*es.Casimir_SP(N)
    WE_scaled = WE*es.Casimir_SP(N)
    
    print("loading flow files")
    fn_bs='pkl_flows_bs/pkl_bs_'+N+'_'+L+'_'+beta+'_'
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


    
    plot_flow = es.avg_flows(bs_flow_symE)

    ax[0].set_ylabel(r'$\mathcal{E}(t)$')
    ax[0].axhline(y=TE_scaled, linestyle='dotted', alpha=0.7)
    color0 = next(plt.gca()._get_lines.prop_cycler)['color']
    lab0=r'$ \\beta='+str(beta)+'$, pl.'
    lab0=r'$N_c='+str(int(N))+'$, $\\beta='+str(beta)+'$, pl.'
    lab0='$\\beta='+str(beta)+'$, pl.'
    ax[0].fill_between(plot_flow['t'], 
        y1=(plot_flow['flow']+plot_flow['err']), 
        y2=(plot_flow['flow']-plot_flow['err']), 
        label=lab0, color=color0, alpha=0.5)
    color1 = next(plt.gca()._get_lines.prop_cycler)['color']
    lab1=r'$N_c='+str(int(N))+'$, $\\beta='+str(beta)+'$, cl.'
    lab1='$\\beta='+str(beta)+'$, cl.'
    plot_flow = es.avg_flows(bs_flow_E)
    ax[0].fill_between(plot_flow['t'], 
        y1=(plot_flow['flow']+plot_flow['err']), 
        y2=(plot_flow['flow']-plot_flow['err']), 
        label=lab1, color=color1, alpha=0.5)
    ax[0].legend(bbox_to_anchor=(1.0,1.0), loc='upper left', frameon=False )
    
    ax[1].set_xlabel(r'$t$')
    ax[1].set_ylabel(r'$\mathcal{W}(t)$')
    ax[1].axhline(y=WE_scaled, linestyle='dotted', alpha=0.7)
    plot_flow = es.avg_flows(w0_flow_symE)
    ax[1].fill_between(plot_flow['t'], 
        y1=(plot_flow['flow']+plot_flow['err']), 
        y2=(plot_flow['flow']-plot_flow['err']), 
        label=lab0, color=color0, alpha=0.5)
    plot_flow = es.avg_flows(w0_flow_E)
    ax[1].fill_between(plot_flow['t'], 
        y1=(plot_flow['flow']+plot_flow['err']), 
        y2=(plot_flow['flow']-plot_flow['err']), 
        label=lab1, color=color1, alpha=0.5)
    ax[1].legend(bbox_to_anchor=(1.0,1.0), loc='upper left', frameon=False )
plt.tight_layout()
plt.savefig(outdir + '/non_scaled_flows_'+str(N)+'_'+str(TE)+'_'+str(WE)+'.pdf')
#plt.show()

