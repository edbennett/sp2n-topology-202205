import os
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.ticker as ticker
from matplotlib.ticker import AutoLocator
import numpy as np
import re
import sys

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern Roman']
plt.rcParams['text.usetex'] = True
plt.rcParams['lines.linewidth'] = 0.3
plt.rcParams['figure.figsize']=[3.54,5.08]
#colours = plt.rcParams['axes.color_cycle']
plt.rcParams['errorbar.capsize'] = 2
plt.rcParams['lines.markersize'] = 2
plt.matplotlib.rc('font', size=9)
#plt.style.use('classic')

outdir = os.environ.get('PLOT_DIR', '.')

patt = re.compile(r"WF_([0-9])_([0-9]+)_([0-9]+.[0-9]+)")
    
@ticker.FuncFormatter
def major_formatter(x, pos):
    return "%d" % x

ngraf=len(sys.argv[1:])
i=0
fig = plt.figure()
gs = gridspec.GridSpec(ngraf, 2, width_ratios=[2.3, 0.6])
gs.update(wspace=0.05, hspace=0.1)

for fname in sys.argv[1:]:
    rawdata = np.genfromtxt(fname, usecols=(0,1,2,3,4,5,6),names=['nconf','t','E','t2E','tsym','t2symE','TC'])
    N, L, beta, = patt.search(fname).groups()
    print(N,L,beta)
    
    tmax=np.max( rawdata['t'] )
    TC = np.compress( rawdata['t'] == tmax, rawdata )
    TC = np.compress( TC['nconf']< 4000, TC)
    maxTC = np.max( np.abs(TC['TC']))+5
    print(tmax)
    

    axs0 = plt.subplot(gs[i,0])
    axs0.yaxis.set_major_formatter(major_formatter)
    axs0.set_ylim( -maxTC, maxTC)
    axs0.set_ylabel(r'$Q_L(t_0)$')
    lab=r'$N_c='+str(int(N))+'$, $\\beta='+str(beta)+'$, $L='+str(L)+'a$'
    axs0.step( TC['nconf'], TC['TC'], label=lab)
    axs0.legend(frameon=False, fontsize=9, loc='upper right')
    axs0.set_xticks([])
    axs1 = plt.subplot(gs[i,1])
    axs1.set_ylim( -maxTC, maxTC)
    axs1.set_xticks([])
    axs1.set_yticklabels([])
    axs1.hist( TC['TC'], histtype='step', orientation='horizontal')
    i=i+1

axs0.xaxis.set_major_locator(AutoLocator())
axs0.set_xlabel(r'Trajectory')
axs1.set_xlabel(r'Count')

gs.tight_layout(fig)
outname=outdir + '/TC_hist_'+N+'.pdf'
#outname='TC_hist_'+N+'_'+str(beta)+'.pdf'
plt.savefig(outname)
plt.show()

#plt.figure()
#plt.xlabel(r'$t$')
#plt.ylabel(r'$t^2 E(sym)$')
#for i in range(len(betalist)):
#    fname = "WF_"+ "%1.2f"%betalist[i] +".dat"
#    rawdata = np.genfromtxt(fname, usecols=(0,1,2,3,4,5,6),names=['nconf','t','E','t2E','tsym','t2symE','TC'])
#    a,c = find_t0(rawdata, TE, deltaT)
#    t0list['beta'][i] = betalist[i]
#    t0list['t'][i] = a
#    t0list['te'][i] = c
#    ts = np.sort(np.unique( rawdata['t'] ))
#    t2symE_max = np.zeros(len(ts))
#    t2symE_min = np.zeros(len(ts))
#    for its in range(len(ts)):
#        data_ts = np.compress( rawdata['t'] == ts[its], rawdata)
#        t2symE_max[its] = np.max( data_ts['t2symE'])
#        t2symE_min[its] = np.min( data_ts['t2symE'])
#    plt.fill_between( ts, y1=t2symE_min, y2=t2symE_max, alpha=0.3, label=r'$\beta = '+str(betalist[i])+'$')
#plt.axhline(y=TE)
#plt.legend(loc=2)
#plt.show()
#
#plt.figure()
#plt.xlabel(r'$\beta$')
#plt.ylabel(r'$\sqrt{8 t_0^{sym}}$')
#plt.title(r'Reference value $t_0^2 E(sym)='+str(TE)+'$')
#plt.errorbar( t0list['beta'], np.sqrt(8.*t0list['t']), yerr = 1/(2.*np.sqrt(8.*t0list['t']))*8.*t0list['te'], linestyle='None', marker="^")
#plt.show()
