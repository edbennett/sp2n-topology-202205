import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import re
import matplotlib.ticker as ticker
from matplotlib import gridspec
from matplotlib.ticker import AutoLocator
import itertools
import lib_topology as es

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern Roman']
plt.rcParams['text.usetex'] = True
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['figure.figsize']=[3.4,3.4]
#colours = plt.rcParams['axes.color_cycle']
plt.rcParams['errorbar.capsize'] = 2
plt.rcParams['lines.markersize'] = 2
plt.matplotlib.rc('font', size=8)
#plt.style.use('classic')

outdir = os.environ.get('PLOT_DIR', '.')


marks = itertools.cycle( ("v" , "s" , "^" ))

faddr=sys.argv[1]
rawdata = np.genfromtxt(faddr,
                        usecols=(0,1,2,3,4,5,6),
                        names=['N', 'beta','Nupd', 't', 'st','tQ', 'stQ'])

print(rawdata)

plt.figure()
plt.xlabel(r'$\sqrt{t_0}/a$')
plt.ylabel(r'$\ln(\tau_Q N_\textrm{sw})$')

plt.yscale('log')
#plt.xlim(0.,2.6)
Ns = np.unique(rawdata['N'])
for iN in Ns:
	data = np.compress( rawdata['N'] == iN, rawdata)
	x = np.sqrt(data['t'])
	y = data['tQ']*data['Nupd']
	ye = data['Nupd']*data['stQ']
	lab=r'$N_c='+str(int(iN))+'$'
	marker1 = next(marks)
	plt.errorbar(x,y,yerr=ye,linestyle='None',alpha=0.6, label=lab, marker=marker1)

plt.legend(loc=5)
plt.tight_layout()
plt.savefig(outdir + '/tauQ_vs_t0_0.5.pdf')
#plt.show()
