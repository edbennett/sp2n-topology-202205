import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import re
import matplotlib.ticker as ticker
from matplotlib import gridspec
from matplotlib.ticker import AutoLocator
from uncertainties import ufloat,umath,unumpy
import itertools
import lib_topology as es
from scipy.optimize import curve_fit, least_squares
from scipy import linalg

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern Roman']
plt.rcParams['text.usetex'] = True
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['figure.figsize']=[5.375,5.375]
#colours = plt.rcParams['axes.color_cycle']
plt.rcParams['errorbar.capsize'] = 2
plt.rcParams['lines.markersize'] = 4
plt.matplotlib.rc('font', size=12)
#plt.style.use('classic')

outdir = os.environ.get('PLOT_DIR', '.')

def f1(x,a,b):
	return a+b*x

def f2(x,a,b,c):
	return a+b*x+c*x**2

def scaling_SUN(n):
	return (n**2-1.)/(4.*n**2)

def scaling_SPN(n):
	return (2.*n+1)/(16.*n)


Nmin1 = 1
Nmin2 = 1

files_SU = { "Athenodorou et al.":"clim_SU_AA_MT.dat",
			"Del Debbio et al.":"clim_SU_DD_EV_HP.dat",
			"Bonati et al.":"clim_SU_CB_MD.dat",
			"Bonanno et al.":"clim_SU_CB_CB_MD.dat",
			"Lucini et al.":"clim_SU_BL_MT.dat"}

files_Sp = { "Bennett et al.":"clim_SP.dat"}


tot_dtype = np.dtype([('label','U100'),('N_c','i'),('d_G','d'),('sc_chi','d'),('sc_chi_err','d')])
chi_total = np.empty(0,dtype=tot_dtype)

plt.figure()
#plt.title(r'$N_{\mathrm{SU}}>'+str(Nmin1)+'$, $N_{\mathrm{Sp}}>'+str(Nmin2)+'$')
plt.xlabel(r'$1/d_G$')
plt.ylabel(r'${\displaystyle \frac{\chi}{\sigma^2} \frac{C_2(F)^2}{d_G}}$')
plt.ylim(0., 0.015)
plt.xlim(-0.05, 0.4)
xr = np.arange(0., 0.65,0.01)

for key,val in files_SU.items():
	fname=val
	chi_SUN_raw = np.genfromtxt(fname, usecols=(0,1,2),names=['N','chi', 'schi'])
	chi_SUN= np.compress( chi_SUN_raw['N']>Nmin1, chi_SUN_raw)
	#xdata = 1/chi_SUN_raw['N']
	xdata = 1./es.d_g_SUN(chi_SUN_raw['N'])
	ydata = chi_SUN_raw['chi']*scaling_SUN(chi_SUN_raw['N'])
	ydata_err = chi_SUN_raw['schi']*scaling_SUN(chi_SUN_raw['N'])
	lab=key
	color1 = next(plt.gca()._get_lines.prop_cycler)['color']
	plt.errorbar(xdata , ydata, yerr=ydata_err, linestyle='None', marker='^', color=color1, label=lab, alpha=0.7)
	
	tmp_total = np.empty(len(xdata),dtype=tot_dtype)
	tmp_total['label']=key
	tmp_total['N_c']=chi_SUN_raw['N']
	tmp_total['d_G']=xdata
	tmp_total['sc_chi']=ydata
	tmp_total['sc_chi_err']=ydata_err
	chi_total = np.append(chi_total, tmp_total)

xdata = chi_total['d_G']
ydata = chi_total['sc_chi']
ydata_err = chi_total['sc_chi_err']
popt, pcov = curve_fit(f1, xdata, ydata, sigma=ydata_err, absolute_sigma=True)
chi2 = np.sum( (ydata-f1(xdata,*popt))**2/ydata_err**2)/( len(ydata)-len(popt))
perr = np.sqrt(np.diag(pcov))
print("Only SU - a+bx")
print(popt,perr, chi2)
print("a = ", '{:.2uS}'.format(ufloat(popt[0],perr[0])))
print("b = ", '{:.2uS}'.format(ufloat(popt[1],perr[1])))
print(" chi2 = ", np.around(chi2,2))
print(" DOF = ", len(ydata)-len(popt))

for key,val in files_Sp.items():
	fname=val
	chi_SPN_raw = np.genfromtxt(fname, usecols=(0,1,2),names=['N','chi', 'schi'])
	chi_SPN = np.compress( chi_SPN_raw['N']>Nmin2, chi_SPN_raw)
	#xdata = 1/(2.*chi_SPN_raw['N'])
	xdata = 1./es.d_g_SPN(chi_SPN_raw['N'])
	ydata = chi_SPN_raw['chi']*scaling_SPN(chi_SPN_raw['N'])
	ydata_err = chi_SPN_raw['schi']*scaling_SPN(chi_SPN_raw['N'])
	lab=key
	color1 = next(plt.gca()._get_lines.prop_cycler)['color']
	plt.errorbar(xdata , ydata, yerr=ydata_err, linestyle='None', marker='o', color=color1, label=lab, alpha=0.7)

	tmp_total = np.empty(len(xdata),dtype=tot_dtype)
	tmp_total['label']=key
	tmp_total['N_c']=chi_SPN_raw['N']
	tmp_total['d_G']=xdata
	tmp_total['sc_chi']=ydata
	tmp_total['sc_chi_err']=ydata_err
	chi_total = np.append(chi_total, tmp_total)
    
	plt.axhline(y=1./(4.*np.pi)**2, linestyle='dotted', alpha=0.7)


print("Sp and SU - a+bx")
xdata = chi_total['d_G']
ydata = chi_total['sc_chi']
ydata_err = chi_total['sc_chi_err']
popt, pcov = curve_fit(f1, xdata, ydata, sigma=ydata_err, absolute_sigma=True)
chi2 = np.sum( (ydata-f1(xdata,*popt))**2/ydata_err**2)/( len(ydata)-len(popt))
perr = np.sqrt(np.diag(pcov))
print(popt,perr, chi2)
print("a = ", '{:.2uS}'.format(ufloat(popt[0],perr[0])))
print("b = ", '{:.2uS}'.format(ufloat(popt[1],perr[1])))
print(" chi2 = ", np.around(chi2,2))
print(" DOF = ", len(ydata)-len(popt))
color1 = next(plt.gca()._get_lines.prop_cycler)['color']
lab=r'$a+b/N_c^2$, $\chi^2/DOF ='+str(np.around(chi2,2))+'$'
plt.plot(xr, f1(xr,*popt), linestyle='dotted',color=color1)
#, label=lab)
plt.errorbar( 0., popt[0], yerr=perr[0], linestyle='None', marker='o', color=color1)

print(" a+ bx+cx^2")
xdata = chi_total['d_G']
ydata = chi_total['sc_chi']
ydata_err = chi_total['sc_chi_err']
popt, pcov = curve_fit(f2, xdata, ydata, sigma=ydata_err, absolute_sigma=True)
chi2 = np.sum( (ydata-f2(xdata,*popt))**2/ydata_err**2)/( len(ydata)-len(popt))
perr = np.sqrt(np.diag(pcov))
print(popt,perr, chi2)
print(" a = ", '{:.3uS}'.format(ufloat(popt[0],perr[0])))
print(" b = ", '{:.3uS}'.format(ufloat(popt[1],perr[1])))
print(" c = ", '{:.3uS}'.format(ufloat(popt[2],perr[2])))
print(" chi2 = ", np.around(chi2,2))
print(" DOF = ", len(ydata)-len(popt))
color1 = next(plt.gca()._get_lines.prop_cycler)['color']
lab=r'$a+b/N_c+c/N_c^2$, $\chi^2/DOF ='+str(np.around(chi2,2))+'$'
plt.plot(xr, f2(xr,*popt), linestyle='dashdot',color=color1)
#, label=lab)
plt.errorbar( 0., popt[0], yerr=perr[0], linestyle='None', marker='o', color=color1)

plt.tight_layout()
plt.legend(loc=2)
plt.savefig(outdir + '/ScaledChi.pdf')
#plt.show()

plt.figure()
#plt.title(r'$N_{\mathrm{SU}}>'+str(Nmin1)+'$, $N_{\mathrm{Sp}}>'+str(Nmin2)+'$')
plt.xlabel(r'$1/N_c$')
plt.ylabel(r'${\displaystyle \frac{\chi}{\sigma^2}}$')
plt.ylim(0., 0.085)
plt.xlim(-0.05, 0.6)
xr = np.arange(0., 0.65,0.01)

for key,val in files_SU.items():
	fname=val
	chi_SUN_raw = np.genfromtxt(fname, usecols=(0,1,2),names=['N','chi', 'schi'])
	chi_SUN= np.compress( chi_SUN_raw['N']>Nmin1, chi_SUN_raw)
	xdata = 1/chi_SUN_raw['N']
	#xdata = 1./es.d_g_SUN(chi_SUN_raw['N'])
	ydata = chi_SUN_raw['chi']
	ydata_err = chi_SUN_raw['schi']
	lab=key
	color1 = next(plt.gca()._get_lines.prop_cycler)['color']
	plt.errorbar(xdata , ydata, yerr=ydata_err, linestyle='None', marker='^', color=color1, label=lab, alpha=0.7)
	
	tmp_total = np.empty(len(xdata),dtype=tot_dtype)
	tmp_total['label']=key
	tmp_total['N_c']=chi_SUN_raw['N']
	tmp_total['d_G']=xdata
	tmp_total['sc_chi']=ydata
	tmp_total['sc_chi_err']=ydata_err
	chi_total = np.append(chi_total, tmp_total)

for key,val in files_Sp.items():
	fname=val
	chi_SPN_raw = np.genfromtxt(fname, usecols=(0,1,2),names=['N','chi', 'schi'])
	chi_SPN = np.compress( chi_SPN_raw['N']>Nmin2, chi_SPN_raw)
	xdata = 1/(2.*chi_SPN_raw['N'])
	#xdata = 1./es.d_g_SPN(chi_SPN_raw['N'])
	ydata = chi_SPN_raw['chi']
	ydata_err = chi_SPN_raw['schi']
	lab=key
	color1 = next(plt.gca()._get_lines.prop_cycler)['color']
	plt.errorbar(xdata , ydata, yerr=ydata_err, linestyle='None', marker='o', color=color1, label=lab, alpha=0.7)

	tmp_total = np.empty(len(xdata),dtype=tot_dtype)
	tmp_total['label']=key
	tmp_total['N_c']=chi_SPN_raw['N']
	tmp_total['d_G']=xdata
	tmp_total['sc_chi']=ydata
	tmp_total['sc_chi_err']=ydata_err
	chi_total = np.append(chi_total, tmp_total)

for i in chi_total:
	print_val = ufloat(i['sc_chi'], i['sc_chi_err'])
	if( i['label'] == "Bennett et al."):
		scaling = np.around(scaling_SPN(int(i['N_c'])),4)
		group="$ Sp("
	else:
		scaling = np.around(scaling_SUN(int(i['N_c'])),4)
		group="$ SU("
	print( group, int(i['N_c']), ")$ &", np.around(i['d_G'],2), "$ & $", 
		i['label']," & $", '{:.2uS}'.format(print_val), "$ & $", scaling, "$ \\\\")

print(chi_total)
plt.tight_layout()
plt.legend(loc=2)
plt.savefig(outdir + '/NONScaledChi.pdf')
#plt.show()
