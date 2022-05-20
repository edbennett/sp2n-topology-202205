import sys
import numpy as np
import matplotlib.pyplot as plt
import re
import matplotlib.ticker as ticker
from matplotlib import gridspec
from matplotlib.ticker import AutoLocator
import itertools
import effstring_module as es
from scipy.optimize import curve_fit
from uncertainties import ufloat,umath

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern Roman']
plt.rcParams['text.usetex'] = True
plt.rcParams['lines.linewidth'] =1 
plt.rcParams['figure.figsize']=[3.74,3.54]
#colours = plt.rcParams['axes.color_cycle']
plt.rcParams['errorbar.capsize'] = 2
plt.rcParams['lines.markersize'] = 2
plt.matplotlib.rc('font', size=8)
#plt.style.use('classic')

marks = itertools.cycle( ("v" , "s" , "^" ))
@ticker.FuncFormatter
def major_formatter(x, pos):
    return "%d" % x

def f(x,a,b):
    return a + b*x

def f2(x,a,b):
    return a + b*x*x

def g(x,a,b,c):
    return a+b*x+c*x**2


faddr=sys.argv[1]

patt = re.compile(r"chi_vs_t0_([0-9]+.[0-9]+)_w0_([0-9]+.[0-9]+)_([a-z]+).dat")
t0, w0, labf, = patt.search(faddr).groups()
print(t0,w0,labf)
chi_SPN = np.genfromtxt(faddr,
                        usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),
                        names=['N','L','nconf','beta','sqrtS','sqrtS_err','TE', 'tscale', 'tscale_err','chi_t', 'schi_t','WE', 'wscale', 'wscale_err','chi_w', 'schi_w'])

#print('$N$ & $L$ & $Nconf$ & $\beta$ & $\mathcal{E}_0$ & $t_0$ & $\chi t_0^2$ & $\mathcal{W}_0$ & $w_0$ & $\chi w_0$ \\\\')
#for i in range(len(chi_SPN)):
#    for c in [ 'N', 'beta']:
#        val=chi_SPN[i][c]
#        print ("$", str(val), '$ &', end="")
#    #val = ufloat(chi_SPN[i]['tscale'], chi_SPN[i]['tscale_err'])
#    #print("$", '{:.2uS}'.format(val), '$', end="")
#    val = ufloat(chi_SPN[i]['chi_t'], chi_SPN[i]['schi_t'])
#    print(" & $", '{:.2uS}'.format(val), '$', end="")
#
#    #c=chi_SPN[i]['WE']
#    #print(" & $", str(np.around(c,4)), '$')
#    #'wscale', 'wscale_err',
#    #val = ufloat(chi_SPN[i]['wscale'], chi_SPN[i]['wscale_err'])
#    #print(" & $", '{:.2uS}'.format(val), '$', end="")
#    #'chi_w', 'schi_w'
#
#    val = ufloat(chi_SPN[i]['chi_w'], chi_SPN[i]['schi_w'])
#    print(" & $", '{:.2uS}'.format(val), '$', end="")
#    print(" \\\\")
#print(chi_SPN)

chi_SUN = np.genfromtxt('TC_SUN_garcia_vera.dat', usecols=(0,1,2,3,4,5,6,7,8),names=['N','T', 'L', 'nconf','beta','scale', 'scale_err','chi', 'schi'])

#print(chi_SPN)
print("#N chi(0) err coeff err chi2")

plt.figure()
plt.xlabel(r'$a^2/t_0$')
plt.ylabel(r'$\chi_L t_0^2$')
#plt.xlim(-0.01, 0.4)
#plt.ylim(0.0, 0.0010)
#plt.ylim(bottom=0.0)

dtclim=np.dtype([('N','i'),('chi','d'),('chi_err','d'),('sl','d'),('sl_err','d'),('chi2','d')])
clim = np.empty(0, dtype=dtclim)

for i in np.unique(chi_SPN['N']):
    pl_data = np.compress( chi_SPN['N'] == i, chi_SPN)
    xdata_t0 = 1./pl_data['tscale']
    ydata_t0 = pl_data['chi_t']*pl_data['tscale']**2/pl_data['L']**4

    ydata_t0_err = np.sqrt( (2.*pl_data['chi_t']*pl_data['tscale']*pl_data['tscale_err'])**2 
                        + (pl_data['schi_t']*pl_data['tscale']**2)**2)/pl_data['L']**4
    xdata_w0 = 1./pl_data['wscale']**2
    ydata_w0 = pl_data['chi_w']*pl_data['wscale']**4/pl_data['L']**4
    ydata_w0_err = np.sqrt( (4.*pl_data['chi_w']*pl_data['wscale']**3*pl_data['wscale_err'])**2
                        + (pl_data['schi_w']*pl_data['wscale']**4)**2)/pl_data['L']**4

    for l in range(len(xdata_t0)):
        val_t0 = 10**3*ufloat(ydata_t0[l], ydata_t0_err[l])
        uval_sqrtS_t0 = ufloat(pl_data['sqrtS'][l]**2.*pl_data['tscale'][l], 
            np.sqrt( (2.*pl_data['sqrtS'][l]*pl_data['sqrtS_err'][l]*pl_data['tscale'][l])**2 + (pl_data['sqrtS'][l]**2*pl_data['tscale_err'][l])**2 ))
        val_w0 = 10**4*ufloat(ydata_w0[l],ydata_w0_err[l])
        uval_sqrtS_w0 = ufloat(pl_data['sqrtS'][l]**2.*pl_data['wscale'][l]**2, 
            np.sqrt( (2.*pl_data['sqrtS'][l]*pl_data['sqrtS_err'][l]*pl_data['wscale'][l]**2)**2 + (pl_data['sqrtS'][l]**2*2.*pl_data['wscale'][l]*pl_data['wscale_err'][l])**2 ))
        print("$", int(pl_data['N'][l]),"$ & $", pl_data['beta'][l],"$ & $", 
                '{:.2uS}'.format(uval_sqrtS_t0), "$ & $", 
                '{:.2uS}'.format(val_t0), "$ & $", 
                '{:.2uS}'.format(uval_sqrtS_w0), "$ & $", 
                '{:.2uS}'.format(val_w0), '$ \\\\')
    
#Continuum extrapolations using the t scale
for i in np.unique(chi_SPN['N']):
    pl_data = np.compress( chi_SPN['N'] == i, chi_SPN)
    marker1 = next(marks)

    xdata = 1./pl_data['tscale']
    ydata = pl_data['chi_t']*pl_data['tscale']**2/pl_data['L']**4

    ydata_err = np.sqrt( (2.*pl_data['chi_t']*pl_data['tscale']*pl_data['tscale_err'])**2 
                        + (pl_data['schi_t']*pl_data['tscale']**2)**2)/pl_data['L']**4


    popt, pcov = curve_fit(f, xdata, ydata, sigma=ydata_err, absolute_sigma=True)
    chi2 = np.sum( (ydata-f(xdata,*popt))**2/ydata_err**2)/( len(ydata)-len(popt))
    perr = np.sqrt(np.diag(pcov))
    lab = r'$N_c= '+str(int(i))+'$, $\\tilde{\mathcal{X}}^2='+str(np.around(chi2,2))+'$'


    color1 = next(plt.gca()._get_lines.prop_cycler)['color']
    xr = np.arange(0., np.max(xdata)*1.15,0.01)
    clim_tmp = np.array([(i,popt[0],perr[0], popt[1], perr[1], chi2)], dtype=dtclim)
    clim = np.append( clim, clim_tmp)
    plt.plot(xr, f(xr,*popt), color=color1, linestyle='dashed')
    #, label=r'$\chi^2/DOF ='+str(np.around(chi2,2))+'$')
    plt.errorbar( 0., popt[0], yerr=perr[0], linestyle='None', marker=marker1, color=color1)
    plt.errorbar(xdata , ydata, yerr=ydata_err, linestyle='None', marker=marker1, color=color1)
    plt.errorbar([np.nan], [np.nan], yerr=[np.nan], color=color1, label=lab, linestyle='dashed')
    #plt.text(0.5, 0.5, 'PRELIMINARY', fontsize=40, color='gray', alpha=0.5,ha='center', va='center', rotation='30')
    
for i in np.unique(chi_SUN['N']):
    pl_data = np.compress( chi_SUN['N'] == i, chi_SUN)
    lab = r'$SU( '+str(int(i))+')$'
    marker1 = next(marks)

    xdata = 1./pl_data['scale']
    ydata = pl_data['chi']
    ydata_err = pl_data['schi']

    popt, pcov = curve_fit(f, xdata, ydata, sigma=ydata_err, absolute_sigma=True)
    chi2 = np.sum( (ydata-f(xdata,*popt))**2/ydata_err**2)/( len(ydata)-2.)
    perr = np.sqrt(np.diag(pcov))

    color1 = next(plt.gca()._get_lines.prop_cycler)['color']
    xr = np.arange(0., 1.,0.05)
    #plt.plot(xr, f(xr,*popt), color=color1)
    #, label=r'$\chi^2/DOF ='+str(np.around(chi2,2))+'$')
    #plt.errorbar( 0., popt[0], yerr=perr[0], linestyle='None', marker=marker1, color=color1)
    #plt.errorbar(xdata , ydata, yerr=ydata_err, linestyle='None', marker=marker1, color=color1)
    #plt.errorbar([np.nan], [np.nan], yerr=[np.nan], color=color1, label=lab)
    #plt.text(0.5, 0.5, 'PRELIMINARY', fontsize=40, color='gray', alpha=0.5,ha='center', va='center', rotation='30')
#plt.legend(loc=1, ncol=2, frameon=False)
plt.legend(bbox_to_anchor=(0.,1.0), loc='lower left', ncol=2, frameon=False)
plt.tight_layout()
plt.savefig('SPN_Topology_contlim_'+str(t0)+'_'+str(w0)+'_'+labf+'.pdf')
plt.show()

for i in clim['N']:
    table_fits = np.compress(clim['N']==i, clim)
    print("N = ", i)
    val = ufloat(table_fits['chi'], table_fits['chi_err'])
    print("     chi t_0^2(a=0 ) = ", '{:.2uS}'.format(val))
    val = ufloat(table_fits['sl'], table_fits['sl_err'])
    print("     coefficient = ", '{:.2uS}'.format(val))
    print("     chi^2 = ", np.around(table_fits['chi2'],2))

for i in clim['N']:
    table_fits = np.compress(clim['N']==i, clim)
    val = ufloat(table_fits['chi'], table_fits['chi_err'])
    c2 = np.around( table_fits['chi2'],2)[0]
    print('${:d}$ & ${:.2uS}$ & ${:f}$ \\\\'.format( i, val, c2 ))

#Continuum extrapolations using the w scale
plt.figure()
plt.xlabel(r'$a^2/w_0^2$')
plt.ylabel(r'$\chi_L w_0^4$')
#plt.xlim(-0.01, 0.4)
#plt.ylim(0.0, 0.0010)
#plt.ylim(bottom=0.0)


clim = np.empty(0, dtype=dtclim)
for i in np.unique(chi_SPN['N']):
    pl_data = np.compress( chi_SPN['N'] == i, chi_SPN)
    marker1 = next(marks)

    xdata = 1./pl_data['wscale']**2
    ydata = pl_data['chi_w']*pl_data['wscale']**4/pl_data['L']**4
    #if( np.max(ydata) > max_y_lim):
    #    max_y_lim = 1.4*np.max(ydata)
    #plt.ylim(top=max_y_lim)
    #if( np.min(ydata) < min_y_lim):
    #    min_y_lim = 0.45*np.max(ydata)
    #plt.ylim(min_y_lim,max_y_lim)
    ydata_err = np.sqrt( (4.*pl_data['chi_w']*pl_data['wscale']**3*pl_data['wscale_err'])**2
                        + (pl_data['schi_w']*pl_data['wscale']**4)**2)/pl_data['L']**4

    popt, pcov = curve_fit(f, xdata, ydata, sigma=ydata_err, absolute_sigma=True)
    chi2 = np.sum( (ydata-f(xdata,*popt))**2/ydata_err**2)/( len(ydata)-len(popt))
    perr = np.sqrt(np.diag(pcov))
    lab = r'$N_c= '+str(int(i))+'$, $\\tilde{\mathcal{X}}^2='+str(np.around(chi2,2))+'$'

    color1 = next(plt.gca()._get_lines.prop_cycler)['color']
    xr = np.arange(0., np.max(xdata)*1.15,0.01)
    clim_tmp = np.array([(i,popt[0],perr[0], popt[1], perr[1], chi2)], dtype=dtclim)
    clim = np.append( clim, clim_tmp)
    plt.plot(xr, f(xr,*popt), color=color1, linestyle='dashed')
    #, label=r'$\chi^2/DOF ='+str(np.around(chi2,2))+'$')
    plt.errorbar( 0., popt[0], yerr=perr[0], linestyle='None', marker=marker1, color=color1)
    plt.errorbar(xdata , ydata, yerr=ydata_err, linestyle='None', marker=marker1, color=color1)
    plt.errorbar([np.nan], [np.nan], yerr=[np.nan], color=color1, label=lab, linestyle='dashed')
    #plt.text(0.5, 0.5, 'PRELIMINARY', fontsize=40, color='gray', alpha=0.5,ha='center', va='center', rotation='30')

for i in clim['N']:
    table_fits = np.compress(clim['N']==i, clim)
    print("N = ", i)
    val = ufloat(table_fits['chi'], table_fits['chi_err'])
    print("     chi w_0(a=0 ) = ", '{:.2uS}'.format(val))
    val = ufloat(table_fits['sl'], table_fits['sl_err'])
    print("     coefficient = ", '{:.2uS}'.format(val))
    print("     chi^2 = ", np.around(table_fits['chi2'],2))

for i in clim['N']:
    table_fits = np.compress(clim['N']==i, clim)
    val = ufloat(table_fits['chi'], table_fits['chi_err'])
    c2 = np.around( table_fits['chi2'],2)[0]
    print('${:d}$ & ${:.2uS}$ & ${:f}$ \\\\'.format( i, val, c2 ))

#plt.legend(loc=1, ncol=2, frameon=False)
plt.legend(bbox_to_anchor=(0.,1.0), loc='lower left', ncol=2, frameon=False)
plt.tight_layout()
plt.savefig('SPN_Topology_contlim_'+str(t0)+'_'+str(w0)+'_'+labf+'_w0.pdf')
plt.show()

clim = np.empty(0, dtype=dtclim)
for i in np.unique(chi_SPN['N']):
    pl_data = np.compress( chi_SPN['N'] == i, chi_SPN)
    marker1 = next(marks)

    xdata = pl_data['sqrtS']**2
    ydata = pl_data['chi_t']/pl_data['sqrtS']**4/pl_data['L']**4
    #if( np.max(ydata) > max_y_lim):
    #    max_y_lim = 1.4*np.max(ydata)
    #plt.ylim(top=max_y_lim)
    #if( np.min(ydata) < min_y_lim):
    #    min_y_lim = 0.45*np.max(ydata)
    #plt.ylim(min_y_lim,max_y_lim)
    ydata_err = np.sqrt( (4.*pl_data['chi_t']/pl_data['sqrtS']**5*pl_data['sqrtS_err'])**2
                        + (pl_data['schi_t']/pl_data['sqrtS']**4)**2)/pl_data['L']**4
    #for l in range(len(xdata)):
    #    val = ufloat(ydata[l],ydata_err[l])
    #    uval_sqrtS = ufloat(pl_data['sqrtS'][l]*pl_data['wscale'][l]**2, 
    #            np.sqrt( (pl_data['sqrtS_err'][l]*pl_data['wscale'][l]**2)**2 + (pl_data['sqrtS'][l]*2.*pl_data['wscale'][l]*pl_data['wscale_err'][l])**2 ))
    #    print(" & $", '{:.2uS}'.format(uval_sqrtS), "$ & $", '{:.2uS}'.format(val), '$ \\\\')

    popt, pcov = curve_fit(f, xdata, ydata, sigma=ydata_err, absolute_sigma=True)
    chi2 = np.sum( (ydata-f(xdata,*popt))**2/ydata_err**2)/( len(ydata)-len(popt))
    perr = np.sqrt(np.diag(pcov))
    lab = r'$N_c= '+str(int(i))+'$, $\\tilde{\mathcal{X}}^2='+str(np.around(chi2,2))+'$'

    color1 = next(plt.gca()._get_lines.prop_cycler)['color']
    xr = np.arange(0., np.max(xdata)*1.15,0.01)
    clim_tmp = np.array([(i,popt[0],perr[0], popt[1], perr[1], chi2)], dtype=dtclim)
    clim = np.append( clim, clim_tmp)
    plt.plot(xr, f(xr,*popt), color=color1, linestyle='dashed')
    #, label=r'$\chi^2/DOF ='+str(np.around(chi2,2))+'$')
    plt.errorbar( 0., popt[0], yerr=perr[0], linestyle='None', marker=marker1, color=color1)
    plt.errorbar(xdata , ydata, yerr=ydata_err, linestyle='None', marker=marker1, color=color1)
    plt.errorbar([np.nan], [np.nan], yerr=[np.nan], color=color1, label=lab, linestyle='dashed')
    #plt.text(0.5, 0.5, 'PRELIMINARY', fontsize=40, color='gray', alpha=0.5,ha='center', va='center', rotation='30')
plt.show()

for i in clim['N']:
    table_fits = np.compress(clim['N']==i, clim)
    print("N = ", i)
    val = ufloat(table_fits['chi'], table_fits['chi_err'])
    print("     chi t_0^2(a=0 ) = ", '{:.2uS}'.format(val))
    val = ufloat(table_fits['sl'], table_fits['sl_err'])
    print("     coefficient = ", '{:.2uS}'.format(val))
    print("     chi^2 = ", np.around(table_fits['chi2'],2))

for i in clim['N']:
    table_fits = np.compress(clim['N']==i, clim)
    val = ufloat(table_fits['chi'], table_fits['chi_err'])
    c2 = np.around( table_fits['chi2'],2)[0]
    print('${:d}$ & ${:.2uS}$ & ${:f}$ \\\\'.format( i, val, c2 ))
#plt.figure()
#plt.xlabel(r'$1/N_c^2$')
#plt.ylabel(r'$\chi_L t_0^2$')
#plt.ylim(0.00, 0.005)
#color1 = next(plt.gca()._get_lines.prop_cycler)['color']
#xr = np.arange(0., 0.3,0.01)
#
#xdata = 1/clim['N']**2
#ydata = clim['chi']
#ydata_err = clim['chi_err']
#plt.errorbar(xdata , ydata, yerr=ydata_err, linestyle='None', marker='^', color=color1)
#
#data=np.compress( clim['N']>1, clim)
#xdata = 1/data['N']**2
#ydata = data['chi']
#ydata_err = data['chi_err']
#popt, pcov = curve_fit(f, xdata, ydata, sigma=ydata_err, absolute_sigma=True)
#chi2 = np.sum( (ydata-f(xdata,*popt))**2/ydata_err**2)/( len(ydata)-2.)
#perr = np.sqrt(np.diag(pcov))[0]
#color1 = next(plt.gca()._get_lines.prop_cycler)['color']
#lab=r'$N_c>0$, $\chi^2/N_\mathrm{d.o.f.} ='+str(np.around(chi2,2))+'$'
#plt.plot(xr, f(xr,*popt), linestyle='dashed',color=color1, label=lab)
#plt.errorbar( 0., popt[0], yerr=perr, linestyle='None', marker='o', color=color1)
#print("a+b/N")
#val = ufloat(popt[0], perr)
#print("     N>1, chi t_0^2 (a=0, Nc=infty) = ", '{:.2uS}'.format(val), "chi2 = ", np.around(chi2,3))
#
#data=np.compress( clim['N']>2, clim)
#xdata = 1/data['N']**2
#ydata = data['chi']
#ydata_err = data['chi_err']
#popt, pcov = curve_fit(f, xdata, ydata, sigma=ydata_err, absolute_sigma=True)
#chi2 = np.sum( (ydata-f(xdata,*popt))**2/ydata_err**2)/( len(ydata)-2.)
#perr = np.sqrt(np.diag(pcov))[0]
#color1 = next(plt.gca()._get_lines.prop_cycler)['color']
#lab=r'$N_c>2$, $\chi^2/N_\mathrm{d.o.f.} ='+str(np.around(chi2,2))+'$'
#plt.plot(xr, f(xr,*popt), linestyle='dashed',color=color1, label=lab)
#plt.errorbar( 0., popt[0], yerr=perr, linestyle='None', marker='o', color=color1)
#print("a+b/N")
#val = ufloat(popt[0], perr)
#print("     N>2, chi t_0^2 (a=0, Nc=infty) = ", '{:.2uS}'.format(val), "chi2 = ", np.around(chi2,3))
#
#plt.errorbar( 0.0, 0.000703, yerr=0.000013, linestyle='None', marker='*', color='black', label=r'$ SU(\infty)$')
#plt.legend(loc=2)
#plt.tight_layout()
#plt.savefig('SPN_Topology_largeN-N_'+str(t0)+'_'+str(w0)+'_'+labf+'.pdf')
##plt.show()
#
#plt.figure()
#plt.xlabel(r'$1/N_c**2$')
#plt.ylabel(r'$\chi t_0^2$')
#plt.ylim(0.00, 0.005)
#color1 = next(plt.gca()._get_lines.prop_cycler)['color']
#xr = np.arange(0., 0.3,0.01)
#
#xdata = 1/clim['N']**2
#ydata = clim['chi']
#ydata_err = clim['chi_err']
#plt.errorbar(xdata , ydata, yerr=ydata_err, linestyle='None', marker='^', color=color1)
#
#data=np.compress( clim['N']>2, clim)
#xdata = 1/data['N']**2
#ydata = data['chi']
#ydata_err = data['chi_err']
#popt, pcov = curve_fit(f, xdata, ydata, sigma=ydata_err, absolute_sigma=True)
#chi2 = np.sum( (ydata-f(xdata,*popt))**2/ydata_err**2)/( len(ydata)-2.)
#perr = np.sqrt(np.diag(pcov))[0]
#color1 = next(plt.gca()._get_lines.prop_cycler)['color']
#lab=r'$a+b/N_c^2$, $N_c>1$, $\chi^2/DOF ='+str(np.around(chi2,2))+'$'
#plt.plot(xr, f(xr,*popt), linestyle='dotted',color=color1, label=lab)
#plt.errorbar( 0., popt[0], yerr=perr, linestyle='None', marker='o', color=color1)
#print("a+b/N^2")
#val = ufloat(popt[0], perr)
#print("     N>2, chi t_0^2 (a=0, Nc=infty) = ", '{:.2uS}'.format(val), "chi2 = ", np.around(chi2,3))
#
#data=np.compress( clim['N']>1, clim)
#xdata = 1/data['N']**2
#ydata = data['chi']
#ydata_err = data['chi_err']
#popt, pcov = curve_fit(f, xdata, ydata, sigma=ydata_err, absolute_sigma=True)
#chi2 = np.sum( (ydata-f(xdata,*popt))**2/ydata_err**2)/( len(ydata)-2.)
#perr = np.sqrt(np.diag(pcov))[0]
#color1 = next(plt.gca()._get_lines.prop_cycler)['color']
#lab=r'$a+b/N_c^2$, $N_c>0$, $\chi^2/DOF ='+str(np.around(chi2,2))+'$'
#plt.plot(xr, f(xr,*popt), linestyle='dotted',color=color1, label=lab)
#plt.errorbar( 0., popt[0], yerr=perr, linestyle='None', marker='o', color=color1)
#print("a+b/N^2")
#val = ufloat(popt[0], perr)
#print("     N>1, chi t_0^2 (a=0, Nc=infty) = ", '{:.2uS}'.format(val), "chi2 = ", np.around(chi2,3))
#
##data=np.compress( clim['N']>1, clim)
##xdata = 1/data['N']**2
##ydata = data['chi']
##ydata_err = data['chi_err']
##xr = np.arange(0., 0.6,0.01)
##popt, pcov = curve_fit(g, xdata, ydata, sigma=ydata_err, absolute_sigma=True)
##chi2 = np.sum( (ydata-g(xdata,*popt))**2/ydata_err**2)/( len(ydata)-3.)
##perr = np.sqrt(np.diag(pcov))[0]
##color1 = next(plt.gca()._get_lines.prop_cycler)['color']
##lab=r'$a+b/N_c+c/N_c^2$, $N_c>0$, $\chi^2/DOF ='+str(np.around(chi2,2))+'$'
##plt.plot(xr, g(xr,*popt), linestyle='dotted',color=color1, label=lab)
##plt.errorbar( 0., popt[0], yerr=perr, linestyle='None', marker='o', color=color1)
##print("a + b/N + c/N^2")
##val = ufloat(popt[0], perr)
##print("     N>1, chi t_0^2 (a=0, Nc=infty) = ", '{:.2uS}'.format(val), "chi2 = ", np.around(chi2,3))
##
#plt.errorbar( 0.0, 0.000703, yerr=0.000013, linestyle='None', marker='*', color='black', label=r'$ SU(\infty)$')
#
#plt.legend(loc=2)
#plt.tight_layout()
##plt.savefig('SPN_Topology_largeN-N2_'+str(t0)+'_'+str(w0)+'_'+labf+'.pdf')
##plt.show()
