import numpy as np
import re
from scipy.interpolate import interp1d
#from scipy.interpolate import interp1d,UnivariateSpline, CubicSpline
#
#
#dt_glueball=np.dtype([("N",'i8'),("L",'i8'),("beta",'d'),('sym','S10'),("m",'d'),("sm",'d'),("chi2", 'd')])
#dt_cont_spectrum=np.dtype( [("N",'i8'), ("sym",'S10'), ("m",'d'), ("sm", 'd'), ("chi2", 'd')] )
#dt_cspectrum_syst=np.dtype( [("N",'i8'), ("sym",'S10'), ("ext",'S10'),("m",'d'), ("sm", 'd'), ("chi2", 'd')] )
#
#effstring_approx = 'NG'
#
#scale_of_spectrum='stringtens'
##scale_of_spectrum='ERPpCp'
##scale_of_spectrum='WF'
#
#FSE_torelons = { '6' : 11, '8' : 7 }
#
#lab_channels_cont = {'0RPpCp' : r'$0^+$',
#                    '0RPpCp+' : r'$0^{+*}$' ,
#                    '0RPmCp' : r'$0^-$' ,
#                    '0RPmCp+' : r'$0^{-*}$' ,
#                    '0RAPpCp' : r'$3^+$' ,
#                    '0RAPmCp' : r'$3^-$' ,
#                    'TRPpCp' : r'$2^+$' ,
#                    'TRPmCp' : r'$2^-$' ,
#                    'TRAPpCp' : r'$1^+$' ,
#                    'TRAPmCp' : r'$1^-$' ,
#                    'ERPpCp' : r'$2^+$' ,
#                    'ERPmCp' : r'$2^-$' }
#
#label_channels = { b'0RPpCp' : r'$A_1^+$',
#                    b'0RPpCp+' : r'$A_1^{+*}$' ,
#                    b'0RPmCp' : r'$A_1^-$' ,
#                    b'0RPmCp+' : r'$A_1^{-*}$' ,
#                    b'0RAPpCp' : r'$A_2^+$' ,
#                    b'0RAPmCp' : r'$A_2^-$' ,
#                    b'TRPpCp' : r'$T_2^+$' ,
#                    b'TRPmCp' : r'$T_2^-$' ,
#                    b'TRAPpCp' : r'$T_1^+$' ,
#                    b'TRAPmCp' : r'$T_1^-$' ,
#                    b'ERPpCp' : r'$E^+$' ,
#                    b'ERPmCp' : r'$E^-$' }
#
#chans = [ [ '0RPpCp+', '0RPpCp'  ], 
#          [ '0RPmCp+', '0RPmCp'  ], 
#          [ '0RAPpCp', '0RAPmCp' ],  
#          [ 'TRPpCp' , 'ERPpCp'  ], 
#          [ 'TRPmCp' , 'ERPmCp'  ], 
#          [ 'TRAPpCp', 'TRAPmCp' ] ]
#
#min_largeN = {  '0RPpCp+' : 2, 
#                '0RPpCp'  : 2, 
#                '0RPmCp'  : 2, 
#                '0RPmCp+' : 2, 
#                '0RAPpCp' : 2,
#                '0RAPmCp' : 2, 
#                'TRPpCp'  : 2, 
#                'TRPmCp'  : 2,
#                'TRAPpCp' : 2, 
#                'TRAPmCp' : 2, 
#                'ERPpCp'  : 2,
#                'ERPmCp'  : 2}
#
#beta_min_2 = {  '0RPpCp+' : 2.2, 
#                '0RPpCp' : 2.2, 
#                '0RPmCp' : 2.2, 
#                '0RPmCp+' : 2.2, 
#                '0RAPpCp': 2.2,
#                '0RAPmCp': 2.2, 
#                'TRPpCp' : 2.2, 
#                'TRPmCp' : 2.2,
#                'TRAPpCp': 2.2, 
#                'TRAPmCp': 2.2, 
#                'ERPpCp' : 2.2,
#                'ERPmCp' : 2.2}
#include_2 = {  '0RPpCp' :  [0,1,2,3,4,5,6], 
#                '0RPpCp+' : [1,2,3,4,5,6], 
#                '0RPmCp' : [0,1,2,3,4,5,6], 
#                '0RPmCp+' : [1,2,3,4,5,6], 
#                '0RAPpCp': [1,2,3,4,5,6],
#                '0RAPmCp': [2,3,4,5], 
#                'TRPpCp' : [0,1,2,3,4,5,6], 
#                'TRPmCp' : [1,2,3,4,5,6],
#                'TRAPpCp': [1,2,3,4,5,6], 
#                'TRAPmCp': [3,4,5,6], 
#                'ERPpCp' : [0,1,2,3,4,5,6],
#                'ERPmCp' : [1,2,3,4,5,6]}
#
#beta_min_4 = {  '0RPpCp'  : 7.6, 
#                '0RPpCp+'  : 7.6, 
#                '0RPmCp'  : 7.6, 
#                '0RPmCp+'  : 7.6, 
#                '0RAPpCp' : 7.6,
#                '0RAPmCp' : 7.6, 
#                'TRPpCp'  : 7.6, 
#                'TRPmCp'  : 7.6,
#                'TRAPpCp' : 7.6, 
#                'TRAPmCp' : 7.6, 
#                'ERPpCp'  : 7.6,
#                'ERPmCp'  : 7.6}
#include_4_new = {  '0RPpCp' :  [0,1,2,3,4,5], 
#                '0RPpCp+' : [1,2,3,4,5], 
#                '0RPmCp' : [0,1,2,3,4,5], 
#                '0RPmCp+' : [1,2,3,4,5], 
#                '0RAPpCp': [0,1,2,3,4,5],
#                '0RAPmCp': [1,2,3,4,5], 
#                'TRPpCp' : [1,2,3,4,5], 
#                'TRPmCp' : [1,2,3,4,5],
#                'TRAPpCp': [1,2,3,4,5], 
#                'TRAPmCp': [1,2,3,4,5], 
#                'ERPpCp' : [1,2,3,4,5],
#                'ERPmCp' : [0,1,2,3,4,5]}
#include_4_old = {  '0RPpCp' :  [0,1,2,3,4,5], 
#                '0RPpCp+' : [0,1,2,3,4,5], 
#                '0RPmCp' : [0,1,2,3,4,5], 
#                '0RPmCp+' : [0,1,2,3,4,5], 
#                '0RAPpCp': [0,1,2,3,4,5],
#                '0RAPmCp': [0,1,2,3,4,5], 
#                'TRPpCp' : [0,1,2,3,4,5], 
#                'TRPmCp' : [0,1,2,3,4,5],
#                'TRAPpCp': [0,1,2,3,4,5], 
#                'TRAPmCp': [0,1,2,3,4,5], 
#                'ERPpCp' : [0,1,2,3,4,5],
#                'ERPmCp' : [0,1,2,3,4,5]}
#
#beta_min_6 = {  '0RPpCp'  : 15.6, 
#                '0RPpCp+'  : 15.6, 
#                '0RPmCp'  : 15.6, 
#                '0RPmCp+'  : 15.6, 
#                '0RAPpCp' : 15.6,
#                '0RAPmCp' : 15.6, 
#                'TRPpCp'  : 15.6, 
#                'TRPmCp'  : 15.6,
#                'TRAPpCp' : 15.6, 
#                'TRAPmCp' : 15.6, 
#                'ERPpCp'  : 15.6,
#                'ERPmCp'  : 15.6}
#include_6 = {  '0RPpCp' :  [0,1,2,3,4,5,6,7,8,9], 
#                '0RPpCp+' : [0,1,2,3,4,5,6,7,8,9], 
#                '0RPmCp' : [0,1,2,3,4,5,6,7,8,9], 
#                '0RPmCp+' : [0,1,2,3,4,5,6,7,8,9], 
#                '0RAPpCp': [0,1,2,3,4,5,6,7,8,9],
#                '0RAPmCp': [0,1,2,3,4,5,6,7,8,9], 
#                'TRPpCp' : [0,1,2,3,4,5,6,7,8,9], 
#                'TRPmCp' : [3,4,5,6,7,8,9],
#                'TRAPpCp': [0,1,2,3,4,5,6,7,8,9], 
#                'TRAPmCp': [0,1,2,3,4,5,6,7,8,9], 
#                'ERPpCp' : [0,1,2,3,4,5,6,7,8,9],
#                'ERPmCp' : [3,4,5,6,7,8,9]}
#
#beta_min_8 = {  '0RPpCp'  : 26.6, 
#                '0RPpCp+'  : 26.4, 
#                '0RPmCp'  : 26.4, 
#                '0RPmCp+'  : 26.4, 
#                '0RAPpCp' : 26.4,
#                '0RAPmCp' : 26.4, 
#                'TRPpCp'  : 26.4, 
#                'TRPmCp'  : 26.4,
#                'TRAPpCp' : 26.4, 
#                'TRAPmCp' : 26.4, 
#                'ERPpCp'  : 26.4,
#                'ERPmCp'  : 26.4}
#include_8 = {  '0RPpCp' :  [0,1,2,3,4,5,6,7,8], 
#                '0RPpCp+' : [0,1,2,3,4,5,6,7,8], 
#                '0RPmCp' : [0,1,2,3,4,5,6,7,8], 
#                '0RPmCp+' : [2,3,4,5,6,7,8], 
#                '0RAPpCp': [4,5,6,7,8],
#                '0RAPmCp': [4,5,6,7,8], 
#                'TRPpCp' : [1,2,3,4,5,6,7,8], 
#                'TRPmCp' : [0,1,2,3,4,5,6,7,8],
#                'TRAPpCp': [1,2,3,4,5,6,7,8], 
#                'TRAPmCp': [1,2,3,4,5,6,7,8], 
#                'ERPpCp' : [1,2,3,4,5,6,7,8],
#                'ERPmCp' : [0,1,2,3,4,5,6,7,8]}
#
#beta_min = {'2' : beta_min_2, 
#            '4' : beta_min_4, 
#            '6' : beta_min_6, 
#            '8' : beta_min_8 }
#
#
#def declare_include_dict(tag):
#    if( tag == 'old' ):
#        include_dict = {'2' : include_2,
#                '4' : include_4_old,
#                '6' : include_6,
#                '8' : include_8 }
#    elif( tag == 'new' ):
#        include_dict = {'2' : include_2,
#                '4' : include_4_new,
#                '6' : include_6,
#                '8' : include_8 }
#    return include_dict
#
########################
##Take length of torelon and string tension,
##and computes the mass of the torelon
#######################
#def m0(L,s):
#    return s*1.*L;
#def mLO(L,s):
#    return m0(L,s) - np.pi/(3.*L);
#def mNLO(L,s):
#    return mLO(L,s) - 0.5*(np.pi/(3.*L))**2*1./(s*L);
#def mNG(L,s):
#    return s*L*np.sqrt(1. - 2*np.pi/(3.*s*L**2))
#
#def torelon_mass(x, approx, s):
#    if( approx == b'linear'):
#        return m0(x,s)
#    if( approx == b'LO'):
#        return mLO(x,s)
#    if( approx == b'NLO'):
#        return mNLO(x,s)
#    if( approx == b'NG'):
#        return mNG(x,s)
#    else:
#        print("Approximation unknown.")
#        exit(1)
#
########################
##Take length of torelon and its mass,
##and computes the string tension
#######################
#def scl(L,m):
#    return m/L
#def errscl(L,m):
#    return 1./L
#def slo(L,m):
#    return m/L + np.pi/(3.*L**2)
#def errslo(L,m):
#    return 1./L
#def snlo(L,m):
#    return (0.5/L**2)* ( m*L + np.pi/3 + np.sqrt( (m*L + np.pi/3. )**2 + (2.*np.pi**2/9.) ))
#def errsnlo(L,m):
#    return (0.5/L**2)*( L + ( m*L*(m*L + np.pi/3. ))/np.sqrt( (m*L + np.pi/3. )**2 + (2.*np.pi**2/9.) ))
#def sNG(L,m):
#    return (0.5/L**2)*( 2.*np.pi/3. + np.sqrt((2.*np.pi/3.)**2 + (2.*L*m)**2)) 
#def errsNG(L,m):
#    return ( (2.*m) /  np.sqrt((2.*np.pi/3.)**2 + (2.*L*m)**2.))
#def string_tension(in_data, approx):
#    sigma = np.empty(len(in_data), dtype=[("N", 'i8'),("L",'i8'),("beta",'d'),('sym','S10'),("s",'d'),("s_err",'d')])
#    sigma['N'] = in_data['N']
#    sigma['L'] = in_data['L']
#    sigma['beta'] = in_data['beta']
#    sigma['sym'] = in_data['sym']
#    L=in_data['L']
#    m=in_data['m']
#    sm=in_data['sm']
#    if( approx == 'linear'):
#        sigma['s'] = scl(L,m)
#        sigma['s_err'] = np.sqrt( (errscl(L,m)*sm)**2 )
#    elif( approx == 'LO'):
#        sigma['s'] = slo(L,m)
#        sigma['s_err'] = np.sqrt( (errslo(L,m)*sm)**2 )
#    elif( approx == 'NLO'):
#        sigma['s'] = snlo(L,m)
#        sigma['s_err'] = np.sqrt( (errsnlo(L,m)*sm)**2 )
#    elif( approx == 'NG'):
#        sigma['s'] = sNG(L,m)
#        sigma['s_err'] = errsNG(L,m)*sm
#    else:
#        print("Approximation unknown.")
#        exit(1)
#    return sigma
#
#
########################
##Load data from file and select specific values
##Takes a variable number of arguments with which one can specify which data to select.
##Example: read_from_file( filename, N=6, sym= '0RPpCp') will select the 0++ channel 
##for Sp(2N=6)
#######################
#def read_from_file( fname, **kwargs):
#    rawdata = np.genfromtxt(fname, names=['N','L','beta','sym','m','sm', 'chi2'], usecols=(0,1,2,3,6,7,10), dtype=('i8','i8','d','S10','d','d','d'))
#    for key,value in kwargs.items():
#        rawdata = np.compress( rawdata[ key ] == value, rawdata)
#    return rawdata
#
#def read_tmin( fname, **kwargs):
#    rawdata = np.genfromtxt(fname, names=['N','L','beta','sym','tmin','tmax','m','sm', 'chi2'], usecols=(0,1,2,3,4,5,6,7,10), dtype=('i8','i8','d','S10','S10','S10','d','d','d'))
#    for key,value in kwargs.items():
#        rawdata = np.compress( rawdata[ key ] == value, rawdata)
#    return rawdata
#
#def read_sigmas( fname, **kwargs):
#    rawdata = np.genfromtxt(fname, names=['N','L','beta','sym','m','sm', 'chi2'], usecols=(0,1,2,3,6,7,10), dtype=('i8','i8','d','S10','d','d','d'))
#    rawdata = np.compress( np.logical_and(rawdata['sym'] == 'SigmaS',rawdata['sym'] == 'SigmaT'), rawdata)
#    for key,value in kwargs.items():
#        rawdata = np.compress( rawdata[ key ] == value, rawdata)
#    return rawdata
#
#
##def extract_s(indata, **kwargs):
##    for key,value in kwargs.items():
##        outdata = np.compress( indata[ key ] == value, indata)
##    return outdata
#
#
#def sqrt_string_tension(in_data, channel, approx):
#    print(channel)
#    if( channel == b'SigmaS'):
#        fmass = np.compress( in_data['sym']==b'SigmaS', in_data)
#    elif( channel == b'SigmaT' ):
#        fmass = np.compress( in_data['sym']==b'SigmaT', in_data)
#    elif( channel == b'both' ):
#        fmass = np.compress( np.logical_or( in_data['sym'] == b'SigmaT', in_data['sym'] == b'SigmaS'), in_data)
#    else:
#        print("Channel unknown, error.")
#        exit(1)
#
#    tmps = string_tension(fmass, approx)
#
#    if( channel == b'SigmaS' or channel == b'SigmaT' ):
#        sqrs = np.empty(len(fmass), dtype=[("L",'i8'),("beta",'d'),("sqrs",'d'),("s_sqrs",'d')])
#        sqrs['L'] = fmass['L']
#        sqrs['beta'] = fmass['beta']
#        sqrs['sqrs'] = np.sqrt(tmps['s'])
#        sqrs['s_sqrs'] = 1./(2.*np.sqrt(tmps['s']))*tmps['s_err']
#    elif( channel == b'both' ):
#        sqrs = np.empty(int(len(fmass)/2), dtype=[("L",'i8'),("beta",'d'),("sqrs",'d'),("s_sqrs",'d')])
#
#        sqrs['L'] = fmass[ fmass['sym'] == b'SigmaS']['L']
#        sqrs['beta'] = fmass[ fmass['sym'] == b'SigmaS']['beta']
#        tmps1 = np.compress( tmps['sym'] == b'SigmaS', tmps)
#        sqrs_s = np.sqrt( tmps1['s'])
#        ssqrs_s = 1/(2.*np.sqrt(tmps1['s']))*tmps1['s_err']
#        tmps1 = np.compress( tmps['sym'] == b'SigmaT', tmps)
#        sqrs_t = np.sqrt( tmps1['s'])
#        ssqrs_t = 1/(2.*np.sqrt(tmps1['s']))*tmps1['s_err']
#        sqrs['sqrs'] = (sqrs_s/ssqrs_s**2 + sqrs_t/ssqrs_t**2)/(1./ssqrs_s**2 + 1./ssqrs_t**2)
#        sqrs['s_sqrs'] =  1./np.sqrt(1./ssqrs_s**2 + 1./ssqrs_t**2)
#    return sqrs
#        
#def weight_string_tension( in_data, approx):
#    fmass = np.compress( np.logical_or( in_data['sym']==b'SigmaT', in_data['sym'] == b'SigmaS'), in_data)
#    tmps_g = string_tension(fmass, approx)
#
#    sqrs = np.empty(0, dtype=[("N", 'i8'),("L",'i8'),("beta",'d'),('sym','S10'),("s",'d'),("s_err",'d')])
#    tmp_iter = zip(fmass['L'], fmass['beta'])
#    for i in tmp_iter:
#        tmps_g = np.compress( np.logical_and( fmass['L'] == i[0], fmass['beta'] == i[1]), tmps)
#
#        tmps_g['beta'] = i[0]
#        tmps_g['L'] = i[1]
#
#        tmps1 = np.compress( tmps_g['sym'] == b'SigmaS', tmps_g)
#        s_s = tmps1['s']
#        ss_s = tmps1['s_err']
#        tmps2 = np.compress( tmps_g['sym'] == b'SigmaT', tmps_g)
#        s_t = tmps2['s']
#        ss_t = tmps2['s_err']
#
#        sqrs['s'] = (s_s/ss_s**2 + s_t/ss_t**2)/(1./ss_s**2 + 1./ss_t**2)
#        sqrs['s_err'] =  1./np.sqrt((1./ss_s**2 + 1./ss_t**2))
#        sqrs = np.append(sqrs, tmps_g)
#    return sqrs
#
#
#def scale_spectrum(in_data, scale_name):
#
#    gmass = np.compress( np.logical_and( in_data['sym'] != b'SigmaS', in_data['sym'] != b'SigmaT'),in_data)
#    
#    if( scale_name == 'stringtens'):
#        scale = sqrt_string_tension(in_data, b'S', effstring_approx)
#        tmp_iter = zip( scale['L'], scale['beta'] )
#        tmp_final = np.empty(0, dtype=dt_glueball)
#        for i in tmp_iter:
#            tmp_gmass = np.compress( np.logical_and( gmass['L'] == i[0], gmass['beta'] == i[1]), gmass) 
#            tmp_scale = np.compress( np.logical_and( scale['L'] == i[0], scale['beta'] == i[1]), scale) 
#
#            tmp_gmass['sm'] =  np.sqrt( ( tmp_gmass['sm']/tmp_scale['sqrs']  )**2 + ( tmp_scale['s_sqrs']*tmp_gmass['m']/tmp_scale['sqrs']**2  )**2 )
#            tmp_gmass['m'] = tmp_gmass['m']/tmp_scale['sqrs']
#            tmp_final = np.append(tmp_final, tmp_gmass)
#        return tmp_final
#    elif( scale_name == 'ERPpCp' ):
#        scale = np.compress( gmass['sym'] == 'ERPpCp' , gmass)
#        tmp_iter = zip( scale['L'], scale['beta'] )
#        tmp_final = np.empty(0, dtype=dt_glueball)
#        for i in tmp_iter:
#            tmp_gmass = np.compress( np.logical_and( gmass['L'] == i[0],gmass['beta'] == i[1]), gmass) 
#            tmp_scale = np.compress( np.logical_and( scale['L'] == i[0],scale['beta'] == i[1]), scale) 
#
#            tmp_gmass['sm'] =  np.sqrt( ( tmp_gmass['sm']/tmp_scale['m']  )**2 + ( tmp_scale['sm']*tmp_gmass['m']/tmp_scale['m']**2  )**2 )
#            tmp_gmass['m'] = tmp_gmass['m']/tmp_scale['m']
#            tmp_final = np.append(tmp_final, tmp_gmass)
#        return tmp_final
#    elif( scale_name == 'WF' ):
#        scale = np.genfromtxt('WF_w0.dat', names=['N','beta','t','st','tsym', 'stsym'], usecols=(0,1,2,3,4,5), dtype=('i8','d','d','d','d','d'))
#        print("scale file:")
#        #print scale
#        tmp_iter = scale['beta']
#        tmp_final = np.empty(0, dtype=dt_glueball)
#        #print tmp_iter
#        for i in tmp_iter:
#            #print i
#            tmp_gmass = np.compress( in_data['beta'] == i, in_data) 
#            tmp_scale_tmp = np.compress( scale['beta'] == i, scale) 
#            #tmp_scale = np.sqrt( 8.* tmp_scale_tmp['tsym'] )
#            tmp_scale = tmp_scale_tmp['tsym'] 
#            #tmp_scale_err = 4.*tmp_scale_tmp['stsym']/tmp_scale
#            tmp_scale_err = tmp_scale_tmp['stsym']
#            #tmp_gmass['sm'] =  np.sqrt( ( tmp_gmass['sm']*tmp_scale  )**2 + ( tmp_scale_err*tmp_gmass['m'] )**2 )
#            tmp_gmass['m'] = tmp_gmass['m']*tmp_scale
#            tmp_gmass['sm'] =  np.sqrt( ( tmp_gmass['sm']*tmp_scale  )**2 + ( tmp_scale_err*tmp_gmass['m'] )**2 )
#            tmp_final = np.append(tmp_final, tmp_gmass)
#        return tmp_final
#    else:
#        print("Scale not understood.")
#        exit(1)
#
#def FSE_data( in_data, channel):
#    if(channel == 'glue'):
#        data = glueball_ratio_extract(in_data, 'stringtens')
#        #print data
#        data_out = np.compress( data['sym'] == 'A1P', data)
#        return data_out
#    elif(channel == 'stringS'):
#        data1 = extract_s(in_data,'S','LO')
#        data2 = extract_s(in_data,'S','NLO')
#        data3 = extract_s(in_data,'S','NG')
#        return data1,data2,data3
#    elif(channel == 'stringT'):
#        data1 = extract_s(in_data,'T','LO')
#        data2 = extract_s(in_data,'T','NLO')
#        data3 = extract_s(in_data,'T','NG')
#        return data1,data2,data3
#    elif(channel == 'string'):
#        data1 = extract_s(in_data,'both','LO')
#        data2 = extract_s(in_data,'both','NLO')
#        data3 = extract_s(in_data,'both','NG')
#        return data1,data2,data3
#    else:
#        print("Error, exiting.")
#        exit(1)
#
#def cas_scal(group, N):
#    if( group == 'Sp2N'):
#        return 4.*(1.*N/2.+1)/(1.*(N+1))
#    if( group == 'SUN'):
#        return 2.*N**2/(1.*(N**2-1.))
#    else:
#        print("Error!")
#        exit(1)
#
#def autocorr(indata):
#
#    L=len(indata)
#    if(L==0):
#        print("Empty data set")
#        autocorr=0.
#        err=0.
#        return autocorr, err
#    avr=0.
#    f2=0.
#    C0=0.
#    for j in range(L):
#        avr = avr + indata[j]
#        f2 = f2 + indata[j]**2
#    f = 2.*avr
#    avr = avr / (1.*L)
#    C0  = f2 / (1.*L) - avr**2
#    tint=0.5
#    val=False
#    for M in range(1,L):
#        avr=0.
#        f=0.
#        f2=0.
#        for i in range(L-M):
#            f2 = f2 + indata[i] * indata[i+M]
#            f = f + indata[i] + indata[i+M]
#            avr = avr + indata[i]
#        avr = avr / (1.*(L-M))
#        Ct  = (f2/(1.*(L-M))) + avr*( avr - f/(1.*(L-M)))
#        rho = Ct/C0
#        tint = tint + rho
#        if( 4.*tint < M):
#            val=True
#            break
#    if(val):
#        autocorr=tint
#        err=np.sqrt( 2.*(2.*M+1)/(1.*L))*tint
#    else:
#        autocorr=0.
#        err=0.
#    return autocorr, err
#
#
def topo_find_alpha(din):
    a_min=0.4
    delta = np.average( ( a_min*din - np.rint(a_min*din) )**2 )
    for ia in np.arange(a_min, 2.0, 0.000025):
        delta_tmp = np.average( ( ia*din - np.rint(ia*din) )**2 )
    #    print(ia, delta_tmp)
        if( delta_tmp < delta):
            a_min = ia
            delta = delta_tmp
    return a_min
#
#
#def load_plaq_data(fname):
#    patt = re.compile(r"WF_([0-9])_([0-9]+)_([0-9]+.[0-9]+)_plaq")
#    N, L, beta, = patt.search(fname).groups()
#    out_data = np.genfromtxt(fname)
#    return N, L, beta, out_data
#
def topo_load_raw_data(fname):
    patt = re.compile(r"WF_([0-9])_([0-9]+)_([0-9]+.[0-9]+)")
    N, L, beta, = patt.search(fname).groups()
    out_data = np.genfromtxt(fname, usecols=(0,1,2,3,4,5,6),
            dtype=[('nconf','int'),('t','f8'),('E','f8'),('t2E','f8'),('Esym','f8'),('t2symE','f8'),('TC','f8')])
    return N, L, beta, out_data
#
#def topo_find_tmax(data_in, L):
#    t_tmp = np.unique(data_in['t'])
#    bool_t_max = 2.*np.sqrt(8*t_tmp) >= int(L)
#    tmax=np.min(t_tmp[bool_t_max])
#    return tmax
#
#def topo_extract_TC_susc(fname):
#    N, L, beta, rawdata = topo_load_raw_data(fname)
#    t_max = topo_find_tmax(rawdata, L)
#    #print(t_max)
#    TC = np.compress( rawdata['t'] == t_max, rawdata['TC'] )
#    a_min = topo_find_alpha(TC)
#    TCdata = np.zeros(len(TC) , dtype=[("nconf",'i8'),("TC",'f8')])
#    TCdata['nconf'] = range(1,len(TC)+1)
#    TCdata['TC'] = np.rint(a_min*TC)
#    return N, L, beta, TCdata
#
def bs_avg_err_TC(din, bin_size=20):
    bin1 = []
    bin2 = []
    for i in range(int(len(din)/bin_size)):
        bin1.append(np.average(din[ i*bin_size : (i+1)*bin_size ]))
        bin2.append(np.average(din[ i*bin_size : (i+1)*bin_size ]**2))

    resampled1 = []
    resampled2 = []
    for j in range(100):
        sam = np.random.randint(0,len(bin1),size=len(bin1))
        avg = np.average( [ bin1[i] for i in sam ] )
        avg2 = np.average( [ bin2[i] for i in sam ])
        resampled1.append( avg )
        resampled2.append( avg2 -avg**2)
        #resampled2.append( avg2)
    #print(N, L, beta)
    return np.average(resampled2), np.std(resampled2)
#
def bstrap( a, nbstrap):
    stat1 = np.zeros(nbstrap)
    stat2 = np.zeros(nbstrap)
    lena = len(a)
    for i in range(nbstrap):
        atmp = np.random.choice( a, lena, replace=True)
        stat1[i] = np.average(atmp)
        stat2[i] = np.std(atmp, ddof=1)
    s1= np.average(stat1)
    e1= np.std(stat1)
    s2= np.average(stat2)
    e2= np.std(stat2)
    return s1, e1, s2, e2
#
#def max_E0(rawdata, obs):
#    t0 = []
#    tmax = np.max( np.unique( rawdata['t'] ) )
#    data = np.compress( rawdata['t'] == tmax , rawdata)
#    emax = np.min( data[obs] )
#    return emax
#
def find_t0(indata, TE):
    t0 = []
    #try:
    for i in np.unique( indata['bs'] ):
        data = np.compress( indata['bs'] == i , indata)
        idx = np.abs( data['flow']-TE ).argmin()
        traj = data[idx-5:idx+5]
        f = interp1d( traj['flow'], traj['t'])
        t0 = np.append(t0,f(TE))
    a,b,c,d =bstrap( t0, 100)
    #except:
    #    a=np.nan
    #    c=np.nan
    return a,c
#
def find_w0(indata, TE):
    t0 = []
    #try:
    for i in np.unique( indata['bs'] ):
        data = np.compress( indata['bs'] == i , indata)
        idx = np.abs( data['flow']-TE ).argmin()
        traj = data[idx-5:idx+5]
        f = interp1d( traj['flow'], traj['t'])
        t0 = np.append(t0,np.sqrt(f(TE)))
    a,b,c,d =bstrap( t0, 100)
    #except:
    #    a=np.nan
    #    c=np.nan
    return a,c
##def find_t0(indata, TE):
##    t0 = []
##    try:
##        for i in np.unique( indata['bs'] ):
##            data = np.compress( indata['bs'] == i , indata)
##            idx = np.abs( data['flow']-TE ).argmin()
##            traj = data[idx-5:idx+5]
##            f = interp1d( traj['flow'], traj['t'])
##            t0 = np.append(t0,f(TE))
##        a,b,c,d =bstrap( t0, 60)
##    except:
##        a=np.nan
##        c=np.nan
##    return a,c
##
##
##def find_w0(indata, TE):
##    t0 = []
##    try:
##        for i in np.unique( indata['bs'] ):
##            data = np.compress( indata['bs'] == i , indata)
##            idx = np.abs( data['flow']-TE ).argmin()
##            traj = data[idx-5:idx+5]
##            f = interp1d( traj['flow'], traj['t'])
##            t0 = np.append(t0, np.sqrt(f(TE)))
##        a,b,c,d =bstrap( t0, 60)
##    except:
##        a=np.nan
##        c=np.nan
##    return a,c
#
rng = np.random.default_rng()
#
def flows(rawdata, N_bs, obs, Ntherm=500):
    Nconf = np.max(np.unique(rawdata['nconf']))
    rawdata['t'] = np.around(rawdata['t'],8)
    ts = np.unique(rawdata['t'])
    delta_t = ts[1]-ts[0]
    bs_flow = np.empty(N_bs*len(ts), dtype=[('bs','i'),('t','f8'),('flow','f8')])
    w0_flow = np.empty(N_bs*(len(ts)-1), dtype=[('bs','i'),('t','f8'),('flow','f8')])
    tmp_bs_flow = np.empty(Nconf*len(ts), dtype=[('nconf','int'),('t','f8'),('E','f8'),('t2E','f8'),('Esym','f8'),('t2symE','f8'),('TC','f8')])
    for i_bs in range(N_bs):
        #print(i_bs)
        sam = rng.integers(Ntherm,Nconf+1,size=Nconf)
        i_c=0
        for s in sam:
            tmp = np.compress(rawdata['nconf'] == s, rawdata)
            for l_t in range(len(ts)):
                try:
                    tmp_bs_flow[i_c*len(ts)+l_t] = tmp[l_t]
                except:
                    print("error at :", s, l_t)
                #print(i_c, l_t, tmp_bs_flow[i_c*len(ts)+l_t], tmp[l_t])
            i_c=i_c+1
        #tmp_bs_flow = np.append(tmp_bs_flow, tmp)
        #print(tmp_bs_flow)
        for i_t in range(len(ts)):
            #print(i_t, ts[i_t])
            t0_tmp = np.around(ts[i_t],8)
            f0tmp = np.compress( tmp_bs_flow['t'] == t0_tmp, tmp_bs_flow)
            f0_avg = np.average( f0tmp[obs] )
            bs_fl_tmp = np.array([(i_bs, t0_tmp, f0_avg)], dtype=[('bs','i'),('t','f8'),('flow','f8')])
            #print(delta_t)
            #bs_flow = np.append( bs_flow, bs_fl_tmp)
            bs_flow[i_bs*len(ts) + i_t] = bs_fl_tmp
            #print(i_bs, i_t, bs_flow[i_bs*len(ts)+i_t], bs_fl_tmp)
            if( i_t < len(ts)-1):
                t1_tmp = np.around(ts[i_t+1],8)
                f1tmp = np.compress( tmp_bs_flow['t'] == t1_tmp, tmp_bs_flow)
                f1_avg = np.average( f1tmp[obs] )
                delta_t = np.around(t1_tmp - t0_tmp,8)
                w0_fl_tmp = np.array([(i_bs, t0_tmp, t0_tmp*(f1_avg-f0_avg)/delta_t)], dtype=[('bs','i'),('t','f8'),('flow','f8')])
                #w0_flow = np.append( w0_flow, w0_fl_tmp)
                w0_flow[i_bs*(len(ts)-1) + i_t] =  w0_fl_tmp
                #print(i_bs, i_t, w0_flow[i_bs*len(ts)+i_t], w0_fl_tmp)
    return bs_flow, w0_flow
##def flows(rawdata, N_bs, obs):
##    Nconf = np.max(np.unique(rawdata['nconf']))
##    rawdata['t'] = np.around(rawdata['t'],8)
##    ts = np.unique(rawdata['t'])
##    delta_t = ts[1]-ts[0]
##    bs_flow = np.empty(0, dtype=[('bs','i'),('t','f8'),('flow','f8')])
##    w0_flow = np.empty(0, dtype=[('bs','i'),('t','f8'),('flow','f8')])
##    for i_bs in range(N_bs):
##        tmp_bs_flow = np.empty(0, dtype=[('nconf','int'),('t','f8'),('E','f8'),('t2E','f8'),('Esym','f8'),('t2symE','f8'),('TC','f8')])
##        sam = np.random.randint(1,Nconf+1,size=Nconf)
##        for s in sam:
##            tmp = np.compress(rawdata['nconf'] == s, rawdata)
##            tmp_bs_flow = np.append(tmp_bs_flow, tmp)
##        for i_t in range(len(ts)-1):
##            t0_tmp = np.around(ts[i_t],8)
##            f0tmp = np.compress( tmp_bs_flow['t'] == t0_tmp, tmp_bs_flow)
##            f0_avg = np.average( f0tmp[obs] )
##            t1_tmp = np.around(ts[i_t+1],8)
##            f1tmp = np.compress( tmp_bs_flow['t'] == t1_tmp, tmp_bs_flow)
##            f1_avg = np.average( f1tmp[obs] )
##            delta_t = t1_tmp - t0_tmp
##            bs_fl_tmp = np.array([(i_bs, t0_tmp, f0_avg)], dtype=[('bs','i'),('t','f8'),('flow','f8')])
##            bs_flow = np.append( bs_flow, bs_fl_tmp)
##            w0_fl_tmp = np.array([(i_bs, t0_tmp, t0_tmp*(f1_avg-f0_avg)/delta_t)], dtype=[('bs','i'),('t','f8'),('flow','f8')])
##            w0_flow = np.append( w0_flow, w0_fl_tmp)
##    return bs_flow, w0_flow
#
#def Casimir_SP(N):
#    return (float(N)+1.)/4.
#
#
def avg_flows(in_flow):
    ts = np.unique(in_flow['t'])
    Lts = len(ts)
    out_flow = np.empty(Lts, dtype=[('t','f8'),('flow','f8'),('err','f8')])
    for i_t in range(Lts):
        out_flow['t'][i_t] = ts[i_t]
        toavg = np.compress(in_flow['t']==ts[i_t], in_flow['flow'])
        out_flow['flow'][i_t] = np.average(toavg)
        out_flow['err'][i_t] = np.std(toavg)
    return out_flow

def Casimir_SP(N):
    return (float(N)+1.)/4.

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

def d_g_SPN(n):
	return n*(2.*n+1.)
