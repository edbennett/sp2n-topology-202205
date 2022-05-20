import sys
import numpy as np
import matplotlib.pyplot as plt
import re
import matplotlib.ticker as ticker
from scipy.interpolate import interp1d,UnivariateSpline, CubicSpline,PchipInterpolator
from matplotlib import gridspec
from matplotlib.ticker import AutoLocator
from uncertainties import ufloat, unumpy
import itertools
import lib_topology as es
import pickle as pkl

for faddr in sys.argv[1:]:
    N,L,beta,rawdata=es.topo_load_raw_data(faddr)
    
    fn_pkl='pkl_flows_bs/pkl_bs_'+N+'_'+L+'_'+beta+'_'

    bs_flow_E, w0_flow_E = es.flows(rawdata, 100, 't2E')

    fn_t_E = open(fn_pkl+'t_E','wb')
    pkl.dump(bs_flow_E, fn_t_E)
    fn_t_E.close()
    
    fn_w_E = open(fn_pkl+'w_E','wb')
    pkl.dump(w0_flow_E, fn_w_E)
    fn_w_E.close()

    bs_flow_symE, w0_flow_symE = es.flows(rawdata, 100, 't2symE')

    fn_t_symE = open(fn_pkl+'t_symE','wb')
    pkl.dump(bs_flow_symE, fn_t_symE)
    fn_t_symE.close()
    
    fn_w_symE = open(fn_pkl+'w_symE','wb')
    pkl.dump(w0_flow_symE, fn_w_symE)
    fn_w_symE.close()
