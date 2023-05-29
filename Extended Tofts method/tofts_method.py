#Packages 
import numpy as np
from scipy.optimize import curve_fit
'''
This script contains two functions, which are need for running ETM.
'''

#Extended Tofts method for analysis of DCE MRI data
def tofts_model(timeline, Kt, ve, vp, cp): 

    nt = len(timeline)
    Ct = np.zeros(nt, dtype='float')
    deltat = timeline[1:]-timeline[0:-1]

    for k in range(1, nt): 
        f1 = cp[0:k]*deltat[0:k]
        f2 = np.exp(-(Kt/ve)*timeline[0:k])
        f = np.convolve(f1, f2, mode='valid')
        w = Kt*f
        Ct[k] = vp*cp[k] + w[-1]
    return Ct

#Run ETM given the concentration of CA in tissue of interest, the timeline and the plasma CA concentration
def run_tofts(timeline, C, cp): 
    
    #Initial values
    Kt = 0.001
    ve = 0.1 
    vp = 0.1

    #curve optimization
    fit_func = lambda timeline, Kt, ve, vp: tofts_model(timeline, Kt, ve,vp, cp)
    popt,pcov,infodict,mes,ier = curve_fit(fit_func, timeline, C, bounds=([-10,0,0], [10,1,1]), full_output=True, maxfev=5000) 

    #Store parameters as dictionary 
    Kt_opt = popt[0]
    ve_opt = popt[1]
    vp_opt = popt[2]
    kep_opt = Kt_opt/ve_opt
    opt_params = {'ktrans':Kt_opt,
                  've': ve_opt,
                  'vp':vp_opt,
                  'kep':kep_opt}
    
    C_model = tofts_model(timeline, Kt_opt, ve_opt, vp_opt, cp)

    return opt_params, C_model
