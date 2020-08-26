"""

Functions to analyze VC test pulses. Functions to measure:

"""

from brian2 import * 
from scipy.optimize import curve_fit

def find_peak_current(t, I, t_start, t_end, dt):
    '''
    Returne the minimum current in the trace between t_start and t_end
    '''
    peak_idx = where((t>t_start)&(t<t_end))
    tr_time = argmin(I[peak_idx[0]])
    I_peak = I[peak_idx[0]][tr_time]
    t_peak = peak_idx[0][tr_time]*dt
    return I_peak, t_peak

def baseline_current(t, I, t_start, t_end):
    base = where((t>t_start)&(t<t_end))
    I_baseline = mean(I[base[0]])
    return I_baseline

def baseline_potential(t, V, t_start, t_end):
    base = where((t>t_start)&(t<t_end))
    V_baseline = mean(V[base[0]])
    return V_baseline

def fit_current_decay(t, I, t_peak, t_end, I_amp_peak, I_plat, dt, plotting = True):
    '''
    Fits an exponential to a capacitive transient.
    Output: tau (RC time constant), peak current
    '''
    def exp_current(t, tau):
        return abs(I_plat - I_amp_peak * exp(-t/tau))
    
    peak_current = abs(I[int(t_peak/dt):int(t_end/dt)]) # pA
    peak_time = t[int(t_peak/dt):int(t_end/dt)]- t[int(t_peak/dt)]  # ms
    #print peak_current
    #print peak_time
    # Tau
    tau_opt = curve_fit(exp_current, peak_time, peak_current)     
    # Extraoplation
    neg_time = [-int(dt)] 
    time_extrapol = hstack((neg_time, peak_time))
    peak_current_extrapol = exp_current(time_extrapol, tau_opt[0])
    neg_time_long = [-2*int(dt), -int(dt)]
    time_full_peak = hstack((neg_time_long, peak_time))
    
    i_peak_extra = max(peak_current_extrapol)
    
    if plotting is True:
        figure('fit')
        plot(peak_time, -peak_current, '-', color='k', label='current')
        plot(peak_time, -exp_current(peak_time, tau_opt[0]), 'r', label='fit tau')
        ylabel('I (pA)')
        xlabel('Time (us)')
        #xlim(-50,1000)
        legend(frameon=False)
        
        figure('extrapol')
        plot(time_full_peak, -abs(I[int((t_peak-2*dt)/dt):int(t_end/dt)]), '-o', color='k', label='current')
        plot(time_extrapol, -peak_current_extrapol, '-', color = 'r', label='extrapolation')
        xlabel('Time (us)')
        ylabel('I (pA)')
        #xlim(-50,1000)
        legend(frameon=False)
        
    return tau_opt[0], i_peak_extra
    
def analyse_test_pulse(t, I, V, dt):
    '''
    Analyse a test pulse:
    Parameters:
    t: time in ms
    I: current in pA
    V: membrane potential in mV
    dt : time step 
    Output:
    tau: time constant in second
    I_baseline: leak current in pA
    I_peak_true: absolute peak current in pA
    Rm : membrane resistance in Mohm
    Rs: series resistance in Mohm
    '''
    
    # Amplitude and time of the peak current
    I_peak, t_peak = find_peak_current(t, I, 20.*ms, 30.*ms, dt) # pA 
    print ('Peak time:', t_peak)
    print ('Measured peak current:', I_peak,'pA')
    plot(t_peak/ms, I_peak, 'ro')

    # Baseline current and amplitude of the peak current        
    I_baseline = baseline_current(t, I, 5.*ms, 20.*ms) # pA
    I_amp_peak = abs(I_peak-I_baseline) # pA        
    print ('Leak current at -80 mV:', I_baseline, 'pA')
    
    # Baseline voltage and amplitude of voltage step
    V_baseline = baseline_potential(t, V, 5.*ms, 20.*ms) # mV
    V_step = abs(V[int(t_peak/dt)]-V_baseline) # true voltage step, in mV
    print ('Voltage step:', V_step, 'mV')

    # Membrane resistance measured from the current during the "plateau"
    I_plat = baseline_current(t, I, 22.*ms, 30.*ms) # pA
    I_amp_plat = abs(I_plat-I_baseline) # pA
    Rm = (V_step/I_amp_plat)*1e3 # Mohm        
    print ('Membrane resistance:', Rm, 'MOhm')

    # Fitting and extrpolation to find the time constant and estimate the true peak current e true peak current
    tau, I_peak_true = fit_current_decay(t, I, t_peak, 29.*ms, I_amp_peak, I_plat, dt, plotting = False) 
    Rs = (V_step/I_peak_true)*Gohm
    Cm = tau*second/Rs        
    print ('Tau RC:', tau * second)
    print ('Rs:', Rs)
    print ('Cm:', Cm[0])
    return I_baseline, tau[0], I_peak_true - I_baseline, Rm, Rs/Mohm, Cm[0]/pF
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    