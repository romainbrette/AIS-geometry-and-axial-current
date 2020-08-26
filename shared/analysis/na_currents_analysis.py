"""

Functions to analyze Na currents recordings at the soma.

"""

from brian2 import * 
from scipy import interpolate
from scipy import signal 
import pyabf
import sys
sys.path.append("/users/sarahgoethals/Documents/repositories/AIS-geometry-and-axonal-Na-current/shared/analysis/")
from vc_test_pulse_analysis import find_peak_current

def remove_passive_component(date, retina, cell, dt, i_data, v_data, dt_tp, i_tp_data, v_tp_data, rec_name, idx_step = None):
    """
    Returns the cut raw traces and the cut traces corrected for the passive component.
    dt: time step for NA currents recording
    i_data: Na currents
    v_data: command potential for Na currents recordings
    dt_tp: time step for test pulses
    i_tp_data: current test pulse
    v_tp_data: command potential test pulse
    
    """
    n_rec = len(i_data)
    print (n_rec)
    t = dt*arange(len(i_data[0])) 
    t_tp = dt_tp*arange(len(i_tp_data[0])) 
    
    vc_tp_amp = min(v_tp_data[0]) - max(v_tp_data[0]) # amplitude of test pulse
    print (vc_tp_amp)
    idx_pulse = where(v_tp_data[0] == min(v_tp_data[0]))[0][0] - 1 # start of test pulse
    len_pulse = int(20.*ms/dt_tp) 
    print ('Test pulse start:', idx_pulse, idx_pulse*dt_tp, dt_tp)
    if idx_step is None:
        idx_step = where(v_data[-1] == max(v_data[-1]))[0][0] - 1 # start of voltage steps
    print ('VC step start:', idx_step, idx_step*dt, dt)
    
    # mean of the TPs
    i_tp_mean = mean(i_tp_data, axis=0) 
    # cut the test pulse trace
    if dt_tp > dt:
        print ('interpolate TP')
        # interpolate the mean because dt_tp = dt/2
        tck = interpolate.splrep(t_tp, i_tp_mean)
        new_t_tp = dt*arange(len(i_tp_data[0]))  #linspace(min(t_tp), max(t_tp), 2*len(t_tp))
        new_i_tp_mean = interpolate.splev(new_t_tp, tck)
        i_tp_cut = new_i_tp_mean[2*idx_pulse+1:2*idx_pulse + 2*len_pulse+1]
        t_tp_cut = new_t_tp[2*idx_pulse+1:2*idx_pulse + 2*len_pulse+1] - new_t_tp[2*idx_pulse+1]

    else:
        i_tp_cut = i_tp_mean[idx_pulse+1:idx_pulse + len_pulse+1]
        t_tp_cut = t_tp[idx_pulse+1:idx_pulse + len_pulse+1] - t_tp[idx_pulse+1]
    

    I_cut = []
    I_corr_passive = []

    figpass = figure('Passive %i,  %s, %i, %s' %(date, retina, cell, rec_name), (10,7))
    cmap = plt.get_cmap('gnuplot')
    cols = [cmap(i) for i in np.linspace(0, 1, n_rec)]

    #subplot(111)
    for i in range(n_rec):
        V_amp =  v_data[i][idx_step + 10] - v_data[i][idx_step-10]
        factor = V_amp/vc_tp_amp
        #print (V_amp, factor)
        # adjusting the length of the rec current to the length of the TP
        if dt_tp > dt:
            i_cut = i_data[i][idx_step:idx_step + 2*len_pulse]
            t_cut = t[idx_step:idx_step + 2*len_pulse] - t[idx_step]
        else:
            i_cut = i_data[i][idx_step:idx_step + len_pulse]
            t_cut = t[idx_step:idx_step + len_pulse] - t[idx_step]
            
        i_corr =  i_cut - i_tp_cut * factor
        I_cut.append(i_cut)
        I_corr_passive.append(i_corr)
    
    for i in range(n_rec):
        plot(t_cut/ms, I_cut[i], '-', color=cols[i])
        plot(t_cut/ms, I_corr_passive[i], '--', color=cols[i])
    plot(t_cut/ms, i_tp_cut, 'g', label='test pulse average')
    plot(t_cut/ms, i_tp_cut * factor, 'b', label='test pulse average')
    legend(frameon=False)
    xlim(0, 10)

    return I_corr_passive, I_cut, t_cut

def p5_subtraction(date, retina, cell, dt, i_data, v_data, rec_name):
    n_rec = len(i_data)
    t = dt*arange(len(i_data[0])) 
    
    pulse1_start = where(v_data[0] > 0)[0][0] - 1 # start of test pulse #25.6*ms
    pulse1_end =  where(v_data[0][pulse1_start+5:] == 0)[0][0] + pulse1_start+5
    pulse_length = pulse1_end - pulse1_start #round(10.*ms/dt) #where(v_data[0][pulse1_start+1:] ==  0)[0][0] + pulse1_start +1 #pulse1_end - pulse1_start #10.*ms
    pulse_amp = 5.*mV
    pulse2_start = where(v_data[0][pulse1_start+pulse_length + 5:] > 0)[0][0] - 1 + pulse1_start+pulse_length + 5
    pulse_interval = pulse2_start - pulse1_start
    pulse3_start = where(v_data[0][pulse2_start+pulse_length + 5:] > 0)[0][0] - 1 + pulse2_start+pulse_length + 5
    pulse4_start = where(v_data[0][pulse3_start+pulse_length + 5:] > 0)[0][0] - 1 + pulse3_start+pulse_length + 5
    pulse5_start = where(v_data[0][pulse4_start+pulse_length + 5:] > 0)[0][0] - 1 + pulse4_start+pulse_length + 5
    pulse_starts = array([pulse1_start, pulse2_start, pulse3_start, pulse4_start, pulse5_start])
    # voltage steps
    steps_start = where(v_data[-1] == max(v_data[-1]))[0][0] - 1 
    print ('Test pulses starts:', pulse_starts*dt)
    print ('Test pulses length:', pulse_length*dt)
    print ('Interval between test pulses:', pulse_interval*dt)
    print ('Steps start:', steps_start*dt)
    
    # figure('P/5 subtraction %i,  %s, %i, %s' %(date, retina, cell, n_rec), figsize=(15,8))
    # cmap = plt.get_cmap('gnuplot')
    # cols = [cmap(i) for i in np.linspace(0, 1, n_rec)]
    
    I_corrected = [] 
    I_cut = []
    for i in range(n_rec):
        # removing baseline current at the beginning of each sweep (before first test pulse)
        i_baseline_steps = mean(i_data[i][steps_start - int(20*ms/dt): steps_start - int(5*ms/dt)])
        #print (i, i_baseline_steps)
        
        # P5 subtraction
        step_amp =  v_data[i][steps_start + 50] - v_data[i][steps_start-50] # amplitude of the voltage step
        n_pulse = int(step_amp/pulse_amp)                                   # how much step pulses fit in the step
        if n_pulse > 5:
            n_pulse = 5
        fraction_last_pulse = step_amp/pulse_amp - n_pulse                  # fraction of the last test pulse that has to be added
        #print (step_amp, n_pulse, fraction_last_pulse)
    
        i_sum_prepulse = 0
        for j in range(n_pulse):
            i_baseline_tp = mean(i_data[i][pulse_starts[j] - int(15*ms/dt):pulse_starts[j] - 10])
            #print (i_baseline_tp)
            i_pulse = i_data[i][pulse_starts[j]: pulse_starts[j] + pulse_length] - i_baseline_tp
            i_sum_prepulse += i_pulse
    
        i_sum_prepulse = i_sum_prepulse + fraction_last_pulse * i_data[i][pulse5_start:\
                        pulse5_start + pulse_length]
        i_step_corrected = i_data[i][steps_start:steps_start + pulse_length] - i_baseline_steps - i_sum_prepulse
        

        I_corrected.append(i_step_corrected)
        I_cut.append(i_data[i][steps_start:steps_start + pulse_length])
        
        # plot(i_data[i][steps_start:steps_start + pulse_length], '-', color=cols[i])
        # plot(i_step_corrected, '--', color=cols[i])
        # plot(i_data[i][pulse1_start:pulse1_start + pulse_length], 'g-')
    
    t_cut = t[steps_start:steps_start + pulse_length] - t[steps_start]

    return I_corrected, I_cut, t_cut
    

def plot_IV(date, retina, cell, dt, i_data, v_data, t_start, rec_name):
    '''
    Returns the IV curve
    '''
    n_rec = len(i_data)
    print (n_rec)
    t = dt*arange(len(i_data[0])) 
    idx_step = where(v_data[-1] == max(v_data[-1]))[0][0] - 1
    
    i_peaks = []
    #i_peaks_amp = []
    vc_peaks = []
    t_peaks = []
    
    figpeaks = figure('Peaks %i,  %s, %i, %s' %(date, retina, cell, rec_name), (8,8))
    cmap = plt.get_cmap('gnuplot')
    cols = [cmap(i) for i in np.linspace(0, 1, n_rec)]
    
    subplot(211)
    for i in range(n_rec):
        i_peak, t_peak = find_peak_current(t, i_data[i], t_start*ms, 10.*ms, dt) # pA
        
        i_peaks.append(i_peak)
        vc_peaks.append(v_data[i][idx_step+10])
        t_peaks.append(t_peak)
        #plotting
        plot(t/ms, i_data[i],  color=cols[i])
        plot(t_peak/ms, i_peak, 'ko')
        #plot(t[t_leak]/ms, i_leak, 'g.')
        ylabel('I (pA)')
        xlabel('Time (ms)')
        #ylim(-1200,500)
    
    plot(t_start*ones(100), linspace(-2000,2000, 100), 'k--')
    
    # Peak axonal current
    dv = vc_peaks[1] - vc_peaks[0]
    didv = array([(i_peaks[j] - i_peaks[j+1])/dv for j in range(len(i_peaks)-1)])
    ddidv = array([(didv[j] - didv[j+1])/dv**2 for j in range(len(didv)-1)])
    
    idx_peak1 = argmax(abs(didv)) + 1
    idx_peak2 = argmax(abs(delete(didv, idx_peak1-1))) + 1
    idx_peak3 = argmax(abs(delete(didv, [idx_peak1-1, idx_peak2-1]))) + 1
    idx_peak4 = argmax(abs(delete(didv, [idx_peak1-1, idx_peak2-1, idx_peak3-1]))) + 1
    idx_peak5 = argmax(abs(delete(didv, [idx_peak1-1, idx_peak2-1, idx_peak3-1, idx_peak4-1]))) + 1
    idx_peak6 = argmax(abs(delete(didv, [idx_peak1-1, idx_peak2-1, idx_peak3-1, idx_peak4-1, idx_peak5-1]))) + 1
    idx_peak7 = argmax(abs(delete(didv, [idx_peak1-1, idx_peak2-1, idx_peak3-1, idx_peak4-1, idx_peak5-1, idx_peak6-1]))) + 1
    
    print (i_peaks[idx_peak1], i_peaks[idx_peak2], i_peaks[idx_peak3])
        
    # if i_peaks[idx_peak] == min(i_peaks): # if the somatic peak current was detected
    if i_peaks[idx_peak1] <= 2*i_peaks[idx_peak2] and i_peaks[idx_peak2] < 3*i_peaks[0]: # and abs(didv[idx_peak1-1]) > abs(delete(didv, idx_peak1-1)[idx_peak2-1]):
        idx_peak = idx_peak2
        if i_peaks[idx_peak2] <= 2*i_peaks[idx_peak3] and i_peaks[idx_peak3] < 3*i_peaks[0]: # and abs(delete(didv, idx_peak1-1)[idx_peak2-1]) > abs(delete(didv, [idx_peak1-1, idx_peak2-1])[idx_peak3-1]):
            idx_peak = idx_peak3
            if i_peaks[idx_peak3] <= 2*i_peaks[idx_peak4] and i_peaks[idx_peak4] < 3*i_peaks[0]: # and abs(delete(didv, [idx_peak1-1, idx_peak2-1])[idx_peak3-1]) > abs(delete(didv, [idx_peak1-1, idx_peak2-1, idx_peak3-1])[idx_peak4-1]):
                idx_peak = idx_peak4
                if i_peaks[idx_peak4] <= 2*i_peaks[idx_peak5]: # and abs(delete(didv, [idx_peak1-1, idx_peak2-1, idx_peak3-1])[idx_peak4-1]) > abs(delete(didv, [idx_peak1-1, idx_peak2-1, idx_peak3-1, idx_peak4-1])[idx_peak5-1]):
                    idx_peak = idx_peak5
                    if i_peaks[idx_peak5] <= 2*i_peaks[idx_peak6]: # and abs(delete(didv, [idx_peak1-1, idx_peak2-1, idx_peak3-1])[idx_peak4-1]) > abs(delete(didv, [idx_peak1-1, idx_peak2-1, idx_peak3-1, idx_peak4-1])[idx_peak5-1]):
                        idx_peak = idx_peak6
                        if i_peaks[idx_peak6] <= 2*i_peaks[idx_peak7]: # and abs(delete(didv, [idx_peak1-1, idx_peak2-1, idx_peak3-1])[idx_peak4-1]) > abs(delete(didv, [idx_peak1-1, idx_peak2-1, idx_peak3-1, idx_peak4-1])[idx_peak5-1]):
                            idx_peak = idx_peak7
    elif i_peaks[idx_peak1] > i_peaks[0]:
        idx_peak = argmax(abs(delete(didv, idx_peak1-1))) + 1
    
    else:
        idx_peak = idx_peak1
    
    print (dv, max(abs(didv)))
    
    
    # mean_subthres_current = mean(i_peaks[:2])
    # print ('Mean subthres current:', mean_subthres_current)
    # #if abs(i_peaks[idx_peak]) < 500.:
    # if abs(mean_subthres_current) < 1: 
    #     #mean_subthres_current = mean(i_peaks[:6])
    #     #print ('Mean subthres current:', mean_subthres_current)
    #     idx_peak = where(abs(array(i_peaks)) > 1000*abs(mean_subthres_current))[0][0]
    # elif abs(mean_subthres_current) < 5: 
    #     #mean_subthres_current = mean(i_peaks[:6])
    #     #print ('Mean subthres current:', mean_subthres_current)
    #     idx_peak = where(abs(array(i_peaks)) > 100*abs(mean_subthres_current))[0][0]
    # elif abs(mean_subthres_current) < 20: 
    #     #mean_subthres_current = mean(i_peaks[:6])
    #     #print ('Mean subthres current:', mean_subthres_current)
    #     idx_peak = where(abs(array(i_peaks)) > 50*abs(mean_subthres_current))[0][0]
    # elif abs(mean_subthres_current) > 400: 
    #     #mean_subthres_current = mean(i_peaks[:6])
    #     #print ('Mean subthres current:', mean_subthres_current)
    #     idx_peak = where(abs(array(i_peaks)) > 1.5*abs(mean_subthres_current))[0][0]
    # else:
    #     idx_peak = where(abs(array(i_peaks)) > 5*abs(mean_subthres_current))[0][0] # 8 works better for a few recs
    
        
    subplot(212)
    plot(vc_peaks, i_peaks, 'o-', color= 'k', label='peak') 
    plot(vc_peaks[idx_peak], i_peaks[idx_peak], 'ro') 
    #plot(vc_peaks, i_peaks_amp, 'o-', color= 'green', label='amplitude') 
    ylabel('I peak (pA)')
    xlabel('V (mV)')
    legend(frameon=False)

    return i_peaks,  vc_peaks, idx_peak, t_peaks

def rs_correction(date, retina, cell, dt, i_data, v_com, f_rs, v_rev, rs, cm):
    '''
    Traynelis 1998 method for linear IV curve, only Rs correction.
    f_rs in the remaining proportion of correction to be applied (1 - % or on-line correction)
    '''

    n_rec = len(i_data)
    
#    # filtering the data
#    sampling_freq = 1./dt
#    nyq_freq = sampling_freq/2 #12.5 Hz
#    b, a = signal.butter(10, 0.25, 'low') # second arg is a fraction of Nyq freq

    I_corr_only_rs = zeros((n_rec, len(i_data[0])))
    
    for i in range(n_rec):
        Ii = i_data[i] #signal.filtfilt(b, a, i_data[i]) * pA
        n_pts = len(Ii)
        
        v_last = v_com[i] - Ii[0] * f_rs * rs
        denom = v_last - v_rev
        
        if denom != 0:
            fracCorr = f_rs * (1-(v_com[i]-v_rev)/denom)
        else:
            fracCorr = 0
                
        for j in range(n_pts):
            v_this = v_com[i] - Ii[j]* f_rs * rs
            #print (v_this, )
            
            if v_this != v_rev:
                fracCorr = f_rs * (1-(v_com[i]-v_rev)/(v_this-v_rev))
                #print (v_this, v_com[i], Ii[j], fracCorr)
            else:
                fracCorr = 0
            
            I_corr_only_rs[i][j-1] = Ii[j-1] * (1-fracCorr)
            
            v_last = v_this 
        #print (fracCorr)
    
    return I_corr_only_rs















