

"""

Peak axonal current analysis. (RGS_Na_current_correction_C_CC)

"""
# import os
import glob2
import pandas as pd
import pyabf
from pandas import ExcelWriter
from pandas import ExcelFile
from vc_test_pulse_analysis import *
from na_currents_analysis import *

### Path to datafiles: load the list of cells used for the analysis
path_to_files = '/Users/sarah/Documents/repositories/AIS-geometry-and-axonal-current/Na currents in RGC/codes submission/data/'
df_cells = pd.read_excel(path_files + 'RGC_electrical_properties.xlsx')
df_pass = pd.read_excel(path_files + 'RGC_passive_properties.xlsx')

first_cell = 20
last_cell = 21 #len(df_cells['Date'])

# Cells identity and age
dates = array(df_cells['Date'])[first_cell:last_cell]
retinas = array(df_cells['Retina'])[first_cell:last_cell]
cells = array(df_cells['Cell'])[first_cell:last_cell]
ages = array(df_cells['Age'])[first_cell:last_cell]
# Recording used to measure axial current
na_recs = array(df_cells['Recording'])[first_cell:last_cell]
# TP with compensation to correct for passive component for cells w/o P5 protocol
tp_corrs = array(df_cells['TP num correction'])[first_cell:last_cell]
# Residual Rs before the recording
rs_residuals = array(df_cells['Residual Rs'])[first_cell:last_cell]
# Holding membrane potential at rest
resting_mp = array(df_cells['V holding (mV)'])[first_cell:last_cell]
# TP w/o compensation used to measure capacitance and series resistance
tp_nums = array(df_cells['TP num passive props'])[first_cell:last_cell]
    
### Path to the data
path_to_data = '/Users/sarah//Documents/Data/Martijn Sierksma/'

selected_dates = []
selected_retinas = []
selected_cells = []
selected_ages = []
selected_tp = []
selected_sweeps = []
t_currents = []
selected_recordings = []

axonal_currents = []
axonal_currents_corr = []
threshold_currents_raw = []

voltage_threshold = []
voltage_threshold_true = []
time_cst_short = []
capacitance_from_tau = []

N = 0 # counting analyzed cells

for date, retina, cell, age, na_rec, tp_corr, rs_res, vh in zip(dates, retinas, cells, ages, na_recs, tp_corrs, rs_residuals, resting_mp): 
    print ('------------------------------')
    print (date, retina, cell)
    
    ### Load Na current recording
    path_to_cell = path_to_data + str(int(date)) + "*/" + '/retina '+ str(retina) +'/cell ' + str(int(cell))
    
    if (date, retina, cell) == (20191115, 'B', 2): # for that cell, a recording of the adaptation protocol is used
        print('adaptation protocol')
        path_to_na_currents = glob2.glob(path_to_cell + '/VC threshold adaptation/' + '*' + str(int(na_rec)).zfill(4) + ".abf")
    else:
        path_to_na_currents = glob2.glob(path_to_cell + '/VC small steps/' + '*' + str(int(na_rec)).zfill(4) + ".abf")
        
    N += 1 
    selected_dates.append(date)
    selected_retinas.append(retina)
    selected_cells.append(cell)
    selected_ages.append(age)
    selected_recordings.append(na_rec)
    
    ### Plotting raw data
    abf = pyabf.ABF(path_to_na_currents[0]) 
    fs = abf.dataRate  * Hz # sampling rate
    dt = 1./fs
    t = dt*arange(len(abf.sweepY)) 
    I = []
    V = []
    n_rec = len(abf.sweepList)
    
    cmap = plt.get_cmap('gnuplot')
    cols = [cmap(i) for i in np.linspace(0, 1, n_rec)]
    
    f1 = figure('%i,  %s, %i' %(date, retina, cell), (12,6))
    ax1 = f1.add_subplot(221)
    ax2 = f1.add_subplot(222)
    
    for sweepNumber in range(n_rec):
        abf.setSweep(sweepNumber)
        I.append(abf.sweepY)
        V.append(abf.sweepC*mV)
        
        ax1.plot(t/ms, abf.sweepY, color=cols[sweepNumber])
        #ax1.set_xlim(55, 70)
        ax1.set_ylabel('I (pA)')
        
        ax2.plot(t/ms, abf.sweepC, color=cols[sweepNumber])
        #xlim(55, 70)
        ax2.set_ylabel('Vc (mV)')
        ax2.set_xlabel('Time (ms)')

    show()
    
    ### Passive component correction and estimation of Tau
    
    if date < 20190610: # no P5 subtraction protocol
        ### Correcting for the passive component of the current
        # Loading the test pulse (w/o compensation ON)
        path_to_test_pulse = glob2.glob(path_to_cell + '/VC test pulses/comp/' + '*' + str(int(tp_corr)).zfill(4) + ".abf") # the TP has Rs compensation ON
        
        abf_tp = pyabf.ABF(path_to_test_pulse[0])
        fs_tp = abf_tp.dataRate  * Hz # sampling rate
        dt_tp = 1./fs_tp
        t_tp = dt_tp*arange(len(abf_tp.sweepY)) 
        I_tp = []
        V_tp = []
        
        ax3 = f1.add_subplot(223)
        
        for sweepNumber in abf_tp.sweepList:
            abf_tp.setSweep(sweepNumber)
            I_tp.append(abf_tp.sweepY)
            V_tp.append(abf_tp.sweepC*mV)
            
            ax3.plot(t_tp/ms, abf_tp.sweepY, color=cols[sweepNumber])
            ax3.set_ylabel('I (pA)')
            
            I_corr_pass, I_cut, t_cut = remove_passive_component(date, retina, cell, dt, I, V, dt_tp, I_tp, V_tp,  str(int(na_rec)).zfill(4))  
    
        ### Tau estimation
        # Loading the test pulse (with compensation ON)
        path_to_test_pulses = glob2.glob(path_to_data + str(int(date)) + "*/" + \
                '/retina '+ str(retina) +'/cell ' + str(int(cell)) + '/VC test pulses/comp/' +'*')
        n_tp_tot = len(path_to_test_pulses) # total number of test pulses w/o Rs compensation
        print ('Number of TP:', n_tp_tot)
                
        # Looping on a cell's test pulses
        for tp in range(1): 
            tp_num = path_to_test_pulses[tp][-8:-4] # the TP number
          
            # Load sweeps
            abf = pyabf.ABF(path_to_test_pulses[tp])
            print ('Number of sweeps:', len(abf.sweepList))
            
            # Measuring passive properties from each individual sweep
            tau_m = []
            
            # Plotting traces
            figure('Capacitance estimation %i,  %s, %i, %s' %(date, retina, cell, tp_num))
            ax1 = subplot(211)
            ax2 = subplot(212)
            
            # Looping on all sweeps 
            for sweepNumber in abf.sweepList:
                abf.setSweep(sweepNumber)
                fs = abf.dataRate  * Hz # sampling rate
                dt = 1./fs
                t = dt*np.arange(len(abf.sweepY))/second
                i = abf.sweepY 
                v = abf.sweepC 
                ax1.plot(t, i)
                ax1.set_xlim(0.055, 0.150)
                
                # test pulse start
                idx_start = int(56.2*ms/dt)
        
                # transient positive peak
                idx_i_max = argmax(i[idx_start+5:idx_start+500]) + idx_start+5
                i_amp_max_peak = i[idx_i_max]
    
                # exponential fitting
                fit_end = idx_i_max + int(1.*ms/dt) 
                i_amp_end = i[fit_end]
                
                max_amp = i_amp_end
                
                ax2.plot(t/ms, i)
                ax2.set_xlim(56, 64)
                ax2.set_xlabel('t (s)')
                
                def exp_current(t, tau):
                    return  (i_amp_max_peak - max_amp) * exp(-t/tau) + max_amp
                
                idx_i_max = idx_i_max + 1
                peak_current = i[idx_i_max: fit_end] # pA
                peak_time = (t[idx_i_max:fit_end]- t[idx_i_max])*1e3  # sec
                
                ax2.plot(t[idx_i_max:fit_end]*1e3, peak_current, 'k-')
    
                tau_opt = curve_fit(exp_current, peak_time, peak_current)
                
                ax2.plot(peak_time + t[idx_i_max]*1e3, exp_current(peak_time, tau_opt[0]), 'r-')
                
                print ('Tau:', tau_opt[0])
                tau_m.append(tau_opt[0])
                
        Tau = median(tau_m)
        time_cst_short.append(Tau)
        
    else: # cells with P5 subtraction protocol
    
        ### Capacitance estimation from P/N pulses 
        tau_m = []
        
        # Plotting traces
        figure('Capacitance estimation %i,  %s, %i' %(date, retina, cell))
        ax1 = subplot(211)
        ax2 = subplot(212)
        
        # Looping on all sweeps 
        for i in range(n_rec):
            ax1.plot(t, I[i])
            ax1.set_xlim(25, 40)
            
            # start and end of test pulse
            if (date, retina, cell) == (20191115, 'B', 2):
                idx_start = int(81.22*ms/dt)
                idx_end = int(100.*ms/dt)
            else:
                idx_start = int(25.6*ms/dt)
                idx_end = int(38.*ms/dt)
            
            i_cut = I[i][idx_start:idx_end]
            t_cut = t[idx_start:idx_end]/ms

            # transient negative large peak
            idx_i_max = argmin(i_cut[5:500]) + 5
            i_amp_max_peak = i_cut[idx_i_max] 

            # exp fitting
            fit_end = idx_i_max + int(.5*ms/dt) 
            i_amp_end = i_cut[fit_end]
            
            max_amp = max(i_cut[idx_i_max:fit_end])
            
            ax2.plot(t_cut, i_cut)
            ax2.set_xlim(idx_start*dt/ms, (idx_start*dt + 4*ms)/ms)
            ax2.set_xlabel('t (s)')
            
            def exp_current(t, tau):
                return  (i_amp_max_peak - max_amp) * exp(-t/tau) + max_amp
            
            idx_i_max = idx_i_max + 1
            peak_current = i_cut[idx_i_max: fit_end] # pA
            peak_time = (t_cut[idx_i_max:fit_end]- t_cut[idx_i_max])  # sec
            
            ax2.plot(t_cut[idx_i_max:fit_end], peak_current, 'k-')

            tau_opt = curve_fit(exp_current, peak_time, peak_current)
            
            ax2.plot(peak_time + t_cut[idx_i_max], \
                     exp_current(peak_time, tau_opt[0]), 'r-')
            print ('Tau m:', tau_opt[0])
            tau_m.append(tau_opt[0])
        
        Tau = median(tau_m)
        time_cst_short.append(Tau)
        
        ### P5 subtraction
        I_corr_pass, I_cut, t_cut = p5_subtraction(date, retina, cell, dt, I, V, rec_name=str(int(na_rec)).zfill(4))
    
    ### Series resistance correction
    
    ### IV curve
    I_peaks,  Vc_peaks, idx_peak_ax_current, t_peaks = plot_IV(date, retina, cell, dt, I_corr_pass, V, 0, str(int(na_rec)).zfill(4))
    
    ###  Correction with Traynelis method
    v_rev = 70. # Na reversal potential
    si = dt/ms
    tauLag = si * 2.
    fc = 1./(2*pi*tauLag)
    filterfactor = 1
     
    # peak of the axonal current
    idx_peak = argmin(I_corr_pass[idx_peak_ax_current][int(t_start*ms/dt):int(9.*ms/dt)]) + int(t_start*ms/dt)
    i_data = I_corr_pass[idx_peak_ax_current][idx_peak-50: idx_peak+20] * 1e-3 #nA 
    
    figure('Rs correction %i,  %s, %i, %0.2f' %(date, retina, cell, Taum), figsize=(12,6))
    
    f_rs = 1 
    f_c = 1
    
    cm = Taum/rs_res # nF
    print ('Cm:', cm)
    
    capacitance_from_tau.append(cm*1e3) #pF
    
    vc_data = Vc_peaks[idx_peak_ax_current]/mV + vh 
    print ('Vh = ', vc_data)
    
    I_cap = zeros(len(i_data))
    I_corr_cap = zeros(len(i_data))
    I_corr_cap_filt = zeros(len(i_data))
    I_corr = zeros(len(i_data))
    
    n_pts = len(i_data)
    
    # first data point
    v_last = vc_data - i_data[0] * rs_res
    denom = v_last - v_rev
    
    if denom != 0:
        fracCorr = f_rs * (1-(vc_data-v_rev)/denom)
    else:
        fracCorr = 0
    
    I_corr[0] = i_data[0] * (1-fracCorr)
    
    # next data points
    for j in range(n_pts):
        v_this = vc_data - i_data[j]*rs_res
    
        if v_this != v_rev:
            fracCorr = f_rs * (1-(vc_data-v_rev)/(v_this-v_rev))
        else:
            fracCorr = 0
        Icap = cm * (v_this-v_last)/si #* 1e-3
        I_cap[j] = Icap
        Icap = Icap * filterfactor
        I_corr_cap_filt[j] = Icap
        I_corr_cap[j-1] = i_data[j-1] - f_c * Icap
        I_corr[j-1] = I_corr_cap[j-1] * (1-fracCorr)
    
        v_last = v_this 
    
    subplot(211)
    plot(i_data, 'k-')
    plot(I_corr, 'r')
    
    subplot(212)
    plot(I_cap, 'b')

    show()
    
    axonal_currents_corr.append(min(I_corr))
    axonal_currents.append(I_peaks[idx_peak_ax_current])
    voltage_threshold.append(Vc_peaks[idx_peak_ax_current-1]/mV + vh)
    
    # # Raw current at threshold (highest non spiking situation)
    # #smoothing
    # n = len(t_cut/ms)
    # i_slide = np.zeros(n)
    # d = 30 # half-window, i.e. number of pixels on each side

    # for j in range(n):
    #     if j < d: # start of the axon, not full window
    #         i_slide[j] = np.mean(I_corr_pass[idx_peak_ax_current-1][0:j+d])
    #     elif j > n-d: # end of the axon, not full window
    #         i_slide[j] = np.mean(I_corr_pass[idx_peak_ax_current-1][j-d:n])
    #     else: 
    #         i_slide[j] = np.mean(I_corr_pass[idx_peak_ax_current-1][j-d:j+d])
        
    # figure('Thres curr %i,  %s, %i' %(date, retina, cell))
    # plot(I_corr_pass[idx_peak_ax_current-1])
    # plot(I_cut[idx_peak_ax_current-1])
    # plot(i_slide, 'k')
    
    # idx_th = argmin(i_slide[50:400]) + 50 #argmin(I_corr_pass[idx_peak_ax_current-1][50:400]) + 50
    # Ie_thres = I_cut[idx_peak_ax_current-1][idx_th]
    
    # plot(idx_th, Ie_thres, 'ko')
    
    # threshold_currents_raw.append(Ie_thres * 1e-3)
    
    # Corrected voltage trace
    idx_step = where(V[-1] == max(V[-1]))[0][0] - 1
    l = len(I_cut[0])
    Vc_th_true = V[idx_peak_ax_current-1][idx_step: idx_step + l]/mV + vh - rs_res * I_cut[idx_peak_ax_current-1] * 1e-3
    idx_vth_max = argmax(Vc_th_true[50:]) + 50

    figure('Corr voltage %i,  %s, %i' %(date, retina, cell))
    plot(V[idx_peak_ax_current-1][idx_step: idx_step + l]/mV + vh)
    plot(Vc_th_true)
    plot(idx_vth_max, Vc_th_true[idx_vth_max], 'ko')
    plot(idx_th,  Vc_th_true[idx_th], 'ro')
    
    voltage_threshold_true.append(Vc_th_true[idx_vth_max])
    


### Write in excel file


# df_select_cells = pd.DataFrame({'Date': selected_dates,
#                   'Retina': selected_retinas,
#                   'Cell': selected_cells,
#                   'Age': selected_ages,
#                   'Recording': selected_recordings,
#                   'Peak axonal current': axonal_currents,
#                   'Peak axonal current corr': axonal_currents_corr,
#                   'Vth': voltage_threshold,
#                   'Vth true': voltage_threshold_true,
#                   'Cm VC residual': capacitance_from_tau, 
#                   'Tau RC': membrane_time_cst_short})

# df_select_cells.to_excel(path_files + "RGC_peak_axonal_current_corrected.xlsx", \
#                           columns=['Date','Retina','Cell','Age','Recording',\
#                                   'Peak axonal current', 'Peak axonal current corr', \
#                                   'Vth','Vth true','Cm VC residual',\
#                                   'Tau RC'
#                                   ])
















        
