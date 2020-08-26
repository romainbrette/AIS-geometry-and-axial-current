
"""

Measurement of the transferred charge from the axonal currents. (RGC_transferred_charge)

"""
import glob2
import pandas as pd
import pyabf
from pandas import ExcelWriter
from pandas import ExcelFile
from vc_test_pulse_analysis import *
from na_currents_analysis import *
from scipy import interpolate

path_files = '/Users/sarahgoethals/Dropbox/Spike initiation/PhD projects/Axonal current and AIS geometry/Data patch/RGC/'

### Path to datafiles: load the list of cells used for the analysis
path_to_files = '/Users/sarah/Documents/repositories/AIS-geometry-and-axonal-current/Na currents in RGC/codes submission/data/'
df_cells = pd.read_excel(path_files + 'RGC_electrical_properties.xlsx')

first_cell = 6
last_cell = 7 #len(df_cells['Date'])

# Cells identity and age
dates = array(df_cells['Date'])[first_cell:last_cell]
retinas = array(df_cells['Retina'])[first_cell:last_cell]
cells = array(df_cells['Cell'])[first_cell:last_cell]
ages = array(df_cells['Age'])[first_cell:last_cell]
# Recording used to measure axial current
na_recs = array(df_cells['Recording'])[first_cell:last_cell]
# Sweep of the peak axonal current
n_sweeps = array(df_cells['Sweep number'])[first_cell:last_cell]
# TP with compensation to correct for passive component for cells w/o P5
tp_corrs = array(df_cells['TP num correction'])[first_cell:last_cell]

### Path to the data
path_to_data = '/Users/sarah//Documents/Data/Martijn Sierksma/'

current_integrals = []
current_integrals_10 = []
current_durations_10 = []
current_durations_50 = []
peak_current_latency = []

N = 0 # counting analyzed cells

for date, retina, cell, age, na_rec, n_sweep, tp_corr in zip(dates, retinas, cells, ages, na_recs, n_sweeps, tp_corrs): 
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
    
    ### Plotting raw data
    abf = pyabf.ABF(path_to_na_currents[0]) 
    fs = abf.dataRate  * Hz # sampling rate
    dt = 1./fs
    t = dt*arange(len(abf.sweepY)) 
    I = []
    V = []
    n_rec = len(abf.sweepList)
        
    f1 = figure('%i,  %s, %i' %(date, retina, cell), (8,6))
    ax1 = f1.add_subplot(211)
    
    for sweepNumber in range(int(n_sweep)-1, int(n_sweep)+3):
        abf.setSweep(sweepNumber)
        I.append(abf.sweepY)
        V.append(abf.sweepC*mV)
        
        ax1.plot(t[int(625.6*ms/dt):int(630.6*ms/dt)]/ms, abf.sweepY[int(625.6*ms/dt):int(630.6*ms/dt)], 'k')
        #ax1.set_xlim(55, 70)
        ax1.set_ylabel('I (pA)')

    show()
    
    ### Passive component correction 
    
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
                
        for sweepNumber in abf_tp.sweepList:
            abf_tp.setSweep(sweepNumber)
            I_tp.append(abf_tp.sweepY)
            V_tp.append(abf_tp.sweepC*mV)
            
        I_corr, I_cut, t_cut = remove_passive_component(date, retina, cell, dt, I, V, dt_tp, I_tp, V_tp,  str(int(na_rec)).zfill(4))  
        
    else: # cells with P5 subtraction protocol
        I_corr, I_cut, t_cut = p5_subtraction(date, retina, cell, dt, I, V, rec_name=str(int(na_rec)).zfill(4))
    
    I_corr = I_corr*pA

    ### Plotting corrected current
    ax2 = f1.add_subplot(212)
    ax2.plot(t_cut/ms, I_corr[1]/pA, 'k')
    ax2.plot(t_cut/ms, I_corr[2]/pA, 'k')
    ax2.plot(t_cut/ms, zeros(len(t_cut)), 'k--')
    #ax2.set_xlim(55, 70)
    ax2.set_ylabel('I (pA)')
    
    ### Integration of the current to estimate the transferred charge
    curr_int = []
    curr_int_10 = []
    curr_duration_50 = []
    curr_duration_10 = []
    peak_lat = []
    for i in range(1, 3): # interating the current just above thres and the next one too
        # interpolating the current
        if date < 20190325: # first experiments are less well leak subtracted
            f = interpolate.interp1d(t_cut[int(0.5*ms/dt):]/ms, I_corr[i][int(0.5*ms/dt):]/pA, fill_value="extrapolate")  
            dt_new = 0.001
            t_new = np.arange(0.5, 10, dt_new)
            i_new = f(t_new)
            ax2.plot(t_new, i_new, 'g-')
        else:
            f = interpolate.interp1d(t_cut/ms, I_corr[i]/pA)
            dt_new = 0.001
            t_new = np.arange(0, 10, dt_new)
            i_new = f(t_new)
            ax2.plot(t_new, i_new, 'g-')
        # detecting current peak
        idx_peak_new = argmin(i_new)
        i_peak_new = i_new[idx_peak_new]
        t_peak_new = t_new[idx_peak_new]
        # current duration
        idx_10 = where(abs(i_new) <= abs(0.1*i_peak_new))[0]
        #ax2.plot(t_new[idx_10], i_new[idx_10], 'r.')
        idx_start_dur_10 = where(idx_10 < idx_peak_new)[0][-1]
        idx_end_dur_10 = where(idx_10 > idx_peak_new)[0][0] 
        ax2.plot(t_new[idx_10][idx_start_dur_10], i_new[idx_10][idx_start_dur_10], 'bo')
        ax2.plot(t_new[idx_10][idx_end_dur_10], i_new[idx_10][idx_end_dur_10], 'bo')
        
        dur_10 = t_new[idx_10][idx_end_dur_10] - t_new[idx_start_dur_10]
        print ('Duration 10:', dur_10, 'ms')
        
        idx_50 = where(abs(i_new) <= abs(0.5*i_peak_new))[0]
        idx_start_dur_50 = where(idx_50 < idx_peak_new)[0][-1]
        idx_end_dur_50 = where(idx_50 > idx_peak_new)[0][0] 
        ax2.plot(t_new[idx_start_dur_50], i_new[idx_start_dur_50], 'bo')
        ax2.plot(t_new[idx_50][idx_end_dur_50], i_new[idx_50][idx_end_dur_50], 'bo')
        
        dur_50 = t_new[idx_50][idx_end_dur_50] - t_new[idx_start_dur_50]
        print ('Duration 10:', dur_50, 'ms')
        
        # total current integration
        idx_peak = argmin(I_corr[i])
        i_peak = I_corr[i][idx_peak]
        idx_end_current = where(I_corr[i][idx_peak:] > 0)[0][0] + idx_peak
    
        area_idx = arange(idx_end_current)
        area_time = t_cut[area_idx]
        area_I = I_corr[i][area_idx]   
        area = 0
        for i in range(len(area_I)-1): 
            area_dt = (area_I[i]*(area_time[i+1]-area_time[i])) + \
                (0.5*((area_I[i+1]-area_I[i])*(area_time[i+1]-area_time[i])))   # the top triangle   
            if area_dt < 0.:
                area += area_dt

        print ('Area:', area, '= ', area/pcoulomb, 'nA ms')
    
        ax2.fill_between(x=area_time/ms, y1=area_I/pA, color='k', alpha=0.3)
        
        # above 10% current integration    
        area_idx_10 = arange(idx_start_dur_10, idx_10[idx_end_dur_10])
        area_time_10 = t_new[area_idx_10]
        area_box_10 = dur_10 * i_new[idx_start_dur_10] * 1e-3
        area_I_10 = i_new[area_idx_10]   
        area_10 = 0
        for i in range(len(area_I_10)-1): 
            area_dt_10 = (area_I_10[i]*(area_time_10[i+1]-area_time_10[i])) + \
                (0.5*((area_I_10[i+1]-area_I_10[i])*(area_time_10[i+1]-area_time_10[i])))   # the top triangle   
            if area_dt_10 < 0.:
                area_10 += area_dt_10
            
        area_10 = (area_10 - area_box_10) * 1e-3
        print ('Area 10:', area_10) 
    
        ax2.fill_between(x=area_time_10, y1=area_I_10, color='k', alpha=0.3)
        
        curr_int.append(area/pcoulomb)
        curr_int_10.append(area_10)
        curr_duration_10.append(dur_10)
        curr_duration_50.append(dur_50)
        peak_lat.append(idx_peak*dt/ms)

    tight_layout()
    
    current_integrals.append(curr_int)
    current_integrals_10.append(curr_int_10)
    current_durations_10.append(curr_duration_10)
    current_durations_50.append(curr_duration_50)
    peak_current_latency.append(peak_lat)
    
# df_select_cells = pd.DataFrame({'Date': dates,
#                   'Retina': retinas,
#                   'Cell': cells,
#                   'Age': ages,
#                   'Charge1': [current_integrals[i][0] for i in range(len(dates))],
#                   'Charge2': [current_integrals[i][1] for i in range(len(dates))],
#                   'Charge1 10': [current_integrals_10[i][0] for i in range(len(dates))],
#                   'Charge2 10': [current_integrals_10[i][1] for i in range(len(dates))],
#                   'Peak latency1': [peak_current_latency[i][0] for i in range(len(dates))],
#                   'Peak latency2': [peak_current_latency[i][1] for i in range(len(dates))],
#                   'Duration1 10': [current_durations_10[i][0] for i in range(len(dates))],
#                   'Duration2 10': [current_durations_10[i][1] for i in range(len(dates))],
#                   'Duration1 50': [current_durations_50[i][0] for i in range(len(dates))],
#                   'Duration2 50': [current_durations_50[i][1] for i in range(len(dates))],
#                   })

# df_select_cells.to_excel(path_files + "RGC_charge_at_SI_1906.xlsx", \
#                           columns=['Date','Retina','Cell','Age',\
#                                     'Charge1', 'Charge2',
#                                     'Charge1 10', 'Charge2 10',
#                                   'Peak latency1', 'Peak latency2',
#                                   'Duration1 10', 'Duration2 10',
#                                   'Duration1 50', 'Duration2 50'])


    
    