


"""

Analysis of adaptation data in cells without P/5 protocol.


"""
from brian2 import *
import glob2
import pandas as pd
import pyabf
from pandas import ExcelWriter
from pandas import ExcelFile
from scipy import signal 
from scipy import stats
from scipy.optimize import curve_fit
from vc_test_pulse_analysis import *
from na_currents_analysis import *

### Loading the cells passive properties and the list of the recording that will be analyzed
path_files = '/Users/sarah/Documents/repositories/AIS-geometry-and-axonal-current/Na currents in RGC/codes submission/data/'

# Load cells passive properties
cells_pp = pd.read_excel(path_to_files + 'RGC_passive_properties.xlsx')

# Load the recordings that will be used for the analysis
df_adapt = pd.read_excel(path_to_files + 'recordings_for_adaptation.xlsx')

# Loading the recording database to have the Rs compensation corresponding to the test pulses that are analysed
df_rec_database = pd.read_excel(path_to_data + 'RGC_recording_database.xlsx')

first_cell = 15
last_cell = 16 #len(cells_pp['Date'])

dates = array(cells_pp['Date'])[first_cell:last_cell]
retinas = array(cells_pp['Retina'])[first_cell:last_cell]
cells = array(cells_pp['Cell'])[first_cell:last_cell]
ages = array(cells_pp['Age'])[first_cell:last_cell]

### Path to the data
path_to_data = '/Volumes/Lolita/Martijn Sierksma/'

# record measures
selected_dates = []
selected_retinas = []
selected_cells = []
selected_ages = []
selected_recs = []
selected_sweeps = []

peak_axonal_currents = []
peak_axonal_currents_corr = []
peak_axonal_currents_raw = []

holding_potentials = []
prepulse_potentials = []
threshold_potentials = []

tp_before = []
tp_after = []
tp_na_rec = []

series_resistance_before = []
series_resistance_after = []
series_resistance_na_rec = []

series_resistance_rec = []
series_resistance_compensation = []
series_resistance_residual = []
compensation_rec = []

capacitance_tau = []
RC_time_cst = []

current_integrals_10 = []
current_durations_10 = []
current_durations_50 = []
peak_current_latency = []

N = 0

### Looping on all cells
for date, retina, cell, age in zip(dates, retinas, cells, ages): 
    if date < 20190611: # only with P5:
        print ('------------------------------')
        print (date, retina, cell)
        
        ### Path to the cell's data
        path_to_cell = path_to_data + str(int(date)) + "*/" + '/retina '+ str(retina) +'/cell ' + str(int(cell))
        
        ### Get the recordings that will be analyzed for this cell
        row_cell = df_adapt[(df_adapt['Date'] == date) & (df_adapt['Retina'] == retina) & (df_adapt['Cell'] == cell)]
        
        if len(row_cell['s1']) > 0:
            
            sweeps = row_cell.values[0][4:]
            
            # delete nans
            sweeps = sweeps[~pd.isnull(sweeps)]
            print (sweeps)
            
            N += 1 # counting recordings used for the analysis
                        
            ### Loop on recordings corresponding to different prepulse potentials
            for rec in sweeps:
                na_rec = int(rec)
                
                selected_dates.append(date)
                selected_retinas.append(retina)
                selected_cells.append(cell)
                selected_ages.append(age)
                selected_recs.append(f"{na_rec:04d}")
                
                ### Extract the series resistance measured during the experiment, 
                ### the series resistance compensation and the holding potential from the recording database.
                ### We make sure to obtain the last measured values before each recording.
                row = df_rec_database[(df_rec_database['Date'] == date) & (df_rec_database['retina'] == retina) & (df_rec_database['cell'] == cell)]
                idx_row = row.index[0]
            
                print ('Analyzed recording:', na_rec)
                
                if  int(array(df_rec_database['first filename'])[idx_row+1][1:]) > na_rec >= int(row['first filename'].values[0][1:]):
                    rs_num = row['first filename'].values[0][1:]
                    rs_rec = row['Rs'].values[0]
                    compensation = row['compensation'].values[0]
                    vh = row['V holding'].values[0]
                    print ('First filename:', row['first filename'].values[0][1:])
                elif int(array(df_rec_database['first filename'])[idx_row+2][1:]) > na_rec >= int(array(df_rec_database['first filename'])[idx_row+1][1:]) :
                    print ('Second filename:', array(df_rec_database['first filename'])[idx_row+1][1:])
                    rs_num = array(df_rec_database['first filename'])[idx_row+1][1:]
                    rs_rec = array(df_rec_database['Rs'])[idx_row+1]
                    compensation = array(df_rec_database['compensation'])[idx_row+1]
                    vh = array(df_rec_database['V holding'])[idx_row+1]
                    if type(compensation) == str:
                        continue
                elif int(array(df_rec_database['first filename'])[idx_row+3][1:]) > na_rec >= int(array(df_rec_database['first filename'])[idx_row+2][1:]) :
                    print ('Third filename:',  array(df_rec_database['first filename'])[idx_row+2][1:])
                    rs_num = array(df_rec_database['first filename'])[idx_row+2][1:]
                    rs_rec = array(df_rec_database['Rs'])[idx_row+2]
                    compensation = array(df_rec_database['compensation'])[idx_row+2]
                    vh = array(df_rec_database['V holding'])[idx_row+2]
                    if type(compensation) == str:
                        continue
                elif int(array(df_rec_database['first filename'])[idx_row+4][1:]) > na_rec >= int(array(df_rec_database['first filename'])[idx_row+3][1:]):
                    print ('Fourth filename:',  array(df_rec_database['first filename'])[idx_row+3][1:])
                    rs_num = array(df_rec_database['first filename'])[idx_row+3][1:]
                    rs_rec = array(df_rec_database['Rs'])[idx_row+3]
                    compensation = array(df_rec_database['compensation'])[idx_row+3]
                    vh = array(df_rec_database['V holding'])[idx_row+3]
                    if type(compensation) == str:
                        continue
                elif int(array(df_rec_database['first filename'])[idx_row+5][1:]) > na_rec >= int(array(df_rec_database['first filename'])[idx_row+4][1:]):
                    print ('Fifth filename:',  array(df_rec_database['first filename'])[idx_row+4][1:])
                    rs_num = array(df_rec_database['first filename'])[idx_row+4][1:]
                    rs_rec = array(df_rec_database['Rs'])[idx_row+4]
                    compensation = array(df_rec_database['compensation'])[idx_row+4]
                    vh = array(df_rec_database['V holding'])[idx_row+4]
                    if type(compensation) == str:
                        continue
                else :
                    print ('Sixth filename:',  array(df_rec_database['first filename'])[idx_row+5][1:])
                    rs_num = array(df_rec_database['first filename'])[idx_row+5][1:]
                    rs_rec = array(df_rec_database['Rs'])[idx_row+5]
                    compensation = array(df_rec_database['compensation'])[idx_row+5]
                    vh = array(df_rec_database['V holding'])[idx_row+5]
                    if type(compensation) == str:
                        continue
                
                series_resistance_compensation.append(compensation)
                series_resistance_rec.append(rs_rec)
                compensation_rec.append(rs_num)
                holding_potentials.append(vh)
                
                ### Extract the series resistance from passive properties measurements,
                ### from the test pulse before and after the recording
                row_pass_props = cells_pp[(cells_pp['Date'] == date) & (cells_pp['Retina'] == retina) & (cells_pp['Cell'] == cell)]
                tp_nums = array([row_pass_props['TP%i num' %i].values[0] for i in range(1, 10)])
                tp_nums = tp_nums[~pd.isnull(tp_nums)]
                idx_tp_na_rec_before = where(tp_nums <= na_rec)[0][-1]
                tp_na_rec_after = where(tp_nums >= na_rec)[0]
                if len(tp_na_rec_after) > 0:
                    idx_tp_na_rec_after = tp_na_rec_after[0]
                    print ('TP after rec:', tp_nums[idx_tp_na_rec_after])
                    rs_after = row_pass_props['Mean Rs%i' %(idx_tp_na_rec_after+1)].values[0]
                    series_resistance_after.append(rs_after)
                    tp_after.append(tp_nums[idx_tp_na_rec_after])
                else:
                    series_resistance_after.append(nan)
                    tp_after.append(nan)
                    
                ### We select the TP num that is the closest to the Na rec
                idx_tp_rs = argmin(abs(tp_nums - na_rec))
                tp_rs = tp_nums[idx_tp_rs]
                
                print ('TP before rec:', tp_nums[idx_tp_na_rec_before])
                
                rs_na_rec = row_pass_props['Mean Rs%i' %(idx_tp_rs+1)].values[0]
                if rs_na_rec != rs_na_rec: # if nan
                    print ('previous TP')
                    idx_tp_rs1 = argmin(abs(delete(tp_nums, idx_tp_rs) - na_rec))
                    tp_rs1 = tp_nums[idx_tp_rs1]
                    print (tp_rs1)
                    rs_na_rec1 = row_pass_props['Mean Rs%i' %(idx_tp_rs1+1)].values[0]
                    tp_rs = tp_rs1
                    rs_na_rec = rs_na_rec1
                    if rs_na_rec1 != rs_na_rec1: # if nan
                        idx_tp_rs2 = argmin(abs(delete(tp_nums, [idx_tp_rs, idx_tp_rs1]) - na_rec))
                        tp_rs2 = tp_nums[idx_tp_rs2]
                        rs_na_rec2 = row_pass_props['Mean Rs%i' %(idx_tp_rs2+1)].values[0]
                        tp_rs = tp_rs2
                        rs_na_rec = rs_na_rec2
                    
                print ('TP Rs:', tp_rs, rs_na_rec)
                
                rs_before = row_pass_props['Mean Rs%i' %(idx_tp_na_rec_before+1)].values[0]
                
                ### Residual series resistance
                rs_residual = rs_na_rec - compensation * rs_rec
                print ('Rs residual:', rs_residual)
                
                series_resistance_before.append(rs_before)
                series_resistance_na_rec.append(rs_na_rec)
                series_resistance_residual.append(rs_residual)
                tp_before.append(tp_nums[idx_tp_na_rec_before])
                tp_na_rec.append(tp_rs)
                
                ### Loading the Na currents
                path_to_na_currents = glob2.glob(path_to_cell + '/VC threshold adaptation/' + '*' + f"{na_rec:04d}" + ".abf")
    
                abf = pyabf.ABF(path_to_na_currents[0]) 
                fs = abf.dataRate  * Hz # sampling rate
                dt = 1./fs
                t = dt*arange(len(abf.sweepY)) 
                I = []
                V = []
                n_rec = len(abf.sweepList)                           
                
                for sweepNumber in range(n_rec):
                    abf.setSweep(sweepNumber)
                    I.append(abf.sweepY)
                    V.append(abf.sweepC*mV)             
                    
                ### RC time constant estimation from a test pulse

                path_to_test_pulses = glob2.glob(path_to_data + str(int(date)) + "*/" + \
                                                 '/retina '+ str(retina) +'/cell ' + str(int(cell)) + '/VC test pulses/comp/' +'*')
                tp_nums_correction = [path_to_test_pulses[i][-8:-4] for i in range(len(path_to_test_pulses))]
                idx_tp_correction = argmin(abs(int(tp_nums_correction[0]) - na_rec))
                tp_corr = tp_nums_correction[idx_tp_correction]
                print ('Tp correction:', tp_corr)
                
                # Loading the test pulse
                abf_tp = pyabf.ABF(path_to_test_pulses[idx_tp_correction])
                fs_tp = abf_tp.dataRate  * Hz # sampling rate
                dt_tp = 1./fs_tp
                t_tp = dt_tp*arange(len(abf_tp.sweepY)) 
                I_tp = []
                V_tp = []
                            
                for sweepNumber in abf_tp.sweepList:
                    abf_tp.setSweep(sweepNumber)
                    I_tp.append(abf_tp.sweepY)
                    V_tp.append(abf_tp.sweepC*mV)
                
                # Correcting the current for the passive component
                I_corr_pass, I_cut, t_cut = remove_passive_component(date, retina, cell, dt, I, V, dt_tp, I_tp, V_tp, str(int(na_rec)).zfill(4))  
            
                ### Tau RC estimation 
                tau_rc = []
                
                figure('Tau Rc estimation %i, %s, %i, %i' %(date, retina, cell, na_rec))
                ax1 = subplot(111)
                
                # Looping on all sweeps 
                for sweepNumber in abf_tp.sweepList:
                    abf_tp.setSweep(sweepNumber)
                    fs = abf_tp.dataRate  * Hz # sampling rate
                    dt_tp = 1./fs
                    t = dt_tp*np.arange(len(abf_tp.sweepY))/second
                    i = abf_tp.sweepY 
                    v = abf_tp.sweepC 
                    ax1.plot(t, i)
                                    
                    idx_start = int(56.2*ms/dt_tp)

                    # transient positive peak
                    idx_i_max = argmax(i[idx_start+5:idx_start+500]) + idx_start+5
                    i_amp_max_peak = i[idx_i_max] 
                    
                    ### exp fitting
                    idx_end = idx_i_max + int(1.*ms/dt_tp) 
                    i_amp_end = i[idx_end]
                    
                    max_amp = i_amp_end 
                    
                    ax1.plot(t/ms, i)
                    ax1.set_xlim(56, 64)
                    ax1.set_xlabel('t (s)')
                    
                    def exp_current(t, tau):
                        return  (i_amp_max_peak - max_amp) * exp(-t/tau) + max_amp
                    
                    idx_i_max = idx_i_max + 1
                    peak_current = i[idx_i_max: idx_end] # pA
                    peak_time = (t[idx_i_max:idx_end]- t[idx_i_max])*1e3  # sec
                    
                    ax1.plot(t[idx_i_max:idx_end]*1e3, peak_current, 'k-')
        
                    tau_opt = curve_fit(exp_current, peak_time, peak_current)
                    
                    ax1.plot(peak_time + t[idx_i_max]*1e3, \
                             exp_current(peak_time, tau_opt[0]), 'r-')
                    print ('Tau m:', tau_opt[0])
                    tau_rc.append(tau_opt[0])
                
                tau_rc_med = median(tau_rc)
                RC_time_cst.append(tau_rc_med)
                
                ### Prepulse potential
                prepulse = where(V[-1] == max(V[-1]))[0][0] - 100
                v0 = V[-1][prepulse]
                print ('Prepulse potential:', v0/mV + vh)
                prepulse_potentials.append(v0/mV + vh)
                
                ### The peak axonal current (just above threshold) is selected by eye
                ### by clicking in the corresping datapoint in the IV curve                
                coords = []
                def onclick(event):
                    global ix, iy
                    ix, iy = event.xdata, event.ydata
                    print ('x = %0.2f, y = %0.2f'%(ix, iy))
                
                    global coords
                    coords.append((ix, iy))
                
                    if len(coords) == 2:
                        fig_curr.canvas.mpl_disconnect(cid)
                        
                    return coords
                
                I_peaks = []
                Vc_peaks = []
                t_peaks = []
                idx_step = where(V[-1] == max(V[-1]))[0][0] - 1
                
                figpeaks = figure('Current traces %i,  %s, %i, %s' %(date, retina, cell, na_rec), (8,8))
                cmap = plt.get_cmap('gnuplot')
                cols = [cmap(i) for i in np.linspace(0, 1, n_rec)]
                
                subplot(211)
                for i in range(n_rec):
                    idx_peak = argmin(I_corr_pass[i])
                    i_peak = I_corr_pass[i][idx_peak]
                    t_peak = t_cut[idx_peak]
                    
                    I_peaks.append(i_peak)
                    Vc_peaks.append(V[i][idx_step+10])
                    t_peaks.append(t_peak)
                    #plotting
                    plot(t_cut/ms, I_corr_pass[i],  color=cols[i])
                    plot(t_peak/ms, i_peak, 'ko')
                    ylabel('I (pA)')
                    xlabel('Time (ms)')

                subplot(212)
                plot(Vc_peaks/mV, I_peaks, 'o-', color= 'k', label='peak') 
                ylabel('I peak (pA)')
                xlabel('V (mV)')
                legend(frameon=False)
                
                cid = figpeaks.canvas.mpl_connect('button_press_event', onclick)
                    
                waitforbuttonpress()
                    
                idx_peak_ax_current = argmin(sqrt((Vc_peaks/mV - coords[0][0])**2 + ((I_peaks - coords[0][1])*1e-3)**2))
                
                peak_axonal_currents.append(I_peaks[idx_peak_ax_current])
                peak_axonal_currents_raw.append(I_cut[idx_peak_ax_current][int(t_peaks[idx_peak_ax_current]/dt)]) # the raw current at peak time
                threshold_potentials.append(Vc_peaks[idx_peak_ax_current-1]/mV + vh)
                
                sweep_peak = idx_peak_ax_current
                selected_sweeps.append(sweep_peak)
                
                ### Measuring the charge and current duration
                # Plotting corrected current
                figure('Charge %i, %s, %i, %i' %(date, retina, cell, na_rec))
                ax1 = subplot(211)
                ax2 = subplot(212)

                ax1.plot(t_cut/ms, I_corr_pass[sweep_peak-1], 'k')
                ax1.plot(t_cut/ms, I_corr_pass[sweep_peak+1], 'k')
                ax1.plot(t_cut/ms, I_corr_pass[sweep_peak+2], 'k')
                ax1.plot(t_cut/ms, zeros(len(t_cut)), 'k--')
                ax1.set_ylabel('I (pA)')
                
                # Integration of the current
                curr_int = []
                curr_int_10 = []
                curr_duration_50 = []
                curr_duration_10 = []
                peak_lat = []
                for i in range(sweep_peak, sweep_peak+2): 
                    if I_peaks[i] > 0.2*I_peaks[idx_peak_ax_current]: # to avoid second current above threshold that has no peak
                        i += 1 # if it is the case we use the third peak above threshold
                    # interpolating the current
                    f = interpolate.interp1d(t_cut/ms, I_corr_pass[i], bounds_error=False, fill_value=0)
                    dt_new = 0.001
                    t_new = np.arange(0, 10, dt_new)
                    i_new = f(t_new)
                    ax2.plot(t_new, i_new, 'g-')
                    # current peak
                    idx_peak_new = argmin(i_new)
                    i_peak_new = i_new[idx_peak_new]
                    t_peak_new = t_new[idx_peak_new]
                    # current duration
                    idx_10 = where(abs(i_new) <= abs(0.1*i_peak_new))[0]
                    start_dur_10 = where(idx_10 < idx_peak_new)[0]
                    if len(start_dur_10) > 0:
                        idx_start_dur_10 = start_dur_10[-1]
                    else:
                        idx_start_dur_10 = 0
                    end_dur_10 = where(idx_10 > idx_peak_new)[0]
                    if len(end_dur_10) > 0:
                        idx_end_dur_10 = end_dur_10[0]
                    else:
                        idx_end_dur_10 = len(idx_10)-1
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
                    
                    idx_peak = argmin(I_corr_pass[i])    
                    
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
                    
                    curr_int_10.append(area_10)
                    curr_duration_10.append(dur_10)
                    curr_duration_50.append(dur_50)
                    peak_lat.append(idx_peak*dt/ms)
            
                tight_layout()
                
                current_integrals_10.append(curr_int_10)
                current_durations_10.append(curr_duration_10)
                current_durations_50.append(curr_duration_50)
                peak_current_latency.append(peak_lat)
                
                ### Correcting the peak axonal current
                # parameters
                v_rev = 70.
                si = dt/ms
                tauLag = si * 2.
                cm = tau_rc_med/rs_residual # nF
                print ('Cm:', cm)
                
                capacitance_tau.append(cm*1e3) #pF
                  
                # data
                idx_peak = argmin(I_corr_pass[idx_peak_ax_current][:int(9.*ms/dt)])
                i_data = I_corr_pass[idx_peak_ax_current][idx_peak - 30:idx_peak+20] * 1e-3 #nA #I_high[idx_ax][int(59.*ms/dt):int(60.*ms/dt)]*1e-3 # nA
                vc_data = Vc_peaks[idx_peak_ax_current]/mV + vh 
                
                figure('Rs correction %i, %s, %i, %i, %0.02f' %(date, retina, cell, na_rec, tau_rc_med), figsize=(12,6))
                
                subplot(211)
                plot(i_data, 'k-')
                
                # correction
                I_cap = zeros(len(i_data))
                I_corr_cap = zeros(len(i_data))
                I_corr = zeros(len(i_data))
                
                n_pts = len(i_data)
                
                # first data point
                v_last = vc_data - i_data[0] * rs_residual
                denom = v_last - v_rev
                
                fracCorr = (1-(vc_data-v_rev)/denom)

                I_corr[0] = i_data[0] * (1-fracCorr)
                
                # next data points
                for j in range(n_pts):
                    v_this = vc_data - i_data[j]*rs_residual
  
                    fracCorr = (1-(vc_data-v_rev)/(v_this-v_rev))

                    Icap = cm * (v_this-v_last)/si #* 1e-3
                    I_cap[j] = Icap
                    I_corr_cap[j-1] = i_data[j-1] - Icap
                    I_corr[j-1] = I_corr_cap[j-1] * (1-fracCorr)
                
                    v_last = v_this 
                
                plot(I_corr, 'r')
                
                subplot(212)
                plot(I_cap, 'b')
            
                show()
                
                peak_axonal_currents_corr.append(min(I_corr))
                print ('-------------------------')
        else:
            print('No adaptation protocol')

# ### Write in excel file

# df_select_cells = pd.DataFrame({'Date': selected_dates,
#                   'Retina': selected_retinas,
#                   'Cell': selected_cells,
#                   'Age': selected_ages,
#                   'Recording': selected_recs,
#                   'Sweep': selected_sweeps,
#                   'Vh': holding_potentials,
#                   'V prepulse': prepulse_potentials,
#                   'TP before': tp_before,
#                   'TP after': tp_after,
#                   'TP Na rec': tp_na_rec,
#                   'Rs before': series_resistance_before,
#                   'Rs after': series_resistance_after,
#                   'Rs Na rec': series_resistance_na_rec,
#                   'Compensation rec': compensation_rec,
#                   'Rs rec': series_resistance_rec,
#                   'Rs compensation': series_resistance_compensation,
#                   'Rs residual': series_resistance_residual,
#                   'Tau RC': RC_time_cst,
#                   'Capacitance': capacitance_tau,
#                   'Peak axonal current': peak_axonal_currents,
#                   'Peak axonal current corrected': peak_axonal_currents_corr,
#                   'Peak axonal current raw': peak_axonal_currents_raw,
#                   'Vth': threshold_potentials,
#                   'Charge1 10': [current_integrals_10[i][0] for i in range(len(selected_dates))],
#                   'Charge2 10': [current_integrals_10[i][1] for i in range(len(selected_dates))],
#                   'Peak latency1': [peak_current_latency[i][0] for i in range(len(selected_dates))],
#                   'Peak latency2': [peak_current_latency[i][1] for i in range(len(selected_dates))],
#                   'Duration1 10': [current_durations_10[i][0] for i in range(len(selected_dates))],
#                   'Duration2 10': [current_durations_10[i][1] for i in range(len(selected_dates))],
#                   'Duration1 50': [current_durations_50[i][0] for i in range(len(selected_dates))],
#                   'Duration2 50': [current_durations_50[i][1] for i in range(len(selected_dates))]
#                   })

# df_select_cells.to_excel(path_to_files + "RGC_adaptation_2905_8.xlsx", \
#                 columns=['Date','Retina','Cell','Age','Recording', 'Sweep',\
#                   'Vh',
#                   'V prepulse',
#                   'TP before',
#                   'TP after',
#                   'TP Na rec',
#                   'Rs before',
#                   'Rs after',
#                   'Rs Na rec',
#                   'Compensation rec',
#                   'Rs rec',
#                   'Rs compensation',
#                   'Rs residual',
#                   'Tau RC',
#                   'Capacitance',
#                   'Peak axonal current',
#                   'Peak axonal current corrected',
#                   'Peak axonal current raw',
#                   'Vth',
#                   'Charge1 10',
#                   'Charge2 10',
#                   'Peak latency1',
#                   'Peak latency2',
#                   'Duration1 10',
#                   'Duration2 10',
#                   'Duration1 50',
#                   'Duration2 50'])
     
                