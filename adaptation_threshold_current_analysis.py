

"""

Threshold current adaptation


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

### Loading the results of analyses
path_files = '/Users/sarah/Documents/repositories/AIS-geometry-and-axial-current/'

# Load the list of the cells that will be used for analysis
cells_pp = pd.read_excel(path_files + 'recordings_for_threshold_current_adaptation.xlsx') # cells already sorted for Rs etc

# # Load the recordings that will be used for the analysis
# df_adapt = pd.read_excel(path_to_files + 'good_recordings_adaptation_all_v2.xlsx')

first_cell = 8
last_cell = 10 #len(cells_pp['Date'])

dates = array(cells_pp['Date'])[first_cell:last_cell]
retinas = array(cells_pp['Retina'])[first_cell:last_cell]
cells = array(cells_pp['Cell'])[first_cell:last_cell]
ages = array(cells_pp['Age'])[first_cell:last_cell]
v_holding = array(cells_pp['V holding (mV)'])[first_cell:last_cell]

#11 Path to the data
path_to_data = '/Users/sarah/Documents/Data/Martijn Sierksma/'

selected_dates = []
selected_retinas = []
selected_cells = []
selected_ages = []
selected_sweeps = []

peak_axonal_currents = []
threshold_current_smoothed = []
holding_potentials = []
threshold_potentials = []

dates_iv = []
retinas_iv = []
cells_iv = []
IV_curves_below_I_all = []
IV_curves_below_V_all = []
v_prepulse_iv_all = []

N = 0

for date, retina, cell, age, vh in zip(dates, retinas, cells, ages, v_holding): 
    if date > 20190611: # only with P5:
        print ('------------------------------')
        print (date, retina, cell)
           
        path_to_cell = path_to_data + str(int(date)) + "*/" + '/retina '+ str(retina) +'/cell ' + str(int(cell))
        
        # get the sweeps number
        row_cell = df_adapt[(df_adapt['Date'] == date) & (df_adapt['Retina'] == retina) & (df_adapt['Cell'] == cell)]
        
        if len(row_cell['s1']) > 0:
            sweeps = row_cell.values[0][4:]
            
            # delete nans
            sweeps = sweeps[~pd.isnull(sweeps)]
            print (sweeps)
            
            N += 1 # counting cells used on the analysis
            
            # recording data for IV curve analysis
            dates_iv.append(date)
            retinas_iv.append(retina)
            cells_iv.append(cell)
            
            IV_curves_below_I = []
            IV_curves_below_V = []
            IV_curve_v0 = []
            v_prepulse_iv = []
            
            # loop on sweeps corresponding to different holding potentials
            for rec in sweeps:
                i_rec = int(rec)
                print (f"{i_rec:04d}")
                
                selected_dates.append(date)
                selected_retinas.append(retina)
                selected_cells.append(cell)
                selected_ages.append(age)
                selected_sweeps.append(f"{i_rec:04d}")
        
                # path to the Na current recordings
                path_to_na_currents = glob2.glob(path_to_cell + '/VC threshold adaptation/' + '*' + f"{i_rec:04d}" + ".abf")
                 
                # Loading and plotting Na currents
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
                
                # Holding potential
                prepulse = where(V[-1] == max(V[-1]))[0][0] - 100
                v0 = V[-1][prepulse]
                print (v0/mV)
                holding_potentials.append(v0/mV)
                
                ### Correcting for the passive component of the current        
                I_corr_pass, I_cut, t_cut = p5_subtraction(date, retina, cell, dt, I, V, f"{int(rec):04d}")

                ### Measuring IV curve: peak axonal current and threshold (by eye)
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
                
                figpeaks = figure('Current traces %i,  %s, %i, %s' %(date, retina, cell, rec), (8,8))
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
                #plot(vc_peaks[idx_peak], i_peaks[idx_peak], 'ro') 
                #plot(vc_peaks, i_peaks_amp, 'o-', color= 'green', label='amplitude') 
                ylabel('I peak (pA)')
                xlabel('V (mV)')
                legend(frameon=False)
                
                cid = figpeaks.canvas.mpl_connect('button_press_event', onclick)
                    
                waitforbuttonpress()
                    
                idx_peak_ax_current = argmin(sqrt((Vc_peaks/mV - coords[0][0])**2 + ((I_peaks - coords[0][1])*1e-3)**2))
                                
                peak_axonal_currents.append(I_peaks[idx_peak_ax_current])
                threshold_potentials.append(Vc_peaks[idx_peak_ax_current-1]/mV + vh)
                
                sweep_peak = idx_peak_ax_current
                # selected_sweeps.append(sweep_peak)
                                
                ### Smoothing traces
                figure('Traces below threshold %i,  %s, %i, %i' %(date, retina, cell, rec))
                
                I_peaks_smoothed = zeros(idx_peak_ax_current)
                I_corr_smoothed = []
    
                cmap = plt.get_cmap('gnuplot')
                cols = [cmap(i) for i in np.linspace(0, 1, idx_peak_ax_current)]
                
                for i in range(idx_peak_ax_current):
                    
                    #smoothing
                    n = len(t_cut/ms)
                    i_slide = np.zeros(n)
                    d = 50 # half-window, i.e. number of pixels on each side
                
                    for j in range(n):
                        if j < d: # start of the axon, not full window
                            i_slide[j] = np.mean(I_corr_pass[i][0:j+d])
                        elif j > n-d: # end of the axon, not full window
                            i_slide[j] = np.mean(I_corr_pass[i][j-d:n])
                        else: 
                            i_slide[j] = np.mean(I_corr_pass[i][j-d:j+d])
                    
                    I_corr_smoothed.append(i_slide)
                    I_peaks_smoothed[i] = min(i_slide)
                    
                    plot(t_cut/ms, I_corr_pass[i], '-', color = cols[i])
                    plot(t_cut/ms, i_slide, 'k')
                
                baseline_peak_current = mean(I_corr_pass[0])
                baseline_peak_current_smoothed = mean(I_corr_smoothed[0])
                
                # plot(t_cut/ms, baseline_peak_current_smoothed * ones(len(t_cut)), '--', color=cols[0])
                plot(t_cut/ms, baseline_peak_current * ones(len(t_cut)), color=cols[0])
                
                Vc_peaks = array( Vc_peaks)*1e3*mV + vh * mV
                I_peaks = array(I_peaks)
                
                ### Threshold current
                
                if len(I_peaks[:idx_peak_ax_current]) > 1. :
                                        
                    I_peaks_below_smoothed = (I_peaks_smoothed - baseline_peak_current_smoothed) * 1e-3 #nA
                    if len(I_peaks[:idx_peak_ax_current]) > 3. :
                        ith_smoothed = min(I_peaks_below_smoothed[idx_peak_ax_current-3:idx_peak_ax_current])
                    else:
                        ith_smoothed = min(I_peaks_below_smoothed[:idx_peak_ax_current])
                    threshold_current_smoothed.append(ith_smoothed)
                    
                    I_peaks_below = (I_peaks - baseline_peak_current) * 1e-3 #nA
                else:
                    threshold_current_smoothed.append(nan)
                
                ### IV curve below threshold
                if len(I_peaks[:idx_peak_ax_current]) > 3. :
                    
                    Vc_peaks_below = Vc_peaks[:idx_peak_ax_current]/mV # mV
                    
                    I_peaks_below = (I_peaks[:idx_peak_ax_current] - baseline_peak_current) * 1e-3 #nA
                    I_peaks_below_smoothed = (I_peaks_smoothed - baseline_peak_current_smoothed) * 1e-3 #nA
    
                    f2 = figure('IV %i,  %s, %i, %i' %(date, retina, cell, rec), (6,5))
                    ax4 = f2.add_subplot(111)
                
                    ax4.plot(Vc_peaks_below, I_peaks_below, 'k-o')
                    ax4.plot(Vc_peaks_below, I_peaks_below_smoothed, 'r-o', label='smoothed')
                    ax4.legend(frameon=False)
                                    
                    IV_curves_below_I.append(I_peaks_below_smoothed)
                    IV_curves_below_V.append(Vc_peaks_below)
                    v_prepulse_iv.append(v0/mV)
                    # v_command.append(Vc_peaks/mV)
                else:
                    IV_curves_below_I.append(nan)
                    IV_curves_below_V.append(nan)
                    v_prepulse_iv.append(nan)
                    # v_command.append(nan)
                
            IV_curves_below_I_all.append(array(IV_curves_below_I))
            IV_curves_below_V_all.append(array(IV_curves_below_V))
            v_prepulse_iv_all.append(array(v_prepulse_iv))
                    
        else:
            print('No adaptation protocol')

# savez('RGC_IV_curves_below_threshold_adaptation_0107_3', \
#       dates_iv, retinas_iv, cells_iv, \
#       IV_curves_below_I_all, IV_curves_below_V_all, v_prepulse_iv_all)

# #Write in excel file

# df_select_cells = pd.DataFrame({'Date': selected_dates,
#                   'Retina': selected_retinas,
#                   'Cell': selected_cells,
#                   'Age': selected_ages,
#                   'Sweep': selected_sweeps,
#                   'V0': holding_potentials,
#                   'Peak current': peak_axonal_currents,
#                   'Threshold current': threshold_current_smoothed,
#                   'Vth': threshold_potentials,
#                   })

# save_path = '/Users/sarah/Dropbox/Spike initiation/PhD projects/Axonal current and AIS geometry/Data patch/RGC/'
# df_select_cells.to_excel(save_path + "RGC_threshold_current_adaptation_0107_3.xlsx", 
#                 columns=['Date','Retina','Cell','Age','Sweep',\
#                           'V0','Peak current', 'Threshold current', 'Vth'])
