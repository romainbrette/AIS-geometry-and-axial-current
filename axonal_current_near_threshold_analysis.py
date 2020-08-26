

"""

Automated measure of Na currents below threshold.


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

### Loading the list of cells used in the analysis
path_files = '/Users/sarah/Documents/repositories/AIS-geometry-and-axonal-current/Na currents in RGC/codes submission/data/'
df_cells = pd.read_excel(path_files + 'RGC_electrical_properties.xlsx')

first_cell = -2
last_cell = len(df_cells['Date'])

dates = array(df_cells['Date'])[:-3][first_cell:last_cell]
retinas = array(df_cells['Retina'])[:-3][first_cell:last_cell]
cells = array(df_cells['Cell'])[:-3][first_cell:last_cell]
ages = array(df_cells['Age'])[:-3][first_cell:last_cell]
v_holding = array(df_cells['V holding (mV)'])[:-3][first_cell:last_cell]

# # Load the number of the good quality recordings
# df_gc = pd.read_excel(path_to_files + 'good_recordings.xlsx')

### Path to the data
path_to_data = '/Users/sarah/Documents/Data/Martijn Sierksma/'

### Keeping track of analyzed cells
selected_dates = []
selected_retinas = []
selected_cells = []

IV_curves_below_I_all = []
IV_curves_below_V_all = []
v_command_all = []

threshold_current_smoothed = []
threshold_current = []

N = 0 #number of analyzed cells

for date, retina, cell, age, vh in zip(dates, retinas, cells, ages, v_holding): 
    # only cells with P5 protocol are analyzed
    if date > 20190611 and age > 9: 
        print ('------------------------------')
        print (date, retina, cell, 'TPs:')
        
        path_to_cell = path_to_data + str(int(date)) + "*/" + '/retina '+ str(retina) +'/cell ' + str(int(cell))
        
        ### Path to the Na current recordings
        row_rec = df_gc[(df_gc['Date'] == date) & (df_gc['Retina'] == retina) & (df_gc['Cell'] == cell)]
        cell_recs = [row_rec['1'].values[0], row_rec['2'].values[0], row_rec['3'].values[0], row_rec['4'].values[0]]
        cell_recs = array(cell_recs)[~isnan(cell_recs)]
        print ('Na rec:', cell_recs)
        
        if len(cell_recs) > 0:
            
            selected_dates.append(date)
            selected_retinas.append(retina)
            selected_cells.append(cell)
                
            IV_curves_below_I = []
            IV_curves_below_V = []
            v_command = []
            thres_curr = []
            thres_curr_smoothed = []
            
            for rec in cell_recs:
                                        
                path_to_na_currents = glob2.glob(path_to_cell + '/VC small steps/20' + '*' + str(int(rec)).zfill(4) + ".abf")
                    
                N += 1 # counting cells used on the analysis
                
                ### Loading and plotting Na currents
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
                    
                ### Correcting for the passive component of the current
                I_corr, I_cut, t_cut_long = p5_subtraction(date, retina, cell, dt, I, V, rec_name=str(int(rec)).zfill(4))
                I_corr_pass = [I_corr[i][:int(20.*ms/dt)] for i in range(n_rec)]
                t_cut = t_cut_long[:int(20.*ms/dt)]
                
                ### IV curve
                I_peaks,  Vc_peaks, idx_peak_ax_current, t_peaks = plot_IV(date, retina, cell, dt, I_corr_pass, V,0, str(int(rec)).zfill(4))
                              
                ### Smoothing current traces below threshold
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
                
                plot(t_cut/ms, baseline_peak_current_smoothed * ones(len(t_cut)), '--', color=cols[0])
                plot(t_cut/ms, baseline_peak_current * ones(len(t_cut)), color=cols[0])
                
                Vc_peaks = array( Vc_peaks)*1e3*mV + vh * mV
                I_peaks = array(I_peaks)
                
                ### Threshold current
                
                if len(I_peaks[:idx_peak_ax_current]) > 1. :                                       
                    I_peaks_below_smoothed = (I_peaks_smoothed - baseline_peak_current_smoothed) * 1e-3 #nA
                    ith_smoothed = min(I_peaks_below_smoothed[idx_peak_ax_current-3:idx_peak_ax_current])
                    thres_curr_smoothed.append(ith_smoothed)
                    
                    I_peaks_below = (I_peaks - baseline_peak_current) * 1e-3 #nA
                else:
                    thres_curr_smoothed.append(nan)
                
                ### IV curve below threshold
                if len(I_peaks[:idx_peak_ax_current]) > 4. :
                    
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
                    v_command.append(Vc_peaks/mV)
                else:
                    IV_curves_below_I.append(nan)
                    IV_curves_below_V.append(nan)
                    v_command.append(nan)
                
            IV_curves_below_I_all.append(array(IV_curves_below_I))
            IV_curves_below_V_all.append(array(IV_curves_below_V))
            v_command_all.append(array(v_command))
            threshold_current_smoothed.append(min(thres_curr_smoothed))
                    
# savez('RGC_IV_curves_below_threshold', selected_dates, selected_retinas, selected_cells, \
#       IV_curves_below_I_all, IV_curves_below_V_all, v_command_all)



# df_select_cells = pd.DataFrame({'Date': selected_dates,
#                   'Retina': selected_retinas,
#                   'Cell': selected_cells,
#                   'Threshold current (nA)': threshold_current_smoothed})

# df_select_cells.to_excel(path_to_files + "RGC_threshold_current_noP8_1606.xlsx", \
#                           columns=['Date','Retina','Cell','Threshold current (nA)'])
        
        











        
