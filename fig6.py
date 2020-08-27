
"""

Figure 6: axial current near threshold.


"""

from brian2 import *
from pandas import ExcelWriter
from pandas import ExcelFile
from scipy import stats
from matplotlib import gridspec
from na_currents_analysis import *
from vc_test_pulse_analysis import *
import params_model_description
import pandas as pd
import seaborn as sns
import glob2
import pyabf

### Figure parameters
rcParams['axes.spines.right'] = False
rcParams['axes.spines.top'] = False

### Model parameters
params = params_model_description

### Loading the results of analyses
path_files = '/Users/sarah/Documents/repositories/AIS-geometry-and-axial-current/'
df_cells = pd.read_excel(path_files + 'RGC_electrical_properties.xlsx')
df_pass = pd.read_excel(path_files + 'RGC_passive_properties.xlsx')

# we select one cell as an exmaple for panel A
first_cell = -5
last_cell = -4 #len(df_cells['Date'])

# cells names and age
dates = array(df_cells['Date'])[first_cell:last_cell]
retinas = array(df_cells['Retina'])[first_cell:last_cell]
cells = array(df_cells['Cell'])[first_cell:last_cell]
ages = array(df_cells['Age'])[first_cell:last_cell]
# recording used to measure axial current
na_recs = array(df_cells['Recording'])[first_cell:last_cell]
# TP with compensation to correct for passive component for cells w/o P5
tp_corrs = array(df_cells['TP num correction'])[first_cell:last_cell]
#t_starts = array(df_cells['t start'])[first_cell:last_cell]
# residual Rs
rs_residuals = array(df_cells['Residual Rs'])[first_cell:last_cell]
# holding membrane potential at rest
resting_mp = array(df_cells['V holding (mV)'])[first_cell:last_cell]
# TP w/o compensation used to measure Cm and Rs
tp_nums = array(df_cells['TP num passive props'])[first_cell:last_cell]

### Path to the data
path_to_data = '/Users/sarah/Documents/Data/Martijn Sierksma/'

### Figure
fig = figure('Threshold and current threshold', figsize=(6,7))

### Panel A: current near threshold: example

for date, retina, cell, age, na_rec, tp_corr, rs_res, vh in zip(dates, retinas, cells, ages, na_recs, tp_corrs, rs_residuals, resting_mp): 
    print ('------------------------------')
    print (date, retina, cell)
    
    ### Load Na current recording
    na_rec = na_rec+1 # to use the same rec as for other panels
    path_to_cell = path_to_data + str(int(date)) + "*" + '/retina '+ str(retina) +'/cell ' + str(int(cell))   
    path_to_na_currents = glob2.glob(path_to_cell + '/VC small steps/' + '*' + str(int(na_rec)).zfill(4) + ".abf")
             
    abf = pyabf.ABF(path_to_na_currents[0]) 
    fs = abf.dataRate  * Hz # sampling rate
    dt = 1./fs
    print('50dt',50*dt)
    t = dt*arange(len(abf.sweepY)) 
    I = []
    V = []
    n_rec = len(abf.sweepList)
    
    cmap = plt.get_cmap('gnuplot')
    cols = [cmap(i) for i in np.linspace(0, 1, n_rec)]

    for sweepNumber in range(n_rec):
        abf.setSweep(sweepNumber)
        I.append(abf.sweepY)
        V.append(abf.sweepC*mV)
    
    ### Correction for passive component
    I_corr_pass, I_cut, t_cut = p5_subtraction(date, retina, cell, dt, I, V, rec_name=str(int(na_rec)).zfill(4))

    ### IV curve
    I_peaks,  Vc_peaks, idx_peak_ax_current, t_peaks = plot_IV(date, retina, cell, dt, I_corr_pass, V, 0, str(int(na_rec)).zfill(4))
    
    # IV surve below threshold: smoothing
    I_peaks_smoothed = zeros(idx_peak_ax_current)
    I_corr_smoothed = []
    
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
        
        I_peaks_smoothed[i] = min(i_slide)
        I_corr_smoothed.append(i_slide)
    
baseline_peak_current_smoothed = mean(I_corr_smoothed[0])

ax1 = fig.add_subplot(321)

i_thres = (I_corr_smoothed[-1] - baseline_peak_current_smoothed)*1e-3
idx_it = argmin(i_thres)

ax1.plot(t_cut/ms, I_corr_pass[idx_peak_ax_current-1]*1e-3, 'gray', label='raw')
ax1.plot(t_cut/ms, i_thres , 'k', label='smoothed')
ax1.plot(linspace(t_cut[0]/ms, t_cut[idx_it]/ms, 10), i_thres[idx_it]*ones(10), '--', color='k', linewidth=1)
ax1.set_yticks(ticks = [0, i_thres[idx_it], -0.1])
ax1.set_yticklabels(['0', '$I_t$', '-0.1']) 
ax1.set_xlim(0, 10)
ax1.set_ylim(-0.15, 0.1)
# ax1.set_xlabel('t (ms)')
ax1.set_ylabel('I (nA)')
ax1.legend(loc='lower right', bbox_to_anchor=(1.05,-0.075), frameon=False, fontsize=8)
ax1.set_xticks([])
sns.despine(bottom=True, ax=ax1)
ax1.plot(linspace(1, 2, 10), -0.147*ones(10), 'k-', linewidth=2)
ax1.text(0.8, -0.17, '1 ms', color='k', fontsize=8)
ax1.annotate("A", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

### Panel C - E - F: example of IV curves below threshold in RGC

# Laoding IV curves data for all cells
data = load('RGC_IV_curves_below_threshold.npz', allow_pickle=True)
dates = data['arr_0']
retinas = data['arr_1']
cells = data['arr_2']
currents = data['arr_3'] # below threshold
voltages = data['arr_4'] # below threshold
v_com = data['arr_5'] # all command potentials

n_cells = len(dates)

# For panel F, we analyze the IV curves below threshold of all cells
# we fit the theoretical prediction to estimate beta
selected_dates = []
selected_retinas = []
selected_cells = []

selected_currents = []
selected_voltages = []

slope_linreg = []
inter_linreg = []
i_threshold = []

# Loop on the selected cells
for i in range(n_cells):
    row = df_cells[(df_cells['Date'] == dates[i]) & (df_cells['Retina'] == retinas[i]) & (df_cells['Cell'] == cells[i])]
    
    if len(row['Date'].values) != 0:
        selected_dates.append(dates[i])
        selected_retinas.append(retinas[i])
        selected_cells.append(cells[i])
        selected_currents.append(currents[i])
        selected_voltages.append(voltages[i])
        
        n_rec = len(currents[i])
        
        if n_rec == 1: # cells we only one currents recording
            slopes = []
            inter = []
            ith = []
            
            fig_lr = figure('Lin reg %i %s %i' %(dates[i], retinas[i], cells[i]), figsize=(5,4))
            subplot(111)
            # fit if possible
            try:
                idx_thres = argmin(currents[i][0][-3:]) + len(currents[i][0][:-3])
                i_thres = currents[i][0][idx_thres] 
                v_thres = voltages[i][0][idx_thres]
                
                slope, intercept, _, _, _ = stats.linregress( \
                        - (currents[i][0][idx_thres-5:idx_thres]/i_thres -1.)**2, \
                        voltages[i][0][idx_thres-5:idx_thres] - v_thres)
                if 0 < slope/2 < 20: 
                    slopes.append(slope)
                    ith.append(i_thres)
                    inter.append(intercept)
                else:
                    slopes.append(nan)
                    ith.append(nan)
                    inter.append(nan)
                    
                title('%i %s %i' %(dates[i], retinas[i], cells[i]))
                ii  = linspace(min(- (currents[i][0]/currents[i][0][-1] -1.)**2) + min(- (currents[i][0]/currents[i][0][-1] -1.)**2), \
                            max(- (currents[i][0]/currents[i][0][-1] -1.)**2) - min(- (currents[i][0]/currents[i][0][-1] -1.)**2), 50)
                plot(ii, intercept+slope*ii, '-', color='k', label=' $k_a$=%.03f mV' %(2*slope))
                plot(- (currents[i][0]/currents[i][0][-1] -1.)**2, voltages[i][0] -  voltages[i][0][-1] , \
                      'o', color='k', label='$I_t = %0.03f$ nA' %i_thres)
                legend(frameon=False)
                
                ylim(-15,5)
                xlim(-2,0.5)
                xlabel('$-0.5(I/I_t-1)^2$')
                ylabel('$V_s-V_s^*$ (mV)')
                tight_layout()
                
            except:
                slopes.append(nan)
                ith.append(nan)
                inter.append(nan)
                pass
            
            slopes.append(nan) # no second rec
            ith.append(nan)
            inter.append(nan)
                
            slope_linreg.append(slopes)
            i_threshold.append(ith)
            inter_linreg.append(inter)
    
        elif n_rec == 2: # cells with two recordings
            fig_lr = figure('Lin reg %i%s%i' %(dates[i], retinas[i], cells[i]), figsize=(9,4))
            
            slopes = []
            ith = []
            inter = []
            
            for j in range(2):
                subplot(1, 2, j+1)
                # fit if possible
                try:
                    idx_thres = argmin(currents[i][j][-3:]) + len(currents[i][j][:-3])
                    i_thres = currents[i][j][idx_thres]
                    v_thres = voltages[i][j][idx_thres]
                    
                    slope, intercept, _, _, _ = stats.linregress( \
                        - (currents[i][j][idx_thres-5:idx_thres]/i_thres -1.)**2, \
                        voltages[i][j][idx_thres-5:idx_thres] - v_thres)
                    
                    if 0 < slope/2 < 20:
                        slopes.append(slope)
                        ith.append(i_thres)
                        inter.append(intercept)
                    else:
                        slopes.append(nan)
                        ith.append(nan)
                        inter.append(nan)
                        
                    ii  = linspace(min(- (currents[i][j]/currents[i][j][-1] -1.)**2) + min(-(currents[i][j]/currents[i][j][-1] -1.)**2), \
                                max(- (currents[i][j]/currents[i][j][-1] -1.)**2) - min(- (currents[i][j]/currents[i][j][-1] -1.)**2), 50)
                    plot(ii, intercept+slope*ii, '-', color='k', label=' $k_a$=%.03f mV' %(2*slope))
                    title('%i %s %i' %(dates[i], retinas[i], cells[i]))
                
                    plot(- (currents[i][j]/currents[i][j][-1] -1.)**2, voltages[i][j] -  voltages[i][j][-1] ,\
                          'o', color='k', label='$I_t = %0.03f$ nA' %i_thres)
                    legend(frameon=False)
                    ylim(-15,5)
                    xlim(-1.5,0.5)
                    xlabel('$-0.5(I/I_t-1)^2$')
                    ylabel('$V_s-V_s^*$ (mV)')
                except:
                    slopes.append(nan)
                    ith.append(nan)
                    inter.append(nan)
                    pass
    
            slope_linreg.append(slopes)
            i_threshold.append(ith)
            inter_linreg.append(inter)
    
            tight_layout()
            
### For the cells with two recordings, we select the one used for the IV curve below threshold analysis

slope_final = []
inter_final = []
i_threshold_final = []
dates_final = []
retinas_final = []
cells_final = []
IV_i_final = []
IV_v_final = []

for i in range(len(slope_linreg)):
    if slope_linreg[i][0] < 1 and slope_linreg[i][1] < 1: 
        print ('too noisy IV curve')
    elif slope_linreg[i][0] != slope_linreg[i][0] and slope_linreg[i][1] != slope_linreg[i][1]:
        print ('removed nans')
    else:
        if slope_linreg[i][0] < 1 :
            idx_max = 1
        elif slope_linreg[i][1] < 1 :
            idx_max = 0
        else:
            idx_max = nanargmin(i_threshold[i])
            
        slope_final.append(slope_linreg[i][idx_max])
        inter_final.append(inter_linreg[i][idx_max])
        i_threshold_final.append(i_threshold[i][idx_max])
        dates_final.append(selected_dates[i])
        retinas_final.append(selected_retinas[i])
        cells_final.append(selected_cells[i])
        IV_i_final.append(selected_currents[i][idx_max])
        IV_v_final.append(selected_voltages[i][idx_max])
        

### We select an exmaple cell for panels C and E 
c = -2

### Panel D: IV curve below threshold in an example cell
ax4 = fig.add_subplot(323)
ax4.plot(selected_voltages[c][1], selected_currents[c][1], '.', color='gray')
ax4.fill_between(linspace(selected_voltages[c][1][-1]-slope_linreg[c][1], selected_voltages[c][1][-1],50), \
                 zeros(50), selected_currents[c][1][-1]*ones(50), color='gray', alpha = 0.2)
ax4.set_xlabel('$V$ (mV)')
ax4.set_ylabel('$I$ (nA)')
ax4.set_ylim(-0.075,0.01)
ax4.set_xlim(-66,-56)
ax4.annotate("C", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

# Panel E: fitting the IV curve below threshold in an exmaple cell
ax5 = fig.add_subplot(325) 
ii  = linspace(-1.2,0.1, 50)
ax5.plot(ii, inter_linreg[c][1] + slope_linreg[c][1] * ii, '-',color='gray')
ax5.plot(- (selected_currents[c][1]/selected_currents[c][1][-1] -1.)**2, selected_voltages[c][1] -  selected_voltages[c][1][-1],'.',color='gray')#, label='$I_t = %0.03f$ nA' %i_thres)
ax5.fill_between(linspace(-1,0,50), -slope_linreg[c][1]*ones(50), zeros(50),\
                 color='gray', alpha = 0.2)
ax5.set_xlabel('$-(I/I_t-1)^2$')
ax5.set_ylabel('$V-V_t$ (mV)')
ax5.set_ylim(-9, 1)
ax5.set_xlim(-1.5, 0.2)
ax5.legend(frameon=False, fontsize=8)
ax5.annotate("E", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

### Panel F: distribution of beta accross cells
ax6 = fig.add_subplot(326)
sns.boxplot(y=slope_final, color='gray', width=0.3, ax=ax6)
sns.swarmplot(y=slope_final,  color='0.2', ax=ax6)
sns.despine(bottom=True, ax=ax6)
ax6.set_ylabel(r'$\beta$ (mV)')
ax6.set_xticks([])
ax6.set_ylim(0, 10)
ax6.annotate("F", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

### Panel B and D: IV curve below threshold in a model

### Model parameters
dt = 0.01*ms
v_h = params.EL # leak potential
start = 10.*um # AIS start position
length = 30.*um # AIS length
ra = (4*params.Ri/(pi*params.axon_diam**2)) # axial resistance per unit length
Ras =  ra * start # axial resistance between soma and AIS start

### Load and plot data from simulations
path_to_data = '/Users/sarah/Documents/repositories/AIS-geometry-and-axial-current/simulations data/fig6/'
dir_name = path_to_data + 'VC dicho APmodel ext AIS x%0.1f L%i' %(start/um, length/um)

Ie = loadtxt(dir_name +'/Steps/I.txt')
Vm = loadtxt(dir_name +'/Steps/V.txt')
Vc = loadtxt(dir_name + '/Steps/Vc.txt')
Im = loadtxt(dir_name+'/Steps/Im.txt')

n_rec = len(Ie)
t = arange(len(Ie[0]))*dt

### Load test pulse data
dir_name_tp = path_to_data + 'TP APmodel ext AIS x%0.1f L%i' %(start/um, length/um)

Ie_tp = loadtxt(dir_name_tp +'/Steps/I.txt')
Ve_tp = loadtxt(dir_name_tp +'/Steps/V.txt')
Vc_tp = loadtxt(dir_name_tp + '/Steps/Vc.txt')

t_tp = arange(len(Ve_tp))*dt
vc_tp_amp = -5
    
### Remove passive response from Na currents recordings
vc = Vc # command potentials

steps_start = int(20.*ms/dt)
steps_end = int(40.*ms/dt)

i_tp_cut = Ie_tp[steps_start:steps_end] - mean(Ie_tp[int(18.*ms/dt):int(19.5*ms/dt)])

Ie_corr = []

for k in range(n_rec):
    V_amp =  vc[k] - v_h/mV
    factor = V_amp/vc_tp_amp
    i_cut = Ie[k][steps_start:steps_end] - mean(Ie[k][int(18.*ms/dt):int(19.5*ms/dt)])
    t_cut = t[steps_start:steps_end]
        
    i_corr =  i_cut - i_tp_cut * factor
    Ie_corr.append(i_corr)
    
### Measure peak axonal current and threshold
spikes = zeros(n_rec)
i_peaks = []
peaks_indexes = []
for j in range(n_rec):
    # peak current
    idx_peak = argmin(Ie_corr[j][int(0.1*ms/dt):int(19.*ms/dt)]) + int(0.1*ms/dt)
    i_peak = Ie_corr[j][idx_peak]
    i_peaks.append(i_peak)
    peaks_indexes.append(idx_peak)
    if i_peak < -1.9:
        spikes[j] = 1
        
vc = array(vc)
i_peaks = array(i_peaks)
I_corr = array(Ie_corr)
spikes = array(spikes)
peaks_indexes = array(peaks_indexes)

idx_spikes = where(spikes == 1)
idx_no_spikes = where(spikes == 0)

vc_below = vc[idx_no_spikes]
i_below = i_peaks[idx_no_spikes]

idx_sort = argsort(vc_below)
IV_i_below = i_below[idx_sort] #- i_below[0]
IV_vc_below = vc_below[idx_sort]

### Panel B: IV curve below threshold in a model
    
ka_model = params.Ka

ii  = linspace(-1.2,0.1, 50)

slope_mod, intercept_mod, _, _, _ = stats.linregress( \
                        - (IV_i_below[-5:]/IV_i_below[-1] -1.)**2, \
                        IV_vc_below[-5:] - IV_vc_below[-1])

ax7 = fig.add_subplot(322)
ax7.plot(IV_vc_below, IV_i_below, '.', color='g')
ax7.fill_between(linspace(IV_vc_below[-1]-slope_mod, IV_vc_below[-1], 50), \
                 zeros(50), IV_i_below[-1]*ones(50), color='g', alpha = 0.2)
ax7.set_xlabel('$V$ (mV)')
ax7.set_ylabel('$I$ (nA)')
ax7.set_ylim(-0.3,0.02)
ax7.set_xlim(-70,-60)
ax7.annotate("B", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

### Panel D: fitting
ax8 = fig.add_subplot(324) 
ax8.plot(ii, intercept_mod + slope_mod * ii, '-',color='g')
ax8.plot(- (IV_i_below/IV_i_below[-1] -1.)**2, IV_vc_below -  IV_vc_below[-1],'.',color='g')#, label='$I_t = %0.03f$ nA' %i_thres)
ax8.fill_between(linspace(-1,0,50), -slope_mod*ones(50), zeros(50),\
                  color='g', alpha = 0.2)
ax8.set_xlabel('$-(I/I_t-1)^2$')
ax8.set_ylabel('$V-V_t$ (mV)')
ax8.set_ylim(-6, 0.5)
ax8.set_xlim(-1.5, 0.1)
ax8.legend(frameon=False, fontsize=8)
ax8.annotate("D", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')
    
fig.tight_layout()

show()

print ('N cells with It measurement:', len(slope_final))

### Saving the figure
# save_path = '/Users/sarah/Documents/repositories/AIS-geometry-and-axonal-current/Na currents in RGC/codes submission/data/'
save_path = '/Users/sarah/Dropbox/Spike initiation/PhD projects/Axonal current and AIS geometry/Paper/Figures/'
# fig.savefig(save_path + "fig6.pdf", bbox_inches='tight')













