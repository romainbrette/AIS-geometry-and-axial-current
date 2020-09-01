

"""

Figure: measurement of axonal current

"""
import glob2
import pandas as pd
import pyabf
import seaborn as sns
import matplotlib.patches as mpatches
from pandas import ExcelWriter
from pandas import ExcelFile
from brian2 import *
from scipy import stats
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from vc_test_pulse_analysis import *
from na_currents_analysis import *
from vc_test_pulse_analysis import *
import params_model_description

### Figure

rcParams['axes.spines.right'] = False
rcParams['axes.spines.top'] = False

fig1 = figure('Ip measurement', figsize=(6,10))

ax7 = fig1.add_subplot(421)
ax8 = fig1.add_subplot(422)
ax1 = fig1.add_subplot(423)
ax2 = fig1.add_subplot(424)
ax3 = fig1.add_subplot(425)
ax4 = fig1.add_subplot(426)
ax5 = fig1.add_subplot(427)
ax6 = fig1.add_subplot(428)

### Panel A and B: axonal currents recording in an exmaple cell
date = 20191115
retina = 'B'
cell = 2

### Loading the data
# path_to_cell = glob2.glob('data/RGC data/' + str(int(date)) + '*' + '/retina '+ str(retina) +'/cell ' + str(int(cell)))[0]
path_to_cell = glob2.glob('/Users/sarah/Documents/Data/Martijn Sierksma/' + str(int(date)) + '*' + '/retina '+ str(retina) +'/cell ' + str(int(cell)))[0]
abf = pyabf.ABF(path_to_cell + '/VC threshold adaptation/2019_11_15_0054.abf')
fs = abf.dataRate  * Hz # sampling rate
dt = 1./fs

cell_name = '%i %s %i' %(date, retina, cell)
vh = -60*mV # holding potential

t = dt*arange(len(abf.sweepY)) 
I = []
V = []

for sweepNumber in abf.sweepList:
    abf.setSweep(sweepNumber)
    I.append(abf.sweepY)
    V.append(abf.sweepC*mV)

### Removing passive component
I_corr_pass, I_cut, t_cut = p5_subtraction(date, retina, cell, dt, I, V, rec_name=str(int(33)).zfill(4))

### IV curves
I_peaks,  Vc_peaks, idx_peak_ax_current,  t_peaks = plot_IV(date, retina, cell, dt, I_corr_pass, V, 0, str(int(6)).zfill(4))
Vc_peaks = array(Vc_peaks/mV)
I_peaks = array(I_peaks)

### Plotting
cmap = plt.get_cmap('binary')
colors = [cmap(i) for i in np.linspace(0, 1, len(abf.sweepList))]

# Currents
ax7.plot(t_cut/ms, I_corr_pass[4] *1e-3, 'k', alpha=0.2) #color= colors[0]) 
ax7.text(0.5, 0.2,'%i'%(Vc_peaks[4] -70), fontsize=8)
ax7.plot(t_cut/ms, I_corr_pass[5] *1e-3,'k') #color= colors[5]) 
ax7.text(2.2, min(I_corr_pass[5])*1e-3,'%i'%(Vc_peaks[5] -70), fontsize=8 )
ax7.plot(t_cut/ms, I_corr_pass[7] *1e-3, 'k') #color= colors[7]) 
ax7.text(0.8, min(I_corr_pass[7])*1e-3,'%i'%(Vc_peaks[7] -70), fontsize=8)
ax7.plot(t_cut/ms, I_corr_pass[9] *1e-3, 'k') #color= colors[9]) 
ax7.text(0.4, min(I_corr_pass[9])*1e-3 -1.2,'%i'%(Vc_peaks[9] -70), fontsize=8)
ax7.plot(t_cut/ms, I_corr_pass[13] *1e-3, 'k', alpha=0.2) #color= colors[13]) 
ax7.text(0.4, (min(I_corr_pass[13])-700)*1e-3,'%i'%(Vc_peaks[13] -70), fontsize=8)
ax7.plot(t_cut/ms, I_corr_pass[17] *1e-3, 'k', alpha=0.2) #color= colors[17]) 
ax7.text(0.4, (min(I_corr_pass[17])-500)*1e-3,'%i'%(Vc_peaks[17] -70),  fontsize=8)
ax7.plot(t_cut/ms, I_corr_pass[-1] *1e-3, 'k', alpha=0.2) #color= colors[-1]) 
ax7.text(0.3, (min(I_corr_pass[-1])-500)*1e-3,'%i'%(Vc_peaks[-1] -70),  fontsize=8)
ax7.set_ylabel('I (nA)')
# ax7.set_xlabel('t (ms)')
ax7.set_xlim(0,3)
ax7.set_ylim(-20, 1)
ax7.set_xticks([])
sns.despine(bottom=True, ax=ax7)
ax7.plot(linspace(2, 2.5, 10), -19.5*ones(10), 'k-', linewidth=2)
ax7.text(1.98, -21.5, '0.5 ms', color='k', fontsize=8)
ax7.annotate("A", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

# IV curve
ax8.plot(Vc_peaks -70, I_peaks *1e-3, '-', color= 'crimson') #, linewidth=2) 
# ax8.annotate('threshold', xy=(-55.5, -0.25), xytext=(-50, -1), fontsize=8,
#             arrowprops=dict(arrowstyle='->', color='k', ls='-'))
# ax8.annotate('peak axonal current', xy=(-55, -6.5), xytext=(-60, -14), fontsize=8,
#             arrowprops=dict(arrowstyle='->', color='k', ls='-'))
# ax8.annotate('somatic currents', xy=(-44.5, -11), xytext=(-50, -19), fontsize=8,
#             arrowprops=dict(arrowstyle='->', color='k', ls='-'))
ax8.set_ylabel('Peak current (nA)')
ax8.set_xlabel('V (mV)')
ax8.set_ylim(-20, 1)
ax8.set_xlim(-61,-36)
ax8.annotate("B", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')


### Panel C and D: currents correction in an example cell

#### Load the list of the cells and the results from the currents analyses that will be used for panels E and F too
df_cells = pd.read_excel('RGC_electrical_properties.xlsx')

# the example cell
first_cell = -9
last_cell = -8 #len(df_cells['Date'])

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

for date, retina, cell, age, na_rec, tp_corr, rs_res, vh in zip(dates, retinas, cells, ages, na_recs, tp_corrs, rs_residuals, resting_mp): 
    print ('------------------------------')
    print (date, retina, cell)
    
    ### Load Na current recording
    path_to_cell = path_to_data + str(int(date)) + "*/" + '/retina '+ str(retina) +'/cell ' + str(int(cell))   

    if (date, retina, cell) == (20191115, 'B', 2):
        print('adaptation protocol')
        path_to_na_currents = glob2.glob(path_to_cell + '/VC threshold adaptation/' + '*' + str(int(na_rec)).zfill(4) + ".abf")
    else:
        path_to_na_currents = glob2.glob(path_to_cell + '/VC small steps/' + '*' + str(int(na_rec)).zfill(4) + ".abf")
             
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
    
    tau_m = []
            
    ### Measuring Tau by fitting an exponential to the response to a test voltage pulse
    for i in range(n_rec):
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

        ### exp fitting
        fit_end = idx_i_max + int(.5*ms/dt) 
        i_amp_end = i_cut[fit_end]
        
        max_amp = max(i_cut[idx_i_max:fit_end])
                
        def exp_current(t, tau):
            return  (i_amp_max_peak - max_amp) * exp(-t/tau) + max_amp
        
        idx_i_max = idx_i_max + 1
        peak_current = i_cut[idx_i_max: fit_end] # pA
        peak_time = (t_cut[idx_i_max:fit_end]- t_cut[idx_i_max])  # sec
        
        tau_opt = curve_fit(exp_current, peak_time, peak_current)
        
        print ('Tau m:', tau_opt[0])
        tau_m.append(tau_opt[0])
    
    Tau = median(tau_m)
    
    ### P5 subtraction to remove the passive component of the currents
    I_corr_pass, I_cut, t_corr = p5_subtraction(date, retina, cell, dt, I, V, rec_name=str(int(na_rec)).zfill(4))
    
    ### Panel C: exponential fit
    ax1.plot(t_cut - idx_start*dt/ms, i_cut*1e-3, 'k', label='data')
    ax1.plot(peak_time + t_cut[idx_i_max] - idx_start*dt/ms, exp_current(peak_time, tau_opt[0])*1e-3, '--',\
             color='crimson', label='exponential fit') #' fit: $\tau$ = %0.02f ms' %Taum)
    ax1.set_xlim(0, 2)
    ax1.set_ylim(-1.5,0.3)
    # ax1.set_xlabel('t (ms)')
    ax1.set_ylabel('I (nA)')
    ax1.legend(frameon=False, loc='center right', fontsize=8)
    ax1.set_xticks([])
    sns.despine(bottom=True, ax=ax1)
    ax1.plot(linspace(1.5,2,10), -1.48*ones(10), 'k-', linewidth=2)
    ax1.text(1.58, -1.65, '0.5 ms', color='k', fontsize=8)
    ax1.annotate("C", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')
         
    ### Series resistance correction
    print ('Off-line series resistance correction:')
    
    ### IV curve
    I_peaks,  Vc_peaks, idx_peak_ax_current, t_peaks = plot_IV(date, retina, cell, dt, I_corr_pass, V, 0, str(int(na_rec)).zfill(4))
    print ('Sweep Ip:', idx_peak_ax_current)
    
    # Traynelis algorithm
    v_rev = 70.
    si = dt/ms
    tauLag = si * 2.
    fc = 1./(2*pi*tauLag)
    filterfactor = 1
      
    idx_peak = argmin(I_corr_pass[idx_peak_ax_current][:int(9.*ms/dt)]) 
    i_data = I_corr_pass[idx_peak_ax_current][0: idx_peak+100] * 1e-3 #nA #I_high[idx_ax][int(59.*ms/dt):int(60.*ms/dt)]*1e-3 # nA
    t_data = t_corr[0: idx_peak+100]
              
    vc_data = Vc_peaks[idx_peak_ax_current]/mV + vh 
    print ('Vh = ', vc_data)
    
    I_cap = zeros(len(i_data))
    I_corr_cap = zeros(len(i_data))
    I_corr = zeros(len(i_data))
    
    n_pts = len(i_data)
    
    # first data point
    v_last = vc_data - i_data[0] * rs_res
    denom = v_last - v_rev
    
    if denom != 0:
        fracCorr = (1-(vc_data-v_rev)/denom)
    else:
        fracCorr = 0
    
    I_corr[0] = i_data[0] * (1-fracCorr)
    
    # next data points
    for j in range(n_pts):
        v_this = vc_data - i_data[j]*rs_res
        if v_this != v_rev:
            fracCorr = (1-(vc_data-v_rev)/(v_this-v_rev))
        else:
            fracCorr = 0
        Icap = Tau * ((i_data[j] - i_data[j-1])/si)  #(v_this-v_last)/si #* 1e-3
        I_cap[j] = Icap
        I_corr_cap[j-1] = i_data[j-1] + Icap
        I_corr[j-1] = I_corr_cap[j-1] * (1-fracCorr)
    
        v_last = v_this 
    
    ### Panel D: uncorrected and corrected axonal current
    ax2.plot(t_data[:-1]/ms, i_data[:-1], 'k-', label='raw')
    ax2.plot(t_data[:-1]/ms, I_corr[:-1], 'crimson', label='corrected')
    ax2.set_ylabel('I (nA)')
    # ax2.set_xlabel('t (ms)')
    ax2.set_ylim(-10,2)
    ax2.set_xlim(0., 2.)
    ax2.legend(loc='lower left', fontsize=8, frameon=False)
    ax2.set_xticks([])
    sns.despine(bottom=True, ax=ax2)
    ax2.plot(linspace(1.5,2,10), -9.87*ones(10), 'k-', linewidth=2)
    ax2.text(1.58, -11., '0.5 ms', color='k', fontsize=8)
    ax2.annotate("D", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

### Panel E and F: effect of series resistance in a model

### Model parameters
params = params_model_description

dt = 0.01*ms
v_h = params.EL # leak potential
start = 5.*um # AIS start position
length = 30.*um # AIS length

# axial resistance between the soma adge and the start/end of the AIS
Ra_start = (4*params.Ri/(pi*params.axon_diam**2)) * (start)
Ra_mid = (4*params.Ri/(pi*params.axon_diam**2)) * (start + 0.5*length)
print('Axial resistance start:', Ra_start/Mohm)
print('Axial resistance mid:', Ra_mid/Mohm)

# the range of series resistance
series_resistances = array([0,1,2,3,4,5])

### Path to the simulations data
path_to_data = 'simulations data/fig2/'

thresholds_command = zeros(len(series_resistances))
thresholds_true = zeros(len(series_resistances))
peak_currents = zeros(len(series_resistances))
threshold_currents = zeros(len(series_resistances))
peak_currents_latency = zeros(len(series_resistances))

peak_currents_corr = zeros(len(series_resistances))
threshold_currents_corr = zeros(len(series_resistances))
peak_currents_corr_max = zeros(len(series_resistances))
threshold_currents_corr_max = zeros(len(series_resistances))

### Voltage threshold and peak axonal current vs series resistance
for i in range(len(series_resistances)):
    rs = series_resistances[i]
    
    print ('Rs:', rs)

    dir_name = path_to_data + 'VC dicho APmodel ext AIS x%0.01f L%0.01f r%0.01f' %(start/um, length/um, rs)

    # Load and plot data
    Ie = loadtxt(dir_name +'/Steps/I.txt') # electrode current
    Vm = loadtxt(dir_name +'/Steps/V.txt') # membrane potential
    Vc = loadtxt(dir_name + '/Steps/Vc.txt') # command potential
    Im = loadtxt(dir_name+'/Steps/Im.txt') # membrane current
    
    n_rec = len(Ie)
    t = arange(len(Ie[0]))*dt
        
    # Load test pulse data
    dir_name_tp = path_to_data + 'Test pulse APmodel ext AIS x%0.01f L%0.01f r%0.01f' %(start/um,length/um,rs)

    Ie_tp = loadtxt(dir_name_tp +'/Steps/I.txt')
    Ve_tp = loadtxt(dir_name_tp +'/Steps/V.txt')
    Vc_tp = loadtxt(dir_name_tp + '/Steps/Vc.txt')
    
    t_tp = arange(len(Ve_tp))*dt
    vc_tp_amp = -5
    
    # estimate the RC time constant from current decay
    idx_start = int(20*ms/dt)
    idx_i_min = argmin(Ie_tp[idx_start:idx_start+100]) + idx_start  #transient negative peak
    idx_end = idx_i_min + int(1.*ms/dt)
    
    i_baseline = mean(Ie_tp[int(18.*ms/dt):int(19.5*ms/dt)])
    i_plateau = mean(Ie_tp[int(30.*ms/dt):int(40.*ms/dt)]) - i_baseline
    
    i_amp_peak = Ie_tp[idx_i_min] #- i_baseline
    i_amp_end = Ie_tp[idx_end] #- i_baseline
    
    def exp_current(t, tau):
        return abs((i_amp_peak - i_amp_end) * exp(-t/tau) + i_amp_end)
 
    peak_current = abs(Ie_tp[idx_i_min:idx_end]) # pA
    peak_time = (t_tp[idx_i_min:idx_end]- t[idx_i_min])/ms  # ms

    # Tau
    tau_opt = curve_fit(exp_current, peak_time, peak_current)   
    tau_RC = tau_opt[0]
    print ('Tau:', tau_RC)
    
    # Remove passive response from Na currents recordings
    vc = Vc 
    
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
        
    # Measure peak axonal current and threshold
    spikes = zeros(n_rec)
    i_peaks = []
    vm_peaks = [] # true membrane potentials
    
    for j in range(n_rec):
        idx_peak = argmin(Ie_corr[j])
        i_peak = Ie_corr[j][idx_peak]
        i_peaks.append(i_peak)
        vm_peak = max(Vm[j][steps_start: steps_end]) #Vm[j][idx_peak + steps_start]
        vm_peaks.append(vm_peak)
        if i_peak < -1.:
            spikes[j] = 1
            
    vc = array(vc)
    vm_peaks = array(vm_peaks)
    i_peaks = array(i_peaks)
    I_corr = array(Ie_corr)
    spikes = array(spikes)
        
    # Current threshold (highest non spiking situation)
    idx_no_spike = where(spikes == 0)[0]
    idx_th = argmin(i_peaks[idx_no_spike])
    Ie_thres = I_corr[idx_no_spike][idx_th]
    
    Vc_thres = vc[idx_no_spike][idx_th]
    Vm_thres = max(Vm[idx_no_spike][idx_th])     #     vm_peaks[idx_no_spike][idx_th]
    
    Ie_thres_cut = Ie[idx_no_spike][idx_th][steps_start+1:steps_end]
    Vm_true_thres = max(Vc_thres - rs * Ie_thres_cut)
    
    # Peak current (first spiking situation)
    idx_spike = where(spikes == 1)[0]
    idx_peak = argmax(i_peaks[idx_spike])
    Ie_peak = I_corr[idx_spike][idx_peak]
    
    thresholds_command[i] = Vc_thres
    thresholds_true[i] = Vm_true_thres #Vm_thres
    threshold_currents[i] = i_peaks[idx_no_spike][idx_th]  
    peak_currents[i] = i_peaks[idx_spike][idx_peak] 
    peak_currents_latency[i] = argmin(Ie_peak)*dt/ms
    
    # Correction for filtering
    v_rev = 70.
    si = dt/ms
    tauLag = si #* 2.
    fc = 1./(2*pi*tauLag)
    filterfactor = 1
        
    # Peak current 
    idx_peak_peak = argmin(Ie_peak)
    print (idx_peak_peak)
    i_data = Ie_peak[:idx_peak_peak+10]
    vc_data = vc[idx_spike][0] 
    print ('Vh = ', vc_data)
    
    I_cap = zeros(len(i_data))
    I_corr_cap = zeros(len(i_data))
    I_corr_rs = zeros(len(i_data))
    
    n_pts = len(i_data)
    
    # first data point
    v_last = vc_data - i_data[0] * rs
    denom = v_last - v_rev
    
    if denom != 0:
        fracCorr = (1-(vc_data-v_rev)/denom)
    else:
        fracCorr = 0
    
    I_corr_rs[0] = i_data[0] * (1-fracCorr)
    
    # next data points
    for j in range(n_pts):
        v_this = vc_data - i_data[j]*rs
    
        if v_this != v_rev:
            fracCorr = (1-(vc_data-v_rev)/(v_this-v_rev))
        else:
            fracCorr = 0
        Icap = tau_RC * ((i_data[j] - i_data[j-1])/si) 
        I_cap[j] = Icap
        I_corr_cap[j-1] = i_data[j-1] + Icap
        I_corr_rs[j-1] = I_corr_cap[j-1] * (1-fracCorr)
                
        v_last = v_this 
    
    peak_currents_corr[i] = min(I_corr_rs) 

### Panel E: uncorrected and corrected peak axonal current vs series resistance
ax3.plot(series_resistances, peak_currents, 'g-', label='raw')
ax3.plot(series_resistances, peak_currents[0] * ones(len(peak_currents)), 'g--')
ax3.plot(series_resistances, peak_currents_corr, 'crimson', label='corrected')
ax3.set_ylabel('$I_p$ (nA)')
ax3.set_xlabel('$R_s$ (M$\Omega$)')
ax3.set_ylim(-15,0)
ax3.set_xlim(0,5)
ax3.legend(frameon=False, fontsize=8)
ax3.annotate("E", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

### Panel F: voltage threshold vs series resistance
ax4.plot(series_resistances, thresholds_command, 'g-', label='command')
ax4.plot(series_resistances, thresholds_true, 'g', alpha=0.5, label='true')
ax4.set_ylabel('$V_t$ (mV)')
ax4.set_xlabel('$R_s$ (M$\Omega$)')
ax4.set_ylim(-75, -45)
ax4.set_xlim(0,5)
ax4.legend(frameon=False, fontsize=8)
ax4.annotate("F", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

### Panel G and H: threshold and peak axonal current vs series reistance in RGC

### Results from axonal currents analyses in all cells
dates = array(df_cells['Date'])[:-3]
retinas = array(df_cells['Retina'])[:-3]
cells = array(df_cells['Cell'])[:-3]
ages = array(df_cells['Age'])[:-3]

voltage_thresholds = array(df_cells['Vth'])[:-3]
voltage_thresholds_true = array(df_cells['Vth true'])[:-3]
v_end = array(df_cells['V end (mV)'])[:-3]

axonal_currents = array(df_cells['Peak axonal current corrected'])[:-3]
series_resistance_residual = array(df_cells['Residual Rs'])[:-3]

n_cells = len(dates)

### Panel G: peak axonal current vs residual series resistance
print ('N panel G:', len(series_resistance_residual))

# Regression
slope_l, intercept_l, r_l, p_l, _ = stats.linregress(series_resistance_residual,axonal_currents) 
ax5.plot(linspace(0.2,4.9,50), intercept_l + slope_l * linspace(0.2,4.9,50), 'k-')

sns.scatterplot(x=series_resistance_residual, y=axonal_currents, color='k', ax=ax5)
ax5.set_ylim(-15,0)
ax5.set_xlim(0,5)
ax5.set_ylabel('$I_p$ (pA)')
ax5.set_xlabel('$R_s$ (M$\Omega$)')
ax5.legend(frameon=False)
ax5.annotate("G", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

### Panel H: voltage threshold vs residual series resistance
# for the threshold analysis, only cells with Vend close to LJP are included
ljp = -11.
idx_discard = where([ljp - 3 <= v_end[i] <= ljp + 3 for i in range(n_cells)])[0]

print ('N panel H:', len(series_resistance_residual[idx_discard]))

# Regression
slope_v, intercept_v, r_v, p_v, _ = stats.linregress(series_resistance_residual[idx_discard], \
                                voltage_thresholds_true[idx_discard]) 

ax6.plot(linspace(0.2,4.8,50), intercept_v + slope_v * linspace(0.2,4.8,50), 'k-')
sns.scatterplot(x=series_resistance_residual[idx_discard], \
         y=voltage_thresholds_true[idx_discard], color='k', ax=ax6)
ax6.set_ylim(-75, -45)
ax6.set_xlim(0,5)
ax6.set_ylabel('$V_t$ (mV)')
ax6.set_xlabel('$R_s$ (M$\Omega$)')
ax6.legend(frameon=False)
ax6.annotate("H", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

fig1.tight_layout()

show()

print ('STATS')
print ('Linregress Rs-Vt:', 'r=', r_v, 'p=', p_v)
print ('Linregress Rs-Ip:', 'r=', r_l, 'p=', p_l)
print ('Pearson Rs-Vt:', stats.pearsonr(series_resistance_residual[idx_discard], voltage_thresholds_true[idx_discard]))
print ('Pearson Rs-Ip:', stats.pearsonr(series_resistance_residual,axonal_currents))

save_path = '/Users/sarah/Dropbox/Spike initiation/PhD projects/Axonal current and AIS geometry/Paper/Figures/'


# save_path = '/Users/sarah/Dropbox/Spike initiation/PhD projects/Axonal current and AIS geometry/Paper/Figures/'

# fig1.savefig(save_path + "fig2.pdf", bbox_inches='tight')














        
