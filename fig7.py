


"""

Figure 7: current and voltage at threshold and AIS geometry.


"""

from brian2 import *
from pandas import ExcelWriter
from pandas import ExcelFile
from scipy import stats
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import pandas as pd
import seaborn as sns
from vc_test_pulse_analysis import *
from na_currents_analysis import *
import params_model_description

rcParams['axes.spines.right'] = False
rcParams['axes.spines.top'] = False

### Loading the results of analyses
df_cells = pd.read_excel('RGC_electrical_properties.xlsx')

dates = array(df_cells['Date'])[:-3]
retinas = array(df_cells['Retina'])[:-3]
cells = array(df_cells['Cell'])[:-3]
ages = array(df_cells['Age'])[:-3]

voltage_thresholds = array(df_cells['Vth'])[:-3]
voltage_thresholds_true = array(df_cells['Vth true'])[:-3]
v_end = array(df_cells['V end (mV)'])[:-3]

series_resistance_residual = array(df_cells['Residual Rs'])[:-3]
ais_starts = array(df_cells['AIS start (um)'])[:-3]
ais_lengths = array(df_cells['AIS length (um)'])[:-3]
diam_starts = array(df_cells['Axon diam start (um)'])[:-3]
ais_mids = ais_starts + 0.5 * ais_lengths

### For the voltage threshold analysis, only cells with Vend close to LJP are included
n_cells = len(dates)
ljp = -11.
idx_discard = where([ljp - 3 <= v_end[i] <= ljp + 3 for i in range(n_cells)])[0]

dates_ljp3 = dates[idx_discard]
retinas_ljp3 = retinas[idx_discard]
cells_ljp3 = cells[idx_discard]

ais_starts_ljp3 = ais_starts[idx_discard]
ais_length_ljp3 = ais_lengths[idx_discard]
ais_mids_ljp3 = ais_mids[idx_discard]
voltage_thresholds_ljp3 = voltage_thresholds[idx_discard]
voltage_thresholds_true_ljp3 = voltage_thresholds_true[idx_discard]

### Removing cells without AIS geometry measurement

dates_plot = dates_ljp3[~isnan(ais_length_ljp3)]
retinas_plot = retinas_ljp3[~isnan(ais_length_ljp3)]
cells_plot = cells_ljp3[~isnan(ais_length_ljp3)]

ais_lengths_plot = ais_length_ljp3[~isnan(ais_length_ljp3)]
ais_mids_plot = ais_mids_ljp3[~isnan(ais_length_ljp3)]
ais_starts_plot = ais_starts_ljp3[~isnan(ais_length_ljp3)]
voltage_thresholds_plot = voltage_thresholds_ljp3[~isnan(ais_length_ljp3)]
voltage_thresholds_true_plot = voltage_thresholds_true_ljp3[~isnan(ais_length_ljp3)]

### Figure

### Voltage threshold

fig = figure('Threshold and current threshold', figsize=(9,6))

print ('Panel A-B-C')
print ('N cells:', len(voltage_thresholds_true_plot))

### Panel A: voltage threshold vs AIS start position
ax1 = fig.add_subplot(231)
# Regression
slope_vd, inter_vd, r_vd, p_vd, _ = stats.linregress(ais_starts_plot,voltage_thresholds_true_plot)

sns.scatterplot(x=ais_starts_plot, y=voltage_thresholds_true_plot, color='k', ax=ax1)
ax1.set_ylim(-60,-45)
ax1.set_xlim(0,20)
ax1.set_xlabel('$\Delta$ ($\mu$m)')
ax1.set_ylabel('$V_t$ (mV)')
ax1.legend(frameon=False, fontsize=8)
ax1.annotate("A", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

### Panel B: voltage threshold vs AIS length
ax2 =  fig.add_subplot(232)

# Regression
slope_vl, inter_vl, r_vl, p_vl, _ = stats.linregress(ais_lengths_plot,voltage_thresholds_true_plot)
sns.scatterplot(x=ais_lengths_plot, y=voltage_thresholds_true_plot, color='k', ax=ax2)
ax2.set_ylim(-60,-45)
ax2.set_xlim(20,50)
ax2.set_xlabel('$L$ ($\mu$m)')
ax2.set_ylabel('$V_t$ (mV)')
ax2.legend(frameon=False, fontsize=8)
ax2.annotate("B", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

# Panel C: voltage threshold vs middle position x length
ax3 =  fig.add_subplot(233)
# regression
slope_vx, inter_vx, r_vx, p_vx, _ = stats.linregress(-log(ais_lengths_plot * ais_mids_plot), \
                                                voltage_thresholds_true_plot)

ax3.set(xscale="log")
ax3.semilogx(linspace(400, 1400, 50), inter_vx - slope_vx * log(linspace(400, 1400, 50)), 'k-')    
sns.scatterplot(x=ais_lengths_plot * ais_mids_plot, y=voltage_thresholds_true_plot, color='k', ax=ax3)
ax3.set_ylim(-60,-45)
ax3.xaxis.set_major_locator(plt.MultipleLocator(1000))
ax3.xaxis.set_minor_locator(plt.MultipleLocator(500))
ax3.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax3.xaxis.set_minor_formatter(FormatStrFormatter('%.0f'))
ax3.set_xlabel('$x_{1/2}L$ ($\mu$m$^2$)')
ax3.set_ylabel('$V_t$ (mV)')
ax3.set_xlim(0, 1500)
ax3.legend(frameon=False, fontsize=8)
ax3.annotate("C", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

### Current at threshold

threshold_currents = array(df_cells['Threshold current (nA)'])[:-3]
thresholds = array(df_cells['Vth'])[:-3]

ais_starts = array(df_cells['AIS start (um)'])[:-3]
ais_lengths = array(df_cells['AIS length (um)'])[:-3]
diam_starts = array(df_cells['Axon diam start (um)'])[:-3]
ais_mids = ais_starts + 0.5 * ais_lengths

residual_rs = array(df_cells['Residual Rs'])[:-3]

# remove cells without threshold current
dates_it = dates[~isnan(threshold_currents)]
retinas_it = retinas[~isnan(threshold_currents)]
cells_it = cells[~isnan(threshold_currents)]

threshold_currents_it = threshold_currents[~isnan(threshold_currents)]
ais_starts_it = ais_starts[~isnan(threshold_currents)]
ais_mids_it = ais_mids[~isnan(threshold_currents)]

# remove cells without AIS geometry measurement
dates_it_plot = dates_it[~isnan(ais_starts_it)]
retinas_it_plot = retinas_it[~isnan(ais_starts_it)]
cells_it_plot = cells_it[~isnan(ais_starts_it)]

threshold_currents_it_plot = threshold_currents_it[~isnan(ais_starts_it)]
ais_starts_it_plot = ais_starts_it[~isnan(ais_starts_it)]
ais_mid_it_plot = ais_mids_it[~isnan(ais_starts_it)]

### Panel D: threshold current vs AIS start position in RGC
print ('Panel D-E')
print ('N cells:', len(threshold_currents_it_plot))

ax4 = fig.add_subplot(234)

# Regression
slope_start, intercept_start, r_start, p_start, _ = stats.linregress(ais_starts_it_plot,\
                                -threshold_currents_it_plot) 

sns.scatterplot(x=ais_starts_it_plot, y=-threshold_currents_it_plot, color='k', ax=ax4)
ax4.plot(linspace(1, 14, 50), intercept_start + slope_start * linspace(1, 14, 50), 'k-') 
ax4.set_xlabel('$\Delta$ ($\mu$m)')
ax4.set_ylabel('$-I_t$ (nA)')
ax4.set_ylim(0,0.3)
ax4.set_xlim(0,20)
ax4.legend(frameon=False, fontsize=8)
ax4.annotate("D", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

### Panel E: threshold current vs AIS middle position in RGC and in a model

### Model parameters
params = params_model_description
ENa = params.ENa
Vs = -55*mV
DeltaV = ENa - Vs

d = 1.*um
dmin = .8*um
dmax = 1.3*um
k = params.Ka
Ri = params.Ri

ais_mid_range = linspace(16,34,100)

ra = Ri * 4 / (pi * d ** 2)
It = -k/(ra*ais_mid_it_plot*um)
It_pred = -k/(ra*ais_mid_range*um)

dt = 0.01*ms
v_h = params.EL

starts = linspace(0, 20, 5)*um 
lengths = [30]*um 

### Loading data from simulations
path_to_data = 'simulations data/fig7/'

peak_currents = zeros((len(lengths), len(starts)))
threshold_currents_sim = zeros((len(lengths), len(starts)))

for l in range(len(lengths)):
    length = lengths[l]

    for i in range(len(starts)):
    
        start = starts[i]
        
        print (start)
        
        dir_name = path_to_data + 'VC dicho APmodel ext AIS x%0.1f L%i' %(start/um, length/um)
    
        # Load and plot data
        Ie = loadtxt(dir_name +'/Steps/I.txt')
        Vm = loadtxt(dir_name +'/Steps/V.txt')
        Vc = loadtxt(dir_name + '/Steps/Vc.txt')
        Im = loadtxt(dir_name+'/Steps/Im.txt')
        
        n_rec = len(Ie)
        t = arange(len(Ie[0]))*dt   
        
        # Load test pulse data
        dir_name_tp = path_to_data + 'TP APmodel ext AIS x%0.1f L%i' %(start/um, length/um)
    
        Ie_tp = loadtxt(dir_name_tp +'/Steps/I.txt')
        Ve_tp = loadtxt(dir_name_tp +'/Steps/V.txt')
        Vc_tp = loadtxt(dir_name_tp + '/Steps/Vc.txt')
        
        t_tp = arange(len(Ve_tp))*dt
        vc_tp_amp = -5
        
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
        
        # Current threshold
        idx_no_spike = where(spikes == 0)[0]
        idx_th = argmin(i_peaks[idx_no_spike])
        Ie_thres = I_corr[idx_no_spike][idx_th]
        
        # Peak current
        idx_spike = where(spikes == 1)[0]
        idx_peak = argmax(i_peaks[idx_spike])
        
        threshold_currents_sim[l, i] = i_peaks[idx_no_spike][idx_th] #min(Ie_thres)
        peak_currents[l, i] = i_peaks[idx_spike][idx_peak]
    
show()

ax7 = fig.add_subplot(235)
idx = 0
ra = (4*params.Ri/(pi*params.axon_diam**2))
### Simulation
ax7.loglog((starts+lengths[idx]/2)/um, -threshold_currents_sim[idx,:], '-', color='g', label='simulation')
### Theory
ax7.loglog((starts+lengths[idx]/2)/um, params.Ka/(ra * (starts+0.5*lengths[idx]))/nA, '--', color='g', label='theory')
### Data
ax7.set(xscale='log', yscale='log')
sns.scatterplot(x=ais_mid_it_plot, y=-threshold_currents_it_plot, color='k', label='data', ax=ax7)
# Regression
slope_ix, intercept_ix, r_ix, p_ix,_ = stats.linregress(ais_mid_it_plot,-threshold_currents_it_plot)
    
ax7.set_ylabel('$-I_t$ (nA)')
ax7.set_xlabel('$x_{1/2}$ ($\mu$m)')
ax7.yaxis.set_major_formatter(FormatStrFormatter('%.01f'))
ax7.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax7.xaxis.set_minor_formatter(FormatStrFormatter('%.0f'))
ax7.set_ylim(0,1)
ax7.set_xlim(15, 35)
ax7.legend(frameon=False, fontsize=8)
ax7.annotate("E", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

fig.tight_layout()

print ('STATS')
print('Panel A')
print ('Linregress Delta-V', 'r=', r_vd, 'p=', p_vd)
print ('Pearson Delta-V', stats.pearsonr(ais_starts_plot, voltage_thresholds_true_plot))

print('Panel B')
print ('Linregress L-V', 'r=', r_vl, 'p=', p_vl)
print ('Pearson L-V', stats.pearsonr(ais_lengths_plot, voltage_thresholds_true_plot))

print('Panel C')
print ('Linregress x-V', 'r=', r_vx, 'p=', p_vx, 'slope=', slope_vx)
print ('Pearson x-V', stats.pearsonr(-log(ais_lengths_plot * ais_mids_plot), \
                                                voltage_thresholds_true_plot))
print('Panel D')
print ('Linregress Delta-It data', 'r=', r_start, 'p=', p_start)
print ('Pearson Delta-It data:', stats.pearsonr(ais_starts_it_plot, -threshold_currents_it_plot))

print('Panel E')
print ('Linrergess x-It data', 'r=', r_ix, 'p=', p_ix)
print ('Pearson x-It data:', stats.pearsonr(ais_mid_it_plot,-threshold_currents_it_plot))

### Saving the figure

save_path = '/Users/sarah/Dropbox/Spike initiation/PhD projects/Axonal current and AIS geometry/Paper/Figures/'
# fig.savefig(save_path + "fig7.pdf", bbox_inches='tight')

# fig.savefig(save_path + "fig7.pdf", bbox_inches='tight')

# fig.savefig(save_path + "fig7.png", dpi=300)


    
  

