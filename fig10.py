


"""

Figure 10: compensation of axial current attenuation.

"""
from brian2 import *
import pandas as pd
import glob2
import pyabf
import statsmodels.api as sm
from pandas import ExcelWriter
from pandas import ExcelFile
from scipy import stats
from scipy import linalg
import seaborn as sns
from matplotlib import gridspec
from trace_analysis import *
from scipy.interpolate import CubicSpline

rcParams['axes.spines.right'] = False
rcParams['axes.spines.top'] = False

### Loading the results of analyses
# Loading the recording database 
df_rec_database = pd.read_excel('RGC_recording_database.xlsx')

# Load the adaptation data
df_cells = pd.read_excel('RGC_adaptation.xlsx')

dates = array(df_cells['Date'])
retinas = array(df_cells['Retina'])
cells = array(df_cells['Cell'])
ages = array(df_cells['Age'])

prepulse_potentials = array(df_cells['V prepulse']) 
holding_potentials = array(df_cells['Vh']) 

axonal_currents = array(df_cells['Peak axonal current corrected'])
threshold_potentials = array(df_cells['Vth']) 
charges = array(df_cells['Charge1 10']) # first peak above threshold because for some recordings, the seocnd is somatic current 
duration10= array(df_cells['Duration1 10']) 
duration50= array(df_cells['Duration1 50']) 

Rs_before = array(df_cells['Rs before']) 
Rs_after = array(df_cells['Rs after']) 
Rs_rec = array(df_cells['Rs Na rec'])  

selected_dates = []
selected_retinas = []
selected_cells = []
selected_ages = []
selected_capa = []
selected_prepulse_pot = []
selected_axonal_currents = []
selected_threshold_potentials = []
selected_charge = []
selected_dur10 = []
selected_dur50 = []

ljp = -11. 

### Removing bad quality recordings
for date, retina, cell, age, v0, vh, rs_before, rs_after, rs_rec, ia, vth, charge, dur10, dur50 in \
                            zip(dates, retinas, cells, ages, \
                                prepulse_potentials, holding_potentials, Rs_before, Rs_after, Rs_rec,\
                                axonal_currents, threshold_potentials,
                                charges, duration10, duration50): 
                                
    if rs_rec < 25 and  rs_after < 1.3 * rs_before:
        selected_dates.append(date)
        selected_retinas.append(retina)
        selected_cells.append(cell)
        selected_ages.append(age)
        selected_prepulse_pot.append(v0)
        selected_axonal_currents.append(ia)
        selected_charge.append(charge)
        selected_dur10.append(dur10)
        selected_dur50.append(dur50)
        
        row = df_rec_database[(df_rec_database['Date'] == date) & (df_rec_database['retina'] == retina) & (df_rec_database['cell'] == cell)]
        v_end = row['Vend'].values[0]
        
        if ljp-3 <= v_end <= ljp+3 or v_end != v_end:
            selected_threshold_potentials.append(vth)
        else:
            selected_threshold_potentials.append(nan)

### Counting cells and sorting measures per cell
N = 0
date_prev = selected_dates[0]
retina_prev = selected_retinas[0]
cell_prev = selected_cells[0]

dates_per_cell = [selected_dates[0]]
retina_per_cell = [selected_retinas[0]]
cell_per_cell = [selected_cells[0]]

v0_per_cell = []
vth_per_cell = []
ia_per_cell = []
charge_per_cell = []
dur10_per_cell = []
dur50_per_cell = []

v0_cell = []
vth_cell = []
ia_cell = []
ch_cell = []
dur10_cell = []
dur50_cell = []

for date, retina, cell, age, v0, ia, vth, charge, dur10, dur50 in zip(selected_dates, selected_retinas, selected_cells, \
                                                selected_ages, selected_prepulse_pot, \
                                                selected_axonal_currents, selected_threshold_potentials,\
                                                selected_charge, selected_dur10, selected_dur50): 

    if date_prev == date and retina_prev == retina and cell_prev == cell:
        v0_cell.append(v0)
        vth_cell.append(vth)
        ia_cell.append(ia)
        ch_cell.append(charge)
        dur10_cell.append(dur10)
        dur50_cell.append(dur50)
    else: 
        N += 1
        dates_per_cell.append(date)
        retina_per_cell.append(retina)
        cell_per_cell.append(cell)
        v0_per_cell.append(v0_cell)
        vth_per_cell.append(vth_cell)
        ia_per_cell.append(ia_cell)
        charge_per_cell.append(ch_cell)
        dur10_per_cell.append(dur10_cell)
        dur50_per_cell.append(dur50_cell)
        
        v0_cell = [v0]
        vth_cell = [vth]
        ia_cell = [ia]
        ch_cell = [charge]
        dur10_cell = [dur10]
        dur50_cell = [dur50]
        
    date_prev = date
    retina_prev = retina
    cell_prev = cell

### Adding the last cell
N += 1
v0_per_cell.append(v0_cell)
vth_per_cell.append(vth_cell)
ia_per_cell.append(ia_cell)
charge_per_cell.append(ch_cell)
dur10_per_cell.append(dur10_cell)
dur50_per_cell.append(dur50_cell)

### Removing doubles
v0_range = linspace(-75, -30, 10)
for i in range(N):
    # Removing the recordings with same v0
    for v0 in v0_range:
        idx_v0 = where(v0_per_cell[i] == v0)[0]
        if len(idx_v0) > 1:
            idx_v0_max = argmin(ia_per_cell[i][idx_v0])
            idx_v0_delete = delete(idx_v0, idx_v0_max)
            v0_per_cell[i] = delete(v0_per_cell[i], idx_v0_delete)
            ia_per_cell[i] = delete(ia_per_cell[i], idx_v0_delete)
            vth_per_cell[i] = delete(vth_per_cell[i], idx_v0_delete)
            charge_per_cell[i] = delete(charge_per_cell[i], idx_v0_delete)
            dur50_per_cell[i] = delete(dur50_per_cell[i], idx_v0_delete)
            dur10_per_cell[i] = delete(dur10_per_cell[i], idx_v0_delete)
        else:
            v0_per_cell[i] = array(v0_per_cell[i])
            ia_per_cell[i] = array(ia_per_cell[i])
            vth_per_cell[i] = array(vth_per_cell[i])
            charge_per_cell[i] = array(charge_per_cell[i])
            dur50_per_cell[i] = array(dur50_per_cell[i])
            dur10_per_cell[i] = array(dur10_per_cell[i])
        
### Current and charge attenuation
current_attenuation = []
charge_attenuation = []

for i in range(N):
    print (dates_per_cell[i], retina_per_cell[i], cell_per_cell[i])
    # attenuation
    idx_60 = where(v0_per_cell[i] == -60.)[0]
    idx_40 = where(v0_per_cell[i] == -40.)[0]
    if len(idx_60) > 0 and len(idx_40) > 0:
        current_attenuation.append(ia_per_cell[i][idx_60[0]]/ia_per_cell[i][idx_40[0]])
        charge_attenuation.append(charge_per_cell[i][idx_60[0]]/charge_per_cell[i][idx_40[0]])
    else:
        current_attenuation.append(nan)
        charge_attenuation.append(nan)
    
### Figure

name1 = "tab20b"
name2 = "tab20c"
name3 = "tab20"
cmap1 = get_cmap(name1)  
cmap2 = get_cmap(name2)  
cmap3 = get_cmap(name3) 
cols = cmap1.colors + cmap2.colors + cmap2.colors 

fig = figure('Compensation', figsize=(9,9))
gs = gridspec.GridSpec(3, 3, width_ratios=[1, 1, 1]) 

ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])
ax3 = fig.add_subplot(gs[2])
ax4 = fig.add_subplot(gs[3])
ax5 = fig.add_subplot(gs[4:6])
ax6 = fig.add_subplot(gs[6])
ax8 = fig.add_subplot(gs[7])
ax9 = fig.add_subplot(gs[8])
# ax9 = fig.add_subplot(gs[11])

### Panel A: peak current vs V0 in an example cell
i = -6
idx_sort = argsort(v0_per_cell[i])

ax1.set_ylabel('$-I_p$')
ax1.set_xlabel('$V_0$ (mV)')
ax1.plot(v0_per_cell[i][idx_sort], -ia_per_cell[i][idx_sort], '-o', color = cols[14])
ax1.set_ylim(0,4)
ax1.set_xlim(-80, -39)
ax1.annotate("A", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

### Panel B: charge vs V0 in example cell
ax2.set_ylabel('-Q (pC)')
ax2.set_xlabel('$V_0$ (mV)')
ax2.plot(v0_per_cell[i][idx_sort], -charge_per_cell[i][idx_sort], '-o', color = cols[14])
ax2.set_ylim(0, 1.5)
ax2.set_xlim(-80, -39)
ax2.annotate("B", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

### Panel C: current curation vs V0 in example cell
ax3.set_ylabel('$t_{50}$ (ms)')
ax3.set_xlabel('$V_0$ (mV)')
ax3.plot(v0_per_cell[i][idx_sort], dur50_per_cell[i][idx_sort], '-o', color = cols[14])
ax3.set_ylim(0.,2)
ax3.set_xlim(-80, -39)
ax3.annotate("C", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

### Panel D: current vs charge attenuation
# ax4.plot(charge_attenuation, current_attenuation, 'k.')
sns.scatterplot(x=charge_attenuation, y=current_attenuation, color = 'k', ax=ax4)
ax4.plot(linspace(0,20,10), linspace(0,20,10), 'k--')
ax4.set_xlim(0, 25)
ax4.set_ylim(0, 25)
ax4.set_ylabel('$I_{60}/I_{40}$')
ax4.set_xlabel('$Q_{60}/Q_{40}$')
ax4.annotate("D", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

### Panel E: current duration vs V0

v0_range = linspace(-75, -35, 9)
mean_duration = []
t50_per_v0 = []

for v0 in v0_range:
    t50_at_v0 = []
    for i in range(len(v0_per_cell)):
        idx_v0 = where(v0_per_cell[i] == v0)[0]
        if len(idx_v0) > 0:
            t50 = dur50_per_cell[i][idx_v0[0]]
            t50_at_v0.append(t50)
        else:
            t50_at_v0.append(nan)
    t50_per_v0.append(t50_at_v0)
    mean_duration.append(nanmean(t50_at_v0))

df = pd.DataFrame({'-75': t50_per_v0[0], 
                   '-70': t50_per_v0[1],
                   '-65': t50_per_v0[2],
                   '-60': t50_per_v0[3],
                   '-55': t50_per_v0[4],
                   '-50': t50_per_v0[5],
                   '-45': t50_per_v0[6],
                   '-40': t50_per_v0[7]})

sns.boxplot(data = df, color='gray', ax=ax5)

sns.set_palette(sns.color_palette(cols))
sns.swarmplot(data=df, ax=ax5)
    
ax5.set_ylabel('$t_{50}$ (ms)')
ax5.set_xlabel('$V_0$ (mV)')
ax5.set_ylim(0,2)
ax5.legend(frameon=False, fontsize=8)
ax5.annotate("E", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

### Panel F: spontaneous activity

day = '20200213'
retina = 'B'
cell = '1'

### Path to recordings
path = '/Users/sarah/Documents/Data/Martijn Sierksma/'
cc_cont_path = glob2.glob(path + '{0}'.format(day) + '*' + '/retina {0}/cell {1}/CC cont/'.format(retina, cell))[0]

abf = pyabf.ABF(cc_cont_path + "2020_02_13_0102.abf".format(day))
fs = abf.dataRate  * Hz # sampling rate
dt = 1./fs

### Spike times 
data =  abf.sweepY
spike_times = find_spikes_at(data[int(110*second/dt):int(112.500*second/dt)], dt, thres=-30) + 110*second
idx_spikes = spike_times/dt

cmap = plt.get_cmap('binary_r')
cols = [cmap(i) for i in np.linspace(0, 1, int(len(spike_times)/1.5))]

### Find the smallest AP
v_peaks = []
for i in range(1, len(spike_times)):
    peak = max(data[int(idx_spikes[i-1]): int(idx_spikes[i])])
    v_peaks.append(peak)
min_peak = argmin(v_peaks) + 1

### Plot spontaneous activity
abf.setSweep(0)
t = dt*np.arange(len(abf.sweepY))
ax6.plot(t/second, abf.sweepY, color='k', linewidth=0.5)
ax6.set_xlim(110.1, 112.4)
ax6.set_ylim(-80, 40)
ax6.set_ylabel('V (mV)')
ax6.plot(linspace(112,112.1,10), -75.*ones(10), 'k-', linewidth=2)
ax6.text(111.85, -85.,'100 ms',color='k', fontsize=8)
ax6.set_xticks([])
sns.despine(bottom=True, ax=ax6)
ax6.annotate("F", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

# Panel G: phase plots
for i in range(int(min_peak/3)+1):
    
    ax6.plot(spike_times[3*i]/second, 35, '|', color=cols[3*i])
    
    idx_spike = int(spike_times[3*i]/dt)
    
    # Measures
    f = data[idx_spike-200:idx_spike+100]
    t_spike = t[idx_spike-200:idx_spike+100]/ms - t[idx_spike-200]/ms
    t_new = (t_spike[:-1] + t_spike[1:])/2
    v = (f[:-1] + f[1:])/2
    dv = (f[1:] - f[:-1])/(dt/ms)
    ddv = (dv[1:] - dv[:-1])/(dt/ms) # shift of dt: add 1 !!! (f[2:] - 2*f[1:-1] + f[:-2])/dt**2 #
    
    ax8.plot(v, dv, color=cols[3*i])
    
ax8.set_ylabel('dV/dt (mV/ms)')
ax8.set_xlabel('V (mV)')
ax8.annotate("G", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

v_onsets = []
v_regeneration = []
for i in range(min_peak):
    idx_spike = int(spike_times[i]/dt)
    
    # Measures
    f = data[idx_spike-200:idx_spike+100]
    t_spike = t[idx_spike-200:idx_spike+100]/ms - t[idx_spike-200]/ms
    t_new = (t_spike[:-1] + t_spike[1:])/2
    v = (f[:-1] + f[1:])/2
    dv = (f[1:] - f[:-1])/(dt/ms)
    ddv = (dv[1:] - dv[:-1])/(dt/ms) # shift of dt: add 1 !!! (f[2:] - 2*f[1:-1] + f[:-2])/dt**2 #
    
    # AP peak
    idx_peak = argmax(v)
    
    # spike onset 
    idx_spike_onset = spike_onsets(v*mV, criterion = 20*volt/second * dt, v_peak = -30.*mV)
    spike_onset = v[idx_spike_onset[0]]
    if spike_onset > -30:
        idx_spike_onset = spike_onsets(v*mV, criterion = 10*volt/second * dt, v_peak = -30.*mV)
        spike_onset = v[idx_spike_onset[0]]
        if spike_onset > -30:
            idx_spike_onset = spike_onsets(v*mV, criterion = 5*volt/second * dt, v_peak = -30.*mV)
            spike_onset = v[idx_spike_onset[0]]
            if spike_onset > -30:
                idx_spike_onset = spike_onsets(v*mV, criterion = 0.5*volt/second * dt, v_peak = -30.*mV)
                spike_onset = v[idx_spike_onset[0]]
    v_onsets.append(spike_onset)
        
for i in range(min_peak-16):
    idx_spike = int(spike_times[i]/dt)
    
    # Measures
    f = data[idx_spike-200:idx_spike+100]
    t_spike = t[idx_spike-200:idx_spike+100]/ms - t[idx_spike-200]/ms
    t_new = (t_spike[:-1] + t_spike[1:])/2
    v = (f[:-1] + f[1:])/2
    dv = (f[1:] - f[:-1])/(dt/ms)
    ddv = (dv[1:] - dv[:-1])/(dt/ms) # shift of dt: add 1 !!! (f[2:] - 2*f[1:-1] + f[:-2])/dt**2 #
    
    # AP peak
    idx_peak = argmax(v)
    
    # spike onset 
    idx_spike_onset = spike_onsets(v*mV, criterion = 20*volt/second * dt, v_peak = -30.*mV)
    spike_onset = v[idx_spike_onset[0]]
    if spike_onset > -30:
        idx_spike_onset = spike_onsets(v*mV, criterion = 10*volt/second * dt, v_peak = -30.*mV)
        spike_onset = v[idx_spike_onset[0]]
        if spike_onset > -30:
            idx_spike_onset = spike_onsets(v*mV, criterion = 5*volt/second * dt, v_peak = -30.*mV)
            spike_onset = v[idx_spike_onset[0]]
            if spike_onset > -30:
                idx_spike_onset = spike_onsets(v*mV, criterion = 0.5*volt/second * dt, v_peak = -30.*mV)
                spike_onset = v[idx_spike_onset[0]]
    
    cs = CubicSpline(v[idx_spike_onset[0]:idx_peak], dv[idx_spike_onset[0]:idx_peak])
    v_new = arange(v[idx_spike_onset[0]], v[idx_peak], 0.1)
    dv_new = cs(v_new)
    
    v = v_new
    dv = dv_new
    ddv = (dv[1:] - dv[:-1])/(dt/ms)
    
    if i > 10:
    
        # smoothing
        #smoothing
        n = len(v)
        i_slide = np.zeros(n)
        d = 60 # half-window, i.e. number of pixels on each side
    
        for j in range(n):
            if j < d: # start of the axon, not full window
                i_slide[j] = np.mean(dv[0:j+d])
            elif j > n-d: # end of the axon, not full window
                i_slide[j] = np.mean(dv[j-d:n])
            else: 
                i_slide[j] = np.mean(dv[j-d:j+d])
                        
        dv = i_slide
        ddv = (dv[1:] - dv[:-1])/(dt/ms)
    
    # somatic regeneration
    
    idx_ax_onset = 0 #idx_spike_onset[0] - 1 # because the function shifts by +1
            
    # global max of dvdt after spike onset
    dvdt_max = argmax(dv[idx_ax_onset:]) + idx_ax_onset
    # global max of the dV^2/dt^2
    ddvdt_max = argmax(ddv[idx_ax_onset:]) + idx_ax_onset
    
    # the global max of dvdt can be in the axonal component:
    # we look for an inflexion pt btw onset and max dvdt:
    # if yes: it is the axonal max, the global max is somatic max
    # if not: the global max is axonal max
    inflexion_before_global_max = where([ddv[i]*ddv[i+1]<0 \
                                         for i in range(idx_ax_onset+1, dvdt_max-2)])[0]
    print(dvdt_max, inflexion_before_global_max + idx_ax_onset+1)
    
    if len(inflexion_before_global_max) < 1: #<= 1: # global max is axonal max
        # the axonal max might not be a local max,
        # so we verifiy that there is no decceleration between spike onset and the max
        if ddvdt_max != idx_ax_onset:
            print('A')
            ddvdt_min = argmin(ddv[idx_ax_onset+1:ddvdt_max+1])+ idx_ax_onset + 1 + 1
        else:
            print('B')
            ddvdt_min = argmin(ddv[idx_ax_onset:ddvdt_max+1])+ idx_ax_onset + 1
        # axonal max
        idx_dvdt_max1 = ddvdt_min
        # somatic max
        idx_dvdt_max2 = dvdt_max
        
        # look for somatic max as next inflexion point
        if len(where([ddv[i]*ddv[i+1]<0  for i in range(dvdt_max+1, idx_peak)])[0]) != 0 : # if another local max after the global max
            print('Global max is axonal max')   
            ddvdt_min = dvdt_max
            extr = where([ddv[i]*ddv[i+1]<0  for i in range(dvdt_max+1, idx_peak)])[0] + dvdt_max + 1 + 1
            dvdt_max = array(extr)[argmax(dv[extr])] 
            # axonal max
            idx_dvdt_max1 = ddvdt_min
            # somatic max
            idx_dvdt_max2 = dvdt_max 
        elif ddvdt_min == ddvdt_max:
            print('C')
            # ddvdt_min = argmin(ddv[ddvdt_max+1:dvdt_max])+ ddvdt_max + 1 + 1
            ddvdt_min = argmin(ddv[idx_ax_onset+1:ddvdt_max])+ ddvdt_max + 1 + 1
            # axonal max
            idx_dvdt_max1 = ddvdt_min
            # somatic max
            idx_dvdt_max2 = dvdt_max 
        elif dv[ddvdt_min] < dv[idx_ax_onset]:
            print('D')
            ddvdt_min = argmin(ddv[idx_ax_onset:dvdt_max+1])+ idx_ax_onset + 1
            # axonal max
            idx_dvdt_max1 = ddvdt_min
            # somatic max
            idx_dvdt_max2 = dvdt_max 
    else:
        print('Global max is somatic max')
        # axonal max
        idx_dvdt_max1 = inflexion_before_global_max[0] + idx_ax_onset + 1 + 1
        # somatic max
        idx_dvdt_max2 = dvdt_max
        
    print(idx_dvdt_max1, idx_dvdt_max2)
        
    # somatic spike onset as the max acceleration between the two local max 
    ddvdt_max_between = argmax(ddv[idx_dvdt_max1:idx_dvdt_max2]) + idx_dvdt_max1 
    idx_som_onset = ddvdt_max_between
            
    somatic_rege = v[idx_som_onset]
    v_regeneration.append(somatic_rege)
    

ax9.plot(arange(0, min_peak), v_onsets, 'k-', label='spike onset')
ax9.plot(arange(0, min_peak-16), v_regeneration, 'k--', label='somatic regeneration')
ax9.set_ylim(-60,0)
# ax9.set_xlim(-60,0)
ax9.set_ylabel('V (mV)')
ax9.set_xlabel('Spike $\#$')
ax9.legend(frameon=False, fontsize=8)
ax9.annotate("H", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

tight_layout()
    
show()


print ('STATS')
print ('N cells:', len(array(charge_attenuation)[~isnan(charge_attenuation)]))
print ('Current attenuation:', nanmean(current_attenuation), '+-', nanstd(current_attenuation))
print ('Charge attenuation:', nanmean(charge_attenuation), '+-', nanstd(charge_attenuation))


save_path = '/Users/sarah/Dropbox/Spike initiation/PhD projects/Axonal current and AIS geometry/Paper/Figures/'

# fig.savefig(save_path + "fig10.pdf", bbox_inches='tight')

# fig.savefig("/Users/sarah/Dropbox/Spike initiation/Thesis/images/fig_rgc_Compensation.pdf", bbox_inches='tight')













