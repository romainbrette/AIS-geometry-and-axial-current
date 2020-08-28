



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
    
### Figure

name1 = "tab20b"
name2 = "tab20c"
name3 = "tab20"
cmap1 = get_cmap(name1)  
cmap2 = get_cmap(name2)  
cmap3 = get_cmap(name3) 
cols = cmap1.colors + cmap2.colors + cmap2.colors 

fig = figure('Cont', figsize=(9,3))
gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 1]) 

# ax1 = fig.add_subplot(gs[0])
# ax2 = fig.add_subplot(gs[1])
# ax3 = fig.add_subplot(gs[2])
# ax4 = fig.add_subplot(gs[3])
# ax5 = fig.add_subplot(gs[4:6])
ax6 = fig.add_subplot(gs[0])
ax8 = fig.add_subplot(gs[1])
ax9 = fig.add_subplot(gs[2])
# ax9 = fig.add_subplot(gs[11])

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
# ax6.set_xlabel('t (s)')
ax6.plot(linspace(112,112.1,10), -75.*ones(10), 'k-', linewidth=2)
ax6.text(111.85, -85.,'100 ms',color='k', fontsize=8)
ax6.set_xticks([])
sns.despine(bottom=True, ax=ax6)
ax6.annotate("F", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

# Panel G: phase plots
# for i in range(int((min_peak-15)/3)+1):
    
#     ax6.plot(spike_times[3*i]/second, 35, '|', color=cols[3*i])
    
#     idx_spike = int(spike_times[3*i]/dt)
    
#     # Measures
#     f = data[idx_spike-200:idx_spike+100]
#     t_spike = t[idx_spike-200:idx_spike+100]/ms - t[idx_spike-200]/ms
    
#     t_new = (t_spike[:-1] + t_spike[1:])/2
#     v = (f[:-1] + f[1:])/2
#     dv = (f[1:] - f[:-1])/(dt/ms)
#     ddv = (dv[1:] - dv[:-1])/(dt/ms) # shift of dt: add 1 !!! (f[2:] - 2*f[1:-1] + f[:-2])/dt**2 #
    
#     ax8.plot(v, dv, color=cols[3*i])
    
#     #AP peak
#     idx_peak = argmax(v)
    
#     # spike onset 
#     idx_spike_onset = spike_onsets(v*mV, criterion = 20*volt/second * dt, v_peak = -30.*mV)
#     spike_onset = v[idx_spike_onset[0]]
#     if spike_onset > -30:
#         idx_spike_onset = spike_onsets(v*mV, criterion = 10*volt/second * dt, v_peak = -30.*mV)
#         spike_onset = v[idx_spike_onset[0]]
#         if spike_onset > -30:
#             idx_spike_onset = spike_onsets(v*mV, criterion = 5*volt/second * dt, v_peak = -30.*mV)
#             spike_onset = v[idx_spike_onset[0]]
#             if spike_onset > -30:
#                 idx_spike_onset = spike_onsets(v*mV, criterion = 0.5*volt/second * dt, v_peak = -30.*mV)
#                 spike_onset = v[idx_spike_onset[0]]
    
#     # interpolation
#     idx_max = argmax(dv)
#     cs = CubicSpline(v[idx_spike_onset[0]:idx_max+1], dv[idx_spike_onset[0]:idx_max+1])
#     v_new = arange(v[idx_spike_onset[0]], v[idx_max+1], 0.1)
#     ax8.plot(v_new, cs(v_new), 'r--')
    
    
ax8.set_ylabel('dV/dt (mV/ms)')
ax8.set_xlabel('V (mV)')
ax8.annotate("G", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')


v_onsets = []
v_regeneration = []
# for i in range(min_peak):
#     idx_spike = int(spike_times[i]/dt)
    
#     # Measures
#     f = data[idx_spike-200:idx_spike+100]
#     t_spike = t[idx_spike-200:idx_spike+100]/ms - t[idx_spike-200]/ms
#     t_new = (t_spike[:-1] + t_spike[1:])/2
#     v = (f[:-1] + f[1:])/2
#     dv = (f[1:] - f[:-1])/(dt/ms)
#     ddv = (dv[1:] - dv[:-1])/(dt/ms) # shift of dt: add 1 !!! (f[2:] - 2*f[1:-1] + f[:-2])/dt**2 #
    
#     # AP peak
#     idx_peak = argmax(v)
    
#     # spike onset 
#     idx_spike_onset = spike_onsets(v*mV, criterion = 20*volt/second * dt, v_peak = -30.*mV)
#     spike_onset = v[idx_spike_onset[0]]
#     if spike_onset > -30:
#         idx_spike_onset = spike_onsets(v*mV, criterion = 10*volt/second * dt, v_peak = -30.*mV)
#         spike_onset = v[idx_spike_onset[0]]
#         if spike_onset > -30:
#             idx_spike_onset = spike_onsets(v*mV, criterion = 5*volt/second * dt, v_peak = -30.*mV)
#             spike_onset = v[idx_spike_onset[0]]
#             if spike_onset > -30:
#                 idx_spike_onset = spike_onsets(v*mV, criterion = 0.5*volt/second * dt, v_peak = -30.*mV)
#                 spike_onset = v[idx_spike_onset[0]]
#     v_onsets.append(spike_onset)
        
for i in range(min_peak-10):
    idx_spike = int(spike_times[i]/dt)
    
    # Measures
    f = data[idx_spike-200:idx_spike+100]
    t_spike = t[idx_spike-200:idx_spike+100]/ms - t[idx_spike-200]/ms
    t_new = (t_spike[:-1] + t_spike[1:])/2
    v = (f[:-1] + f[1:])/2
    dv = (f[1:] - f[:-1])/(dt/ms)
    ddv = (dv[1:] - dv[:-1])/(dt/ms) # shift of dt: add 1 !!! (f[2:] - 2*f[1:-1] + f[:-2])/dt**2 #
    
    ax8.plot(v, dv, color=cols[i])
    
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
    
    # if i < 10:
    # interpolation
    # idx_max = argmax(dv)
    cs = CubicSpline(v[idx_spike_onset[0]:idx_peak], dv[idx_spike_onset[0]:idx_peak])
    v_new = arange(v[idx_spike_onset[0]], v[idx_peak], 0.1)
    dv_new = cs(v_new)
    ax8.plot(v_new, dv_new, 'r--')
    
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
                
        ax8.plot(v, i_slide, 'g--')
        
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
    

ax9.plot(arange(0, min_peak-16), v_onsets, 'k-', label='spike onset')
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

save_path = '/Users/sarah/Dropbox/Spike initiation/PhD projects/Axonal current and AIS geometry/Paper/Figures/'

# fig.savefig(save_path + "fig_Current_attenuation_CCcont.pdf", bbox_inches='tight')

# fig.savefig("/Users/sarah/Dropbox/Spike initiation/Thesis/images/fig_rgc_Compensation.pdf", bbox_inches='tight')













