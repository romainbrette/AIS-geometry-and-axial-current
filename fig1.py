

"""

Figure 1: Spontaneous APs of retinal ganglion cells.

"""

from brian2 import *
import glob2
import pandas as pd
import pyabf
import seaborn as sns
from scipy import interpolate
from scipy import stats
from pandas import ExcelWriter
from pandas import ExcelFile
from matplotlib import gridspec
from trace_analysis import *

### Figure parameters
rcParams['axes.spines.right'] = False
rcParams['axes.spines.top'] = False

name1 = "tab20b"
name2 = "tab20c"
cmap1 = get_cmap(name1)  
cmap2 = get_cmap(name2)  
cols = cmap1.colors + cmap2.colors 

fig = figure('AP analysis', (6,7))

gs = gridspec.GridSpec(3, 4, width_ratios=[1, 1, 1, 1]) 
ax1 = fig.add_subplot(gs[0, 0:2])
ax2 = fig.add_subplot(gs[0, 2:4])
ax3 = fig.add_subplot(gs[1, 0:2])
ax4 = fig.add_subplot(gs[1, 2:4])
ax5 = fig.add_subplot(gs[2, 0:2])
ax7 = fig.add_subplot(gs[2, 2:4])

##### PANEL A-D: example of a spontaneous AP of a RGC

### Load the list of cells used for the analysis
df_cells = pd.read_excel('RGC_electrical_properties.xlsx')
idx_cell = -4 # the cell shown in exmaple in the figure

date = array(df_cells['Date'])[idx_cell]
retina = array(df_cells['Retina'])[idx_cell]
cell = array(df_cells['Cell'])[idx_cell]

### Path to the data
# path_to_data = 'data/RGC data/'
path_to_data = '/Users/sarah/Documents/Data/Martijn Sierksma/'
path_to_cell = path_to_data + str(int(date)) + "*" + '/retina '+ str(retina) +'/cell ' + str(int(cell))
path_to_cc_cont = glob2.glob(path_to_cell + '/CC cont/' + '*' + ".abf")
print (path_to_cc_cont)

### Loading and analysing the data
abf = pyabf.ABF(path_to_cc_cont[0]) 
V = abf.sweepY * mV
fs = abf.dataRate  * Hz # sampling rate
dt = 1./fs
t = dt*arange(len(V)) 

### Analyze AP shape
# First spike
spike_times = find_spikes_at(V, dt/ms, -20.*mV)/(dt/ms)
idx_spike1 = int(spike_times[0] )

# Interpolating V and computing dV and d2V
f = V[idx_spike1-70:idx_spike1+10]
t_spike = t[idx_spike1-70:idx_spike1+10]
t_new = (t_spike[:-1] + t_spike[1:])/2
v = (f[:-1] + f[1:])/2
dv = (f[1:] - f[:-1])/dt
ddv = (dv[1:] - dv[:-1])/dt # shift of dt: add 1 !!!

# AP peak
idx_peak = argmax(v)

# Spike onset 
spike_onset = spike_onsets(v, criterion = 20*volt/second * dt, v_peak = -30.*mV)

print ('Onset:', spike_onset)
idx_ax_onset = spike_onset[0] - 1 # because the function shifts by +1

# Global max of dvdt after spike onset
dvdt_max = argmax(dv[idx_ax_onset:]) + idx_ax_onset
# Global max of the dV^2/dt^2
ddvdt_max = argmax(ddv[idx_ax_onset:]) + idx_ax_onset

inflexion_before_global_max = where([ddv[i]*ddv[i+1]<0 \
                                     for i in range(idx_ax_onset+1, dvdt_max-2)])[0]
# Axonal max
idx_dvdt_max1 = inflexion_before_global_max[0] + idx_ax_onset + 1 + 1
# Somatic max
idx_dvdt_max2 = dvdt_max
    
print(idx_dvdt_max1, idx_dvdt_max2)
    
# Somatic regeneration as the max acceleration between the two local max 
ddvdt_max_between = argmax(ddv[idx_dvdt_max1:idx_dvdt_max2]) + idx_dvdt_max1 
idx_som_onset = ddvdt_max_between
            
t_dvdt_max1 = t_new[idx_dvdt_max1]/ms
t_dvdt_max2 = t_new[idx_dvdt_max2]/ms
dv_dvdt_max1 = dv[idx_dvdt_max1]/(mV/ms)
dv_dvdt_max2 = dv[idx_dvdt_max2]/(mV/ms)
v_dvdt_max1 = v[idx_dvdt_max1]/mV
v_dvdt_max2 = v[idx_dvdt_max2]/mV

t_ax_onset = t_new[idx_ax_onset]/ms
v_ax_onset = v[idx_ax_onset]/mV
dvdt_ax_onset = dv[idx_ax_onset]/(mV/ms)

t_som_onset = t_new[idx_som_onset]/ms
v_som_onset = v[idx_som_onset]/mV
dvdt_som_onset = dv[idx_som_onset]/(mV/ms)
               
### Plotting
# Panel A: V vs t
ax1.plot(t_new/ms - t_ax_onset, v/mV, 'k-')
ax1.plot(t_ax_onset- t_ax_onset, v_ax_onset, 'o', c=cols[2], label='spike onset', markerfacecolor='none')
ax1.plot(t_dvdt_max1 - t_ax_onset, v_dvdt_max1, 'o',  c=cols[10], label='first max dV/dt', markerfacecolor='none')
ax1.plot(t_dvdt_max2 - t_ax_onset, v_dvdt_max2, 'o',  c=cols[6], label='second max dV/dt', markerfacecolor='none')
ax1.plot(t_som_onset - t_ax_onset, v_som_onset, 'o',  c=cols[14], label='som. regeneration', markerfacecolor='none')
ax1.legend(frameon=False, fontsize=6)
ax1.set_xlim(-0.5, 1.5)
ax1.set_ylim(-60, 35)
ax1.set_ylabel('$V$ (mV)')
# ax1.set_xlabel('$t$ (ms)')
ax1.set_xticks([])
sns.despine(bottom=True, ax=ax1)
ax1.plot(linspace(-0.25,.25,10), -59.*ones(10), 'k-', linewidth=2)
ax1.text(-0.15, -67,'0.5 ms', color='k', fontsize=8)
ax1.annotate("A", xy=(0,1.1), xycoords="axes fraction",
                xytext=(5,-5), textcoords="offset points",
                ha="left", va="top",
                fontsize=12, weight='bold')

# Panel B: dV vs t
ax2.plot(t_new/ms- t_ax_onset, dv, 'k-')
ax2.plot(t_ax_onset- t_ax_onset, dvdt_ax_onset, 'o', c=cols[2], markerfacecolor='none')
ax2.plot(t_dvdt_max2- t_ax_onset, dv_dvdt_max2, 'o', c=cols[6], markerfacecolor='none')
ax2.plot(t_dvdt_max1- t_ax_onset, dv_dvdt_max1, 'o', c=cols[10], markerfacecolor='none')
ax2.plot(t_som_onset- t_ax_onset, dvdt_som_onset, 'o', c=cols[14], markerfacecolor='none')
ax2.set_ylabel('$dV/dt$ (mV/ms)')
# ax2.set_xlabel('$t$ (ms)')
ax2.set_xlim(-0.5, 1.5)
ax2.set_ylim(-100, 350)
ax2.set_xticks([])
sns.despine(bottom=True, ax=ax2)
ax2.plot(linspace(-0.25,.25,10), -95.*ones(10), 'k-', linewidth=2)
ax2.text(-0.15, -135,'0.5 ms', color='k', fontsize=8)
ax2.annotate("B", xy=(0,1.1), xycoords="axes fraction",
                xytext=(5,-5), textcoords="offset points",
                ha="left", va="top",
                fontsize=12, weight='bold')

# Panel C: dV vs V   
ax3.plot(v/mV, dv, 'k-')
ax3.plot(v_ax_onset*ones(50), linspace(-100, 320,50), 'k--', alpha=0.3)
ax3.plot(v_som_onset*ones(50), linspace(-100, 320,50), 'k--', alpha=0.3)
ax3.annotate('', xy=(v_ax_onset, 250), xytext=(v_som_onset, 250),
        arrowprops=dict(arrowstyle='<->', color='k', ls='-', alpha=0.3))
ax3.text((v_ax_onset + v_som_onset)/2 - 5, 275, "$\Delta V$", alpha=0.3)
ax3.plot(v_ax_onset, dvdt_ax_onset, 'o', c=cols[2], markerfacecolor='none')
ax3.plot(v_dvdt_max2, dv_dvdt_max2, 'o', c=cols[6], markerfacecolor='none')
ax3.plot(v_dvdt_max1, dv_dvdt_max1, 'o', c=cols[10], markerfacecolor='none')
ax3.plot(v_som_onset, dvdt_som_onset, 'o', c=cols[14], markerfacecolor='none')
ax3.set_ylabel('$dV/dt$ (mV/ms)')
ax3.set_xlabel('$V$ (mV)')
ax3.set_ylim(-100, 350)
ax3.annotate("C", xy=(0,1.1), xycoords="axes fraction",
                xytext=(5,-5), textcoords="offset points",
                ha="left", va="top",
                fontsize=12, weight='bold')

# Panel D: d2V vs V
ax4.plot(t_new[:-1]/ms- t_ax_onset, ddv/(mV/ms**2), 'k-')
ax4.plot(t_ax_onset- t_ax_onset, ddv[idx_ax_onset]/(mV/ms**2), 'o', c=cols[2], markerfacecolor='none')
ax4.plot(t_dvdt_max2- t_ax_onset, ddv[idx_dvdt_max2]/(mV/ms**2), 'o', c=cols[6], markerfacecolor='none')
ax4.plot(t_dvdt_max1- t_ax_onset, ddv[idx_dvdt_max1]/(mV/ms**2), 'o', c=cols[10], markerfacecolor='none')
ax4.plot(t_som_onset- t_ax_onset, ddv[idx_som_onset]/(mV/ms**2), 'o', c=cols[14], markerfacecolor='none')
ax4.set_ylabel('$d^2V/dt^2$ (mV/ms$^2$)')
# ax4.set_xlabel('$t$ (ms)')
ax4.set_xlim(-0.5, 1.5)
ax4.set_ylim(-2700, 2100)
ax4.set_xticks([])
sns.despine(bottom=True, ax=ax4)
ax4.plot(linspace(-0.25,.25,10), -2600.*ones(10), 'k-', linewidth=2)
ax4.text(-0.15, -3000, '0.5 ms', color='k', fontsize=8)
ax4.annotate("D", xy=(0,1.1), xycoords="axes fraction",
                xytext=(5,-5), textcoords="offset points",
                ha="left", va="top",
                fontsize=12, weight='bold')

##### Panel E-F: statistics over the neuron population

### Load the results of the AP shape analysis 
df_cells = pd.read_excel('RGC_action_potential.xlsx')

ljp = -11. # liquid junction potential

dates = array(df_cells['Date'])
retinas = array(df_cells['Retina'])
cells = array(df_cells['Cell'])
ages = array(df_cells['Age'])
vend = array(df_cells['Vend'])

### If the reference potential changed too much during teh recording, the cell is discarded
idx = where([ljp-3 <= vend[i] <= ljp+3 for i in range(len(vend))])

ap_onsets = array(df_cells['AP onset'])
dvdt_max1 = array(df_cells['dvdt max1']) 
dvdt_max2 = array(df_cells['dvdt max2'])
v_max1 = array(df_cells['v max1'])
v_max2 = array(df_cells['v max2'])
dvdt_regeneration = array(df_cells['dvdt somatic onset'])
v_regeneration = array(df_cells['v somatic onset'])

### Statistics
print('Spike onset:', mean(ap_onsets[idx]), std(ap_onsets[idx]))
print('Spike regeneration:', mean(v_regeneration[idx]), std(v_regeneration[idx]))
print('Delta V', mean(v_regeneration[idx] - ap_onsets[idx]), std(v_regeneration[idx] - ap_onsets[idx]))
              
### Plotting

# Panel E: somatic regeneration vs spike onset
slope, intercept, r, p, _ = stats.linregress(ap_onsets[idx], v_regeneration[idx]) 

sns.scatterplot(x=ap_onsets[idx], y=v_regeneration[idx], color='0.2', ax=ax5)
ax5.plot(linspace(-60,-5, 10), linspace(-60, -5, 10), 'k--')
ax5.set_ylim(-60, 0)
ax5.set_xlim(-60, 0)
ax5.set_ylabel('Somatic regeneration (mV)')
ax5.set_xlabel('Spike onset (mV)')
ax5.annotate("E", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                     ha="left", va="top",
                     fontsize=12, weight='bold')

# Panel F: difference between the spike onset and the somatic regeneration
sns.boxplot(y=v_regeneration[idx] - ap_onsets[idx], color='gray', width=0.3, ax=ax7)
sns.swarmplot(y=v_regeneration[idx] - ap_onsets[idx],  color='0.2', ax=ax7)
sns.despine(bottom=True, ax=ax7)
ax7.set_ylabel('$\Delta V$ (mV)')
ax7.set_ylim(20, 40)
ax7.set_xticks([])
ax7.annotate("F", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

tight_layout()

show()

print ('N cells:', len(ap_onsets))
print ('N cells with constraint on Vend:', len(ap_onsets[idx]))


### Saving the figure
save_path = '/Users/sarah/Dropbox/Spike initiation/PhD projects/Axonal current and AIS geometry/Paper/Figures/'
# fig.savefig(save_path + "fig1.pdf", bbox_inches='tight')

# fig.savefig("fig1.pdf", bbox_inches='tight')

# fig.savefig(save_path + "fig1.png", dpi=300)



# ### Figure thesis

# fig_bis = figure('AP thesis', (6,2.5))

# ax1 = fig_bis.add_subplot(121)
# ax1.plot(t_new/ms - t_ax_onset, v/mV, 'k-')
# ax1.plot(t_ax_onset- t_ax_onset, v_ax_onset, 'o', c=cols[2], label='spike onset', markerfacecolor='none')
# ax1.plot(t_dvdt_max1 - t_ax_onset, v_dvdt_max1, 'o',  c=cols[10], label='first max dV/dt', markerfacecolor='none')
# # ax1.plot(t_dvdt_max2 - t_ax_onset, v_dvdt_max2, 'o',  c=cols[6], label='second max dV/dt', markerfacecolor='none')
# ax1.plot(t_som_onset - t_ax_onset, v_som_onset, 'o',  c=cols[14], label='som. regeneration', markerfacecolor='none')
# ax1.legend(frameon=False, fontsize=6)
# ax1.set_xlim(-0.5, 1.5)
# ax1.set_ylim(-60, 35)
# ax1.set_ylabel('$V$ (mV)')
# # ax1.set_xlabel('$t$ (ms)')
# ax1.set_xticks([])
# sns.despine(bottom=True, ax=ax1)
# ax1.plot(linspace(-0.25,.25,10), -59.*ones(10), 'k-', linewidth=2)
# ax1.text(-0.15, -67,'0.5 ms', color='k', fontsize=8)



# ax3 = fig_bis.add_subplot(122)
# ax3.plot(v/mV, dv, 'k-')
# ax3.plot(v_ax_onset*ones(50), linspace(-100, 320,50), 'k--', alpha=0.3)
# ax3.plot(v_som_onset*ones(50), linspace(-100, 320,50), 'k--', alpha=0.3)
# ax3.annotate('', xy=(v_ax_onset, 250), xytext=(v_som_onset, 250),
#         arrowprops=dict(arrowstyle='<->', color='k', ls='-', alpha=0.3))
# ax3.text((v_ax_onset + v_som_onset)/2 - 5, 275, "$\Delta V$", alpha=0.3)
# ax3.plot(v_ax_onset, dvdt_ax_onset, 'o', c=cols[2], markerfacecolor='none')
# # ax3.plot(v_dvdt_max2, dv_dvdt_max2, 'o', c=cols[6], markerfacecolor='none')
# ax3.plot(v_dvdt_max1, dv_dvdt_max1, 'o', c=cols[10], markerfacecolor='none')
# ax3.plot(v_som_onset, dvdt_som_onset, 'o', c=cols[14], markerfacecolor='none')
# ax3.set_ylabel('$dV/dt$ (mV/ms)')
# ax3.set_xlabel('$V$ (mV)')
# ax3.set_ylim(-100, 350)


# tight_layout()


# save_path = '/Users/sarah/Dropbox/Spike initiation/Thesis/images/'
# fig_bis.savefig(save_path + "ppt_RGC_AP_with_dvmax.pdf", bbox_inches='tight')






