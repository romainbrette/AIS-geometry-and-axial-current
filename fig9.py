

"""

Figure 9: threshold current adaptation.

"""

from brian2 import *
import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile
from scipy.optimize import curve_fit
from scipy.stats import linregress
from matplotlib import gridspec, cm
import seaborn as sns

rcParams['axes.spines.right'] = False
rcParams['axes.spines.top'] = False

### Loading the results of analyses
path_files = '/Users/sarah/Documents/repositories/AIS-geometry-and-axial-current/'

### Threshold current at different V0 (from RGC_adaptation_threshold_current.py)
df_cells = pd.read_excel(path_files + 'RGC_threshold_current_adaptation.xlsx')

dates = array(df_cells['Date'])
retinas = array(df_cells['Retina'])
cells = array(df_cells['Cell'])
ages = array(df_cells['Age'])
recordings = array(df_cells['Sweep'])
holding_potentials = array(df_cells['V0']) 
axonal_currents = array(df_cells['Threshold current'])
threshold_potentials = array(df_cells['Vth']) 

#### Load Rs during adaptation protocol
df_rs = pd.read_excel(path_files + 'RGC_adaptation.xlsx')

#### Counting cells and sorting measures per cell
N = 0

# initialisation
date_prev = dates[0]
retina_prev = retinas[0]
cell_prev = cells[0]

dates_per_cell = [dates[0]]
retina_per_cell = [retinas[0]]
cell_per_cell = [cells[0]]
v0_per_cell = []
vth_per_cell = []
ith_per_cell = []
ia_per_cell = []

v0_cell = []
vth_cell = []
ith_cell = []
ia_cell = []

for date, retina, cell, age, rec, v0, ith, vth in zip(dates, retinas, cells, \
                                        ages, recordings, holding_potentials, \
                                        axonal_currents, threshold_potentials): 
    
    row_rs = df_rs[(df_rs['Date'] == date) & (df_rs['Retina'] == retina) &\
                   ( df_rs['Cell'] == cell) & (df_rs['Recording'] == rec)]
    if len(row_rs) > 0:
        vh = row_rs['Vh'].values[0]
        rs_rec = row_rs['Rs Na rec'].values[0]
        rs_before = row_rs['Rs before'].values[0]
        rs_after = row_rs['Rs after'].values[0]
        ia = row_rs['Peak axonal current corrected'].values[0]
    else:
        continue
    if rs_rec < 25 and  rs_after < 1.3 * rs_before and age > 9:
        if date_prev == date and retina_prev == retina and cell_prev == cell:
            v0_cell.append(v0 + vh)
            vth_cell.append(vth + vh)
            ia_cell.append(ia)
            ith_cell.append(ith)
        else: 
            N += 1
            dates_per_cell.append(date)
            retina_per_cell.append(retina)
            cell_per_cell.append(cell)
            v0_per_cell.append(v0_cell)
            vth_per_cell.append(vth_cell)
            ia_per_cell.append(ia_cell)
            ith_per_cell.append(ith_cell)
            
            v0_cell = [v0 + vh]
            vth_cell = [vth + vh]
            ia_cell = [ia]
            ith_cell = [ith]
        
        date_prev = date
        retina_prev = retina
        cell_prev = cell
    else:
        continue

### Adding the last cell
N += 1
v0_per_cell.append(v0_cell)
vth_per_cell.append(vth_cell)
ia_per_cell.append(ia_cell)
ith_per_cell.append(ith_cell)
    
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
            ith_per_cell[i] = delete(ith_per_cell[i], idx_v0_delete)
        else:
            v0_per_cell[i] = array(v0_per_cell[i])
            ia_per_cell[i] = array(ia_per_cell[i])
            vth_per_cell[i] = array(vth_per_cell[i])
            ith_per_cell[i] = array(ith_per_cell[i])
    
### Threshold current attenuation from -60 to -40
threshold_current_attenuation = []
peak_current_attenuation = []

for i in range(len(dates_per_cell)):
    idx_60 = where(v0_per_cell[i] == -60.)[0]
    idx_40 = where(v0_per_cell[i] == -40.)[0]
    if len(idx_60) > 0 and len(idx_40) > 0:
        threshold_current_attenuation.append(ith_per_cell[i][idx_60[0]]/ith_per_cell[i][idx_40[0]])
        peak_current_attenuation.append(ia_per_cell[i][idx_60[0]]/ia_per_cell[i][idx_40[0]])
    else:
        threshold_current_attenuation.append(nan)
        peak_current_attenuation.append(nan)
        
### IV curves at different V0 (from RGC_adaptation_threshold_current.py) for one exmaple cell
data = load('RGC_IV_curves_below_threshold_adaptation.npz', allow_pickle=True)

dates_iv = data['arr_0']
retinas_iv = data['arr_1']
cells_iv = data['arr_2']
currents_iv= data['arr_3'] 
voltages_iv = data['arr_4']
v_prepulse_iv = data['arr_5']

# n_cells = len(dates)
    
name1 = "tab20b"
name2 = "tab20c"
name3 = "tab20"
cmap1 = get_cmap(name1)  
cmap2 = get_cmap(name2)  
cmap3 = get_cmap(name3) 
cols = cmap1.colors + cmap2.colors + cmap3.colors

### Figure
fig = figure('Threshold current adaptation', figsize=(9, 3))
gs = gridspec.GridSpec(1, 3, width_ratios=[2, 2, 1]) 

### Panel A: example of IV curves below threshold at different V0 in one cell
ax1 = fig.add_subplot(gs[0])
m = 1

idx_sort = argsort(v_prepulse_iv[m])
IV_v0s = v_prepulse_iv[m][idx_sort] - 70
IV_is = currents_iv[m][idx_sort]
IV_vs = voltages_iv[m][idx_sort]

n_sweeps = len(IV_is)

for j in range(n_sweeps):
    if IV_v0s[j] == IV_v0s[j]:
        ax1.plot(IV_vs[j], IV_is[j], '-', c=cols[j+1])
        idx_it = argmin(IV_is[j][-3:]) + len(IV_is[j][:-3])
        ax1.plot(IV_vs[j][idx_it], IV_is[j][idx_it], 'o', c=cols[j+1])

ax1.set_xticks(ticks = [ -60, IV_vs[0][-1], -50, -40])
ax1.set_xticklabels(['-60', '$V_t$', '-50', '-40']) 
ax1.set_yticks(ticks = [ 0.05, 0, -0.05, IV_is[0][-1], -0.1, -0.15])
ax1.set_yticklabels(['0.050', '0', '-0.050', '$I_t$', '-0.100', '-0.150']) 
ax1.plot(linspace(-62, IV_vs[0][-1], 10), IV_is[0][-1]*ones(10), '--', color = 'k', linewidth=1)
ax1.plot(IV_vs[0][-1]*ones(10), linspace(-0.15, IV_is[0][-1], 10),  '--', color = 'k', linewidth=1)
ax1.set_ylim(-0.15,0.0)
ax1.set_xlim(-61, -35)
ax1.set_ylabel('$I$ (nA)')
ax1.set_xlabel('$V$ (mV)')
ax1.annotate("A", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

### Panel B: It vs V0 in an example cell
ax2 = fig.add_subplot(gs[1])
m = -9
idx_sort = argsort(v0_per_cell[m])
for i in range(len(v0_per_cell[m][:-1])): # tor emove the -35 point for which we have no IV curve
    ax2.plot(array(v0_per_cell[m])[idx_sort][i], array(ith_per_cell[m])[idx_sort][i], 'o', color=cols[i+1])
ax2.set_ylim(-0.15, 0)
ax2.set_xlim(-75, -35)
ax2.set_ylabel('$I_t$ (nA)')
ax2.set_xlabel('$V_0$ (mV)')
ax2.annotate("B", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

### Panel C: threshold current attenuation
ax3 = fig.add_subplot(gs[2])

sns.boxplot(y=threshold_current_attenuation, color='gray', ax=ax3)
sns.swarmplot(y=threshold_current_attenuation,  color='0.2', ax=ax3)
sns.despine(bottom=True, ax=ax3)
ax3.set_ylabel('$I_t^{60}/I_t^{40}$')
ax3.set_xticks([])
ax3.set_ylim(0, 5)
ax3.annotate("C", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

tight_layout() 

print ('Stats IT attenuation:', nanmean(threshold_current_attenuation), nanstd(threshold_current_attenuation))

### Saving the figure
save_path = '/Users/sarah/Documents/repositories/AIS-geometry-and-axonal-current/Na currents in RGC/codes submission/data/'

# fig.savefig(save_path + "fig9.pdf", bbox_inches='tight')




