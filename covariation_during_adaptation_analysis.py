


"""

Analysis of voltage threshold and peak axonal current covariation.

OK

"""

from brian2 import *
import pandas as pd
import seaborn as sns
from pandas import ExcelWriter
from pandas import ExcelFile
from scipy.optimize import curve_fit
from scipy.stats import linregress
from matplotlib import gridspec

ljp = -11.

### Loading the results of analyses
df_cells = pd.read_excel('RGC_adaptation.xlsx')

# Loading the recording database 
df_rec_database = pd.read_excel('RGC_recording_database.xlsx')

dates = array(df_cells['Date'])
retinas = array(df_cells['Retina'])
cells = array(df_cells['Cell'])
ages = array(df_cells['Age'])
prepulse_potentials = array(df_cells['V prepulse']) 
holding_potentials = array(df_cells['Vh']) 
axonal_currents = array(df_cells['Peak axonal current corrected']) #* 1e-3
threshold_potentials = array(df_cells['Vth']) 

Rs_before = array(df_cells['Rs before']) 
Rs_after = array(df_cells['Rs after'])
Rs_rec = array(df_cells['Rs Na rec'])  

selected_dates = []
selected_retinas = []
selected_cells = []
selected_ages = []
selected_prepulse_pot = []
selected_axonal_currents = []
selected_threshold_potentials = []

for date, retina, cell, age, v0, vh, rs_before, rs_after, rs_rec, ia, vth in zip(dates, retinas, cells, ages, \
                                prepulse_potentials, holding_potentials, Rs_before, Rs_after, Rs_rec,\
                                    axonal_currents, threshold_potentials): 
    row = df_rec_database[(df_rec_database['Date'] == date) & (df_rec_database['retina'] == retina) & (df_rec_database['cell'] == cell)]
    v_end = row['Vend'].values[0]
    idx_row = row.index[0]

    if ljp-3 <= v_end <= ljp+3 or v_end != v_end:
        if rs_rec < 25 and  rs_after < 1.3 * rs_before:
            selected_dates.append(date)
            selected_retinas.append(retina)
            selected_cells.append(cell)
            selected_ages.append(age)
            selected_prepulse_pot.append(v0)
            selected_axonal_currents.append(ia)
            selected_threshold_potentials.append(vth)

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

v0_cell = []
vth_cell = []
ia_cell = []

for date, retina, cell, age, v0, ia, vth in zip(selected_dates, selected_retinas, selected_cells, \
                                                selected_ages, selected_prepulse_pot, \
                                                selected_axonal_currents, selected_threshold_potentials): 

    if date_prev == date and retina_prev == retina and cell_prev == cell:
        v0_cell.append(v0)
        vth_cell.append(vth)
        ia_cell.append(ia)
    else: 
        N += 1
        dates_per_cell.append(date)
        retina_per_cell.append(retina)
        cell_per_cell.append(cell)
        v0_per_cell.append(v0_cell)
        vth_per_cell.append(vth_cell)
        ia_per_cell.append(ia_cell)
        
        v0_cell = [v0]
        vth_cell = [vth]
        ia_cell = [ia]
    
    date_prev = date
    retina_prev = retina
    cell_prev = cell

### Adding the last cell
N += 1
v0_per_cell.append(v0_cell)
vth_per_cell.append(vth_cell)
ia_per_cell.append(ia_cell)
    
thresholds_all = []
axonal_currents_all = []
axonal_currents_norm_all = []
prepulse_pot_all = []

for i in range(N):
    # Removing the recordings with same v0
    if v0_per_cell[i][0] == -60. and v0_per_cell[i][-1] == -60.:
        v0_per_cell[i] = delete(v0_per_cell[i], -1)
        ia_per_cell[i] = delete(ia_per_cell[i], -1)
        vth_per_cell[i] = delete(vth_per_cell[i], -1)
    else:
        v0_per_cell[i] = v0_per_cell[i]
        ia_per_cell[i] = ia_per_cell[i]
        vth_per_cell[i] = vth_per_cell[i]
        
    thresholds_all += list(vth_per_cell[i])
    axonal_currents_all += list(ia_per_cell[i])
    axonal_currents_norm_all += list(ia_per_cell[i]/min(ia_per_cell[i]))
    prepulse_pot_all += list(v0_per_cell[i])
    
### Threshold and peak axonal current covariation

slope_per_cell = []
intercept_per_cell = []
r_per_cell = []

name1 = "tab20b"
name2 = "tab20c"
name3 = "tab20"
cmap1 = get_cmap(name1)  
cmap2 = get_cmap(name2)  
cmap3 = get_cmap(name3) 
cols = cmap1.colors + cmap2.colors + cmap3.colors

fig_vth_ia = figure('Threshold and peak current', figsize=(10,4))
gs = gridspec.GridSpec(1, 3, width_ratios=[2, 2, 1]) 
ax1 = fig_vth_ia.add_subplot(gs[0])
ax2 = fig_vth_ia.add_subplot(gs[1])

for i in range(N):
    if len(ia_per_cell[i]) > 5:
        slope, intercept, r_value, p_value, std_err = linregress(-log(-array(ia_per_cell[i])), vth_per_cell[i])
        slope_per_cell.append(slope/2)
        r_per_cell.append(r_value)
        intercept_per_cell.append(intercept)
        ax1.semilogx(-array(ia_per_cell[i]), vth_per_cell[i], 'o', color=cols[i])
    else:
        slope_per_cell.append(nan)
        r_per_cell.append(nan)
        intercept_per_cell.append(nan)
        
ax1.set_ylabel('$V_t$ (mV)')
ax1.set_xlabel('$-I_p$ (nA)')
ax1.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax1.set_ylim(-65,-25)
ax1.set_xlim(0,12)
ax1.legend(frameon=False, fontsize=6)

m = 8
iii = linspace(min(-array(ia_per_cell[m])), max(-array(ia_per_cell[m])), 100)
ax2.semilogx(iii, intercept_per_cell[m] - 2*slope_per_cell[m]*log(iii), '-', color=cols[m], \
             alpha=0.5, label='$k_a$= %0.02f mV, r= %0.02f' %(slope_per_cell[m], r_per_cell[m]))
ax2.semilogx(-array(ia_per_cell[m]), vth_per_cell[m], 'o', color=cols[m])
ax2.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax2.set_ylim(-65, -25)
ax2.set_xlabel('$-I_p$ (nA)')
ax2.legend(frameon=False)

ax3 = fig_vth_ia.add_subplot(gs[2])
sns.boxplot(y=slope_per_cell, color='gray')
sns.swarmplot(y=slope_per_cell,  color='0.2')
ax3.set_ylabel('$k_a$ (mV)')
ax3.set_xticks([])
ax3.set_ylim(0, 6)

tight_layout() 

# savez('RGC_adaptation_covariation', ia_per_cell, vth_per_cell, v0_per_cell)








