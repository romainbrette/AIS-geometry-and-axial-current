


"""

Analysis of threshold adaptation in RGC.

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

# Load the axonal currents
df_cells = pd.read_excel('RGC_adaptation.xlsx')

# Loading the recording database 
df_rec_database = pd.read_excel('RGC_recording_database.xlsx')

dates = array(df_cells['Date'])
retinas = array(df_cells['Retina'])
cells = array(df_cells['Cell'])
ages = array(df_cells['Age'])
prepulse_potentials = array(df_cells['V prepulse']) 
holding_potentials = array(df_cells['Vh']) 
axonal_currents = array(df_cells['Peak axonal current corrected']) * 1e-3
threshold_potentials = array(df_cells['Vth']) 
rs_during_recording = array(df_cells['Rs Na rec'])

### Counting cells, removing the ones with Vend too far from LJP and sorting measures per cell
N = 0
date_prev = dates[0]
retina_prev = retinas[0]
cell_prev = cells[0]

dates_per_cell = [dates[0]]
retina_per_cell = [retinas[0]]
cell_per_cell = [cells[0]]
v0_per_cell = []
vth_per_cell = []
ia_per_cell = []

v0_cell = []
vth_cell = []
ia_cell = []

for date, retina, cell, age, v0, vh, rs_rec, ia, vth in zip(dates, retinas, cells, ages, \
                    prepulse_potentials, holding_potentials, rs_during_recording, axonal_currents, threshold_potentials): 
    row = df_rec_database[(df_rec_database['Date'] == date) & (df_rec_database['retina'] == retina) & (df_rec_database['cell'] == cell)]
    v_end = row['Vend'].values[0]
    idx_row = row.index[0]
        
    if ljp-3 <= v_end <= ljp+3 or v_end != v_end:
        if rs_rec < 25. :
            if date_prev == date and retina_prev == retina and cell_prev == cell: # same cell
                v0_cell.append(v0)
                vth_cell.append(vth)
                ia_cell.append(ia)
            else: # next cell
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
    
### Removing doubles
for i in range(N):
    # Removing the recordings with same v0
    #print (ia_per_cell[i])
    if v0_per_cell[i][0] == -60. and v0_per_cell[i][-1] == -60.:
        v0_per_cell[i] = delete(v0_per_cell[i], -1)
        ia_per_cell[i] = delete(ia_per_cell[i], -1)
        vth_per_cell[i] = delete(vth_per_cell[i], -1)
    else:
        v0_per_cell[i] = v0_per_cell[i]
        ia_per_cell[i] = ia_per_cell[i]
        vth_per_cell[i] = vth_per_cell[i]
    
### Threshold adaptation

name1 = "tab20b"
name2 = "tab20c"
name3 = "tab20"
cmap1 = get_cmap(name1)  
cmap2 = get_cmap(name2)  
cmap3 = get_cmap(name3) 
cols = cmap1.colors + cmap2.colors + cmap3.colors
    
# Fitting each cells's data
f = lambda v0, c, Ka, Ki, Vh : c + Ka*log(1.+exp((v0-Vh)/Ki))
x = linspace(-75., -30.,50)

cst_cells = zeros(N)
ka_cells = zeros(N)
ki_cells = zeros(N)
vh_cells = zeros(N)
kh_cells = zeros(N)

for i in range(len(dates_per_cell)): 
    
    if len(v0_per_cell[i]) > 5: # and (dates_per_cell[i], retina_per_cell[i], cell_per_cell[i]) != (20191121, 'A',1):
        fig = figure('Adaptation %i,  %s, %i' %(dates_per_cell[i], retina_per_cell[i], cell_per_cell[i]), (4,3.5)) 
        ax2 = fig.add_subplot(111)
        ax2.set_title('%i,  %s, %i' %(dates_per_cell[i], retina_per_cell[i], cell_per_cell[i]))
        
        popt,_ = curve_fit(f, v0_per_cell[i], vth_per_cell[i], p0=[-55., 5., 5., -55.],\
                           bounds = ([-70, 0, 0, -70], [-40, 10, 10,-40]))
            

        cst, ka, ki, vh = popt
        ax2.plot(x, f(x, cst, ka, ki, vh), color=cols[i], label='fit: C=%0.1f, Vh=%0.1f,  ka=%0.1f, ki=%0.1f' %(cst, vh, ka, ki))
                
        cst_cells[i] = cst
        ka_cells[i] = ka
        ki_cells[i] = ki
        vh_cells[i] = vh            

        ax2.plot(v0_per_cell[i], vth_per_cell[i], 'o', color = cols[i])
        ax2.plot(linspace(-80, -25,100), linspace(-80, -25,100), 'k--')
        # ax2.plot(v0_per_cell[i], result.best_fit, 'r.')
        ax2.set_ylim(-80, -30)
        ax2.set_xlim(-80, -30)
        ax2.set_ylabel('Threshold (mV)')
        ax2.set_xlabel('Holding potential (mV)')
        ax2.legend(frameon=False, fontsize=8)
        tight_layout()
    else:
        cst_cells[i] =nan
        ka_cells[i] =nan
        ki_cells[i] =nan
        vh_cells[i] =nan
            

# savez('RGC_threshold_adaptation_test', v0_per_cell[m], vth_per_cell[m], vh_cells, cst_cells, \
#                           ka_cells, ki_cells, outliers, dates_per_cell, retina_per_cell, cell_per_cell)


# ### Save fit results in DF

# df_select_cells = pd.DataFrame({'Date': dates_per_cell,
#                   'Retina': retina_per_cell,
#                   'Cell': cell_per_cell,
#                   'k': ka_cells,
#                   'Vi':vh_cells
#                   })

# save_path = '/Users/sarahgoethals/Dropbox/Spike initiation/PhD projects/Axonal current and AIS geometry/Data patch/RGC/'
# df_select_cells.to_excel(save_path + "RGC_Vi_from_adaptation_test.xlsx", \
#         columns=['Date','Retina','Cell','k','Vi'])




