


"""

Figure 8: adaptation in RGC.


"""

from brian2 import *
import pandas as pd
import seaborn as sns
import params_model_description
from pandas import ExcelWriter
from pandas import ExcelFile
from scipy.optimize import curve_fit
from scipy.stats import linregress
from matplotlib import gridspec

rcParams['axes.spines.right'] = False
rcParams['axes.spines.top'] = False

params = params_model_description

### Loading the results of analyses

# Load the axonal currents
df_cells = pd.read_excel('RGC_adaptation.xlsx')

# Load the Vi from fitting threshold adaptation
df_fit = pd.read_excel('RGC_Vi_from_adaptation.xlsx')

dates = array(df_cells['Date'])
retinas = array(df_cells['Retina'])
cells = array(df_cells['Cell'])
ages = array(df_cells['Age'])

prepulse_potentials = array(df_cells['V prepulse']) 
holding_potentials = array(df_cells['Vh']) 

axonal_currents = array(df_cells['Peak axonal current']) * 1e-3
axonal_currents_corr = array(df_cells['Peak axonal current corrected']) #* 1e-3
threshold_potentials = array(df_cells['Vth']) 
charges = array(df_cells['Charge1 10']) # first peak above threshold because for some recordings, the seocnd is somatic current 
duration50= array(df_cells['Duration1 10']) 

Rs_before = array(df_cells['Rs before']) 
Rs_after = array(df_cells['Rs after'])
Rs_rec = array(df_cells['Rs Na rec'])  

selected_dates = []
selected_retinas = []
selected_cells = []
selected_ages = []
selected_prepulse_pot = []
selected_axonal_currents = []
selected_axonal_currents_corr = []
selected_threshold_potentials = []
selected_charge = []
selected_dur50 = []

### Removing the recordings for which the series resistance is larger than 25 Mohm 
### and the cells for which the series resistance changed too much during the adaptation protocol
for date, retina, cell, age, v0, vh, rs_before, rs_after, rs_rec, ia, ia_corr, vth, charge, dur50 in zip(dates, retinas, cells, ages, \
                                prepulse_potentials, holding_potentials, Rs_before, Rs_after, Rs_rec,\
                                    axonal_currents, axonal_currents_corr, threshold_potentials,\
                                    charges, duration50): 
    if rs_rec < 25 and rs_after < 1.3 * rs_before:
        selected_dates.append(date)
        selected_retinas.append(retina)
        selected_cells.append(cell)
        selected_ages.append(age)
        selected_prepulse_pot.append(v0)
        selected_axonal_currents.append(ia)
        selected_axonal_currents_corr.append(ia_corr)
        selected_threshold_potentials.append(vth)
        selected_charge.append(charge)
        selected_dur50.append(dur50)

### Counting cells and sorting measures of threshold and peak current per cell
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
ia_corr_per_cell = []
charge_per_cell = []
dur50_per_cell = []

v0_cell = []
vth_cell = []
ia_cell = []
ia_corr_cell = []
ch_cell = []
dur50_cell = []

for date, retina, cell, age, v0, ia, ia_corr, vth, charge, dur50 in zip(selected_dates, selected_retinas, selected_cells, \
                                                selected_ages, selected_prepulse_pot, \
                                                selected_axonal_currents, selected_axonal_currents_corr,\
                                                selected_threshold_potentials,\
                                                selected_charge, selected_dur50): 

    if date_prev == date and retina_prev == retina and cell_prev == cell:
        v0_cell.append(v0)
        vth_cell.append(vth)
        ia_cell.append(ia)
        ia_corr_cell.append(ia_corr)
        ch_cell.append(charge)
        dur50_cell.append(dur50)
    else: 
        N += 1
        dates_per_cell.append(date)
        retina_per_cell.append(retina)
        cell_per_cell.append(cell)
        v0_per_cell.append(array(v0_cell))
        vth_per_cell.append(array(vth_cell))
        ia_per_cell.append(array(ia_cell))
        ia_corr_per_cell.append(array(ia_corr_cell))
        charge_per_cell.append(ch_cell)
        dur50_per_cell.append(dur50_cell)
        
        v0_cell = [v0]
        vth_cell = [vth]
        ia_cell = [ia]
        ia_corr_cell = [ia_corr]
        ch_cell = [charge]
        dur50_cell = [dur50]
    
    date_prev = date
    retina_prev = retina
    cell_prev = cell

### Adding the last cell
N += 1
v0_per_cell.append(array(v0_cell))
vth_per_cell.append(array(vth_cell))
ia_per_cell.append(array(ia_cell))
ia_corr_per_cell.append(array(ia_corr_cell))
charge_per_cell.append(ch_cell)
dur50_per_cell.append(dur50_cell)
    
thresholds_all = []
axonal_currents_all = []
axonal_currents_norm_all = []
prepulse_pot_all = []

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
            ia_corr_per_cell[i] = delete(ia_corr_per_cell[i], idx_v0_delete)
            vth_per_cell[i] = delete(vth_per_cell[i], idx_v0_delete)
            charge_per_cell[i] = delete(charge_per_cell[i], idx_v0_delete)
            dur50_per_cell[i] = delete(dur50_per_cell[i], idx_v0_delete)
        else:
            v0_per_cell[i] = array(v0_per_cell[i])
            ia_per_cell[i] = array(ia_per_cell[i])
            ia_corr_per_cell[i] = array(ia_corr_per_cell[i])
            vth_per_cell[i] = array(vth_per_cell[i])
            charge_per_cell[i] = array(charge_per_cell[i])
            dur50_per_cell[i] = array(dur50_per_cell[i])
        
    thresholds_all += list(vth_per_cell[i])
    axonal_currents_all += list(ia_per_cell[i])
    axonal_currents_norm_all += list(ia_per_cell[i]/min(ia_per_cell[i]))
    prepulse_pot_all += list(v0_per_cell[i])
        
### Measuring current and charge attenuation 
current_attenuation = []
current_attenuation_Vi = []
charge_attenuation = []
current_max = []
Vi_cells = []
Vi_thres_cells = []

### the function that is fit to the relationship between the peak axonal current and V0
vv = linspace(-80, -30,100)
f_abs = lambda v0, ki, vi, c, b : c + b/sqrt(1+exp((v0 - vi)/ki))

for i in range(N):
    if len(ia_per_cell[i]) > 5: # cells with more than 4 data points
        print (dates_per_cell[i], retina_per_cell[i], cell_per_cell[i])
        ### current and charge attenuation
        idx_60 = where(v0_per_cell[i] == -60.)[0]
        idx_40 = where(v0_per_cell[i] == -40.)[0]
        if len(idx_60) > 0 and len(idx_40) > 0:
            current_attenuation.append(ia_corr_per_cell[i][idx_60[0]]/ia_corr_per_cell[i][idx_40[0]])
            charge_attenuation.append(charge_per_cell[i][idx_60[0]]/charge_per_cell[i][idx_40[0]])
        else:
            current_attenuation.append(nan)
            charge_attenuation.append(nan)
            
        ### Vi from threshold adaptation
        row_fit = df_fit[(df_fit['Date'] == dates_per_cell[i]) & (df_fit['Retina'] == retina_per_cell[i])\
                      & (df_fit['Cell'] == cell_per_cell[i]) ]
            
        ### Fitting the current adaptation curve
        try:
            popt,_ = curve_fit(f_abs, v0_per_cell[i], abs(ia_corr_per_cell[i]), p0=[5., -55., 2, 1])
            k, V, cst, cst2 = popt
            v0_half = V + k*log((cst2/(max(abs(ia_corr_per_cell[i]))/sqrt(2) - cst))**2 - 1)
            Vi_cells.append(v0_half)
            
            # Current attenuation at Vi (from fitted function)
            i_max = max(abs(ia_corr_per_cell[i]))
            current_max.append(i_max)
            
            if len(row_fit) > 0:
                k_cell = row_fit['k'].values[0]
                vi_cell = row_fit['Vi'].values[0]
                print ('Vi from thres adapt:', vi_cell)
                
                Vi_thres_cells.append(vi_cell)
                i_at_vi = f_abs(vi_cell, k, V, cst, cst2)
                current_attenuation_Vi.append(i_at_vi/i_max)
            else:
                Vi_thres_cells.append(nan)
                current_attenuation_Vi.append(nan)
            
        except:
            Vi_cells.append(nan)
            current_attenuation_Vi.append(nan)
            current_max.append(nan)
            pass        
    else:
        current_attenuation.append(nan)
        charge_attenuation.append(nan)
        current_max.append(nan)
        Vi_cells.append(nan)
        current_attenuation_Vi.append(nan)
        
### Figure

name1 = "tab20b"
name2 = "tab20c"
name3 = "tab20"
cmap1 = get_cmap(name1)  
cmap2 = get_cmap(name2)  
cmap3 = get_cmap(name3) 
cols = cmap1.colors + cmap2.colors + cmap3.colors

fig = figure('Current adaptation', figsize=(9, 7.5))
gs = gridspec.GridSpec(3, 5, width_ratios=[3, 1, 1, 1, 1]) 
#threshold
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[0, 2])
ax4 = fig.add_subplot(gs[0, 3])
ax4bis = fig.add_subplot(gs[0, 4])
#current
ax5 = fig.add_subplot(gs[1, 0])
ax6 = fig.add_subplot(gs[1, 1])
ax7 = fig.add_subplot(gs[1, 3:5])
ax8 = fig.add_subplot(gs[1, 2])
#covariation
ax9 = fig.add_subplot(gs[2, 0])
ax10 = fig.add_subplot(gs[2, 1:3])
ax11 = fig.add_subplot(gs[2, 3:5])

### Threshold adaptation

### loading the results from the fits to the threshold adaptation data
data = load('RGC_threshold_adaptation.npz')
thres_v0 = data['arr_0']
thres_vth = data['arr_1']
thres_vi = data['arr_2']
thres_vmin = data['arr_3']
thres_ka = data['arr_4']
thres_ki = data['arr_5']
outliers = data['arr_6']
dates_vi = data['arr_7']
retinas_vi = data['arr_8']
cells_vi = data['arr_9']

print (' STATS THRESHOLD ADAPTATION')
print('Vi:', nanmean(thres_vi), '+-', nanstd(thres_vi))
print('Vi-Vmin:', nanmean(thres_vi-thres_vmin), '+-', nanstd(thres_vi-thres_vmin))
print('ka:', nanmean(thres_ka), '+-', nanstd(thres_ka))
print('ki/ka:', nanmean(thres_ki/thres_ka), '+-', nanstd(thres_ki/thres_ka))

### Panel A: example of threshold adaptation for one cell
m=-1
f = lambda v0, c, K, Vh : c + K * log(1.+exp((v0-Vh)/K))
x = linspace(-75., -30.,50)

popt,_ = curve_fit(f, thres_v0, thres_vth, p0=[-55., 5., -55.],\
                   bounds = ([-70, 1, -70], [-40, 9, -40]))
cst, ka, vh = popt
ax1.plot(x, f(x, cst, ka, vh), color=cols[m])
ax1.plot(thres_v0, thres_vth, 'o', color = cols[m])
ax1.plot(linspace(-80, -30, 100), linspace(-80, -30, 100), 'k--')
ax1.plot(thres_vi[m] * ones((10)), linspace(-76, f(thres_vi[m], cst, ka, vh), 10), '--', color='k', alpha=0.3)
ax1.set_ylim(-76, -30)
ax1.set_xlim(-76, -30)
ax1.set_yticks(ticks = [-70, -60, thres_vmin[m], -50, -40, -30])
ax1.set_yticklabels(['-70',  '-60', '$V_{min}$','-50', '-40', '-30']) 
ax1.set_xticks(ticks = [ -70, -60, thres_vi[m], -50, -40, -30])
ax1.set_xticklabels(['-70', '', '$V_i$', '-50', '-40', '-30']) 
ax1.set_ylabel('$V_{t}$ (mV)')
ax1.set_xlabel('$V_0$ (mV)')
ax1.legend(frameon=False, fontsize=8)
ax1.annotate("A", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')


### Panel B: stats of Vi
sns.boxplot(y=thres_vi, color='gray', ax=ax2)
sns.swarmplot(y=delete(thres_vi, outliers),  color=cols[m], ax=ax2)
sns.swarmplot(y=thres_vi[outliers], color=cols[m], ax=ax2)
sns.despine(bottom=True, ax=ax2)
ax2.set_ylabel('$V_i$ (mV)')
ax2.set_xticks([])
ax2.set_ylim(-65, -45)
ax2.annotate("B", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

### Panel C: stats of Vmin
sns.boxplot(y=thres_vi - thres_vmin, color='gray', ax=ax3)
sns.swarmplot(y=delete((thres_vi - thres_vmin), outliers),  color=cols[m], ax=ax3)
sns.swarmplot(y=(thres_vi - thres_vmin)[outliers],  color=cols[m],  ax=ax3)
sns.despine(bottom=True, ax=ax3)
ax3.set_ylabel('$V_i - V_{min}$ (mV)')
ax3.set_xticks([])
ax3.set_ylim(-10, 10)
ax3.annotate("C", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

### Panel D: stats of ki
sns.boxplot(y=thres_ka, color='gray', ax=ax4)
sns.swarmplot(y=delete(thres_ka, outliers),  color=cols[m], ax=ax4)
sns.swarmplot(y=thres_ka[outliers],  color=cols[m], ax=ax4)
sns.despine(bottom=True, ax=ax4)
ax4.set_ylabel('$k_a$ (mV)')
ax4.set_xticks([])
ax4.set_ylim(-1, 9)
ax4.annotate("D", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

### Panel E: stats of ki/ka
sns.boxplot(y=thres_ki/thres_ka, color='gray', ax=ax4bis)
sns.swarmplot(y=delete(thres_ki/thres_ka, outliers),  color=cols[m], ax=ax4bis)
sns.swarmplot(y=(thres_ki/thres_ka)[outliers],  color=cols[m],  ax=ax4bis)
sns.despine(bottom=True, ax=ax4bis)
ax4bis.set_ylabel('$k_i/k_a$')
ax4bis.set_xticks([])
ax4bis.set_ylim(0, 2)
ax4bis.annotate("E", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')


### Current adaptation 

print ('STATS CURRENT ADAPTATION')
print('I60/I40:', nanmean(current_attenuation), '+-', nanstd(current_attenuation))
print('lambda:', nanmean(current_attenuation_Vi), '+-', nanstd(current_attenuation_Vi))
print('Vi*:', nanmean(Vi_cells), '+-', nanstd(Vi_cells))

### Panel F: example of current adaptation for one cell
m = 2 
popt,_ = curve_fit(f_abs, v0_per_cell[m], abs(ia_corr_per_cell[m]), p0=[5., -60., 2, 1])
k, V, cst, cst2 = popt
v0_half = V + k*log((cst2/(max(abs(ia_corr_per_cell[m]))/sqrt(2) - cst))**2 - 1)
ax5.plot(vv, f_abs(vv, k, V, cst, cst2), color = cols[15])
ax5.plot(V * ones((10)), linspace(0, f_abs(V,k,V,cst,cst2), 10), '--', color='k', alpha=0.3)
ax5.plot(v0_per_cell[m], -ia_corr_per_cell[m], 'o', color = cols[15])
ax5.set_xlim(-76, -30)
ax5.set_ylim(0,10)
ax5.set_xticks(ticks = [ -70, -60, V, -50, -40, -30])
ax5.set_xticklabels(['-70', '', '$V_i^*$','-50', '-40', '-30']) 
ax5.set_ylabel('$-I_{p}$ (nA)')
ax5.set_xlabel('$V_0$ (mV)')
ax5.legend(frameon=False, fontsize=8)
ax5.annotate("F", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

### Panel G: stats of Imax/Imin
sns.boxplot(y=current_attenuation, color='gray', ax=ax6)
sns.swarmplot(y=current_attenuation,  color=cols[15], ax=ax6)
sns.despine(bottom=True, ax=ax6)
ax6.set_ylabel('$I_{60}/I_{40}$')
ax6.set_xticks([])
ax6.set_ylim(0, 25)
ax6.annotate("G", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

### Panel H: current attenuation
# c = linspace(0,8,9)
sns.boxplot(y=current_attenuation_Vi, color='gray', ax=ax8)
sns.swarmplot(y=current_attenuation_Vi, color=cols[15], ax=ax8)
sns.despine(bottom=True, ax=ax8)
ax8.set_ylabel('$\lambda$')
ax8.set_xticks([])
ax8.set_ylim(0, 1)
ax8.annotate("H", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

# Panel I: Vi* vs Vi
thres_vi_corr = []
Vi_cells_corr = []

for i in range(len(Vi_cells)):
    for j in range(len(thres_vi)):
        if (dates_per_cell[i], retina_per_cell[i], cell_per_cell[i]) == (dates_vi[j], retinas_vi[j], cells_vi[j]):
            Vi_cells_corr.append(Vi_cells[i])
            thres_vi_corr.append(thres_vi[j])

# slope_vi, intercept_vi, r_vi, p_vi, _ = linregress(array(thres_vi_corr)[~isnan(Vi_cells_corr)], array(Vi_cells_corr)[~isnan(Vi_cells_corr)])
# print ('Vi - Vi* correlation:', 'r=', r_vi, 'p=', p_vi) 

sns.scatterplot(x=thres_vi_corr, y=Vi_cells_corr, color='0.2', ax=ax7)
ax7.plot(linspace(-65,-46,10), linspace(-65,-46,10), 'k--')
ax7.set_ylim(-65,-45)
ax7.set_xlim(-65,-45)
ax7.set_ylabel('$V_i^*$ (mV)')
ax7.set_xlabel('$V_i$ (mV)')
ax7.annotate("I", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

### Covariation 

data_cov = load('RGC_adaptation_covariation.npz', allow_pickle=True)
cov_ia = data_cov['arr_0']
cov_vth = data_cov['arr_1']
cov_v0 = data_cov['arr_2']

slope_per_cell = []
intercept_per_cell = []
r_per_cell = []

N = len(cov_ia)

for i in range(N):
    if len(ia_per_cell[i]) > 5:
        slope, intercept, r_value, p_value, std_err = linregress(-log(-array(cov_ia[i])), cov_vth[i])
        slope_per_cell.append(slope/2)
        r_per_cell.append(r_value)
        intercept_per_cell.append(intercept)
        ax9.semilogx(-array(cov_ia[i]), cov_vth[i], 'o', color=cols[i])#, label='%i %s %i' %(dates_per_cell[i], retina_per_cell[i], cell_per_cell[i]))
    else:
        slope_per_cell.append(nan)
        r_per_cell.append(nan) 
        intercept_per_cell.append(nan)

### Panel J: voltage threshold vs peak axonal current        
ax9.set_ylabel('$V_t$ (mV)')
ax9.set_xlabel('$-I_p$ (nA)')
ax9.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax9.set_ylim(-65,-25)
ax9.set_xlim(0,12)
ax9.annotate("J", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

### Panel K: voltage threshold vs peak axonal current  in an example cell
m = 7
iii = linspace(min(-array(cov_ia[m])), max(-array(cov_ia[m])), 100)
ax10.semilogx(iii, intercept_per_cell[m] - 2*slope_per_cell[m]*log(iii), '-', color=cols[m],alpha=0.5)
ax10.semilogx(-array(cov_ia[m]), cov_vth[m], 'o', color=cols[m])
ax10.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax10.set_ylim(-65, -25)
ax10.set_xlabel('$-I_p$ (nA)')
ax10.legend(frameon=False)
ax10.annotate("K", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

### Panel L: estimation of ka from fit
sns.boxplot(y=slope_per_cell, color='gray', width=0.3, ax=ax11)
sns.swarmplot(y=slope_per_cell,  color='0.2', ax=ax11)
sns.despine(bottom=True, ax=ax11)
ax11.set_ylabel('$k_a$ (mV)')
ax11.set_xticks([])
ax11.set_ylim(0, 6)
ax11.annotate("L", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

fig.tight_layout()

print('STATS COVARIATION')
print('ka:', nanmean(slope_per_cell),'+-', nanstd(slope_per_cell))
print('ka median:', nanmedian(slope_per_cell))

### Saving the figure
save_path = '/Users/sarah/Dropbox/Spike initiation/PhD projects/Axonal current and AIS geometry/Paper/Figures/'
# fig.savefig(save_path + "fig8.pdf", bbox_inches='tight')

# fig.savefig(save_path + "fig8.pdf", bbox_inches='tight')

# fig.savefig(save_path + "fig8.png", dpi=300)
