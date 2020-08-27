

"""

Figure 3: Transmission of the axial current to the soma.

"""
from brian2 import *
import pandas as pd
import seaborn as sns
from pandas import ExcelWriter
from pandas import ExcelFile
from scipy import stats

### Figure parameters
rcParams['axes.spines.right'] = False
rcParams['axes.spines.top'] = False

fig = figure('Capacitive current', figsize=(6,5))
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)

### Laoding the results of analyses
path_files = '/Users/sarah/Documents/repositories/AIS-geometry-and-axial-current/'
df_cells = pd.read_excel(path_files + 'RGC_electrical_properties.xlsx')

dates = array(df_cells['Date'])[:-3]
retinas = array(df_cells['Retina'])[:-3]
cells = array(df_cells['Cell'])[:-3]
ages = array(df_cells['Age'])[:-3]

axonal_currents_corr = array(df_cells['Peak axonal current corrected'])[:-3]
capacitance_cc = array(df_cells['Cm CC'])[:-3]
dvdt_peak1 = array(df_cells['dvdt peak1'])[:-3]

charge = array(df_cells['Charge1 10'])[:-3]
duration50 = array(df_cells['Duration1 50'])[:-3]

# remove nans
capa_cc_nonan = capacitance_cc[~isnan(capacitance_cc)]
charge_nonan = charge[~isnan(capacitance_cc)]
dvdt_peak_nonan = dvdt_peak1[~isnan(capacitance_cc)]
axonal_current_nonan = axonal_currents_corr[~isnan(capacitance_cc)]
duration_nonan = duration50[~isnan(capacitance_cc)]

### Statistics
print ('Mean Ip:', mean(axonal_currents_corr), '+-', std(axonal_currents_corr))
print ('Mean dV/dt1:', nanmean(dvdt_peak1), '+-', nanstd(dvdt_peak1))
print ('Mean Q:', mean(charge), '+-', std(charge))
print ('Mean C:', mean(capacitance_cc[~isnan(capacitance_cc)]), '+-', std(capacitance_cc[~isnan(capacitance_cc)]))
print ('Mean t50:', mean(duration50), '+-', std(duration50))

### Panel A: axonal current vs capacitive current 
print ('PANEL A')

sns.scatterplot(x=capa_cc_nonan[~isnan(dvdt_peak_nonan)]*dvdt_peak_nonan[~isnan(dvdt_peak_nonan)] *1e-3, \
            y=-axonal_current_nonan[~isnan(dvdt_peak_nonan)], color='k', ax=ax1)
ax1.plot(linspace(0,20,50), linspace(0,20,50), 'k--')
ax1.set_xlabel('C dV/dt (nA)')
ax1.set_ylabel('$-I_p$ (nA)')
ax1.set_ylim(0, 20)
ax1.set_xlim(0, 20)
ax1.annotate("A", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

# Panel B: axonal current vs capacitance
print ('PANEL B')
print ('N cells:', len(capa_cc_nonan))

# Regression
slope_ic_o, intercept_ic_o, r_ic_o, p_ic_o, _ = stats.linregress(capa_cc_nonan, -axonal_current_nonan)
print ('Regression Ip - C:', 'r=%0.02f' %r_ic_o, 'p=%0.02f' %p_ic_o)
print ('Pearson Ip - C :', stats.pearsonr(capa_cc_nonan, -axonal_current_nonan))

# Linear regression
M = capa_cc_nonan[:, np.newaxis]
a, residuals, _, _ = lstsq(M, -axonal_current_nonan)
print ('Slope best linear fit:', a[0]*1e3, 'mV/ms')

ax2.plot(linspace(0,100,50), a[0] * linspace(0,100,50), 'k-')
sns.scatterplot(x=capa_cc_nonan, y=-axonal_current_nonan, color='k', ax=ax2)
ax2.set_xlabel('C (pF)')
ax2.set_ylabel('$-I_p$ (nA)')
ax2.set_ylim(0, 12.5)
ax2.set_xlim(0, 100)
ax2.annotate("B", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')


# Panel C: charge vs capacitance
print ('PANEL C')
print ('N cells:', len(capa_cc_nonan))

# Regression
slope_qc, intercept_qc, r_qc, p_qc, _ = stats.linregress(capa_cc_nonan, -charge_nonan)
print ('Regression Q - C:', 'r=', r_qc, 'p=', p_qc)
print ('Pearson Q - C :', stats.pearsonr(capa_cc_nonan, -charge_nonan))

# Linear regression
M = capa_cc_nonan[:, np.newaxis]
a, residuals, _, _ = lstsq(M, -charge_nonan)
print ('Slope best linear fit:', a[0]*1e3)

ax3.plot(linspace(0,100,50), a[0] * linspace(0,100,50), 'k-')
sns.scatterplot(x=capa_cc_nonan, y=-charge_nonan, color='k', ax=ax3)
ax3.set_xlabel('C (pF)')
ax3.set_ylabel('-Q (pC)')
ax3.set_ylim(0, 3.5)
ax3.set_xlim(0, 100)
ax3.annotate("C", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

# Panel D: current duration vs capacitance
print ('PANEL D')
print ('N cells:', len(capa_cc_nonan))

# Regression
slope_dc, intercept_dc, r_dc, p_dc, _ = stats.linregress(capa_cc_nonan, duration_nonan)
print ('Regression capacitance - duration50:', 'r=%0.02f'%r_dc, 'p=%0.02f'%p_dc)
print ('Pearson t50 - C :', stats.pearsonr(capa_cc_nonan, duration_nonan))

# Outlier
outliers_t50 = where(duration_nonan > 0.6)[0]
print ('Mean T50 without outlier:', mean(delete(duration_nonan, outliers_t50)), \
        '+-', std(delete(duration_nonan, outliers_t50)))


ax4.plot(linspace(20,90,50), intercept_dc + slope_dc * linspace(20,90,50), 'k-')
sns.scatterplot(x=capa_cc_nonan, y=duration_nonan, color='k', ax=ax4)
ax4.set_xlabel('C (pF)')
ax4.set_ylabel('$t_{50}$ (ms)')
ax4.set_ylim(0, 1)
ax4.set_xlim(0, 100)
ax4.annotate("D", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

tight_layout()

### Saving the figure
save_path = '/Users/sarah/Documents/repositories/AIS-geometry-and-axial-current/'
# fig.savefig("/Users/sarah/Dropbox/Spike initiation/Thesis/images/fig3.pdf", bbox_inches='tight')










