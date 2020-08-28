
'''

Figure 5: predictions of axial current with resistive coupling theory.

'''
from brian2 import *

from pandas import ExcelWriter
from pandas import ExcelFile
from scipy.stats import linregress
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib import gridspec
import pandas as pd
import matplotlib.patches as mpatches
import seaborn as sns
import params_model_description
import matplotlib.image as mpimg
import matplotlib.patches as patches

### Figure parameters
rcParams['axes.spines.right'] = False
rcParams['axes.spines.top'] = False

### Model parameters
params = params_model_description

### Loading the results of analyses
df_cells = pd.read_excel('RGC_electrical_properties.xlsx')

dates = array(df_cells['Date'])[:-3]
retinas = array(df_cells['Retina'])[:-3]
cells = array(df_cells['Cell'])[:-3]
ages = array(df_cells['Age'])[:-3]

axonal_currents = array(df_cells['Peak axonal current corrected'])[:-3] 
ais_starts = array(df_cells['AIS start (um)'])[:-3]
ais_lengths = array(df_cells['AIS length (um)'])[:-3]

# remove nans
dates = dates[~isnan(ais_lengths)]
Ip = axonal_currents[~isnan(ais_lengths)] * nA
Delta = ais_starts[~isnan(ais_lengths)] * um
L = ais_lengths[~isnan(ais_lengths)] * um

Delta_range = linspace(min(Delta),max(Delta),100)

### Figure

fig = figure('Ip and geometry, theory', figsize=(6, 5))

ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)

# Panel A: spatial profile of Vm
# loading simulations data
data_v = load('model_SI_voltage_profile.npz')
v_60 = data_v['arr_0'][1]
x_axon = linspace(0,100, 100) #data_v['arr_1']

# simulations parameters
ENa = params.ENa
Vs = -55*mV
DeltaV = ENa - Vs
d = 1.*um
dmin = .8*um
dmax = 1.2*um
k = params.Ka

Ri = params.Ri

ax1.fill_between(linspace(5, 35,50), -80, 80*ones(50), color='k', alpha = 0.1)
ax1.plot(x_axon, v_60, 'g-', alpha=0.5)

v0 = min(v_60)
ax1.plot(linspace(0,5,50), ((ENa/mV -v0)/5)*linspace(0,5,50)+v0, 'g-')
ax1.plot(linspace(5,100,50), ENa/mV*ones(50), 'g-')
# slope
slope = (v_60[5]-v_60[0])/5
delta = sqrt(d/(4*Ri*5000.*siemens/meter**2))
ax1.plot(linspace(0,5+delta/um,50), slope*linspace(0,5+delta/um,50)+v0, 'g--')
# delta
ax1.annotate('', xy=(4, 76), xytext=(6+delta/um, 76),
            arrowprops=dict(arrowstyle='<->', color='k', ls='-'))
ax1.text(5 + (delta/um)/2 - 1, 81, "$\delta$")
ax1.annotate('', xy=(0, 76), xytext=(6, 76),
            arrowprops=dict(arrowstyle='<->', color='k', ls='-'))
ax1.text(1.5, 81, "$\Delta$")
ax1.set_ylabel('V (mV)')
ax1.set_xlabel('distance along axon ($\mu$m)')
ax1.set_xlim(0,50)
ax1.set_ylim(-80,80)
ax1.annotate("A", xy=(-0.1,1.2), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

# Panel B: peaka xonal current vs AIS start position in the model
# loading simulations data
data = load('model_SI_peak_current_ext_AIS_g5000.npz')
Delta_sim = data['arr_0']*um # AIS start positions
L_sim = data['arr_1']*um # AIS length
Ip_sim = data['arr_3']*nA # peak axonal currents

# theoretical prediction for the peak axonal current
def I(g, d, Delta=Delta):
    ra = Ri * 4 / (pi * d ** 2)
    x = (4*Ri*g/d)**.5
    #L=1e6*um
    #u = tanh(x*L)
    u=1
    return -DeltaV/ra*1/(Delta+1/(x*u))

g_sim = 5000.*siemens/meter**2
Delta_range_sim = linspace(0, 20, 100)*um
ax2.plot(Delta_sim/um, Ip_sim/nA, '.-', color='g', label='simulation')
ax2.plot(Delta_range_sim/um, I(g_sim,d,Delta_range_sim)/nA, 'g--', label='theory')
ax2.set_ylabel('$I_p$ (nA)')
ax2.set_xlabel('$\Delta$ ($\mu$m)')
ax2.set_ylim(-15,0)
ax2.set_xlim(0,20)
ax2.legend(frameon=False, fontsize=8)
ax2.annotate("B", xy=(-0.1,1.2), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

# Panel C: error on I
g0=linspace(300,10000,1000)*siemens/meter**2
error = array([mean((Ip-I(u,d))**2) for u in g0])

ax3.plot(g0, error/(nA**2) , 'k-')
ax3.set_ylim(0, 30)
ax3.set_xlim(0,10000)
ax3.set_xlabel('g (S/m$^2$)')
ax3.set_ylabel('Error (nA$^2$)')
ax3.annotate("C", xy=(-0.1,1.2), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

g = g0[argmin(error)] #5000*siemens/meter**2
print("g = {} S/m^2".format(g/(siemens/meter**2)))
print("Erreur = {} nA".format(min(error)**.5/nA))


# Panel D: Ip vs Delta in data

Id = I(g,d, Delta_range)
Imin = I(g,dmin,Delta_range)
Imax = I(g,dmax,Delta_range)

ax4.plot(Delta_range/um, Id/nA,'g', alpha=0.5)
ax4.plot(Delta_range/um, Imin/nA,'g', alpha=0.25)
ax4.plot(Delta_range/um, Imax/nA,'g', alpha=1)
ax4.plot(Delta/um, Ip/nA,'k.')
ax4.set_ylabel('$I_p$ (nA)')
ax4.set_xlabel('$\Delta$ ($\mu$m)')
ax4.set_ylim(-15,0)
ax4.set_xlim(0,20)
ax4.legend(frameon=False)
ax4.annotate("D", xy=(-0.1,1.2), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

tight_layout()

### Saving the figure

save_path = '/Users/sarah/Dropbox/Spike initiation/PhD projects/Axonal current and AIS geometry/Paper/Figures/'
# fig.savefig(save_path + "fig5.pdf", bbox_inches='tight')

# fig.savefig(save_path + "fig_Peak_axonal_current_theory.pdf", bbox_inches='tight')




