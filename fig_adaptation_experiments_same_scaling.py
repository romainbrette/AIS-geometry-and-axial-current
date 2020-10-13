#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 10:30:19 2020

@author: sarah
"""


"""

Figure: measurement of axonal current

"""
import glob2
import pandas as pd
import pyabf
import seaborn as sns
import matplotlib.patches as mpatches
from pandas import ExcelWriter
from pandas import ExcelFile
from brian2 import *
from scipy import stats
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from vc_test_pulse_analysis import *
from na_currents_analysis import *
from vc_test_pulse_analysis import *
import params_model_description

### Figure

rcParams['axes.spines.right'] = False
rcParams['axes.spines.top'] = False

fig1 = figure('Ip adapta', figsize=(6,3))

ax7 = fig1.add_subplot(121)
ax8 = fig1.add_subplot(122)

### Panel A and B: axonal currents recording in an exmaple cell
date = 20191114
retina = 'C'
cell = 1
cell_name = '%i %s %i' %(date, retina, cell)

### Loading the data
# path_to_cell = glob2.glob('data/RGC data/' + str(int(date)) + '*' + '/retina '+ str(retina) +'/cell ' + str(int(cell)))[0]
path_to_cell = glob2.glob('/Users/sarah/Documents/Data/Martijn Sierksma/' + str(int(date)) + '*' + '/retina '+ str(retina) +'/cell ' + str(int(cell)))[0]

### -60 mV
abf60 = pyabf.ABF(path_to_cell + '/VC threshold adaptation/2019_11_14_0044.abf')
fs60 = abf60.dataRate  * Hz # sampling rate
dt60 = 1./fs60

t = dt60*arange(len(abf60.sweepY)) 
I = []
V = []

for sweepNumber in abf60.sweepList:
    abf60.setSweep(sweepNumber)
    I.append(abf60.sweepY)
    V.append(abf60.sweepC*mV)

### Removing passive component
I_corr_pass, I_cut, t_cut = p5_subtraction(date, retina, cell, dt60, I, V, rec_name=str(int(33)).zfill(4))

### IV curves
I_peaks,  Vc_peaks, idx_peak_ax_current,  t_peaks = plot_IV(date, retina, cell, dt60, I_corr_pass, V, 0, str(int(6)).zfill(4))
Vc_peaks = array(Vc_peaks/mV)
I_peaks = array(I_peaks)

### Plotting
cmap = plt.get_cmap('binary')
colors = [cmap(i) for i in np.linspace(0, 1, len(abf60.sweepList))]

# Currents
ax7.plot(t_cut/ms, I_corr_pass[4] *1e-3, 'k', alpha=0.2) #color= colors[0]) 
ax7.text(0.5, 0.2,'%i'%(Vc_peaks[4] -70), fontsize=8)
ax7.plot(t_cut/ms, I_corr_pass[5] *1e-3,'k') #color= colors[5]) 
ax7.text(4.2, min(I_corr_pass[5])*1e-3,'%i'%(Vc_peaks[5] -70), fontsize=8 )
ax7.plot(t_cut/ms, I_corr_pass[7] *1e-3, 'k') #color= colors[7]) 
ax7.text(1.5, min(I_corr_pass[7])*1e-3,'%i'%(Vc_peaks[7] -70), fontsize=8)
#ax7.plot(t_cut/ms, I_corr_pass[9] *1e-3, 'k') #color= colors[9]) 
#ax7.text(0.4, min(I_corr_pass[9])*1e-3 -1.2,'%i'%(Vc_peaks[9] -70), fontsize=8)
#ax7.plot(t_cut/ms, I_corr_pass[13] *1e-3, 'k', alpha=0.2) #color= colors[13]) 
#ax7.text(0.4, (min(I_corr_pass[13])-700)*1e-3,'%i'%(Vc_peaks[13] -70), fontsize=8)
ax7.plot(t_cut/ms, I_corr_pass[17] *1e-3, 'k', alpha=0.2) #color= colors[17]) 
ax7.text(0.5, (min(I_corr_pass[17])-500)*1e-3,'%i'%(Vc_peaks[17] -70),  fontsize=8)
ax7.plot(t_cut/ms, I_corr_pass[-1] *1e-3, 'k', alpha=0.2) #color= colors[-1]) 
ax7.text(1, (min(I_corr_pass[-1])-500)*1e-3,'%i'%(Vc_peaks[-1] -70),  fontsize=8)
ax7.set_ylabel('I (nA)')
# ax7.set_xlabel('t (ms)')
ax7.set_xlim(0,5)
ax7.set_ylim(-15, 1)
ax7.set_xticks([])
sns.despine(bottom=True, ax=ax7)
ax7.plot(linspace(3.4, 4, 10), -14.5*ones(10), 'k-', linewidth=2)
ax7.text(3.3, -15.5, '0.5 ms', color='k', fontsize=8)
ax7.annotate("A", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

### -60 mV
abf40 = pyabf.ABF(path_to_cell + '/VC threshold adaptation/2019_11_14_0046.abf')
fs40 = abf40.dataRate  * Hz # sampling rate
dt40 = 1./fs40

t = dt40*arange(len(abf40.sweepY)) 
I = []
V = []

for sweepNumber in abf40.sweepList:
    abf40.setSweep(sweepNumber)
    I.append(abf40.sweepY)
    V.append(abf40.sweepC*mV)

### Removing passive component
I_corr_pass, I_cut, t_cut = p5_subtraction(date, retina, cell, dt40, I, V, rec_name=str(int(33)).zfill(4))

### IV curves
I_peaks,  Vc_peaks, idx_peak_ax_current,  t_peaks = plot_IV(date, retina, cell, dt40, I_corr_pass, V, 0, str(int(6)).zfill(4))
Vc_peaks = array(Vc_peaks/mV)
I_peaks = array(I_peaks)

### Plotting
cmap = plt.get_cmap('binary')
colors = [cmap(i) for i in np.linspace(0, 1, len(abf40.sweepList))]

# Currents
ax8.plot(t_cut/ms, I_corr_pass[0] *1e-3, 'k', alpha=0.2) #color= colors[0]) 
ax8.text(0.5, 0.2,'%i'%(Vc_peaks[0] -70), fontsize=8)
ax8.plot(t_cut/ms, I_corr_pass[3] *1e-3, 'k') #color= colors[0]) 
ax8.text(3, -0.9,'%i'%(Vc_peaks[3] -70), fontsize=8)
ax8.plot(t_cut/ms, I_corr_pass[9] *1e-3, 'k', alpha=0.2) #color= colors[9]) 
ax8.text(0.8, min(I_corr_pass[9]-500)*1e-3,'%i'%(Vc_peaks[9] -70), fontsize=8)
# ax8.plot(t_cut/ms, I_corr_pass[13] *1e-3, 'k', alpha=0.2) #color= colors[13]) 
# ax8.text(0.4, (min(I_corr_pass[13])-700)*1e-3,'%i'%(Vc_peaks[13] -70), fontsize=8)
ax8.plot(t_cut/ms, I_corr_pass[-1] *1e-3, 'k', alpha=0.2) #color= colors[-1]) 
ax8.text(1, (min(I_corr_pass[-1]-700))*1e-3,'%i'%(Vc_peaks[-1] -70),  fontsize=8)
ax8.set_ylabel('I (nA)')
# ax7.set_xlabel('t (ms)')
ax8.set_xlim(0,5)
ax8.set_ylim(-15, 1)
ax8.set_xticks([])
sns.despine(bottom=True, ax=ax8)
ax8.plot(linspace(3.5, 4, 10), -14.5*ones(10), 'k-', linewidth=2)
ax8.text(3.3, -15.5, '0.5 ms', color='k', fontsize=8)
ax8.annotate("B", xy=(0,1.1), xycoords="axes fraction",
                    xytext=(5,-5), textcoords="offset points",
                    ha="left", va="top",
                    fontsize=12, weight='bold')

tight_layout()

### Saving the figure
save_path = '/Users/sarah/Dropbox/Spike initiation/PhD projects/Axonal current and AIS geometry/Paper/Figures/'

# fig1.savefig(save_path + "fig_adaptation_experiments_same_scaling.pdf", bbox_inches='tight')

# fig1.savefig(save_path + "fig_adaptation_experiments_same_scaling.png", dpi=300)