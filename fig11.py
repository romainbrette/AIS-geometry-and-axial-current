


"""

Figure 11: threoretical estimate of axial current at threshold compared to simulations.

"""

from brian2 import * 
import glob
import pandas as pd
from vc_test_pulse_analysis import *
from na_currents_analysis import *
import params_model_description

rcParams['axes.spines.right'] = False
rcParams['axes.spines.top'] = False

### Model parameters
params = params_model_description

dt = 0.01*ms
v_h = params.EL

starts = linspace(0, 20, 5)*um 
lengths = linspace(10,30,3)*um

ra = (4*params.Ri/(pi*params.axon_diam**2))
Ras =  ra * starts
print('Axial resistance:', Ras/Mohm)

### Path to simulation results
path_to_data = "/users/sarah/Documents/Models simulations/current at threshold/"

threshold_currents = zeros((len(lengths), len(starts)))
peak_currents_latency = zeros((len(lengths), len(starts)))

IV_i = []
IV_vc = []

for j in range(len(lengths)):
    length = lengths[j]
    for i in range(len(starts)):

        start = starts[i]
        
        print (length, start)
        
        dir_name = path_to_data + 'VC dicho SImodel ext AIS x%0.1f L%i' %(start/um, length/um)
    
        # Load and plot data
        Ie = loadtxt(dir_name +'/Steps/I.txt')
        Vm = loadtxt(dir_name +'/Steps/V.txt')
        Vc = loadtxt(dir_name + '/Steps/Vc.txt')
        Im = loadtxt(dir_name+'/Steps/Im.txt')
        
        n_rec = len(Ie)
        t = arange(len(Ie[0]))*dt
        
        # Load test pulse data
        dir_name_tp = path_to_data + 'TP SImodel ext AIS x%0.1f L%i' %(start/um, length/um)
    
        Ie_tp = loadtxt(dir_name_tp +'/Steps/I.txt')
        Ve_tp = loadtxt(dir_name_tp +'/Steps/V.txt')
        Vc_tp = loadtxt(dir_name_tp + '/Steps/Vc.txt')
        
        t_tp = arange(len(Ve_tp))*dt
        vc_tp_amp = -5
            
        # Remove passive response from Na currents recordings
        vc = Vc #[Vc[i][int(30.*ms/dt)] for i in range(n_rec)] # command potentials
        
        steps_start = int(20.*ms/dt)
        steps_end = int(170.*ms/dt)
        
        i_tp_cut = Ie_tp[steps_start:steps_end] - mean(Ie_tp[int(18.*ms/dt):int(19.5*ms/dt)])
        
        Ie_corr = []
                
        for k in range(n_rec):
            V_amp =  vc[k] - v_h/mV
            factor = V_amp/vc_tp_amp
            #print (V_amp, factor)
            i_cut = Ie[k][steps_start:steps_end] - mean(Ie[k][int(18.*ms/dt):int(19.5*ms/dt)])
            t_cut = t[steps_start:steps_end]
                
            i_corr =  i_cut - i_tp_cut * factor
            Ie_corr.append(i_corr)
            
        # Measure peak axonal current and threshold
        spikes = zeros(n_rec)
        i_peaks = []
        for l in range(n_rec):
            # peak current
            idx_peak = argmin(Ie_corr[l][int(0.1*ms/dt):int(150*ms/dt)]) + int(0.1*ms/dt)
            i_peak = Ie_corr[l][idx_peak]
            #i_peak = mean(Ie_corr[j][int(120*ms/dt):int(150*ms/dt)])
            i_peaks.append(i_peak)
            if i_peak < -1.:
                spikes[l] = 1
                                
        vc = array(vc)
        i_peaks = array(i_peaks)
        I_corr = array(Ie_corr)
        spikes = array(spikes)
        
        idx_sort = argsort(vc)
        IV_i.append(i_peaks[idx_sort])
        IV_vc.append(vc[idx_sort])
        
        # Current threshold
        idx_no_spike = where(spikes == 0)[0]
        Ie_thres = I_corr[idx_no_spike][-1]
                
        threshold_currents[j,i] = i_peaks[idx_no_spike][-1] 
    
show()

### Predictions

f2 = figure('Thres current', figsize=(4, 3))

ax = subplot(111)
ax.loglog(linspace(5,35,50), (params.Ka/(ra*linspace(5,35,50)*um))/nA, 'g--', label='theory')

ax.loglog((starts[:4]+lengths[0]/2)/um, -threshold_currents[0,:][:4], 'go', label='simulation, L=%i $\mu$m' %(lengths[0]/um))
ax.loglog((starts[2:]+lengths[2]/2)/um, -threshold_currents[2,:][2:], 'go', label='simulation, L=%i $\mu$m' %(lengths[2]/um), alpha=0.5)
ax.set_xlabel('$\Delta+L/2$ ($\mu$m)')
ax.set_ylabel('$-I_t$ (nA)')
ax.set_ylim(0.1, 1)
ax.set_xlim(1,50)
ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax.xaxis.set_minor_formatter(FormatStrFormatter(''))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.set_minor_formatter(FormatStrFormatter(''))
ax.legend(frameon=False)

tight_layout()

show()

### Saving the figure
save_path = '/Users/sarah/Documents/repositories/AIS-geometry-and-axonal-current/Na currents in RGC/codes submission/data/'
# f2.savefig(save_path + "fig11.pdf", bbox_inches='tight')


    
        