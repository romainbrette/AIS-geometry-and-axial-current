

"""

Peak axonal current in a simple RGC model.

Extended AIS.

SI model.

Dicho protocol.

"""

from brian2 import * 
import glob
import pandas as pd
import params_model_description
from vc_test_pulse_analysis import *
from na_currents_analysis import *

### Model parameters
params = params_model_description
dt = 0.01*ms
v_h = params.EL # leak potential
length = 30.*um # AIS length
starts = linspace(0, 20, 5)*um # AIS start position
gna_density = 4300.*siemens/meter**2 # Na conductance density
gna_tot = gna_density * length * params.axon_diam * pi # total Na conductance
ra = (4*params.Ri/(pi*params.axon_diam**2)) # axial resistance per unit length

### Load data from simulations
path_to_data = 'simulations data/fig5/'

### Measure peak currents
peak_currents = zeros( len(starts))

for i in range(len(starts)):
    start = starts[i]    
    print (start)
    
    dir_name = path_to_data + 'VC dicho SImodel ext AIS x%0.1f L%i g5000' %(start/um, length/um)

    # Load and plot axonal currents
    Ie = loadtxt(dir_name +'/Steps/I.txt')
    Vm = loadtxt(dir_name +'/Steps/V.txt')
    Vc = loadtxt(dir_name + '/Steps/Vc.txt')
    Im = loadtxt(dir_name+'/Steps/Im.txt')
    
    n_rec = len(Ie)
    t = arange(len(Ie[0]))*dt
    
    cmap = plt.get_cmap('gnuplot')
    cols = [cmap(i) for i in np.linspace(0, 1, n_rec)]
        
    # Load test pulse data
    dir_name_tp = path_to_data + 'TP SImodel x%0.1f L%0.01f' %(start/um, length/um)

    Ie_tp = loadtxt(dir_name_tp +'/Steps/I.txt')
    Ve_tp = loadtxt(dir_name_tp +'/Steps/V.txt')
    Vc_tp = loadtxt(dir_name_tp + '/Steps/Vc.txt')
    
    t_tp = arange(len(Ve_tp))*dt
    vc_tp_amp = -5
    
    vc = Vc #[Vc[i][int(30.*ms/dt)] for i in range(n_rec)] # command potentials
    
    steps_start = int(20.*ms/dt)
    steps_end = int(170.*ms/dt)
    
    i_tp_cut = Ie_tp[steps_start:steps_end] - mean(Ie_tp[int(18.*ms/dt):int(19.5*ms/dt)])
    
    Ie_corr = []
    
    # Correction for the passive component of the currents
    fig_corr = figure('Correction x%0.1f L%i' %(start/um, length/um), (10,12))
    
    for k in range(n_rec):
        V_amp =  vc[k] - v_h/mV
        factor = V_amp/vc_tp_amp
        i_cut = Ie[k][steps_start:steps_end] - mean(Ie[k][int(18.*ms/dt):int(19.5*ms/dt)])
        t_cut = t[steps_start:steps_end]
            
        i_corr =  i_cut - i_tp_cut * factor
        Ie_corr.append(i_corr)
        
        plot(t_cut/ms, Ie[k][steps_start:steps_end], color=cols[k])
        plot(t_cut/ms, i_corr, '--', color=cols[k])
        plot(t_cut/ms, i_tp_cut, 'g-')
        ylim(-15,15)
        
    # Measure peak axonal current
    spikes = zeros(n_rec)
    i_peaks = []
    peaks_indexes = []
    for j in range(n_rec):
        idx_peak = argmin(Ie_corr[j][int(0.1*ms/dt):int(149.*ms/dt)]) + int(0.1*ms/dt)
        i_peak = Ie_corr[j][idx_peak]
        i_peaks.append(i_peak)
        peaks_indexes.append(idx_peak)
        if i_peak < -3.:
            spikes[j] = 1
    
    i_peaks = array(i_peaks)
    
    # Peak current
    idx_spike = where(spikes == 1)[0]
    idx_peak = argmax(i_peaks[idx_spike])
    peak_currents[i] = i_peaks[idx_spike][idx_peak]
    
show()

### Save the measured peak currents in an npz file
# savez('model_SI_peak_current_ext_AIS_g5000_test.npz', starts/um, length/um, gna_tot, peak_currents)

### Plotting results
f3 = figure('It pred', figsize=(5,4))

subplot(111)
plot(starts/um, peak_currents, 'k-')
ylabel('$I_p$ (nA)')
xlabel('$\Delta$ ($\mu$m)')
legend(frameon=False)

tight_layout()

show()