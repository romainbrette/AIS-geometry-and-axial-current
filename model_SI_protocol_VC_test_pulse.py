


"""
Test pulse.

"""

import os
from brian2 import *
import params_model_description, params_simple_model
from model_spike_initiation import model_spike_initiation

path_save_data = "/users/sarahgoethals/Documents/Thesis/Models simulations/peak axonal current/"

params = params_model_description

defaultclock.dt = 0.01*ms
dt = defaultclock.dt

V0 = params.EL

length = 30.*um
starts = linspace(0, 20, 5)*um 
gna_density = 4300.*siemens/meter**2
gna_tot = gna_density * length * params.axon_diam * pi
    
for start in starts:
    do_experiment = not os.path.exists('Steps')
    
    if do_experiment:
        neuron = model_spike_initiation(params=params, Na_start=start, Na_end=start+length, gna_tot=gna_tot)
        path = path_save_data + 'TP SImodel x%0.01f L%0.01f' % (start/um, length/um)
            
        # Make a data folder
        if not os.path.exists('data'):
            os.mkdir('data')
        os.mkdir(path)
        
        os.mkdir(path+'/Steps')
        I = []
        Im = []
        V = []
        Vcom = []
        
        M = StateMonitor(neuron, ('v','I_VC', 'Im'), record = 0)
        
        figure('TP  x%0.01f' % (start/um), figsize=(10,6))
        
        # VC protocol
        neuron.V_VC[0] = V0
        neuron.VC_on[0] = 1
        run(20*ms)
        neuron.V_VC[0] = V0-5*mV
        neuron.VC_on[0] = 1
        run(150*ms)
        neuron.V_VC[0] = V0
        neuron.VC_on[0] = 1
        run(20*ms)
        
        
        subplot(211)
        plot(M.t/ms, M.v[0]/mV) 
        ylabel('Voltage (mV)')
        xlabel('Time (ms)')
        
        subplot(212)
        plot(M.t/ms, M.I_VC[0]) 
        xlabel('Time (ms)')
        ylabel('Electrode current (nA)')
        
        tight_layout()
        
        I.append(M.I_VC[0])
        Im.append(M.Im[0])
        V.append(M.v[0])
            
        Vcom.append(V0-5*mV)
        
        # Save data
        savetxt(path+'/Steps/I.txt',array(I)/nA)
        savetxt(path+'/Steps/Im.txt',array(Im)/nA)
        savetxt(path+'/Steps/V.txt',array(V)/mV)
        savetxt(path+'/Steps/Vc.txt',array(Vcom)/mV)





















