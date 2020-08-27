

"""
TP without clampy in model.

"""

import os
from brian2 import *
import params_model_description, params_simple_model
from model_Na_Kv1 import *
from model_Na_Kv1_with_Rs import *

path_save_data = '/Users/sarah/Documents/repositories/AIS-geometry-and-axial-current/simulations data/fig2/'

params = params_model_description

defaultclock.dt = 0.01*ms
dt = defaultclock.dt

#series_resistance = False
#rs = 1.*Mohm
V0 = params.EL

start = 5.*um
length = 30.*um 

series_resistance= linspace(0, 5, 6)*Mohm

for rs in series_resistance:
    do_experiment = not os.path.exists('Steps')
    
    if do_experiment:
        if rs/Mohm == 0:
            print ('no Rs')
            neuron = model_Na_Kv1(params=params,resting_vm=V0, Na_start=start, \
                                  Na_end=start+length, density=False, gna_tot=700.*nS)
            path = path_save_data + 'Test pulse APmodel ext AIS x%0.01f L%0.01f r0.0' \
                % (start/um, length/um)
        else:
            print ('Rs:', rs)
            neuron = model_Na_Kv1_with_Rs(params=params,resting_vm=V0, Na_start=start, \
                                          Na_end=start+length, density=False, Rs=rs, gna_tot=700.*nS)
            path = path_save_data + 'Test pulse APmodel ext AIS x%0.01f L%0.01f r%0.01f' \
                % (start/um, length/um, rs/Mohm)
      
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
            
        figure('TP', figsize=(10,6))
        
        # VC protocol
        neuron.V_VC[0] = V0
        neuron.VC_on[0] = 1
        run(20*ms)
        neuron.V_VC[0] = V0-5*mV
        neuron.VC_on[0] = 1
        run(20*ms)
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
            
        # Save data
        savetxt(path+'/Steps/I.txt',array(I)/nA)
        savetxt(path+'/Steps/V.txt',array(V)/mV)
        savetxt(path+'/Steps/Vc.txt',array([(V0-5*mV)/mV]))
    
    



















