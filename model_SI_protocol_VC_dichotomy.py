
"""
Dichotomy method for precise threshold measurement.

SI model.

"""

import os
from brian2 import *
import params_model_description, params_simple_model
from model_Na_Kv1 import *
from model_Na_Kv1_with_Rs import *
from model_spike_initiation import model_spike_initiation, model_spike_initiation_with_Rs

path_save_data = "simulations data/fig5/"
# path_save_data = "simulations data/fig11/"

params = params_model_description

defaultclock.dt = 0.01*ms
dt = defaultclock.dt

V0 = params.EL

length = 30.*um
starts = linspace(0, 20, 5)*um 
gna_density = 5000.*siemens/meter**2
gna_tot = gna_density * length * params.axon_diam * pi

for start in starts:
    do_experiment = not os.path.exists('Steps')
    
    if do_experiment:
        neuron = model_spike_initiation(params=params, Na_start=start, Na_end=start+length, gna_tot=gna_tot)
        path = path_save_data + 'VC dicho SImodel ext AIS x%0.01f L30 g5000' % (start/um)
        
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
        
        store()
        
        figure('Dichotomy x=%i' %(start/um))

        ampli_min = V0
        ampli_current = V0 + 30.*mV
        ampli_max = V0 + 60.*mV
        spike = False
        
        n_it = 0
                                      
        while True:
            print (n_it, ampli_current/mV)
            
            restore()
            
            # VC protocol
            neuron.V_VC[0] = V0
            neuron.VC_on[0] = 1
            run(20*ms)
            neuron.V_VC[0] = ampli_current
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
            Vcom.append(ampli_current)
            
            # Measuring the peak axonal current
            i_max = max(abs(M.I_VC[0][int(20.25 * ms / dt):int(169 * ms / dt)]))
            print ('i=', i_max)
            i_threshold = 1.5*nA
            
            if n_it > 51:
                print ('too much iterations')
                break
            if i_max >= i_threshold and abs(ampli_current - ampli_min) <= 0.05*mV and spike is False:
                print (' stop ')
                break
            if i_max <= i_threshold:
                ampli_min = ampli_current
                spike = False
            else:
                ampli_max = ampli_current
                spike = True
                    
            ampli_current = 0.5*ampli_max + 0.5*ampli_min
            
            n_it += 1
        
        # Save data
        savetxt(path+'/Steps/I.txt',array(I)/nA)
        savetxt(path+'/Steps/Im.txt',array(I)/nA)
        savetxt(path+'/Steps/V.txt',array(V)/mV)
        savetxt(path+'/Steps/Vc.txt',array(Vcom)/mV)
            
        
        
    
        





