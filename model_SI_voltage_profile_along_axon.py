
"""

Voltage profile along axon.

Spike initiation model with extended AIS.

"""

from brian2 import *
import params_model_description, params_simple_model
from model_spike_initiation import *

defaultclock.dt = 0.01*ms
dt = defaultclock.dt

### Model parameters
params = params_model_description
V0 = params.EL # leak reversal potential
length = 30.*um # AIS length
start = 5.*um # AIS start position
end = start+length # AIS end position
gna_tot = 5000*(siemens/metre**2) * params.axon_diam * length * pi # total Na conductance

### Neuron model
neuron = model_spike_initiation(params=params, Na_start=start, \
                      Na_end=start+length, gna_tot=gna_tot)

I = []
Im = []
Im_ais = []
INa = []
V = []
V_ais = []
Vcom = []

### Running simulation: voltage clamp below and above threshold
M = StateMonitor(neuron, ('v','I_VC', 'Im'), record = 0)
M_AIS = StateMonitor(neuron, ('v', 'Im', 'INa'), record = neuron.morphology.axon[:100.*um])

store()

fig = figure('V profile', (6,3))

amplis = linspace(V0+5*mV, V0+7.445*mV, 2)

v_profile_at_max = []
x = arange(0, (start+end)/um, 1)                                
 
for ampli in amplis:
    print (ampli/mV)
    
    restore()
    
    # VC protocol
    neuron.V_VC[0] = V0
    neuron.VC_on[0] = 1
    run(20*ms)
    neuron.V_VC[0] = ampli
    neuron.VC_on[0] = 1
    run(20*ms)
    neuron.V_VC[0] = V0
    neuron.VC_on[0] = 1
    run(20*ms)
     
    v_profile_at_max.append(M_AIS.v[:, argmax(M_AIS.v[int(end/um)])]/mV)
    
    subplot(111)
    plot(M_AIS.v[:, argmax(M_AIS.v[int(end/um)])]/mV)  
    ylabel('V (nA)')
    xlabel('x')

tight_layout()
    
show()

### Saving data
# savez('/Users/sarah/Documents/repositories/AIS-geometry-and-axonal-current/Na currents in RGC/codes submission/data/SI_model_voltage_profile',\
#       v_profile_at_max, x)

       

        
    
        





