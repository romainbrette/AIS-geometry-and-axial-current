
'''

A simple model to understand the effect of the series resistance on Na current recordings.

The pipette-electrode is represented as a dendrite

'''
from brian2 import *

defaultclock.dt = 0.01*ms

# General parameters
Ri = 100. * ohm * cm

# Measured CGC parameters
Cm_tot = 3.* pF 
Ce_tot = 0.*pF #5.*pF

# Neuron
EL = -90.*mV
soma_diam = 5.*um
soma_surface = pi * soma_diam**2
axon_len = 1000.*um
axon_diam = 0.5*um
Cm = Cm_tot/soma_surface
gL_neuron = 10. * (siemens / meter ** 2) # to obtain a realistic input resistance at the soma

# Electrode
elec_diam = 1.*um

def neuron_electrode_model(Rs, i_inj = 0):
    elec_len = Rs * pi * elec_diam**2 /(4 * Ri) # for a target series resistance
    print 'Electrode length:', elec_len
    
    # Morphology: 
    # 0: soma, 1: elec, from 2 to 1001: axon
    morpho = Soma(soma_diam)
    morpho.electrode = Cylinder(diameter=elec_diam, length=elec_len, n=1 ) #int(elec_len/um))
    morpho.axon = Cylinder(diameter=axon_diam, length=axon_len, n=1000)
    
    elec_surface = morpho.electrode.area
    print 'Electrode surface:', elec_surface, 'Soma surface:', soma_surface
    Ce = Ce_tot/elec_surface # F/m**2
    print 'Cm:', Cm/(uF/cm**2), 'Ce:', Ce/(uF/cm**2)
    Cm_total = hstack((Cm, Ce * ones(morpho.electrode.n), Cm * ones(morpho.axon.n))) * farad/meter**2
     
    # Equations
    eqs = '''
    Im = gL*(EL-v) : amp/meter**2
    I : amp (point current)
    gL : siemens / meter ** 2
    '''
    #Current-clamp and voltage-clamp stimulation
    gclamp = 100*usiemens
    
    eqs += '''
    I_CC : amp (point current)
    I_VC = gclamp*VC_on*(V_VC - v) : amp (point current)
    V_VC : volt
    VC_on : 1 # if 1, voltage-clamp is on
    '''
    
    neuron = SpatialNeuron(morphology=morpho, model=eqs, Cm=Cm_total, Ri=Ri, method="exponential_euler")
    
    neuron.v = EL
    neuron.gL = gL_neuron
    
    # Current injection in the proximal axon to mimic axonal currents at SI
    injection_pt = 10 

    # Record current and membrane potential at soma, electrode and injection point
    M_soma = StateMonitor(neuron, ('v', 'Im'), record = [0])
    M_elec = StateMonitor(neuron, ('v', 'Im'), record = morpho.electrode[-1])
    M_inj = StateMonitor(neuron, ('v', 'Im'), record = morpho.axon[injection_pt])

    Ra = (4*Ri)/(pi*axon_diam**2) * injection_pt*(axon_len/morpho.axon.n)
    
    # Run VC simulation with current injected in the axon
    # voltage-clamp at electrode
    neuron.electrode.V_VC[-1] = EL
    neuron.electrode.VC_on[-1] = 1
    run(20*ms, report='text')
    # voltage-clamp at electrode
    neuron.electrode.V_VC[-1] = EL + 10.*mV
    neuron.electrode.VC_on[-1] = 1
    run(1*ms)
    # voltage-clamp and current injection in the axon
    neuron.electrode.V_VC[-1] = EL + 10.*mV
    neuron.electrode.VC_on[-1] = 1
    neuron.axon.I[injection_pt] = i_inj
    run(1*ms, report='text')
    # voltage-clamp at electrode
    neuron.electrode.V_VC[-1] = EL + 10.*mV
    neuron.electrode.VC_on[-1] = 1
    neuron.I = 0. * amp
    run(18*ms)
    # voltage-clamp at electrode
    neuron.electrode.V_VC[-1] = EL
    neuron.electrode.VC_on[-1] = 1
    run(20*ms)
    
    return M_elec.t/ms, M_elec.v[0]/mV, M_soma.v[0]/mV, M_inj.v[0]/mV, M_soma.Im[0], M_elec.Im[0], M_inj.Im[0], Ra, elec_surface, soma_surface








