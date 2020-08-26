
"""
A simple model of spike initiation.

"""

from brian2 import *

__all__ = ['model_spike_initiation', 'model_spike_initiation_with_Rs']

def model_spike_initiation(params, Na_start = 5.*um, Na_end = 35.*um, gna_tot = 350*nS, morpho = None):
    '''
    If the morphology is not given, then it is created to a default.

    neuron.I_CC is a current-clamp stimulation.
    neuron.V_VC is a voltage-clamp stimulation.
    neuron.VC_on is a voltage-clamp switch (0/1).
    neuron.I_VC is the voltage-clamp current.
    '''

    ### Passive parameters
    EL = params.EL 
    Cm = params.Cm 
    gL = params.gL 
    Ri = params.Ri 

    if morpho is None:
        dend_diam = 6. * um
        dend_length = 1000. * um
        axon_diam = 1. * um
        axon_length = 500. * um
        soma_diameter = 30. * um
        morpho = Soma(diameter=soma_diameter)
        dendrite = Cylinder(diameter=dend_diam, length=dend_length, n=500)
        axon = Cylinder(diameter=axon_diam, length=axon_length, n=500)
        morpho.dendrite = dendrite
        morpho.axon = axon

    # Na channels distribution
    AIS_length = Na_end-Na_start

    # Na channels parameters
    ENa = params.ENa 
    Gna = gna_tot
    
    # Channels kinetics
    Va = params. Va 
    Ka = params.Ka 
    Taum_max = params.Taum_max 
    
    # Equations
    eqs='''
    Im = gL*(EL - v) + gNa*m*(ENa - v) : amp/meter**2
    INa = gNa*m*(ENa - v) : amp/meter**2
    dm/dt = (minf - m) / Taum_max: 1  # simplified Na channel
    minf = 1 / (1 + exp((Va - v) / Ka)) : 1
    gNa : siemens/meter**2
    
    '''
    # Current-clamp and voltage-clamp stimulation
    gclamp = 1000*usiemens

    eqs += '''
    I_CC : amp (point current)
    I_VC = gclamp*VC_on*(V_VC - v) : amp (point current)
    V_VC : volt
    VC_on : 1 # if 1, voltage-clamp is on
    
    I_hyp : amp (point current)
    
    '''
    neuron = SpatialNeuron(morphology=morpho, model=eqs, Cm=Cm, Ri=Ri,
                                    namespace=dict(gL=gL, EL=EL, ENa=ENa, 
                                                   Ka=Ka, Va=Va, Taum_max=Taum_max,                 
                                                   gclamp = gclamp),
                           method="exponential_euler")
    
    # Initial segment                           
    if Na_start == Na_end:
        initial_segment = morpho.axon[Na_start]
        neuron.gNa[initial_segment] = Gna/neuron.area[initial_segment]
    else:
        initial_segment = morpho.axon[Na_start:Na_end]
        neuron.gNa[initial_segment] =  Gna/(pi*AIS_length*axon_diam) 
    
    # Initialisation
    neuron.v = EL

    return neuron

def model_spike_initiation_with_Rs(params, Na_start = 5.*um, Na_end = 35.*um, gna_tot = 350*nS, Rs=0, morpho = None):
    '''
    If the morphology is not given, then it is created to a default.

    neuron.I_CC is a current-clamp stimulation.
    neuron.V_VC is a voltage-clamp stimulation.
    neuron.VC_on is a voltage-clamp switch (0/1).
    neuron.I_VC is the voltage-clamp current.
    '''

    ### Passive parameters
    EL = params.EL 
    Cm = params.Cm 
    gL = params.gL 
    Ri = params.Ri 

    if morpho is None:
        dend_diam = 6. * um
        dend_length = 1000. * um
        axon_diam = 1. * um
        axon_length = 500. * um
        soma_diameter = 30. * um
        morpho = Soma(diameter=soma_diameter)
        dendrite = Cylinder(diameter=dend_diam, length=dend_length, n=500)
        axon = Cylinder(diameter=axon_diam, length=axon_length, n=500)
        morpho.dendrite = dendrite
        morpho.axon = axon

    # Na channels distribution
    AIS_length = Na_end-Na_start

    # Na channels parameters
    ENa = params.ENa 
    Gna = gna_tot
    
    # Channels kinetics
    Va = params. Va 
    Ka = params.Ka 
    Taum_max = params.Taum_max 
    
    # Equations
    eqs='''
    Im = gL*(EL - v) + gNa*m*(ENa - v) : amp/meter**2
    INa = gNa*m*(ENa - v) : amp/meter**2
    dm/dt = (minf - m) / Taum_max: 1  # simplified Na channel
    minf = 1 / (1 + exp((Va - v) / Ka)) : 1
    gNa : siemens/meter**2
    
    '''
    # Current-clamp and voltage-clamp stimulation
    gclamp = 1000*usiemens

    eqs += '''
    I_CC : amp (point current)
    I_VC = (V_VC-v)/Rs * VC_on : amp (point current)
    V_VC : volt
    VC_on : 1 # if 1, voltage-clamp is on
    
    I_hyp : amp (point current)
    
    '''
    neuron = SpatialNeuron(morphology=morpho, model=eqs, Cm=Cm, Ri=Ri,
                                    namespace=dict(gL=gL, EL=EL, ENa=ENa, Rs=Rs,
                                                   Ka=Ka, Va=Va, Taum_max=Taum_max,                 
                                                   gclamp = gclamp),
                           method="exponential_euler")
    
    # Initial segment                           
    if Na_start == Na_end:
        initial_segment = morpho.axon[Na_start]
        neuron.gNa[initial_segment] = Gna/neuron.area[initial_segment]
    else:
        initial_segment = morpho.axon[Na_start:Na_end]
        neuron.gNa[initial_segment] =  Gna/(pi*AIS_length*axon_diam) 
    
    # Initialisation
    neuron.v = EL

    return neuron
    
    

        
    
    