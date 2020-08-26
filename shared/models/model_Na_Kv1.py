

'''
Biophysical model of an AP.

This is a model with a spherical isopotential soma, a large dendrite and an unmyelinated axon. 
The AIS is located in the proximal axon and has a spatial extend or not (if Na_start = Na_end). 

The AIS contains a higher density of Na and K channels than the soma, the dendrite and the rest of the axon.

'''
from brian2 import *

__all__ = ['model_Na_Kv1']

def model_Na_Kv1(params,  resting_vm = -75.*mV, Na_start = 5.*um, Na_end = 35.*um, \
                 density = False, gna_tot = 350.*nS, gk_tot = 150.*nS, gL_soma = 1./(15000. * ohm * cm**2) , \
                 gna_density = 4000.* (siemens / meter ** 2), gk_density = 1500.* (siemens / meter ** 2), morpho = None):
    '''
    If the morphology is not given, then it is created to a default. If density is False, the total number of Na and Kv channels
    in the AIS is fixed. 
    
    params: a file that contains all the passive and channels parameters,
    
    neuron.I_CC is a current-clamp stimulation.
    neuron.V_VC is a voltage-clamp stimulation.
    neuron.VC_on is a voltage-clamp switch (0/1).
    neuron.I_VC is the voltage-clamp current.
    '''

    ### Passive parameters
    EL = resting_vm 
    Cm = params.Cm 
    gL = params.gL 
    Ri = params.Ri 

    if morpho is None:
        dend_diam = 6.*um 
        dend_length = 1000.*um 
        axon_diam = 1. * um
        axon_length = 500. * um 
        soma_diameter = 30. * um
        morpho = Soma(diameter=soma_diameter)
        dendrite = Cylinder(diameter=dend_diam, length=dend_length, n=500)
        axon = Cylinder(diameter=axon_diam, length=axon_length, n=500)
        morpho.dendrite = dendrite
        morpho.axon = axon
 
    ### Na channels distribution
    AIS_length = Na_end-Na_start

    ### AIS

    # Na channels parameters
    ENa = params.ENa 
    Gna = gna_tot 
    
    # K channels parameters
    EK = params.EK 
    Gk = gk_tot
    
    Va = params.Va 
    Ka = params.Ka 
    Taum_max = params.Taum_max 
    Vh = params.Vh 
    Kh = params.Kh 
    Tauh_max = params.Tauh_max

    # K+:
    Vn = params.Vn 
    Kn = params.Kn 
    Taun_max = params.Taun_max 

    ### Soma

    # Na channels parameters
    gna_soma = params.gna_soma 
    gna_dend = params.gna_dend 

    # K channels parameters
    gk_soma = params.gk_soma 
    gk_dend = params.gk_dend 

    ## Channels kinetics
    # Na+:
    Va_soma = params.Va_soma 
    Vh_soma = params.Vh_soma 

    # Equations
    eqs = '''
    Im = (gL*(EL-v) + gNa*m*h*(ENa-v) + gK*n**8*(EK-v) ) : amp/meter**2
    INa = gNa*m*h*(ENa-v) : amp/meter**2
    IK = gK*n**8*(EK-v) : amp/meter**2

    dm/dt = alpham*(1-m) - betam*m : 1
    dh/dt = alphah*(1-h) - betah*h : 1
    dn/dt = alphan*(1-n) - betan*n : 1

    alpham = (1/ka)*(v-va) / (1-exp(-(v-va)/ka)) /(2*taum_max) : Hz
    betam = -(1/ka)*(v-va) / (1-exp((v-va)/ka)) /(2*taum_max) : Hz

    alphah = -(1/kh)*(v-vh) / (1-exp((v-vh)/kh)) /(2*tauh_max) : Hz
    betah = (1/kh)*(v-vh) / (1-exp(-(v-vh)/kh)) /(2*tauh_max) : Hz

    alphan = (1/kn)*(v-vn) / (1-exp(-(v-vn)/kn)) /(2*taun_max) : Hz
    betan = -(1/kn)*(v-vn) / (1-exp((v-vn)/kn)) /(2*taun_max): Hz
    
    gL: siemens/meter**2
    gNa : siemens/meter**2
    gK : siemens/meter**2
    va : volt
    vh : volt
    vn : volt
    ka : volt
    kh : volt
    kn : volt
    taum_max : second
    tauh_max : second
    taun_max : second

    I : amp (point current)
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
                                    namespace=dict(EL=EL, ENa=ENa, EK=EK,
                                                   Ka=Ka, Va=Va, Taum_max=Taum_max,
                                                   Kh=Kh, Vh=Vh, Tauh_max=Tauh_max,
                                                   Kn=Kn, Vn=Vn, Taun_max=Taun_max,
                                                   Va_soma=Va_soma, 
                                                   Vh_soma=Vh_soma,                                                  
                                                   gclamp = gclamp),
                                                   method="exponential_euler")

    # Parameters of the soma are assigned in the entire neuron for integration purpose, 
    # but gNa and gK at the soma and AIS are then modified
    neuron.va = Va_soma
    neuron.ka = Ka
    neuron.taum_max = Taum_max

    neuron.vh = Vh_soma
    neuron.kh = Kh
    neuron.tauh_max = Tauh_max

    neuron.vn = Vn
    neuron.kn = Kn
    neuron.taun_max = Taun_max

    #Dendrites and axon
    neuron.gNa = gna_dend
    neuron.gK = gk_dend
    neuron.gL = gL

    # Soma
    neuron.gNa[0] = gna_soma
    neuron.gK[0] = gk_soma
    neuron.gL[0] = gL_soma
    
    # Initial segment
    if Na_start == Na_end: # point AIS
        initial_segment = morpho.axon[Na_start]
        if density == False:
            neuron.gNa[initial_segment] = Gna/neuron.area[initial_segment]
            neuron.gK[initial_segment] = Gk/neuron.area[initial_segment]
        else:
            neuron.gNa[initial_segment] = gna_density
            neuron.gK[initial_segment] = gk_density
    else: # extended AIS
        initial_segment = morpho.axon[Na_start:Na_end] 
        if density == False:
            neuron.gNa[initial_segment] =  Gna/(pi*AIS_length*morpho.axon.diameter[0])
            neuron.gK[initial_segment] =  Gk/(pi*AIS_length*morpho.axon.diameter[0])
        else:
            neuron.gNa[initial_segment] = gna_density
            neuron.gK[initial_segment] = gk_density
    
    neuron.va[initial_segment] = Va
    neuron.vh[initial_segment] = Vh
    
    # Initialisation
    neuron.v = EL
    neuron.h = 1

    return neuron
