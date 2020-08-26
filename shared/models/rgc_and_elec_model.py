

'''

A simple biophysical model of a large cell with an electrode.

'''
from brian2 import *

__all__ = ['rgc_model', 'rgc_and_series_resistance_model']

def rgc_model(params, resting_vm = -75.*mV, Na_start = 5.*um, Na_end = 35.*um, density = False, gna_tot = 350.*nS, gk_tot = 150.*nS, gna_density = 4000.* (siemens / meter ** 2), morpho = None, i_inj = 0):
    '''
    A simple model of a big cell with a large dendrite and a thin axon.
    
    If the morphology is not given, then it is created to a default.

    neuron.I_CC is a current-clamp stimulation.
    neuron.V_VC is a voltage-clamp stimulation.
    neuron.VC_on is a voltage-clamp switch (0/1).
    neuron.I_VC is the voltage-clamp current.
    '''
    
    # Passive parameters
    EL = resting_vm 
    Ri = params.Ri 
    
    # Morphology: 
    # Neuron
    soma_diam = 25.*um
    dend_diam = 6.*um 
    dend_length = 500.*um 
    axon_diam = 1. * um
    axon_length = 1000. * um
    
    # Morpho: soma at 0, electrode, dendrite from compartment 1 to 501, axon from compartment 501 to 1501
    morpho = Soma(soma_diam)
    morpho.dendrite = Cylinder(diameter=dend_diam, length=dend_length, n=20)
    morpho.axon = Cylinder(diameter=axon_diam, length=axon_length, n=1000)
    
    # Neuron
    soma_area = morpho.area[0]
    Cm = params.Cm_tot/soma_area
    gL = params.gL #1./(100*Mohm * soma_area)
    print (gL)
    
    
    ### Na channels distribution
    AIS_length = Na_end-Na_start

    ### AIS

    # Na channels parameters
    ENa = params.ENa #70. * mV
    Gna = gna_tot #params.Gna #370.* nS  # NEW
    
    # K channels parameters
    EK = params.EK #-90. * mV
    Gk = gk_tot
    
    Va = params.Va #-35. * mV  
    Ka = params.Ka #5. * mV  
    Taum_max = params.Taum_max # factor * 0.15 * ms 
    Vh = params.Vh #-65. * mV  
    Kh = params.Kh #5. * mV  
    Tauh_max = params.Tauh_max# factor * 5. * ms  

    # K+:
    Vn = params.Vn #
    Kn = params.Kn #18. * mV  
    Taun_max = params.Taun_max #1. * ms  

    ### Soma

    # Na channels parameters
    gna_soma = params.gna_soma #700. * (siemens / meter ** 2) 
    gna_dend = params.gna_dend #50. * (siemens / meter ** 2) 

    # K channels parameters
    gk_soma = params.gk_soma #500. * (siemens / meter ** 2)
    gk_dend = params.gk_dend #50. * (siemens / meter ** 2) 

    ## Channels kinetics
    # Na+:
    Va_soma = params.Va_soma #-30.*mV #-29 * mV  
    Vh_soma = params.Vh_soma #-60.*mV #-59 * mV 

    # Equations
    eqs = '''
    Im = (gL*(EL-v) + gNa*m*h*(ENa-v) + gK*n**8*(EK-v)) : amp/meter**2
    INa = gNa*m*h*(ENa-v) : amp/meter**2
    IK = gK*n**8*(EK-v) : amp/meter**2

    dm/dt = alpham*(1-m) - betam*m : 1
    dh/dt = alphah*(1-h) - betah*h : 1
    dn/dt = alphan*(1-n) - betan*n : 1

    alpham = (1/ka)*(v-va) / (1-exp(-(v-va)/ka)) /(2*taum_max) : Hz
    betam = (1/ka)*(-v+va) / (1-exp((v-va)/ka)) /(2*taum_max) : Hz

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
    gclamp = 100*usiemens
    
    eqs += '''
    I_CC : amp (point current)
    I_VC = gclamp*VC_on*(V_VC - v) : amp (point current)
    V_VC : volt
    VC_on : 1 # if 1, voltage-clamp is on
    
    I_hyp : amp (point current)
    '''
    
    neuron = SpatialNeuron(morphology=morpho, model=eqs, Cm=Cm, Ri=Ri,
                                    namespace=dict(gL=gL, EL=EL, ENa=ENa, EK=EK,
                                                   Ka=Ka, Va=Va, Taum_max=Taum_max,
                                                   Kh=Kh, Vh=Vh, Tauh_max=Tauh_max,
                                                   Kn=Kn, Vn=Vn, Taun_max=Taun_max,
                                                   Va_soma=Va_soma, 
                                                   Vh_soma=Vh_soma,                                                   
                                                   gclamp = gclamp
                                                   ),
                           method="exponential_euler")

    # Parameters of the soma are assigned in the entire neuron for integration purpose, 
    # but gNa, gK=0 outside the soma and AIS
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
    #neuron.gL = gL_true

    # Soma
    neuron.gNa[0] = gna_soma
    neuron.gK[0] = gk_soma
    
    # Initial segment
    
    if Na_start == Na_end:
        initial_segment = morpho.axon[Na_start]
        if density == False:
            neuron.gNa[initial_segment] = Gna/neuron.area[initial_segment]
            neuron.gK[initial_segment] = Gk/neuron.area[initial_segment]
        else:
            neuron.gNa[initial_segment] = gna_density
            neuron.gK[initial_segment] = params.gk_dens
    else:
        initial_segment = morpho.axon[Na_start:Na_end]
        if density == False:
            neuron.gNa[initial_segment] =  Gna/(pi*AIS_length*morpho.axon.diameter[0])
            neuron.gK[initial_segment] =  Gk/(pi*AIS_length*morpho.axon.diameter[0])
        else:
            neuron.gNa[initial_segment] = gna_density
            neuron.gK[initial_segment] = params.gk_dens
    
    neuron.va[initial_segment] = Va
    neuron.vh[initial_segment] = Vh
    
    
    
    # Initialisation
    neuron.v = EL
    neuron.h = 1

    return neuron

def rgc_and_series_resistance_model(params,  resting_vm = -75.*mV, Na_start = 5.*um, Na_end = 35.*um, density = False, gna_tot = 350.*nS, gk_tot = 150.*nS, gna_density = 4000.* (siemens / meter ** 2), morpho = None, i_inj = 0, Rs = 10*mohm):
    '''
    If the morphology is not given, then it is created to a default.

    neuron.I_CC is a current-clamp stimulation.
    neuron.V_VC is a voltage-clamp stimulation.
    neuron.VC_on is a voltage-clamp switch (0/1).
    neuron.I_VC is the voltage-clamp current.
    '''
    
    # Passive parameters
    EL = resting_vm
    Ri = params.Ri

    # Morphology: 
    # Neuron
    soma_diam = 25.*um
    dend_diam = 6.*um 
    dend_length = 500.*um 
    axon_diam = 1. * um
    axon_length = 1000. * um
    
    # Morpho: soma at 0, electrode, dendrite from compartment 1 to 501, axon from compartment 501 to 1501
    morpho = Soma(soma_diam)
    morpho.dendrite = Cylinder(diameter=dend_diam, length=dend_length, n=20)
    morpho.axon = Cylinder(diameter=axon_diam, length=axon_length, n=1000)
    
    # Neuron
    soma_area = morpho.area[0]
    Cm = params.Cm_tot/soma_area
    gL_true = params.gL #1./(100*Mohm * soma_area)
    print (gL_true)
    
    # Electrode
    Ce = params.Ce_tot
    
    ### Na channels distribution
    AIS_length = Na_end-Na_start

    ### AIS

    # Na channels parameters
    ENa = params.ENa #70. * mV
    Gna = gna_tot #params.Gna #370.* nS  # NEW
    
    # K channels parameters
    EK = params.EK #-90. * mV
    Gk = gk_tot
    
    Va = params.Va #-35. * mV  
    Ka = params.Ka #5. * mV  
    Taum_max = params.Taum_max # factor * 0.15 * ms 
    Vh = params.Vh #-65. * mV  
    Kh = params.Kh #5. * mV  
    Tauh_max = params.Tauh_max# factor * 5. * ms  

    # K+:
    Vn = params.Vn #
    Kn = params.Kn #18. * mV  
    Taun_max = params.Taun_max #1. * ms  

    ### Soma

    # Na channels parameters
    gna_soma = params.gna_soma #700. * (siemens / meter ** 2) 
    gna_dend = params.gna_dend #50. * (siemens / meter ** 2) 

    # K channels parameters
    gk_soma = params.gk_soma #500. * (siemens / meter ** 2)
    gk_dend = params.gk_dend #50. * (siemens / meter ** 2) 

    ## Channels kinetics
    # Na+:
    Va_soma = params.Va_soma #-30.*mV #-29 * mV  
    Vh_soma = params.Vh_soma #-60.*mV #-59 * mV 

    # Equations
    eqs = '''
    Im = (gL*(EL-v) + gNa*m*h*(ENa-v) + gK*n**8*(EK-v)) : amp/meter**2
    INa = gNa*m*h*(ENa-v) : amp/meter**2
    IK = gK*n**8*(EK-v) : amp/meter**2
    
    ve : volt
    I_VC = (ve-v)/Rs * VC_on : amp (point current)
    VC_on : 1 # if 1, voltage-clamp is on 

    
    dm/dt = alpham*(1-m) - betam*m : 1
    dh/dt = alphah*(1-h) - betah*h : 1
    dn/dt = alphan*(1-n) - betan*n : 1

    alpham = (1/ka)*(v-va) / (1-exp(-(v-va)/ka)) /(2*taum_max) : Hz
    betam = (1/ka)*(-v+va) / (1-exp((v-va)/ka)) /(2*taum_max) : Hz

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
    #gclamp = 100*usiemens

    #eqs += '''
    
    
    
    #I_CC : amp (point current)
    #I_VC = gclamp*VC_on*(V_VC - ve) : amp (point current)
    #V_VC : volt
    #VC_on : 1 # if 1, voltage-clamp is on
    
    #I_hyp : amp (point current)
    #'''

    neuron = SpatialNeuron(morphology=morpho, model=eqs, Cm=Cm, Ri=Ri,
                                    namespace=dict(EL=EL, ENa=ENa, EK=EK, Ce = Ce, soma_area=soma_area,
                                                   Ka=Ka, Va=Va, Taum_max=Taum_max,
                                                   Kh=Kh, Vh=Vh, Tauh_max=Tauh_max,
                                                   Kn=Kn, Vn=Vn, Taun_max=Taun_max,
                                                   Va_soma=Va_soma, 
                                                   Vh_soma=Vh_soma, 
                                                   Rs=Rs
                                                   #gclamp = gclamp
                                                   ),
                           method="exponential_euler")

    # Parameters of the soma are assigned in the entire neuron for integration purpose, 
    # but gNa, gK=0 outside the soma and AIS
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
    neuron.gL = gL_true

    # Soma
    neuron.gNa[0] = gna_soma
    neuron.gK[0] = gk_soma
    
    # Initial segment
    
    if Na_start == Na_end:
        initial_segment = morpho.axon[Na_start]
        if density == False:
            neuron.gNa[initial_segment] = Gna/neuron.area[initial_segment]
            neuron.gK[initial_segment] = Gk/neuron.area[initial_segment]
        else:
            neuron.gNa[initial_segment] = gna_density
            neuron.gK[initial_segment] = params.gk_dens
    else:
        initial_segment = morpho.axon[Na_start:Na_end]
        if density == False:
            neuron.gNa[initial_segment] =  Gna/(pi*AIS_length*morpho.axon.diameter[0])
            neuron.gK[initial_segment] =  Gk/(pi*AIS_length*morpho.axon.diameter[0])
        else:
            neuron.gNa[initial_segment] = gna_density
            neuron.gK[initial_segment] = params.gk_dens
    
    neuron.va[initial_segment] = Va
    neuron.vh[initial_segment] = Vh
    
    # Initialisation
    neuron.v = EL
    #neuron.ve = EL
    neuron.h = 1
    neuron.VC_on = 0

    return neuron