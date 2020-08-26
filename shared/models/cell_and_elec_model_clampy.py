'''
Model with an electrode. 
'''
import sys
sys.path.append("/Users/sarah/Documents/repositories/clampy/") #clampy/clampy/")
from clampy import *
from brian2 import *
from clampy.brianmodels import *

__all__ = ['AxonalInitiationModel_with_Rs', 'CGCandElectrodeModel', 'CGCandSeriesResistanceModel', 'RGCandSeriesResistanceModel']

class AxonalInitiationModel_with_Rs(SpatialBrianExperiment):
    '''
    A spatially extended model with a dendrite and axon, including an initial segment.
    
    With electrode.
    '''
    def __init__(self, params, resting_vm = -70.*mV, Na_start = 5.*um, Na_end = 35.*um, \
                 density = False, i_inj = 0*amp, gclamp = 10000*usiemens, dt = 0.02*ms, Rs=0):
         
        # Passive parameters
        EL = resting_vm
        Ri = params.Ri
    
        # Morphology: 
        # Neuron
        soma_diam = 30.*um
        dend_diam = 6.*um 
        dend_length = 1000.*um 
        axon_diam = 1. * um
        axon_length = 1000. * um
        
        # Morpho: soma at 0, electrode, dendrite from compartment 1 to 501, axon from compartment 501 to 1501
        morpho = Soma(soma_diam)
        morpho.dendrite = Cylinder(diameter=dend_diam, length=dend_length, n=20)
        morpho.axon = Cylinder(diameter=axon_diam, length=axon_length, n=1000)
        
        # Neuron
        Cm = params.Cm 
        gL_true = params.gL 
                
        ### Na channels distribution
        AIS_length = Na_end-Na_start
    
        ### AIS
    
        # Na channels parameters
        ENa = params.ENa #70. * mV
        Gna = params.Gna #370.* nS  # NEW
        
        # K channels parameters
        EK = params.EK #-90. * mV
        Gk = params.Gk
        
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
        
        I_VC = (Vcommand(t-t_start)-v)/Rs * VC_on : amp (point current)
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
        
        V = v : volt
        t_start : second
        gclamp : siemens
        
        I : amp (point current)
    
        '''

        Board.__init__(self)

        self.dt = dt
        self.gclamp = gclamp
        self.eqs = Equations(eqs)        
        
        neuron = SpatialNeuron(morphology=morpho, model=eqs, Cm=Cm, Ri=Ri,
                               method="exponential_euler",
                               namespace = dict(EL=EL, ENa=ENa, EK=EK,  Rs=Rs,
                                                       Ka=Ka, Va=Va, Taum_max=Taum_max,
                                                       Kh=Kh, Vh=Vh, Tauh_max=Tauh_max,
                                                       Kn=Kn, Vn=Vn, Taun_max=Taun_max,
                                                       Va_soma=Va_soma, 
                                                       Vh_soma=Vh_soma,
                                                       gclamp=gclamp))
    
        self.neuron = neuron
            
        # Parameters of the soma are assigned in the entire neuron for integration purpose, 
        # but gNa, gK=0 outside the soma and AIS
        self.neuron.va = Va_soma
        self.neuron.ka = Ka
        self.neuron.taum_max = Taum_max
    
        self.neuron.vh = Vh_soma
        self.neuron.kh = Kh
        self.neuron.tauh_max = Tauh_max
    
        self.neuron.vn = Vn
        self.neuron.kn = Kn
        self.neuron.taun_max = Taun_max
    
        #Dendrites and axon
        self.neuron.gNa = gna_dend
        self.neuron.gK = gk_dend
        self.neuron.gL = gL_true
    
        # Soma
        self.neuron.gNa[0] = gna_soma
        self.neuron.gK[0] = gk_soma
        
        # Initial segment
        
        if Na_start == Na_end:
            initial_segment = morpho.axon[Na_start]
            if density == False:
                self.neuron.gNa[initial_segment] = Gna/neuron.area[initial_segment]
                self.neuron.gK[initial_segment] = Gk/neuron.area[initial_segment]
            else:
                self.neuron.gNa[initial_segment] = params.gna_dens
                self.neuron.gK[initial_segment] = params.gk_dens
        else:
            initial_segment = morpho.axon[Na_start:Na_end]
            if density == False:
                self.neuron.gNa[initial_segment] =  Gna/(pi*AIS_length*morpho.axon.diameter[0])
                self.neuron.gK[initial_segment] =  Gk/(pi*AIS_length*morpho.axon.diameter[0])
            else:
                self.neuron.gNa[initial_segment] = params.gna_dens
                self.neuron.gK[initial_segment] = params.gk_dens
        
        self.neuron.va[initial_segment] = Va
        self.neuron.vh[initial_segment] = Vh
        
        self.neuron.gclamp = 10000*usiemens
        
        # Initialisation
        self.neuron.v = EL
        self.neuron.h = 1
        self.neuron.VC_on = 0
        self.neuron.VC_on[0] = 1
        
        # Set up the experiment      
        self.network = Network(self.neuron)
        self.configure_board()
        self.is_voltage_clamp = False # Initially in current clamp  
        



class CGCandElectrodeModel(BrianExperiment):
    '''
    A spatially extended model with a dendrite and axon, including an initial segment.
    '''
    def __init__(self, params, resting_vm = -90.*mV, Na_start = 5.*um, Na_end = 35.*um, \
                 density = False, i_inj = 0*amp, Rs = 10.*Mohm,  gclamp = 1*usiemens, dt = 0.02*ms):
        
        
        # Passive parameters
        EL = resting_vm
        Ri = params.Ri
    
        # Morphology: 
        # Neuron
        soma_diam = 5.*um
        dend_diam = 2.*um 
        dend_length = 20.*um 
        axon_diam = .2 * um
        axon_length = 1000. * um
        
        # Morpho: soma at 0, electrode, dendrite from compartment 1 to 501, axon from compartment 501 to 1501
        morpho = Soma(soma_diam)
        morpho.dendrite = Cylinder(diameter=dend_diam, length=dend_length, n=20)
        morpho.axon = Cylinder(diameter=axon_diam, length=axon_length, n=1000)
        
        # Neuron
        soma_area = morpho.area[0]
        Cm = params.Cm_tot/soma_area
        gL_true = 1./(1.5*Gohm * soma_area)
        
        # Electrode
        Ce = params.Ce_tot
        
        ### Na channels distribution
        AIS_length = Na_end-Na_start
    
        ### AIS
    
        # Na channels parameters
        ENa = params.ENa #70. * mV
        Gna = params.Gna #370.* nS  # NEW
        
        # K channels parameters
        EK = params.EK #-90. * mV
        Gk = params.Gk
        
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
        
        Ie = (ve-v)/Rs : amp (point current)
        dve/dt = (I - Ie)/Ce : volt
        
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
        #ve : volt
    
        '''
    
        # Current clamp or voltage-clamp       
        eqs += '''
        V = v : volt
        I = CC_switch*Icommand(t-t_start) + Iclamp : amp (point current)
        Iclamp = gclamp*(Vcommand(t-t_start)-ve) : amp
        CC_switch : 1
        gclamp : siemens
        t_start : second
        '''
        
        Board.__init__(self)

        self.dt = dt
        self.gclamp = gclamp
        self.eqs = Equations(eqs)        
        
        neuron = SpatialNeuron(morphology=morpho, model=eqs, Cm=Cm, Ri=Ri,
                               method="exponential_euler",
                               namespace = dict(EL=EL, ENa=ENa, EK=EK, Ce = Ce, Rs=Rs, soma_area=soma_area,
                                                       Ka=Ka, Va=Va, Taum_max=Taum_max,
                                                       Kh=Kh, Vh=Vh, Tauh_max=Tauh_max,
                                                       Kn=Kn, Vn=Vn, Taun_max=Taun_max,
                                                       Va_soma=Va_soma, 
                                                       Vh_soma=Vh_soma))
    
        neuron.CC_switch[0] = 1
        self.neuron = neuron
            
        # Parameters of the soma are assigned in the entire neuron for integration purpose, 
        # but gNa, gK=0 outside the soma and AIS
        self.neuron.va = Va_soma
        self.neuron.ka = Ka
        self.neuron.taum_max = Taum_max
    
        self.neuron.vh = Vh_soma
        self.neuron.kh = Kh
        self.neuron.tauh_max = Tauh_max
    
        self.neuron.vn = Vn
        self.neuron.kn = Kn
        self.neuron.taun_max = Taun_max
    
        #Dendrites and axon
        self.neuron.gNa = gna_dend
        self.neuron.gK = gk_dend
        self.neuron.gL = gL_true
    
        # Soma
        self.neuron.gNa[0] = gna_soma
        self.neuron.gK[0] = gk_soma
        
        # Initial segment
        
        if Na_start == Na_end:
            initial_segment = morpho.axon[Na_start]
            if density == False:
                self.neuron.gNa[initial_segment] = Gna/neuron.area[initial_segment]
                self.neuron.gK[initial_segment] = Gk/neuron.area[initial_segment]
            else:
                self.neuron.gNa[initial_segment] = params.gna_dens
                self.neuron.gK[initial_segment] = params.gk_dens
        else:
            initial_segment = morpho.axon[Na_start:Na_end]
            if density == False:
                self.neuron.gNa[initial_segment] =  Gna/(pi*AIS_length*morpho.axon.diameter[0])
                self.neuron.gK[initial_segment] =  Gk/(pi*AIS_length*morpho.axon.diameter[0])
            else:
                self.neuron.gNa[initial_segment] = params.gna_dens
                self.neuron.gK[initial_segment] = params.gk_dens
        
        self.neuron.va[initial_segment] = Va
        self.neuron.vh[initial_segment] = Vh
        
        # Initialisation
        self.neuron.v = EL
        self.neuron.ve = EL
        self.neuron.h = 1
    
        # Set up the experiment      
        self.network = Network(self.neuron)
        self.configure_board()
        self.is_voltage_clamp = False # Initially in current clamp

class CGCandSeriesResistanceModel(BrianExperiment):
    '''
    A spatially extended model with a dendrite and axon, including an initial segment.
    '''
    def __init__(self, params, resting_vm = -90.*mV, Na_start = 5.*um, Na_end = 35.*um, \
                 density = False, i_inj = 0*amp, Rs = 10.*Mohm, Cm_tot = 3.*pF, gclamp = 1*usiemens, dt = 0.02*ms):
         
        # Passive parameters
        EL = resting_vm
        Ri = params.Ri
    
        # Morphology: 
        # Neuron
        soma_diam = 5.*um
        dend_diam = 1.*um 
        dend_length = 20.*um 
        axon_diam = .2 * um
        axon_length = 1000. * um
        
        # Morpho: soma at 0, electrode, dendrite from compartment 1 to 501, axon from compartment 501 to 1501
        morpho = Soma(soma_diam)
        morpho.dendrite = Cylinder(diameter=dend_diam, length=dend_length, n=20)
        morpho.axon = Cylinder(diameter=axon_diam, length=axon_length, n=1000)
        
        # Neuron
        soma_area = morpho.area[0]
        print (soma_area)
        Cm = Cm_tot/soma_area
        gL_true = 1./(1.5*Gohm * soma_area)
                
        ### Na channels distribution
        AIS_length = Na_end-Na_start
    
        ### AIS
    
        # Na channels parameters
        ENa = params.ENa #70. * mV
        Gna = params.Gna #370.* nS  # NEW
        
        # K channels parameters
        EK = params.EK #-90. * mV
        Gk = params.Gk
        
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
        
        I_VC = (Vcommand(t-t_start)-v)/Rs * VC_on : amp (point current)
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
        
        V = v : volt
        t_start : second
        gclamp : siemens
        
        I : amp (point current)
    
        '''

        Board.__init__(self)

        self.dt = dt
        self.gclamp = gclamp
        self.eqs = Equations(eqs)        
        
        neuron = SpatialNeuron(morphology=morpho, model=eqs, Cm=Cm, Ri=Ri,
                               method="exponential_euler",
                               namespace = dict(EL=EL, ENa=ENa, EK=EK, Rs=Rs, soma_area=soma_area,
                                                       Ka=Ka, Va=Va, Taum_max=Taum_max,
                                                       Kh=Kh, Vh=Vh, Tauh_max=Tauh_max,
                                                       Kn=Kn, Vn=Vn, Taun_max=Taun_max,
                                                       Va_soma=Va_soma, 
                                                       Vh_soma=Vh_soma,
                                                       gclamp=gclamp))
    
        self.neuron = neuron
            
        # Parameters of the soma are assigned in the entire neuron for integration purpose, 
        # but gNa, gK=0 outside the soma and AIS
        self.neuron.va = Va_soma
        self.neuron.ka = Ka
        self.neuron.taum_max = Taum_max
    
        self.neuron.vh = Vh_soma
        self.neuron.kh = Kh
        self.neuron.tauh_max = Tauh_max
    
        self.neuron.vn = Vn
        self.neuron.kn = Kn
        self.neuron.taun_max = Taun_max
    
        #Dendrites and axon
        self.neuron.gNa = gna_dend
        self.neuron.gK = gk_dend
        self.neuron.gL = gL_true
    
        # Soma
        self.neuron.gNa[0] = gna_soma
        self.neuron.gK[0] = gk_soma
        
        # Initial segment
        
        if Na_start == Na_end:
            initial_segment = morpho.axon[Na_start]
            if density == False:
                self.neuron.gNa[initial_segment] = Gna/neuron.area[initial_segment]
                self.neuron.gK[initial_segment] = Gk/neuron.area[initial_segment]
            else:
                self.neuron.gNa[initial_segment] = params.gna_dens
                self.neuron.gK[initial_segment] = params.gk_dens
        else:
            initial_segment = morpho.axon[Na_start:Na_end]
            if density == False:
                self.neuron.gNa[initial_segment] =  Gna/(pi*AIS_length*morpho.axon.diameter[0])
                self.neuron.gK[initial_segment] =  Gk/(pi*AIS_length*morpho.axon.diameter[0])
            else:
                self.neuron.gNa[initial_segment] = params.gna_dens
                self.neuron.gK[initial_segment] = params.gk_dens
        
        self.neuron.va[initial_segment] = Va
        self.neuron.vh[initial_segment] = Vh
        
        # Initialisation
        self.neuron.v = EL
        self.neuron.h = 1
        self.neuron.VC_on = 0
        self.neuron.VC_on[0] = 1
        
        # Set up the experiment      
        self.network = Network(self.neuron)
        self.configure_board()
        self.is_voltage_clamp = False # Initially in current clamp  
        
class CGCandElectrodeModel(BrianExperiment):
    '''
    A spatially extended model with a dendrite and axon, including an initial segment.
    '''
    def __init__(self, params, resting_vm = -90.*mV, Na_start = 5.*um, Na_end = 35.*um, \
                 density = False, i_inj = 0*amp, Rs = 10.*Mohm,  gclamp = 1*usiemens, dt = 0.02*ms):
        
        
        # Passive parameters
        EL = resting_vm
        Ri = params.Ri
    
        # Morphology: 
        # Neuron
        soma_diam = 5.*um
        dend_diam = 2.*um 
        dend_length = 20.*um 
        axon_diam = .2 * um
        axon_length = 1000. * um
        
        # Morpho: soma at 0, electrode, dendrite from compartment 1 to 501, axon from compartment 501 to 1501
        morpho = Soma(soma_diam)
        morpho.dendrite = Cylinder(diameter=dend_diam, length=dend_length, n=20)
        morpho.axon = Cylinder(diameter=axon_diam, length=axon_length, n=1000)
        
        # Neuron
        soma_area = morpho.area[0]
        Cm = params.Cm_tot/soma_area
        gL_true = 1./(1.5*Gohm * soma_area)
        
        # Electrode
        Ce = params.Ce_tot
        
        ### Na channels distribution
        AIS_length = Na_end-Na_start
    
        ### AIS
    
        # Na channels parameters
        ENa = params.ENa #70. * mV
        Gna = params.Gna #370.* nS  # NEW
        
        # K channels parameters
        EK = params.EK #-90. * mV
        Gk = params.Gk
        
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
        
        Ie = (ve-v)/Rs : amp (point current)
        dve/dt = (I - Ie)/Ce : volt
        
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
        #ve : volt
    
        '''
    
        # Current clamp or voltage-clamp       
        eqs += '''
        V = v : volt
        I = CC_switch*Icommand(t-t_start) + Iclamp : amp (point current)
        Iclamp = gclamp*(Vcommand(t-t_start)-ve) : amp
        CC_switch : 1
        gclamp : siemens
        t_start : second
        '''
        
        Board.__init__(self)

        self.dt = dt
        self.gclamp = gclamp
        self.eqs = Equations(eqs)        
        
        neuron = SpatialNeuron(morphology=morpho, model=eqs, Cm=Cm, Ri=Ri,
                               method="exponential_euler",
                               namespace = dict(EL=EL, ENa=ENa, EK=EK, Ce = Ce, Rs=Rs, soma_area=soma_area,
                                                       Ka=Ka, Va=Va, Taum_max=Taum_max,
                                                       Kh=Kh, Vh=Vh, Tauh_max=Tauh_max,
                                                       Kn=Kn, Vn=Vn, Taun_max=Taun_max,
                                                       Va_soma=Va_soma, 
                                                       Vh_soma=Vh_soma))
    
        neuron.CC_switch[0] = 1
        self.neuron = neuron
            
        # Parameters of the soma are assigned in the entire neuron for integration purpose, 
        # but gNa, gK=0 outside the soma and AIS
        self.neuron.va = Va_soma
        self.neuron.ka = Ka
        self.neuron.taum_max = Taum_max
    
        self.neuron.vh = Vh_soma
        self.neuron.kh = Kh
        self.neuron.tauh_max = Tauh_max
    
        self.neuron.vn = Vn
        self.neuron.kn = Kn
        self.neuron.taun_max = Taun_max
    
        #Dendrites and axon
        self.neuron.gNa = gna_dend
        self.neuron.gK = gk_dend
        self.neuron.gL = gL_true
    
        # Soma
        self.neuron.gNa[0] = gna_soma
        self.neuron.gK[0] = gk_soma
        
        # Initial segment
        
        if Na_start == Na_end:
            initial_segment = morpho.axon[Na_start]
            if density == False:
                self.neuron.gNa[initial_segment] = Gna/neuron.area[initial_segment]
                self.neuron.gK[initial_segment] = Gk/neuron.area[initial_segment]
            else:
                self.neuron.gNa[initial_segment] = params.gna_dens
                self.neuron.gK[initial_segment] = params.gk_dens
        else:
            initial_segment = morpho.axon[Na_start:Na_end]
            if density == False:
                self.neuron.gNa[initial_segment] =  Gna/(pi*AIS_length*morpho.axon.diameter[0])
                self.neuron.gK[initial_segment] =  Gk/(pi*AIS_length*morpho.axon.diameter[0])
            else:
                self.neuron.gNa[initial_segment] = params.gna_dens
                self.neuron.gK[initial_segment] = params.gk_dens
        
        self.neuron.va[initial_segment] = Va
        self.neuron.vh[initial_segment] = Vh
        
        # Initialisation
        self.neuron.v = EL
        self.neuron.ve = EL
        self.neuron.h = 1
    
        # Set up the experiment      
        self.network = Network(self.neuron)
        self.configure_board()
        self.is_voltage_clamp = False # Initially in current clamp

class RGCandSeriesResistanceModel(BrianExperiment):
    '''
    A spatially extended model with a dendrite and axon, including an initial segment.
    '''
    def __init__(self, params, resting_vm = -70.*mV, Na_start = 5.*um, Na_end = 35.*um, \
                 density = False, i_inj = 0*amp, Rs = 10.*Mohm, Cm_tot = 50.*pF, gclamp = 1*usiemens, dt = 0.02*ms):
         
        # Passive parameters
        EL = resting_vm
        Ri = params.Ri
    
        # Morphology: 
        # Neuron
        soma_diam = 30.*um
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
        print (soma_area)
        Cm = Cm_tot/soma_area
        gL_true = 10. * (siemens / meter ** 2) # to obtain a realistic input resistance
        #1./(.1*Gohm * soma_area)
        #print gL_true
                
        ### Na channels distribution
        AIS_length = Na_end-Na_start
    
        ### AIS
    
        # Na channels parameters
        ENa = params.ENa #70. * mV
        Gna = params.Gna #370.* nS  # NEW
        
        # K channels parameters
        EK = params.EK #-90. * mV
        Gk = params.Gk
        
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
        
        I_VC = (Vcommand(t-t_start)-v)/Rs * VC_on : amp (point current)
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
        
        V = v : volt
        t_start : second
        gclamp : siemens
        
        I : amp (point current)
    
        '''

        Board.__init__(self)

        self.dt = dt
        self.gclamp = gclamp
        self.eqs = Equations(eqs)        
        
        neuron = SpatialNeuron(morphology=morpho, model=eqs, Cm=Cm, Ri=Ri,
                               method="exponential_euler",
                               namespace = dict(EL=EL, ENa=ENa, EK=EK, Rs=Rs, soma_area=soma_area,
                                                       Ka=Ka, Va=Va, Taum_max=Taum_max,
                                                       Kh=Kh, Vh=Vh, Tauh_max=Tauh_max,
                                                       Kn=Kn, Vn=Vn, Taun_max=Taun_max,
                                                       Va_soma=Va_soma, 
                                                       Vh_soma=Vh_soma,
                                                       gclamp=gclamp))
    
        self.neuron = neuron
            
        # Parameters of the soma are assigned in the entire neuron for integration purpose, 
        # but gNa, gK=0 outside the soma and AIS
        self.neuron.va = Va_soma
        self.neuron.ka = Ka
        self.neuron.taum_max = Taum_max
    
        self.neuron.vh = Vh_soma
        self.neuron.kh = Kh
        self.neuron.tauh_max = Tauh_max
    
        self.neuron.vn = Vn
        self.neuron.kn = Kn
        self.neuron.taun_max = Taun_max
    
        #Dendrites and axon
        self.neuron.gNa = gna_dend
        self.neuron.gK = gk_dend
        self.neuron.gL = gL_true
    
        # Soma
        self.neuron.gNa[0] = gna_soma
        self.neuron.gK[0] = gk_soma
        
        # Initial segment
        
        if Na_start == Na_end:
            initial_segment = morpho.axon[Na_start]
            if density == False:
                self.neuron.gNa[initial_segment] = Gna/neuron.area[initial_segment]
                self.neuron.gK[initial_segment] = Gk/neuron.area[initial_segment]
            else:
                self.neuron.gNa[initial_segment] = params.gna_dens
                self.neuron.gK[initial_segment] = params.gk_dens
        else:
            initial_segment = morpho.axon[Na_start:Na_end]
            if density == False:
                self.neuron.gNa[initial_segment] =  Gna/(pi*AIS_length*morpho.axon.diameter[0])
                self.neuron.gK[initial_segment] =  Gk/(pi*AIS_length*morpho.axon.diameter[0])
            else:
                self.neuron.gNa[initial_segment] = params.gna_dens
                self.neuron.gK[initial_segment] = params.gk_dens
        
        self.neuron.va[initial_segment] = Va
        self.neuron.vh[initial_segment] = Vh
        
        # Initialisation
        self.neuron.v = EL
        self.neuron.h = 1
        self.neuron.VC_on = 0
        self.neuron.VC_on[0] = 1
        
        # Set up the experiment      
        self.network = Network(self.neuron)
        self.configure_board()
        self.is_voltage_clamp = False # Initially in current clamp  
        
        
        