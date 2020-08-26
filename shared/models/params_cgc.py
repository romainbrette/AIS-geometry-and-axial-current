#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Parameters describing a cerebellar granule cell.
"""

from brian2 import *

### Morphology: soma at 0, dendrite from compartment 1 to 501, axon from compartment 501 to 1501
morpho = Soma(5. * um)
dend_diam = 2.*um 
dend_length = 20.*um 
axon_diam = .2 * um
axon_length = 1000. * um
dendrite = Cylinder(diameter=dend_diam, length=dend_length, n=20)
axon = Cylinder(diameter=axon_diam, length=axon_length, n=1000)
morpho.dendrite = dendrite
morpho.axon = axon

# Passive parameters
EL = -90. * mV
Cm = 0.9* uF / cm ** 2
gL = 10.* (siemens / meter ** 2)
Ri = 100. * ohm * cm

### AIS
# Na channels parameters
ENa = 70. * mV
gna_dens = 3500. * (siemens / meter ** 2)
Gna = 400.* nS

# K channels parameters
EK = -90. * mV
gk_dens = 1000. * (siemens / meter ** 2)
Gk = 150.*nS

# Channels kinetics
T = 33. #33. 
factor = (1. / 2.5) ** ((T - 23.) / 10.) #(1 / 2.8) ** ((T - 23.) / 10.)

# Na+
Va = -35. * mV  
Ka = 5. * mV  
Taum_max = factor * 0.15 * ms
Vh = -65.*mV 
Kh = 5. * mV  
Tauh_max =  factor * 5. * ms  

# K+:
Vn = -73. * mV  
Kn = 18. * mV  
Taun_max = 1. * ms 

### Soma

# Na channels parameters
gna_soma = 350. * (siemens / meter ** 2) #  350. * (siemens / meter ** 2)  #
gna_dend = 20. * (siemens / meter ** 2) 

# K channels parameters
gk_soma = 100. * (siemens / meter ** 2) # 150. * (siemens / meter ** 2)  
gk_dend = 20. * (siemens / meter ** 2) 

## Channels kinetics
# Na+:
Va_soma = -30.*mV #-29 * mV  
Vh_soma = -60.*mV  # # #-59 * mV  

# prediction
y0 = 1.19968

