#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Parameters describing a cerebellar granule cell and a patch clamp pipette.
"""

from brian2 import *

# Passive parameters
gL = 5.* (siemens / meter ** 2)
Ri = 100. * ohm * cm

# Measured CGC parameters
Cm_tot = 3.*pF #3.* pF # cell capacitance
Ce_tot = 5.*pF/100 # electrode capacitance

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
gna_soma = 1000. * (siemens / meter ** 2) #  to have a true spike at the soma
gna_dend = 20. * (siemens / meter ** 2) 

# K channels parameters
gk_soma = 500. * (siemens / meter ** 2) # 150. * (siemens / meter ** 2)  
gk_dend = 20. * (siemens / meter ** 2) 

## Channels kinetics
# Na+:
Va_soma = -30.*mV #-29 * mV  
Vh_soma = -60.*mV  # # #-59 * mV  

# prediction
y0 = 1.19968

