#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Functions to measure diverse features of the models:
    - resting membrane potential
    - threshold for spike initiation
"""
from brian2 import *
#from shared.analysis import *

__all__ = ['calculate_resting_state', 'measure_threshold', 'measure_rheobase', 'measure_threshold_in_vc', 'measure_rheobase_peak', 'spike_initiation_site']

defaultclock.dt = 0.005*ms
#defaultclock.dt = 0.005*ms

def calculate_resting_state(neuron, Na_start=5.*um, Na_end=30.*um, hyperpol_current = False):
 
    M = StateMonitor(neuron, ('v'), record = 0)
    M_AIS = StateMonitor(neuron, ('v'), record = neuron.morphology.axon[Na_end]) 
    
    if hyperpol_current == True:
        ais_idx = neuron.morphology.n + neuron.morphology.dendrite.n + int((Na_end*2)/um)
        neuron.I_CC[ais_idx] = -0.2*nA
        run(150*ms)
    else: 
        run(150.*ms)
    
    return M.t/ms, M.v[0]/mV, M_AIS.t/ms, M_AIS.v[0]/mV


def spike_initiation_site(params, neuron, resting_vm=-75.*mV, ais_start=5.*um, ais_end=30.*um, pulse_length = 2.*ms, i_pulse = 1.5*nA, i_inj = 0):
    
    M = StateMonitor(neuron, ('v','I_VC'), record = 0)
    M_AIS = StateMonitor(neuron, ('v', 'm'), record = neuron.morphology.axon[ais_start:ais_end]) 
        
    #inj_idx = neuron.morphology.n + neuron.morphology.dendrite.n + int((ais_end*2)/um) 

    neuron.V_VC[0] = resting_vm 
    neuron.VC_on[0] = 1
    #neuron.I_CC[inj_idx] = i_inj
    run(20*ms)
    neuron.VC_on[0] = 0
    neuron.I_CC[0] = i_pulse
    #neuron.I_CC[inj_idx] = i_inj
    run(pulse_length)
    neuron.I_CC[0] = 0*amp 
    #neuron.I_CC[inj_idx] = i_inj
    run(20*ms)

    #find the site of spike initiation
    print len(M_AIS.m)
    time_of_initiation = np.where(M_AIS.m[int((ais_end-ais_start)/um)-1, :] >= 0.5)[0][0] 
    print time_of_initiation

    # Site of spike initiation: look for the location of the first occuring spike 
    si_idx = argmax(M_AIS.v[:int((ais_end-ais_start)/um), time_of_initiation]) 
    print 'SI site:', si_idx*um + ais_start
    print M_AIS.v[si_idx-1, time_of_initiation], M_AIS.v[si_idx, time_of_initiation], M_AIS.v[si_idx+1, time_of_initiation]
    
    return si_idx*um + ais_start


def measure_rheobase(params, neuron, resting_vm=-75.*mV, ais_start=5.*um, ais_end=30.*um, pulse_length = 2.*ms, i_inj = 0, si_site = False, latency = 20.*ms, plot_v = True, color='k'):
    print ais_start, ais_end
    
    if ais_start > ais_end:
        print 'break'
        i_current = nan
        i_previous = nan
    else:      
        #peak = 20.*mV
        i_max = 3.*nA
        i_min = 0.*nA
        i_current = 1.5*nA
        i_previous = i_current
        spike = False
        
        M = StateMonitor(neuron, ('v'), record = 0)
        
        if si_site:
            M_AIS = StateMonitor(neuron, ('v', 'm'), record = neuron.morphology.axon[si_site])
        else:
            M_AIS = StateMonitor(neuron, ('v', 'm'), record = neuron.morphology.axon[ais_end])
            
        store()
        
        inj_idx = neuron.morphology.n + neuron.morphology.dendrite.n + int(ais_end/um) 
       
        while True:
            #print i_min, i_current, i_max
            print i_previous
            
            restore()
            
            neuron.ve[0] = resting_vm 
            neuron.VC_on[0] = 1
            neuron.I[inj_idx] = i_inj
            run(latency)
            neuron.VC_on[0] = 0
            neuron.I[0] = i_current
            neuron.I[inj_idx] = i_inj
            run(pulse_length)
            neuron.I[0] = 0*amp 
            neuron.I[inj_idx] = i_inj
            run(20*ms)
            
            figure('measure rheobase')
            plot(M.t/ms, M.v[0]/mV, color=color)
            plot(M_AIS.t/ms, M_AIS.v[0]/mV, '--', color=color)
            
            m_max_ais = max(M_AIS.m[0][int(latency/defaultclock.dt):]) #int(22.*ms/defaultclock.dt)])
            #v_max = max(M_AIS.v[0][int(20.*ms/defaultclock.dt):])

            if m_max_ais >= 0.5 and abs(i_current - i_min) <= 0.5*pA and spike == False :
                print 'stop'
                plot(M.t/ms, M.v[0]/mV, 'r', label='V soma rheobase')
                plot(M_AIS.t/ms, M_AIS.v[0]/mV, '--r', label='V AIS rheobase')
                legend()
                break
            if m_max_ais <= 0.5:
                i_min = i_current
                spike = False
            else: 
                i_max = i_current
                spike = True
            
            i_previous = i_current
            i_current = 0.5*i_max + 0.5*i_min
                
    return i_previous, i_current

def measure_threshold(params, neuron, resting_vm=-75.*mV, ais_start=5.*um, ais_end=30.*um, i_pulse = 0.5*nA, pulse_length = 2.*ms, i_inj = 0, latency = 20.*ms, si_site = False, plot_v = True, color = 'k'):

    if ais_start > ais_end:
        print 'break'
        th_rheo = nan
    else:   
        M = StateMonitor(neuron, ('v'), record = 0)
        
        if si_site:
            M_AIS = StateMonitor(neuron, ('v', 'm'), record = neuron.morphology.axon[si_site])
        else:
            M_AIS = StateMonitor(neuron, ('v', 'm'), record = neuron.morphology.axon[ais_end])
        
        inj_idx = neuron.morphology.n + neuron.morphology.dendrite.n + int(ais_end/um)
        
        neuron.ve[0] = resting_vm #params.EL #-80*mV
        neuron.VC_on[0] = 1
        neuron.I[inj_idx] = i_inj
        run(latency)
        neuron.VC_on[0] = 0
        neuron.I[inj_idx] = i_inj
        neuron.I[0] = i_pulse
        run(pulse_length)
        neuron.I[inj_idx] = i_inj
        neuron.I[0] = 0*amp 
        run(20*ms)
        
        th_rheo_soma = max(M.v[0][int(latency/defaultclock.dt):]) #int(22.*ms/defaultclock.dt)])
        th_rheo_ais = max(M_AIS.v[0][int(latency/defaultclock.dt):]) 
                     
    return th_rheo_soma, th_rheo_ais, M.v[0] #th_s_dV, th_s_m005, th_s_m05, th_a, th_rheo

def measure_rheobase_peak(params, neuron, resting_vm=-75.*mV, ais_start=5.*um, ais_end=30.*um, i_inj = 0, plot_v = True, color='k'):
    print ais_start, ais_end
    
    if ais_start > ais_end:
        print 'break'
        i_current = nan
        i_previous = nan
    else:      
        peak = 0.*mV
        i_max = 3.*nA
        i_min = 0.*nA
        i_current = 1.5*nA
        i_previous = i_current
        spike = False
        
        M = StateMonitor(neuron, ('v','I_VC'), record = 0)
        M_AIS = StateMonitor(neuron, ('v', 'm'), record = neuron.morphology.axon[ais_end]) # -0.5*um]) # because at axon[ais_end]: no channels
        store()
        
        inj_idx = neuron.morphology.n + neuron.morphology.dendrite.n + int((ais_end*2)/um) 
        #print inj_idx
        
        while True:
            #print i_min, i_current, i_max
            print i_previous
            
            restore()
            
            neuron.V_VC[0] = resting_vm #params.EL #-75.*mV
            neuron.VC_on[0] = 1
            neuron.I_CC[inj_idx] = i_inj
            run(20*ms)
            neuron.VC_on[0] = 0
            neuron.I_CC[0] = i_current
            neuron.I_CC[inj_idx] = i_inj
            run(2.*ms)
            neuron.I_CC[0] = 0*amp 
            neuron.I_CC[inj_idx] = i_inj
            run(10*ms)
            
            figure('rheobase')
            plot(M.t/ms, M.v[0]/mV, color=color)
            plot(M_AIS.t/ms, M_AIS.v[0]/mV, '--', color=color)
            
            v_max = max(M_AIS.v[0][int(20.*ms/defaultclock.dt):])

            if v_max >= peak and abs(i_current - i_min) <= .5*pA and spike == False :
                print 'stop'
                plot(M.t/ms, M.v[0]/mV, 'r', label='V soma rheobase')
                plot(M_AIS.t/ms, M_AIS.v[0]/mV, '--r', label='V AIS rheobase')
                legend()
                break
            if v_max < peak:
                i_min = i_current
                spike = False
            else: 
                i_max = i_current
                spike = True
            
            i_previous = i_current
            i_current = 0.5*i_max + 0.5*i_min
                
    return i_previous, i_current

def measure_axonal_threshold(params, neuron, resting_vm=-75.*mV, ais_start=5.*um, ais_end=30.*um, i_pulse = 0.5*nA, hyperpol_current = False, plot_v = True, color = 'k'):

    if ais_start > ais_end:
        print 'break'
        th_rheo = nan
    else:   
        M = StateMonitor(neuron, ('v','I_VC'), record = 0)
        M_AIS = StateMonitor(neuron, ('v', 'm'), record = neuron.morphology.axon[ais_end]) 
        
        #end_idx = neuron.morphology.n + neuron.morphology.dendrite.n + int((ais_end*2)/um) 
        inj_idx = neuron.morphology.n + neuron.morphology.dendrite.n + int((ais_end*2)/um) 
        
        if hyperpol_current == True:
            i_inj = -0.1*nA
        else:
            i_inj = 0.*nA
        
        neuron.V_VC[0] = resting_vm #params.EL #-80*mV
        neuron.VC_on[0] = 1
        neuron.I_CC[inj_idx] = i_inj
        run(20*ms)
        neuron.VC_on[0] = 0
        neuron.I_CC[inj_idx] = i_inj
        neuron.I_CC[0] = i_pulse
        run(2.*ms)
        neuron.I_CC[inj_idx] = i_inj
        neuron.I_CC[0] = 0*amp 
        run(20*ms)
    
        i_spike_a = spike_onsets(M_AIS.v[0], criterion = 20*volt/second * defaultclock.dt, v_peak=-20 * mV)
        th_a = M_AIS.v[0][i_spike_a[0]]/mV
    
        return th_a
    
    
def measure_threshold_in_vc(params, neuron, ais_start=5.*um, ais_end=30.*um, hyperpol_current = False, plot_v = True, color='k'):
    
    """
    Measures the threshold as the lowest somatic membrane potential that triggers a spike. 
    """
    print ais_start, ais_end
    
    if ais_start > ais_end:
        print 'break'
        v_current = nan
    else:
        
        peak = 0.*mV
        v_max = -25.*mV
        v_min = -75.*mV
        v_current = -50.*mV
        #v_previous = v_current
        spike = False
        
        M = StateMonitor(neuron, ('v','I_VC'), record = 0)
        M_AIS = StateMonitor(neuron, ('v', 'm'), record = neuron.morphology.axon[ais_end]) 
        
        store()
        
        while True:
            print v_min, v_current, v_max
            #print v_previous
            restore()
            neuron.V_VC[0] = params.EL #-75.*mV
            neuron.VC_on[0] = 1
            run(20*ms)
            neuron.V_VC[0] = v_current
            neuron.VC_on[0] = 1
            run(20.*ms)
            neuron.VC_on[0] = 0
            run(20*ms)
            
            #figure('Threshold VC')
            #plot(M.t/ms, M.v[0]/mV, color=color)
            #plot(M_AIS.t/ms, M_AIS.v[0]/mV, '--', color=color)
    
            v_max_ais = max(M_AIS.v[0][int(20.*ms/defaultclock.dt):int(40.*ms/defaultclock.dt)])
            #m_max_ais = max(M_AIS.m[0][int(20.*ms/defaultclock.dt):int(40.*ms/defaultclock.dt)])
            
            if v_max_ais >= peak and abs(v_current - v_min) <= .1*mV and spike == False :
            #if 0.095 <= m_max_ais <= 0.105: 
                print 'threshold:', v_current
                break
            if v_max_ais < peak:
            #if m_max_ais <= 0.5:
                v_min = v_current
                spike = False
            else: 
                v_max = v_current
                spike = True
            
            #v_previous = v_current
            v_current = 0.5*v_max + 0.5*v_min
        
    return v_current #v_previous

### With Si site

#def measure_threshold(params, neuron, resting_vm=-75.*mV, ais_start=5.*um, ais_end=30.*um, i_pulse = 0.5*nA, pulse_length = 2.*ms, i_inj = 0, plot_v = True, color = 'k'):
#
#    if ais_start > ais_end:
#        print 'break'
#        th_rheo = nan
#    else:   
#        M = StateMonitor(neuron, ('v','I_VC'), record = 0)
#        M_AIS = StateMonitor(neuron, ('v', 'm'), record = neuron.morphology.axon[ais_end]) 
#        
#        #end_idx = neuron.morphology.n + neuron.morphology.dendrite.n + int((ais_end*2)/um) 
#        inj_idx = neuron.morphology.n + neuron.morphology.dendrite.n + int((ais_end*2)/um)
#        #inj_idx = neuron.morphology.n + neuron.morphology.dendrite.n + int((ais_start+ais_end)/um) # injection in the AIS middle
#        
#        neuron.V_VC[0] = resting_vm #params.EL #-80*mV
#        neuron.VC_on[0] = 1
#        neuron.I_CC[inj_idx] = i_inj
#        run(20*ms)
#        neuron.VC_on[0] = 0
#        neuron.I_CC[inj_idx] = i_inj
#        neuron.I_CC[0] = i_pulse
#        run(pulse_length)
#        neuron.I_CC[inj_idx] = i_inj
#        neuron.I_CC[0] = 0*amp 
#        run(20*ms)
#        
#        #if plot_v ==True:
#        #    figure('Vm')
#        #    plot(M.t/ms, M.v[0]/mV, color=color)
#        #    plot(M_AIS.t/ms, M_AIS.v[0]/mV, '--', color = color)     
#        
#        # soma
#        #i_spike_dV = spike_onsets(M.v[0], criterion = 20*volt/second * defaultclock.dt, v_peak=-20 * mV)
#        #i_spike_m005 = [k for k in range(len(M_AIS.t/ms)) if M_AIS.m[0][k] >= 0.1][0]
#        #i_spike_m05 = [k for k in range(len(M_AIS.t/ms)) if M_AIS.m[0][k] >= 0.5][0]
#        #i_spike_rheo = argmax(M.v[0][int(20.*ms/defaultclock.dt):int(22.*ms/defaultclock.dt)+1]/mV) + 2000
#        #thres_rheo = M.v[0][i_spike_rheo] /mV
#        
#        ## AIS
#        #i_spike_a = spike_onsets(M_AIS.v[0], criterion = 20*volt/second * defaultclock.dt, v_peak=-20 * mV)
#        
#        
#        #th_s_m005 = 0. #M.v[0][i_spike_m005]/mV
#        #th_s_m05 = 0. #M.v[0][i_spike_m05]/mV
#        #th_a = 0. #M_AIS.v[0][i_spike_a[0]]/mV
#        th_rheo_soma = max(M.v[0][int(20.*ms/defaultclock.dt):]) #int(22.*ms/defaultclock.dt)])
#        th_rheo_ais = max(M_AIS.v[0][int(20.*ms/defaultclock.dt):]) 
#        
#        #print "Threshold m>0.05:",M.v[0][i_spike_m005]
#        #print "Threshold m>0.5:",M.v[0][i_spike_m05]
#            
#        #if size(i_spike_dV) == 0:
#        #    th_s_dV = nan
#        #else:
#        #    th_s_dV = M.v[0][i_spike_dV[0]]/mV
#        #    print "Threshold dV:",M.v[0][i_spike_dV[0]]
#            
#        #if plot_v ==True:  
#            #plot(M.t[int(22.*ms/defaultclock.dt)+1]/ms,th_rheo/mV,'xr')
#            #if size(i_spike_dV) != 0:
#            #    plot(M.t[i_spike_dV[0]]/ms,th_s_dV,'^r')
#                
#            #plot(M.t[i_spike_m005]/ms,th_s_m005,'.g')
#            #plot(M.t[i_spike_m05]/ms,th_s_m05,'og')
#            #plot(M_AIS.t[i_spike_a]/ms, M_AIS.v[0][i_spike_a]/mV,'^b')
#            #show()
#            
#        
#    return th_rheo_soma, th_rheo_ais, M.v[0] #th_s_dV, th_s_m005, th_s_m05, th_a, th_rheo
#    
    
    
    
    
    
    
    
    
    
    
    
    

