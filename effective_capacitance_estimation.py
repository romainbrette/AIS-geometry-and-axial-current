

"""

Exponential fit to voltage decay measured in current-clamp to estimate the membrane time constant tau_m.

OK

"""
from brian2 import *
import os
import glob2
import pandas as pd
import pyabf
from pandas import ExcelWriter
from pandas import ExcelFile
from trace_analysis import *
from scipy.optimize import curve_fit
import statistics

rcParams['axes.spines.right'] = False
rcParams['axes.spines.top'] = False

### Laod list of cells used for the analysis
df_cells = pd.read_excel('RGC_electrical_properties.xlsx')

first_cell = -5
last_cell = -3 #    len(df_cells['Date'])
dates = array(df_cells['Date'])[first_cell:last_cell]
retinas = array(df_cells['Retina'])[first_cell:last_cell]
cells = array(df_cells['Cell'])[first_cell:last_cell]
ages = array(df_cells['Age'])[first_cell:last_cell]
          
# Path to the data
path_to_data = 'data/RGC data/'

N = 0

selected_dates = []
selected_retinas = []
selected_cells = []
selected_ages = []
selected_sweeps = []
delta_v_2ms = []
membrane_time_cst = []
capacitances = []


for date, retina, cell, age in zip(dates, retinas, cells, ages): 
    print (date, retina, cell)
    
    path_to_cell = path_to_data + str(int(date)) + "*/" + '/retina '+ str(retina) +'/cell ' + str(int(cell))

    # path to the steps CC recordings
    path_to_cc_steps = glob2.glob(path_to_cell + '/CC steps/' + '*' + ".abf")
    print (path_to_cc_steps)
    
    if len(path_to_cc_steps) > 0:
        
        # Loading and plotting the data
        abf = pyabf.ABF(path_to_cc_steps[0]) 
        fs = abf.dataRate  * Hz # sampling rate
        dt = 1./fs
        t = dt*arange(len(abf.sweepY)) 
        n_rec = len(abf.sweepList)
        
        cell_time_cst = []
        cell_cm = []
        
        if n_rec > 5:
            selected_dates.append(date)
            selected_retinas.append(retina)
            selected_cells.append(cell)
        
            f1 = figure('%i,  %s, %i' %(date, retina, cell), (12,9))
            ax2 = f1.add_subplot(211)
            ax2.set_xlim(0, 1500)
            ax2.set_ylabel('Vc (mV)')
            ax3 = f1.add_subplot(212)
            #ax3.set_xlim(0, 200)
            ax3.set_ylabel('Vc (mV)')
            ax3.set_xlabel('Time (ms)')
            
            I = []
            V = []
            V_baseline = zeros(n_rec)
            V_pulse = zeros(n_rec)
            I_pulse = zeros(n_rec)
        
            for sweepNumber in range(n_rec):
                abf.setSweep(sweepNumber)
                V.append(abf.sweepY)
                I.append(abf.sweepC)

                I_pulse[sweepNumber] = I[sweepNumber][int(250.*ms/dt)]
                V_baseline[sweepNumber] = mean(V[sweepNumber][int(50*ms/dt):int(125*ms/dt)]) #mV
                V_pulse[sweepNumber] = mean(V[sweepNumber][int(800*ms/dt):int(1100*ms/dt)]) #mV
                
                ax2.plot(t/ms, abf.sweepY, 'k')
            
            idx_pulse_start = where(I[-1] == max(I[-1]))[0][0] + int(.25*ms/dt)

            # exponential fit to the voltage decay
            colors = ['green', 'blue', 'orange','magenta','cyan', 'yellow','purple']
            for i in range(0, 7):
                
                voltage_decay = V[i][idx_pulse_start:idx_pulse_start+int(1*ms/dt)] # mV
                decay_time = (t[idx_pulse_start:idx_pulse_start+int(1*ms/dt)]- t[idx_pulse_start])/ms  # ms
                i_pulse = I[i][idx_pulse_start+int(1*ms/dt)]
                
                ax2.plot(decay_time + t[idx_pulse_start]/ms, voltage_decay, 'r-')#, color = colors[i])
                ax2.set_xlim(150, 1200)
                
                v_start = V[i][idx_pulse_start]
                
                def exp_voltage(t, tau1, V1):
                    return v_start - (v_start - V1) * (1.-exp(-t/tau1)) 
                
                # Tau
                p0 = [2., -60.]
                tau_opt = curve_fit(exp_voltage, decay_time, voltage_decay, p0, bounds=([0, -90], [10, -50]))
                print ('Taum:', tau_opt[0])
                tau_m = tau_opt[0][0]
                cm = tau_opt[0][0] / ((tau_opt[0][1]-v_start)/(i_pulse))
                print ('Cm:', cm)
                
                cell_time_cst.append(tau_m)
                cell_cm.append(cm)
                
                ax3.plot(decay_time, voltage_decay, 'k-')
                ax3.plot(decay_time, exp_voltage(decay_time, tau_opt[0][0], tau_opt[0][1]), '-', color = colors[i], 
                                                 label= '$t_m$ = %0.2f ms, $C_m = %0.2f$ pF' %(tau_m, cm))
                ax3.legend(frameon=False)
                
            time_cst = statistics.median(cell_time_cst)
            capacitance = statistics.median(cell_cm)
            
            print ('Cm final:', capacitance)
            
            ax2.set_title('Membrane time constant: %0.2f ms' %time_cst)
            membrane_time_cst.append(time_cst)
            capacitances.append(capacitance)
             
            tight_layout()

show()

# ### Write in excel file

# df_select_cells = pd.DataFrame({'Date': selected_dates,
#                   'Retina': selected_retinas,
#                   'Cell': selected_cells,
#                   'Tau m': membrane_time_cst,
#                   'Cm': capacitances
#                   })

# df_select_cells.to_excel("RGC_capacitance_test.xlsx", columns=['Date','Retina','Cell','Tau m', 'Cm'])









