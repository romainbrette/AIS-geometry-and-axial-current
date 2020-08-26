
"""

Analysis of spontaneous APs of retinal ganglion cells.

"""

from brian2 import *
import glob2
import pandas as pd
import pyabf
from scipy import interpolate
from pandas import ExcelWriter
from pandas import ExcelFile
from trace_analysis import *

rcParams['axes.spines.right'] = False
rcParams['axes.spines.top'] = False

### Path to datafiles: load the list of cells used for the analysis
path_to_files = '/Users/sarah/Documents/repositories/AIS-geometry-and-axonal-current/Na currents in RGC/codes submission/data/'

df_cells = pd.read_excel(path_to_files + 'RGC_electrical_properties.xlsx')
first_cell = 0
last_cell = 3 #len(df_cells['Date']) 

dates = array(df_cells['Date'])[first_cell:last_cell]
retinas = array(df_cells['Retina'])[first_cell:last_cell]
cells = array(df_cells['Cell'])[first_cell:last_cell]
ages = array(df_cells['Age'])[first_cell:last_cell]
v_ends = array(df_cells['V end (mV)'])[first_cell:last_cell]

### Path to the data
path_to_data = '/Users/sarah//Documents/Data/Martijn Sierksma/'

### Loading and analysing the data
N = 0

selected_dates = []
selected_retinas = []
selected_cells = []
selected_ages = []
selected_sweeps = []

dvdt_axonal_peak = []
v_axonal_peak = []
dvdt_max1 = []
dvdt_max2 = []
v_max1 = []
v_max2 = []
dvdt_somatic_onset = []
v_somatic_onset = []

ap_onsets = []
rise_times = []
spike_widths = []
ap_peaks = []
ap_amplitude = []

for date, retina, cell, age in zip(dates, retinas, cells, ages): 
    print ('----------------')
    print (date, retina, cell)

    ### Path to spontaneous activity data
    path_to_cell = path_to_data + str(int(date)) + "*/" + '/retina '+ str(retina) +'/cell ' + str(int(cell))
    path_to_cc_cont = glob2.glob(path_to_cell + '/CC cont/' + '*' + ".abf")
    print (path_to_cc_cont)
    
    if len(path_to_cc_cont) > 0:
        abf = pyabf.ABF(path_to_cc_cont[0]) 
        V = abf.sweepY * mV
        fs = abf.dataRate  * Hz # sampling rate
        dt = 1./fs
        t = dt*arange(len(V)) 
        
        ### Analysis of AP shape (first spontaneous spike)
        # First spike
        spike_times = find_spikes_at(V, dt/ms, -20.*mV)/(dt/ms)
        if len(spike_times) > 0:
            idx_spike1 = int(spike_times[0] )
            print ('Spike at:', idx_spike1)
        else: # no spike in CC cont
            print ('No spike in CC cont recording')
            continue
        
        # AP shape analysis
        if len(V[idx_spike1-100:idx_spike1]) == 0:
            pass
        elif mean(V[idx_spike1-100:idx_spike1])/mV > 20.:
            pass
        else:
            # Add the cell to the list of selected cells
            selected_dates.append(date)
            selected_retinas.append(retina)
            selected_cells.append(cell)
            selected_ages.append(age)
            
            # Interpolating V and computing dV and d2V
            f = V[idx_spike1-70:idx_spike1+30]
            t_spike = t[idx_spike1-70:idx_spike1+30]
            t_new = (t_spike[:-1] + t_spike[1:])/2
            v = (f[:-1] + f[1:])/2
            dv = (f[1:] - f[:-1])/dt
            ddv = (dv[1:] - dv[:-1])/dt # shift of dt: add 1 !!! (f[2:] - 2*f[1:-1] + f[:-2])/dt**2 #
            
            # AP peak
            idx_peak = argmax(v)
            
            # Spike onset 
            spike_onset = spike_onsets(v, criterion = 20*volt/second * dt, v_peak = -30.*mV)

            if len(spike_onset) > 0:
                if v[spike_onset[0]] > 0:
                    spike_onset = spike_onsets(v/mV, v_peak = -30.)
            else:
                dvdt_max1.append(nan)
                dvdt_max2.append(nan)
                v_max1.append(nan)
                v_max2.append(nan)
                ap_onsets.append(nan)
                dvdt_somatic_onset.append(nan)
                v_somatic_onset.append(nan)
                continue

            print ('Onset:', spike_onset)
            idx_ax_onset = spike_onset[0] - 1 # because the function shifts by +1
            
            # Global max of dvdt after spike onset
            dvdt_max = argmax(dv[idx_ax_onset:]) + idx_ax_onset
            # Global max of the dV^2/dt^2
            ddvdt_max = argmax(ddv[idx_ax_onset:]) + idx_ax_onset
            
            # The global max of dvdt can be in the axonal component:
            # we look for an inflexion point between the onset and the max dvdt:
            # if yes: it is the axonal max, the global max is somatic max
            # if not: the global max is axonal max
            inflexion_before_global_max = where([ddv[i]*ddv[i+1]<0 for i in range(idx_ax_onset+1, dvdt_max-2)])[0]
            print(dvdt_max, inflexion_before_global_max + idx_ax_onset+1)
            
            if len(inflexion_before_global_max) < 1:  # the global max is the axonal max
                # the axonal max might not be a local max,
                # so we verifiy that there is no decceleration between spike onset and the max
                if ddvdt_max != idx_ax_onset:
                    print('A')
                    ddvdt_min = argmin(ddv[idx_ax_onset+1:ddvdt_max+1])+ idx_ax_onset + 1 + 1
                else:
                    print('B')
                    ddvdt_min = argmin(ddv[idx_ax_onset:ddvdt_max+1])+ idx_ax_onset + 1
                # axonal max
                idx_dvdt_max1 = ddvdt_min
                # somatic max
                idx_dvdt_max2 = dvdt_max
                
                # we look for the somatic max as the next inflexion point
                if len(where([ddv[i]*ddv[i+1]<0  for i in range(dvdt_max+1, idx_peak)])[0]) != 0 : # if another local max after the global max
                    print('Global max is axonal max')   
                    ddvdt_min = dvdt_max
                    extr = where([ddv[i]*ddv[i+1]<0  for i in range(dvdt_max+1, idx_peak)])[0] + dvdt_max + 1 + 1
                    dvdt_max = array(extr)[argmax(dv[extr])] 
                    # axonal max
                    idx_dvdt_max1 = ddvdt_min
                    # somatic max
                    idx_dvdt_max2 = dvdt_max 
                elif ddvdt_min == ddvdt_max:
                    print('C')
                    ddvdt_min = argmin(ddv[idx_ax_onset+1:ddvdt_max])+ ddvdt_max + 1 + 1
                    # axonal max
                    idx_dvdt_max1 = ddvdt_min
                    # somatic max
                    idx_dvdt_max2 = dvdt_max 
                elif dv[ddvdt_min] < dv[idx_ax_onset]:
                    print('D')
                    ddvdt_min = argmin(ddv[idx_ax_onset:dvdt_max+1])+ idx_ax_onset + 1
                    # axonal max
                    idx_dvdt_max1 = ddvdt_min
                    # somatic max
                    idx_dvdt_max2 = dvdt_max 
            else: # the global max is the somatic max
                print('Global max is somatic max')
                # axonal max
                idx_dvdt_max1 = inflexion_before_global_max[0] + idx_ax_onset + 1 + 1
                # somatic max
                idx_dvdt_max2 = dvdt_max
                
            print(idx_dvdt_max1, idx_dvdt_max2)
                
            # Somatic regeneration as the max acceleration between the two local max 
            ddvdt_max_between = argmax(ddv[idx_dvdt_max1:idx_dvdt_max2]) + idx_dvdt_max1 
            idx_som_onset = ddvdt_max_between
            print (ddv[idx_dvdt_max1:idx_dvdt_max2])
            
            t_dvdt_max1 = t_new[idx_dvdt_max1]/ms
            t_dvdt_max2 = t_new[idx_dvdt_max2]/ms
            dv_dvdt_max1 = dv[idx_dvdt_max1]/(mV/ms)
            dv_dvdt_max2 = dv[idx_dvdt_max2]/(mV/ms)
            v_dvdt_max1 = v[idx_dvdt_max1]/mV
            v_dvdt_max2 = v[idx_dvdt_max2]/mV
            
            t_ax_onset = t_new[idx_ax_onset]/ms
            v_ax_onset = v[idx_ax_onset]/mV
            dvdt_ax_onset = dv[idx_ax_onset]/(mV/ms)
            
            t_som_onset = t_new[idx_som_onset]/ms
            v_som_onset = v[idx_som_onset]/mV
            dvdt_som_onset = dv[idx_som_onset]/(mV/ms)
            
            dvdt_max1.append(dv_dvdt_max1)
            dvdt_max2.append(dv_dvdt_max2)
            v_max1.append(v_dvdt_max1)
            v_max2.append(v_dvdt_max2)
            ap_onsets.append(v_ax_onset)
            dvdt_somatic_onset.append(dvdt_som_onset)
            v_somatic_onset.append(v_som_onset)

            N +=1
               
            ### Plotting
            f1 = figure('AP %i,  %s, %i' %(date, retina, cell), (10,8))
            
            # V vs t
            subplot(221)
            plot(t_new/ms, v/mV, 'k-')
            plot(t_ax_onset, v_ax_onset, 'bo', label='spike onset')
            plot(t_dvdt_max1, v_dvdt_max1, 'ro', label='first max dV/dt')
            plot(t_dvdt_max2, v_dvdt_max2, 'go', label='second max dV/dt')
            plot(t_som_onset, v_som_onset, 'yo', label='somatic regeneration')
            legend(frameon=False)
            ylabel('V (mV)')
            xlabel('t (ms)')
            
            # dV vs t
            subplot(222)
            plot(t_new/ms, dv, 'k-')
            plot(t_ax_onset, dvdt_ax_onset, 'bo')
            plot(t_dvdt_max2, dv_dvdt_max2, 'go')
            plot(t_dvdt_max1, dv_dvdt_max1, 'ro')
            plot(t_som_onset, dvdt_som_onset, 'yo')
            ylabel('dV/dt (mV/ms)')
            xlabel('t (ms)')
            
            # dV vs V (phase plot)
            subplot(223)
            plot(v/mV, dv, 'k-')
            plot(v_ax_onset, dvdt_ax_onset, 'bo')
            plot(v_dvdt_max2, dv_dvdt_max2, 'go')
            plot(v_dvdt_max1, dv_dvdt_max1, 'ro')
            plot(v_som_onset, dvdt_som_onset, 'yo')
            ylabel('dV/dt (mV/ms)')
            xlabel('V (mV)')
                
            # d2V vs t
            subplot(224)
            plot(t_new[:-1]/ms, ddv, 'k-')
            plot(t_ax_onset, ddv[idx_ax_onset], 'bo')
            plot(t_dvdt_max2, ddv[idx_dvdt_max2], 'go')
            plot(t_dvdt_max1, ddv[idx_dvdt_max1], 'ro')
            plot(t_som_onset, ddv[idx_som_onset], 'yo')
            ylabel('d2V/dt2 (mV/ms2)')
            xlabel('t (ms)')
                                             
tight_layout()

show()



# ### Write the results in an excel file

# df_select_cells = pd.DataFrame({'Date': selected_dates,
#                   'Retina': selected_retinas,
#                   'Cell': selected_cells,
#                   'Age': selected_ages,
#                   'AP onset': ap_onsets,
#                   'dvdt max1': dvdt_max1,
#                   'dvdt max2': dvdt_max2,
#                   'v max1': v_max1,
#                   'v max2': v_max2,
#                   'dvdt somatic onset': dvdt_somatic_onset,
#                   'v somatic onset': v_somatic_onset
#                   })

# save_path = '/Users/sarahgoethals/Dropbox/Spike initiation/PhD projects/Axonal current and AIS geometry/Data patch/RGC/'
# df_select_cells.to_excel(save_path + "RGC_action_potential.xlsx", \
#                 columns=['Date','Retina','Cell','Age','AP onset',\
#                           'dvdt max1', 'dvdt max2', 'v max1', 'v max2', \
#                            'dvdt somatic onset', 'v somatic onset'])
















