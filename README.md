# AIS geometry and axial current

### Requirements

In addition to standard scientific packages (scipy, matplotlib), this code requires:
* Seaborn
* [Brian 2](http://briansimulator.org)

Electrophysiological data from recorded in retina ganglion cells
(DOI: ) must be
[downloaded]() into the folder
`shared/data/`.

### Measure of AIS geometry

...

## The axial current at spike initiation

In this part we show ...

* `fig1.py`: spontaneous APs of retinal ganglion cells.\
We measure the spike onset, axonal and somatic max of dV/dt and somatic regeneration across the cell population.\
List of cells used for the population analysis (panels E and F): `RGC_electrical_properties.xlsx`\
Analysis of the raw recordings: `action_potential_analysis.py`\
Results of the analysis used to plot the figure: `RGC_action_potential.xlsx`\

* `fig2.py`: recording the axial current.\
We illustrate how we measure the axial currents, how we correct the peak axial current for the effect of the series resistance and how the peak axial current and the voltage threshold depend on the series resistance (in a model and in RGC).\
Model (panel E and F): \
Neuron model (with series resistance): `shared/models/model_Na_Kv1_with_Rs.py`\
voltage clamp protocol to record axial currents for different series resistance: `model_AP_protocol_VC_dichotomy_with_series_resistance.py`\
voltage clamp protocol to record a test pulse: `model_AP_protocol_VC_test_pulse_with_series_resistance.py`\
simulation results are stored in `simulations data/fig2/`\
RGC data (panels G and H):\
List of cells used for the analysis: `RGC_electrical_properties.xlsx`\
Measure of the peak axonal current and voltage threshold from axial current recordings: `axonal_current_analysis.py`\
Results of the analysis used to plot panels G and H: `RGC_electrical_properties.xlsx`\

* `fig3.py`: transmission of the axial current to the soma.\
We compare the peak axial current with the capacitive current at the soma and look at the correlation of the peak axonal current, the charge transferred to the soma and the current duration with the effective capacitance (see Methods).\
List of cells used for the analysis: `RGC_electrical_properties.xlsx`\
Capacitance estimation from the response to current pulses: ?? TODO
Measure of the peak axonal current: `axonal_current_analysis.py`\
Measure of the charge transferred to the soma: `charge_measurement.py`\
Results of the analysis used to plot the figure: `RGC_electrical_properties.xlsx`

* `fig4.py`: geometry of the AIS.\
We illustrate the diversity of the AIS start position and length in our cell population.\
Measure of the AIS geometry: ??\
Results of the images analysis used to plot the figure: `RGC_electrical_properties.xlsx`

* `fig5.py`: predictions of axial current with resistive coupling theory.\
We compare our measures of the peak axonal current in RGC with measures in a model and with theoretical prediction.\
Panel A: voltage profile in the spike initiation model.\
Neuron model: `shared/models/model_spike_initiation.py`\
Run the model: `model_SI_voltage_profile_along_axon.py`\
Results stored in: `model_SI_voltage_profile.npz`\
Panel B: peak axonal current vs AIS start position in the spike initiation model.\
Neuron model: `shared/models/model_spike_initiation.py`\
Run the model: `model_SI_protocol_VC_dichotomy.py` and `model_SI_protocol_VC_test_pulse.py`\
Simulations results analyzed in: `model_SI_peak_current_vs_delta.py`
Results stored in: `model_SI_peak_current_ext_AIS_g5000.npz`\
Panel D: peak axonal current vs AIS start position in RGC compared to the theoretical prediction
Measure of the peak axonal current: `axonal_current_analysis.py`\
Results of the analysis used to plot the panel: `RGC_electrical_properties.xlsx`

## The threshold axial current

In this part, we examinate the sodium current at threshold and how it depends on the AIS geometry.

* `fig6.py`: axial current near threshold.\
We show how the axial current at threshold is measured and that the Na current just below threshold increases very steeply with voltage. We compare the RGC data with the IV curve below threshold in the biophysical model.\
Panel C-E-F: RGC data.\
Measure of the IV curve below threshold: `axonal_current_near_threshold_analysis.py`\
IV curves stored in: `RGC_IV_curves_below_threshold.npz`\
Panel B-D: model.\
Neuron model: `shared/models/model_Na_Kv1.py`\
Run the model: `model_SI_protocol_VC_dichotomy.py` and `model_SI_protocol_VC_test_pulse.py`\
Results stored in: `simulations data/fig6/`\

* `fig7.py`: threshold vs AIS geometry.\
We examinate the dependance of the voltage threshold and the axial current near threshold on the AIS geometry in RGC, and compare with the theoretical prediction and simulations in a biophysicla model.\
Panel A-B-C: voltage threshold\
Measure of the voltage threshold: `axonal_current_analysis.py`\
Results of the analysis used to plot the panel: `RGC_electrical_properties.xlsx`
Panel D-E:
Measure of the axial current near threshold: `axonal_current_near_threshold_analysis.py`\
Results of the analysis used to plot the panel: `RGC_electrical_properties.xlsx`\
Neuron model: `shared/models/model_Na_Kv1.py`\
Run the model: `model_SI_protocol_VC_dichotomy.py` and `model_SI_protocol_VC_test_pulse.py`\
Results stored in: `simulations data/fig7/`\

## Adaptation of the axial current

* `fig8.py`:

* `fig9.py`:

* `fig10.py`:

