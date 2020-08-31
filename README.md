# AIS geometry and axial current

### Requirements

In addition to standard scientific packages (scipy, matplotlib), this code requires:
* Seaborn
* [Brian 2](http://briansimulator.org)

Electrophysiological data recorded in retina ganglion cells and confocal images
(DOI: 10.5281/zenodo.4005628) must be
[downloaded](https://zenodo.org/record/4005629#.X0zMTi3pNBw) into the folder
`/data/`.

### Measure of AIS geometry

To measure the AIS start and end position, we first trace the axon with Vaa3D and save its coordinate in a SWC file. The AIS start position and length are measured from the traced axon and the confocal images og ankyrin G labellings with the script: `axon_tracing_analysis.py`. The measures are stored in `RGC_electrical_properties.xlsx`.

## The axial current at spike initiation

### Action potentials of retinal ganglion cells

* `fig1.py`: we measure the spike onset, axonal and somatic max of dV/dt and somatic regeneration across the cell population.

Analysis of the raw recordings: `action_potential_analysis.py`\
Results of the analysis used to plot the figure: `RGC_action_potential.xlsx`

### Measuring the axial current

* `fig2.py`: we illustrate how we measure the axial currents, how we correct the peak axial current for the effect of the series resistance and how the peak axial current and the voltage threshold depend on the series resistance (in a model and in RGC).

Panel E and F: model \
Neuron model (with series resistance): `shared/models/model_Na_Kv1_with_Rs.py`\
Voltage clamp protocol to record axial currents for different series resistance: `model_AP_protocol_VC_dichotomy_with_series_resistance.py`\
Voltage clamp protocol to record a test pulse: `model_AP_protocol_VC_test_pulse_with_series_resistance.py`\
Simulation results are stored in `simulations data/fig2/`

Panels G and H: RGC data\
Measure of the peak axonal current and voltage threshold from axial current recordings: `axonal_current_analysis.py`\
Results of the analysis used to plot panels G and H: `RGC_electrical_properties.xlsx`

### Transmission of the axial current to the soma

* `fig3.py`: we compare the peak axial current with the capacitive current at the soma and look at the correlation of the peak axonal current, the charge transferred to the soma and the current duration with the effective capacitance (see Methods).

Capacitance estimation from the response to current pulses:  `effective_capacitance_estimation.py`\
Measure of the peak axonal current: `axonal_current_analysis.py`\
Measure of the charge transferred to the soma and the current duration: `charge_measurement.py`\
Results of the analysis used to plot the figure: `RGC_electrical_properties.xlsx`

### Nav conductance density

* `fig4.py`: we illustrate the diversity of the AIS start position and length in our cell population.

* `fig5.py`: predictions of axial current with resistive coupling theory. We compare our measures of the peak axonal current in RGC with measures in a model and with theoretical prediction.

Panel A: voltage profile in the spike initiation model.\
Neuron model: `shared/models/model_spike_initiation.py`\
Run the model: `model_SI_voltage_profile_along_axon.py`\
Results stored in: `model_SI_voltage_profile.npz`

Panel B: peak axonal current vs AIS start position in the spike initiation model.\
Neuron model: `shared/models/model_spike_initiation.py`\
Run the model: `model_SI_protocol_VC_dichotomy.py` and `model_SI_protocol_VC_test_pulse.py`\
Simulations results analyzed in: `model_SI_peak_current_vs_delta.py`
Results stored in: `model_SI_peak_current_ext_AIS_g5000.npz`

Panel D: peak axonal current vs AIS start position in RGC compared to the theoretical prediction\
Measure of the peak axonal current: `axonal_current_analysis.py`\
Results of the analysis used to plot the panel: `RGC_electrical_properties.xlsx`

## The threshold axial current

### Variation of axial current near threshold

* `fig6.py`: we show how the axial current at threshold is measured and that the Na current just below threshold increases very steeply with voltage. We compare the RGC data with the IV curve below threshold in the biophysical model.

Panel C-E-F: RGC data.\
Measure of the IV curve below threshold: `axonal_current_near_threshold_analysis.py`\
IV curves stored in: `RGC_IV_curves_below_threshold.npz`

Panel B-D: model.\
Neuron model: `shared/models/model_Na_Kv1.py`\
Run the model: `model_SI_protocol_VC_dichotomy.py` and `model_SI_protocol_VC_test_pulse.py`\
Results stored in: `simulations data/fig6/`

### Threshold vs. AIS geometry

* `fig7.py`: we examinate the dependance of the voltage threshold and the axial current near threshold on the AIS geometry in RGC, and compare with the theoretical prediction and simulations in a biophysicla model.

Panel A-B-C: voltage threshold\
Measure of the voltage threshold: `axonal_current_analysis.py`\
Results of the analysis used to plot the panel: `RGC_electrical_properties.xlsx`

Panel D-E:\
Measure of the axial current near threshold: `axonal_current_near_threshold_analysis.py`\
Results of the analysis used to plot the panel: `RGC_electrical_properties.xlsx`\
Neuron model: `shared/models/model_Na_Kv1.py`\
Run the model: `model_SI_protocol_VC_dichotomy.py` and `model_SI_protocol_VC_test_pulse.py`\
Results stored in: `simulations data/fig7/`

## Adaptation of the axial current

### Properties of adaptation

* `fig8.py`: we illustrate adaptation of the voltage threshold and of the axial current at spike initiation, and how they co-variate.

Measure of the voltage threshold and peak axial current during adaptation: `adaptation_analysis_with_PN.py` for the cells with P/5 protocol and `adaptation_analysis_no_PN.py` for the other cells (before 20190611)\
Results are stored in: `RGC_adaptation.xlsx`

Panel A-B-C-D-E: voltage threshold adaptation\
Individual Vth vs V0 curves are fitted (see Methods) in a separate script to estimate different parameters: `threshold_adaptation_analysis.py`\
Results are stored in: `RGC_Vi_from_adaptation.xlsx`

Panel J-K-L: covariation\
Individual Vth vs Ip curves are fitted in a separate script: `covariation_during_adaptation_analysis.py`\
Curves are stored in: `RGC_adaptation_covariation.npz`

* `fig9.py`: adaptation of the axial current at threshold

Measure of the axial current near threshold and IV curve below threshold during adaptation: `adaptation_threshold_current_analysis.py`\
Currents are stored in: `RGC_threshold_current_adaptation.xlsx`\
IV curves are stored in: `RGC_IV_curves_below_threshold_adaptation.npz`

### Compensation of axial current attenuation

* `fig10.py`: we compare the attenuation of the peak axial current with that of the transferred charge and current duration.

Panel A-B-C-D-E: axial current, charge and current duration adaptation \
Measure of the peak axial current, charge and current duration during adaptation: `adaptation_analysis_with_PN.py` for the cells with P/5 protocol and `adaptation_analysis_no_PN.py` for the other cells (before 20190611)\
Results are stored in: `RGC_adaptation.xlsx`

## Methods

* `fig11.py`: Theoretical estimation of axial current at threshold compared to simulations

Neuron model: `shared/models/model_spike_initiation.py`\
Run the model: `model_SI_protocol_VC_dichotomy.py` and `model_SI_protocol_VC_test_pulse.py`\
Results stored in: `simulations data/fig11/`


