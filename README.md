# AIS geometry and axial current

### Requirements

In addition to standard scientific packages (scipy, matplotlib), this code requires:
* Seaborn
* [Brian 2](http://briansimulator.org)

Electrophysiological data from recorded in retina ganglion cells
(DOI: ) must be
[downloaded]() into the folder
`shared/data/`.

## The axial current at spike initiation

* `fig1.py`: spontaneous APs of retinal ganglion cells.
We measure the spike onset, axonal and somatic max of dV/dt and somatic regeneration across the cell population.
List of cells used for the population analysis (panels E and F): `RGC_electrical_properties.xlsx`
Analysis of the raw recordings: `action_potential_analysis.py`
Results of the analysis used to plot the figure: `RGC_action_potential.xlsx`

* `fig2.py`: recording the axial current.
We illustrate how we measure the axial currents, how we correct the peak axial current for the effect of the series resistance and how the peak axial current and the voltage threshold depend on the series resistance (in a model and in RGC).
Model (panel E and F): 
voltage clamp protocol to record axial currents for different series resistance: `model_AP_protocol_VC_dichotomy_with_series_resistance.py`
voltage clamp protocol to record a test pulse: `model_AP_protocol_VC_test_pulse_with_series_resistance.py`
simulation results are stored in `simulations data/fig2/`
RGC data (panels G and H):
List of cells used for the analysis: `RGC_electrical_properties.xlsx`
Measure of the peak axonal current and voltage threshold from axial current recordings: `axonal_current_analysis.py`

## The threshold axial current

In this part, we show ...

## Adaptation of the axial current

