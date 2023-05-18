This repo contains the plotting scripts used to make the figures in the nonlinear spectroscopy paper:

### Requirements
These scripts were gnerated with 
* Python 3.10.6
* Numpy 1.23.2
* Matplotlib 3.6.3

### Instructions
To generate the figures for the paper, run the following scripts:
* plot_2DES_models.py
* plot_2DES_solvent_1.py
* plot_2DES_solvent_2.py
* plot_SS_vs_time.py
* plot_transient_abs.py

The resulting figures will be stored in the `png` directory

### Customization
The style used for all figures is controlled by `style.mplstyle`, but some scripts use internal rcParams for that particular plot.

Global settings, including the root directory where all data is loaded from, is controlled by `global_settings.py`, which is imported by all scripts.
