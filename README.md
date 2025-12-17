# Quality Factor & Frequency Shift Fitting (CPW Resonators)
> Tools to fit internal quality factor $Q_{int}$ and frequency shift $\Delta f_0/f_0$ for CPW resonators, visualize loss channels, and persist best-fit parameters as JSON for iterative runs.

- Works from simple CSVs (examples included)

- Produces publication-ready PNGs

- Saves best parameters as JSON (used as the next run’s initial guess)


[![Build](https://img.shields.io/badge/build-passing-brightgreen)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)]()
[![Docs](https://img.shields.io/badge/docs-available-informational)]()

<img src="docs/demo.gif" alt="Demo" width="600"/>

# Contents (what each file does)

- [cpwfit/__init__.py](cpwfit/__init__.py) — marks cpwfit as a Python package, can expose a small public API (e.g., __version__) so import cpwfit works cleanly
- [cpwfit/models/__init__.py](cpwfit/models/__init__.py) — gathers and re-exports functions/classes from individual model scripts so users can write from cpwfit.models import ... without knowing file paths
- [cpwfit/models/qint_allcurves_semilog.py](cpwfit/models/qint_allcurves_semilog.py) — fits all curves simultaneously with 7 free parameters, minimizing total residuals; weights −160 dBm points 5×; plots on a semilog axis to emphasize low-power behavior; saves PNG and best-params JSON
- [cpwfit/models/qint_allcurves_lin.py](cpwfit/models/qint_allcurves_lin.py) — same as above (same physics and −160 dBm 5× weighting) but plots on a linear axis; saves PNG and best-params JSON
- [cpwfit/models/qint_loss_channels.py](cpwfit/models/qint_loss_channels.py) — takes a fitted parameter set (from JSON) and for a chosen power plots loss-channel composition (TLS, qp, residual); saves PNG
- [cpwfit/models/qint_singlecurve.py](cpwfit/models/qint_singlecurve.py) — fits one curve independently, returning parameters that yield the smallest residuals between the datapoints and model; saves PNG and best-params JSON
- [cpwfit/models/freqshift_loss_channels.py](cpwfit/models/freqshift_loss_channels.py) — fits Δf/f datapoints; outputs two PNGs: one with loss-channel overlays and one without for a clean view; saves best-params JSON
- [examples/csv_data/](examples/csv_data/) — example CSVs that define the expected column formats
- [examples/plots/](examples/plots/) — example PNGs produced by the scripts and the example best-params JSONs

# Equations used
##Quality 

The equations and fitting models used in this work are shown below. If you want to use different fitting equations also adjust the code and parameters accordingly. 

### Total $Q_{int}$ 

$$\frac{1}{Q_\mathrm{int}}=\frac{1}{Q_\mathrm{TLS}}+\frac{1}{Q_\mathrm{QP}}+\frac{1}{Q_\mathrm{other}}$$

### TLS term

$$Q_{\mathrm{TLS}}(\bar{n},T)
= Q_{\mathrm{TLS},0}
  \frac{\sqrt{1+
    \left(\dfrac{\bar{n}^{\beta_{2}}}{DT^{\beta_{1}}}\right)
    \tanh\left(\frac{\hbar\omega}{2 k_{\mathrm{B}} T}\right)}}
       {\tanh\left(\frac{\hbar\omega}{2 k_{\mathrm{B}} T}\right)}$$
       
K. D. Crowley et al., “Disentangling Losses in Tantalum Superconducting Circuits,” Physical Review X 13, 041005 (2023).

### QP term

$$Q_{\mathrm{QP}}(T)
= A_{\mathrm{QP}}
  \frac{e^{\Delta_{0}/(k_{\mathrm B} T)}}
       {\sinh\left(\frac{\hbar\omega}{2 k_{\mathrm B} T}\right)
        K_{0}\left(\frac{\hbar\omega}{2 k_{\mathrm B} T}\right)}$$

K. D. Crowley et al., “Disentangling Losses in Tantalum Superconducting Circuits,” Physical Review X 13, 041005 (2023).
        
### TLS-induced frequency shift
$$
\left(\frac{\delta f(T)}{f_0}\right)_{\mathrm{TLS}}
= \frac{1}{\pi Q_{\mathrm{TLS},0}}
\{Re}\left[
  \psi\left(\tfrac{1}{2}+i\frac{\hbar\omega}{2\pi k_{\mathrm B}T}\right)
-\ln\left(\frac{\hbar\omega}{2\pi k_{\mathrm B}T}\right)
\right]
$$

J. Gao, “The Physics of Superconducting Microwave Resonators,” PhD thesis, Caltech (2008).

### Quasiparticle-induced shift
$$
\left(\frac{\delta f(T)}{f_0}\right)_{\mathrm{QP}} = -\frac{1}{2}\frac{\Delta L}{L} = -\tfrac{1}{2}\alpha f_{0} \frac{\Delta L_{l,kin}}{L_{l,kin}}
$$

M. Scigliuzzo et al., “Phononic Loss in Superconducting Resonators on Piezoelectric Substrates,” New Journal of Physics 22, 053027 (2020).

### Total shift
$$
\frac{\delta f(T)}{f_0} =
\left(\frac{\delta f(T)}{f_0}\right)_{\mathrm{TLS}} +
\left(\frac{\delta f(T)}{f_0}\right)_{\mathrm{QP}}
$$

 ### **Note on methods**
> The $\boldsymbol{\Delta f_0/f_0}$ model and the $Q_{\mathrm{int}}$ model currently use different fitting strategies (**stochastic mean-absolute error search** vs. **bounded non-linear least squares on relative error**). This will likely be unified/adjusted in future versions.


# Quickstart

## 1) Installation (one time)

Install the virtual environment to isolate this project’s Python packages from your system so dependencies don’t conflict and installs stay reproducible.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
```
## 2) Check correct data formats (CSV)

Import your own CSV files for the Qint and Δf/f. Match the example files for column names, delimiters, and layout.

Example: 

```bash
Temperature (mK);Power (dBm);Qint
100;-140;1.64e6
200;-140;1.53e6
300;-140;2.44e6
400;-140;2.71e6
500;-140;2.92e6
600;-140;3.52e6
700;-140;3.95e6
800;-140;4.31e6
900;-140;4.42e6
1000;-140;4.38e6
1100;-140;4.28e6
1200;-140;3.464415e6
1300;-140;3.042468e6
1400;-140;2.657917e6
1500;-140;1.974064e6
1600;-140;1.306835e6
1700;-140;861599
```
## Usage

### 1) Fit all Qint curves (semilog)
```bash
python -m cpwfit.models.qint_allcurves_semilog \
  --csv examples/csv_data/qint_all.csv \
  --out examples/plots/qint_semilog.png \
  --save-best examples/plots/qint_best.json
```
### 2) Fit a single curve
```bash
python -m cpwfit.models.qint_singlecurve \
  --csv examples/csv_data/qint_single.csv \
  --out examples/plots/qint_single.png \
  --save-best examples/plots/qint_single_best.json
```
### 3) Fit Frequency-shift
```bash
python -m cpwfit.models.freqshift_loss_channels \
  --csv examples/csv_data/freqshift.csv \
  --out examples/plots/freqshift.png \
  --save-best examples/plots/freqshift_best.json
```
### 4) Plot loss-channel composition (from saved params)
```bash
python -m cpwfit.models.qint_loss_channels \
  --best examples/plots/qint_best.json \
  --power -160 \
  --out examples/plots/qint_channels.png
```

## Configuration knobs for your adjustment 
Adjust these to your own devices and needs:
- Resonance frequency $f_0$: set to your device.
- $\bar n$: adjust photon-number according to your own data/calculation (if you don't have the value, remove as a constant and add as parameter in the script)
- Bounds & initial guesses: Adjust accordingly.
- Dataset weights: e.g., weight −160 dBm ×5 for TLS sensitivity.






