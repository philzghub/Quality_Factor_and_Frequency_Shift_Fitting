# CPW Resonator Fitting (Qi & Frequency Shift)

This repo contains:
- **Equations** used in the models (`docs/README_equations.md`)
- **Example graphs** (`examples/`)
- **Code**:
  - `qi7params.py` — Qi model (7 free params) with linear & semilog plotting
  - `freq_shift.py` — Δf/f₀ model (off-resonant TLS)
  - `single_curve.py` — single-trace circle fit → Qi, Qe, f₀
  - `loss_channels.py` — composite loss channel model

## Quick start
```bash
pip install -e .
python examples/make_example_plots.py


## docs/README_equations.md (put the math here first; you can refine later)
```markdown
# Equations used

## Dissipative (near-resonant TLS) contribution to 1/Q
\[
\frac{1}{Q_\mathrm{TLS}} = F\,p\,\delta_0 \, s(P, T)
\]
with a saturation function \(s(P,T)\) decreasing with field/power and temperature.

## Dispersive (off-resonant TLS) frequency shift
\[
\frac{\Delta f}{f_0} = \mathcal{D}_\mathrm{TLS}(T) \quad (\text{weak power dependence})
\]

## Composite internal loss
\[
\frac{1}{Q_\mathrm{int}} = \frac{1}{Q_\mathrm{TLS}} + \frac{1}{Q_\mathrm{qp}} + \frac{1}{Q_\mathrm{res}}
\]

*(Add exact functional forms + parameter definitions as you finalize. Keep symbols consistent with your thesis.)*
