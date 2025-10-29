# filepath: /Users/philip/fit_all.py/Nbonly/Nbonly_SinglePhoton/losschanels_160dBm_fixedparams.py
import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from scipy.special import kv as K_v

# --- Plot style ---
plt.rcParams.update({
    "font.size": 17,
    "axes.labelsize": 21,
    "legend.fontsize": 17,
    "xtick.labelsize": 17,
    "ytick.labelsize": 17,
})

# --- Physical constants ---
h = 6.62607015e-34
hbar = h / (2 * np.pi)
kB = 1.380649e-23
f_res = 5.61989e9
NBAR_CONST = 1.35032  # fixed nbar (NOT taken from JSON)

# --- Models (using NBAR_CONST) ---
def Q_TLS(T, Q_TLS_0, beta1, beta2, D):
    omega = 2 * np.pi * f_res
    arg = hbar * omega / (2 * kB * T)
    tanh_arg = np.tanh(arg)
    core = 1.0 + (NBAR_CONST * beta2 / (D * T**beta1)) * tanh_arg
    return Q_TLS_0 * (np.sqrt(core) / tanh_arg)

def Q_QP(T, A_QP, Tc):
    omega = 2 * np.pi * f_res
    Delta = 1.76 * kB * Tc
    arg = hbar * omega / (2 * kB * T)
    denom = np.sinh(arg) * K_v(0, arg)
    numer = np.exp(Delta / (kB * T))
    return A_QP * numer / denom

def Q_total(T, Q_TLS_0, beta1, beta2, D, A_QP, Q_other, Tc):
    tls = Q_TLS(T, Q_TLS_0, beta1, beta2, D)
    qp  = Q_QP(T, A_QP, Tc)
    return 1.0 / (1.0 / tls + 1.0 / qp + 1.0 / Q_other)

# --- I/O helpers ---
def load_data(csv_path):
    df = pd.read_csv(csv_path, delimiter=";")
    df.columns = df.columns.str.strip()
    tcol = next((c for c in df.columns if "temp" in c.lower()), None)
    qcol = next((c for c in df.columns if "qint" in c.lower()), None)
    if tcol is None:
        raise KeyError("No temperature column found (expects 'temp' in name).")
    if qcol is None:
        raise KeyError("No Qint column found (expects 'Qint' in name).")
    df[qcol] = df[qcol].astype(str).str.replace(",", ".").astype(float)
    T = df[tcol].to_numpy() / 1000.0  # mK -> K
    Q = df[qcol].to_numpy()
    return T, Q

def load_fixed_params(fixed_json_path, expected, required_keys):
    with open(fixed_json_path, "r") as f:
        fixed = json.load(f)
    if "params" not in fixed:
        raise KeyError("Key 'params' missing in fixed params JSON.")
    p = fixed["params"]
    missing = [k for k in required_keys if k not in p]
    if missing:
        raise KeyError(f"Missing parameter keys in JSON: {missing}")
    # Use JSON if present; override to EXPECTED if any difference
    params = {k: float(p.get(k, expected[k])) for k in required_keys}
    mismatch = [k for k in required_keys if not np.isclose(params[k], expected[k], rtol=0, atol=0)]
    if mismatch:
        print(f"Overriding to provided weighted params for: {', '.join(mismatch)}")
        params.update(expected)
    return params

def crosspoint_upper_y(Q_tls_curve, Q_qp_curve, margin=1.30):
    # Find TLS–QP crossing (or closest approach), compute y-upper (×1e7 units) with margin
    diff = Q_tls_curve - Q_qp_curve
    sign = np.sign(diff)
    idx = np.where(sign[:-1] * sign[1:] <= 0)[0]
    if len(idx) > 0 and np.isfinite(diff[idx[0]]) and np.isfinite(diff[idx[0] + 1]):
        i = idx[0]
        return float(margin * (0.5 * (Q_tls_curve[i] + Q_qp_curve[i]) / 1e7))
    # Fallback: closest approach on grid
    i_min = int(np.argmin(np.abs(diff)))
    return float(margin * (0.5 * (Q_tls_curve[i_min] + Q_qp_curve[i_min]) / 1e7))

def set_y_and_report(ax, ymin, ymax):
    ax.set_ylim(ymin, ymax)
    lo, hi = ax.get_ylim()
    print(f"Y-limits (×1e7 units): {lo:.3g} — {hi:.3g}")

# --- Paths ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))

# CSV source
csv_path = os.path.join(REPO_ROOT, "examples", "csv_data", "qint", "160dBm.csv")

# Fixed-params JSON: try linear first, then semilog, then legacy name
plots_qint_dir = os.path.join(REPO_ROOT, "examples", "plots", "qint")
json_candidates = [
    os.path.join(plots_qint_dir, "qint_lin.json"),
    os.path.join(plots_qint_dir, "qint_log.json"),
    os.path.join(plots_qint_dir, "fixed_shared_params_results_weighted.json"),
]
fixed_json_path = next((p for p in json_candidates if os.path.exists(p)), None)

if not os.path.exists(csv_path):
    raise FileNotFoundError(f"CSV not found: {csv_path}")
if fixed_json_path is None:
    raise FileNotFoundError(
        "Fixed-params JSON not found. Looked for: " + ", ".join(json_candidates)
    )

# --- Load data and params ---
T_data, Qint_data = load_data(csv_path)

REQUIRED = ["Q_TLS_0", "beta1", "beta2", "D", "Q_other", "A_QP", "Tc"]
EXPECTED = {
    "Q_TLS_0": 702034.5791713039,
    "beta1":   0.5633765494406937,
    "beta2":   0.31438135746391943,
    "D":       1.0771171131840873,
    "Q_other": 42738370.852098696,
    "A_QP":    450.4365833076074,
    "Tc":      6.0333925474881,
}
params = load_fixed_params(fixed_json_path, EXPECTED, REQUIRED)

Q_TLS_0 = params["Q_TLS_0"]; beta1 = params["beta1"]; beta2 = params["beta2"]
D = params["D"]; Q_other = params["Q_other"]; A_QP = params["A_QP"]; Tc = params["Tc"]

print(f"Using weighted params: Q_TLS_0={Q_TLS_0:.6g}, beta1={beta1:.6g}, "
      f"beta2={beta2:.6g}, D={D:.6g}, Q_other={Q_other:.6g}, "
      f"A_QP={A_QP:.6g}, Tc={Tc:.6g}, nbar={NBAR_CONST:.6g}")

# --- Compute curves (nbar constant) ---
T_fit = np.linspace(T_data.min(), T_data.max(), 400)
Q_total_fit = Q_total(T_fit, Q_TLS_0, beta1, beta2, D, A_QP, Q_other, Tc)
Q_tls_curve = Q_TLS(T_fit, Q_TLS_0, beta1, beta2, D)
Q_qp_curve  = Q_QP(T_fit, A_QP, Tc)
Q_other_curve = np.full_like(T_fit, Q_other)

# Pre-scale for plotting (×1e7 units)
y_data  = Qint_data   / 1e7
y_total = Q_total_fit / 1e7
y_tls   = Q_tls_curve / 1e7
y_qp    = Q_qp_curve  / 1e7
y_other = Q_other_curve / 1e7

# Determine y-upper from TLS–QP crosspoint, add 30% headroom
upper = max(1e-9, crosspoint_upper_y(Q_tls_curve, Q_qp_curve, margin=1.30))

# --- Plot ---
fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(T_data * 1000, y_data,  'o', label="Data",  color="#2F6C9E")
ax.plot(T_fit  * 1000, y_total, '-', label="Fit",   color="#8FC2E8", linewidth=2.2)
ax.plot(T_fit  * 1000, y_tls,   '--', label="$Q_{TLS}$",   color="green",  linewidth=2.2)
ax.plot(T_fit  * 1000, y_qp,    '-.', label="$Q_{QP}$",    color="orange", linewidth=2.2)
ax.plot(T_fit  * 1000, y_other, ':',  label="$Q_{other}$", color="grey",   linewidth=2.2)

ax.set_xlabel("Temperature (mK)")
ax.set_ylabel(r"$Q_{\mathrm{int}}$ ($\times 10^7$)")
set_y_and_report(ax, 0, upper)

ax.set_xlim(0, 1700)
ax.xaxis.set_major_locator(MultipleLocator(400))
ax.xaxis.set_minor_locator(MultipleLocator(200))
ax.grid(True, which="major", alpha=0.35)
ax.grid(True, which="minor", alpha=0.30, linestyle="-")
fig.tight_layout()
ax.legend()
plt.grid(alpha=0.3)

# --- Save ---
out_dir = os.path.join(REPO_ROOT, "examples", "plots", "qint")
os.makedirs(out_dir, exist_ok=True)
final_path = os.path.join(out_dir, "160dBm_losschannels.png")
tmp_path = os.path.join(out_dir, "160dBm_losschannels.tmp.png")

plt.savefig(tmp_path, dpi=300, format="png")
plt.close(fig)
os.replace(tmp_path, final_path)
print(f"Saved: {final_path}")