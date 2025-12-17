import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import kv as K_v
import os
import json

# --- Constants ---
h = 6.62607015e-34
hbar = h / (2 * np.pi)
kB = 1.380649e-23
f_res = 5.60157e9  # Adjust to your resonator frequency (Hz)
NBAR_CONST = 1.35032  # fixed nbar (NOT taken from JSON) for 160 dBm input power, ADJUST AS NEEDED!!!

# --- Paths (repo-relative) ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
CSV_PATH = os.path.join(REPO_ROOT, "examples", "csv_data", "qint", "160dBm.csv")
PLOT_DIR = os.path.join(REPO_ROOT, "examples", "plots", "qint")
JSON_PATH = os.path.join(PLOT_DIR, "qint_singlecurve.json")
os.makedirs(PLOT_DIR, exist_ok=True)

# --- Data Loading ---
def load_data(filename):
    df = pd.read_csv(filename, delimiter=';')
    df.columns = df.columns.str.strip()
    # Temperature column: prefer names containing 'temp', fallback to exact German header
    t_candidates = [c for c in df.columns if 'temp' in c.lower()]
    t_col = t_candidates[0] if t_candidates else ('Temperatur (mK)' if 'Temperatur (mK)' in df.columns else None)
    if t_col is None:
        raise KeyError("No temperature column found. Available columns: " + str(df.columns.tolist()))
    # Qint column
    q_candidates = [c for c in df.columns if 'qint' in c.lower()]
    if not q_candidates:
        raise KeyError("No Qint column found. Available columns: " + str(df.columns.tolist()))
    q_col = q_candidates[0]
    # Normalize decimal comma and types
    T_data = df[t_col].astype(str).str.replace(',', '.').astype(float).to_numpy() / 1000.0
    Qint_data = df[q_col].astype(str).str.replace(',', '.').astype(float).to_numpy()
    return T_data, Qint_data

# --- Model Functions ---
def Q_TLS(T, Q_TLS_0, beta1, beta2, D):
    omega = 2 * np.pi * f_res
    arg = hbar * omega / (2 * kB * T)
    numerator = np.sqrt(1 + (NBAR_CONST ** beta2 / (D * T**beta1)) * np.tanh(arg))
    denominator = np.tanh(arg)
    return Q_TLS_0 * (numerator / denominator)

def Q_QP(T, A_QP, Tc):
    omega = 2 * np.pi * f_res
    Delta = 1.76 * kB * Tc
    arg = hbar * omega / (2 * kB * T)
    denominator = np.sinh(arg) * K_v(0, arg)
    numerator = np.exp(Delta / (kB * T))
    return A_QP * numerator / denominator

def Q_total(T, Q_TLS_0, beta1, beta2, D, A_QP, Q_other, Tc):
    return 1 / (
        1 / Q_TLS(T, Q_TLS_0, beta1, beta2, D) +
        1 / Q_QP(T, A_QP, Tc) +
        1 / Q_other
    )

# --- Fitting Functions ---
def fit_best_params(T, Qint, initial_guess, bounds, attempts=1000, random_factor=0.1):
    best_params = initial_guess
    best_mean_diff = np.inf
    for i in range(attempts):
        p0 = np.random.uniform(bounds[0], bounds[1])
        try:
            popt, _ = curve_fit(Q_total, T, Qint, p0=p0, bounds=bounds, maxfev=10000)
            Qfit = Q_total(T, *popt)
            mean_diff = np.mean(np.abs(Qint - Qfit))
            print(f"Attempt {i+1}/{attempts}: Mean difference = {mean_diff:.3e}")
            if mean_diff < best_mean_diff:
                best_mean_diff = mean_diff
                best_params = popt
        except Exception as e:
            print(f"Attempt {i+1}/{attempts}: Fit failed ({e})")
            continue
    print("\nBest parameters for this round:")
    param_names = [
        "Q_TLS_0", "beta1", "beta2", "D", "A_QP", "Q_other", "Tc"
    ]
    for name, val in zip(param_names, best_params):
        print(f"  {name}: {val:.6g}")
    print(f"  Mean difference: {best_mean_diff:.3e}\n")
    return best_params, best_mean_diff

def iterative_fit(T, Qint, initial_guess, bounds, attempts=1000, random_factor=0.1):
    best_params = np.array(initial_guess)
    best_mean_diff = np.mean(np.abs(Qint - Q_total(T, *best_params)))
    print(f"Starting parameters mean difference: {best_mean_diff:.3e}")
    round_counter = 1
    while True:
        params, mean_diff = fit_best_params(T, Qint, best_params, bounds, attempts, random_factor)
        print(f"Round {round_counter}: Mean difference = {mean_diff:.3e}")
        if mean_diff < best_mean_diff - 1e-3:
            best_params = params
            best_mean_diff = mean_diff
            round_counter += 1
        else:
            print("Parameters did not improve. Stopping and saving previous best parameters.")
            return best_params, best_mean_diff

# Load data
if not os.path.exists(CSV_PATH):
    raise FileNotFoundError(f"CSV not found: {CSV_PATH}")
T_data, Qint_data = load_data(CSV_PATH)

# Run iterative fit
# Initialize from previous JSON if available
if os.path.exists(JSON_PATH):
    with open(JSON_PATH, "r") as f:
        loaded = json.load(f).get("params", [])
        initial_guess = np.array(loaded[:7]) if len(loaded) >= 7 else np.array([1e6, 0.1, 0.1, 500, 500, 1e7, 5.5])
else:
    initial_guess = np.array([1e6, 0.1, 0.1, 500, 500, 1e7, 5.5])

bounds = (initial_guess * 0.1, initial_guess * 10.0)

best_params, best_mean_diff = iterative_fit(
    T_data, Qint_data, initial_guess, bounds, attempts=1000, random_factor=0.1
)
print("Best parameters found:", best_params)
print("Best mean difference:", best_mean_diff)

# Save best parameters to JSON
with open(JSON_PATH, "w") as f:
    json.dump({"params": best_params.tolist(), "mean_diff": float(best_mean_diff)}, f, indent=2)
print(f"Best parameters saved to: {JSON_PATH}")

# Plot results for best parameters only
T_fit = np.linspace(float(np.min(T_data)), float(np.max(T_data)), 300)
Q_fit_curve = Q_total(T_fit, *best_params)
plt.figure(figsize=(8, 6))
plt.plot(T_data * 1000, Qint_data, 'o', label='Data')
plt.plot(T_fit * 1000, Q_fit_curve, '-', label='Best Fit', color='red')
plt.title("Iterative Optimized Fit")
plt.xlabel("Temperature (mK)")
plt.ylabel("Qint")
plt.legend()
plt.grid(True)
plt.tight_layout()
save_path = os.path.join(PLOT_DIR, "qint_singlecurve.png")
plt.savefig(save_path, dpi=300)
plt.close()
print(f"Plot saved to: {save_path}")