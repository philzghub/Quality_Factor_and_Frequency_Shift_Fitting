import os, sys, json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.special import digamma
from matplotlib.ticker import MultipleLocator

# ---- Styling ----
plt.rcParams.update({
    "font.size": 17, 
    "axes.labelsize": 21, 
    "legend.fontsize": 17,
    "xtick.labelsize": 17, 
    "ytick.labelsize": 17
})

# ---- Constants ----
kB = 1.380649e-23
h  = 6.62607015e-34
pi = np.pi
f0 = 5619901372.0  # Adjust to your resonator frequency here (Hz)
C  = h * f0 / (2 * pi * kB)  # convenience: h f0 / (2π kB)

# ---- Model ----
def tls_shift(T, Q0_tls):
    x = C / np.maximum(T, 1e-12)
    term = np.real(digamma(0.5 + 1j * x)) - np.log(C / np.maximum(T, 1e-12))
    return term / (Q0_tls * pi)

def qp_shift(T, Delta_K, alpha):
    y = Delta_K / np.maximum(T, 1e-12)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        return -(alpha / 2.0) * (y / np.sinh(y))

def total_shift(T, Q0_tls, Delta_K, alpha):
    return tls_shift(T, Q0_tls) + qp_shift(T, Delta_K, alpha)

# ---- Data I/O ----
def load_data(path):
    df = pd.read_csv(path, delimiter=';')
    df.columns = [c.strip() for c in df.columns]
    t_col = next((c for c in df.columns if 'temp'  in c.lower()), None)
    s_col = next((c for c in df.columns if 'shift' in c.lower()), None)
    if t_col is None or s_col is None:
        raise RuntimeError("Need 'Temperature' and 'Frequency Shift' columns.")
    TmK  = df[t_col].astype(str).str.replace(',', '.').astype(float).to_numpy()
    frac = df[s_col].astype(str).str.replace(',', '.').astype(float).to_numpy()
    return TmK / 1000.0, frac  # K, dimensionless

# ---- Fit utils ----
def mean_abs_diff(y_true, y_pred):
    return np.mean(np.abs(y_true - y_pred))

def _logu(lo, hi):
    return 10 ** np.random.uniform(np.log10(lo), np.log10(hi))

def random_params(bounds):
    (Qmin, Dmin, Amin), (Qmax, Dmax, Amax) = bounds
    return (_logu(Qmin, Qmax), _logu(Dmin, Dmax), _logu(Amin, Amax))

def vary_params(p):
    return tuple(np.array(p) * (10 ** np.random.uniform(-1, 1, 3)))

def load_previous(json_path):
    if not os.path.exists(json_path):
        return None, np.inf
    try:
        with open(json_path, "r") as f:
            saved = json.load(f)
        sp = saved.get("params")
        sl = saved.get("loss", np.inf)
        if sp and np.isfinite(sl):
            print(f"Loaded previous best: loss={sl:.3e} params={tuple(sp)}")
            return tuple(sp), float(sl)
    except Exception:
        pass
    return None, np.inf

def save_best(json_path, p, loss, attempt, threshold):
    with open(json_path, "w") as f:
        json.dump({
            "params": list(p),
            "loss": float(loss),
            "attempt": int(attempt),
            "param_names": ["Q0_TLS", "Delta_over_kB(K)", "alpha"],
            "threshold": float(threshold)
        }, f, indent=2)

# ---- Core fit (stochastic search + restarts) ----
def fit(T, y, bounds, json_path, max_attempts=200000, threshold=1e-8, restart_every=500):
    best_p, best_loss = load_previous(json_path)

    (Qmin, Dmin, Amin), (Qmax, Dmax, Amax) = bounds
    for attempt in range(1, max_attempts + 1):
        cand = random_params(bounds) if (best_p is None or attempt % restart_every == 0) else vary_params(best_p)
        Q0, Dk, a = cand
        Q0 = np.clip(Q0, Qmin, Qmax)
        Dk = np.clip(Dk, Dmin, Dmax)
        a  = np.clip(a,  Amin, Amax)
        cand = (Q0, Dk, a)

        loss = mean_abs_diff(y, total_shift(T, *cand))
        improved = loss < best_loss
        print(f"Attempt {attempt}/{max_attempts} loss={loss:.6e} best={best_loss:.6e}{' *improved*' if improved else ''}")
        if improved:
            best_p, best_loss = cand, loss
            save_best(json_path, best_p, best_loss, attempt, threshold)
        if best_loss <= threshold:
            print(f"Threshold reached: {best_loss:.3e} <= {threshold:.3e}")
            break

    return best_p, best_loss

# ---- Plotting ----
def _axes_common(ax):
    ax.set_xlim(0, 1700)
    ax.xaxis.set_major_locator(MultipleLocator(400))
    ax.xaxis.set_minor_locator(MultipleLocator(200))
    ax.grid(True, which='major', alpha=0.35)
    ax.grid(True, which='minor', alpha=0.20, linestyle='-')
    ax.set_ylim(-4.5, 1.2)                   # fixed requested bounds (ppm)
    ax.ticklabel_format(style='plain', axis='y')

def plot_total(T, frac_ppm, T_plot, fit_ppm, out_path):
    data_blue = "#2F6C9E"
    fit_blue  = "#8FC2E8"
    plt.figure(figsize=(7, 5))
    plt.plot(T * 1000, frac_ppm, 'o', label="Data",
             color=data_blue, markeredgecolor=data_blue, zorder=4)
    plt.plot(T_plot * 1000, fit_ppm, '-', label="Total Fit",
             color=fit_blue, linewidth=2.2, zorder=5)
    plt.xlabel("Temperature (mK)")
    plt.ylabel("δf/f0 (ppm)")
    ax = plt.gca()
    _axes_common(ax)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Saved plot: {out_path}")

def plot_components(T, frac_ppm, T_plot, tls_ppm, qp_ppm, fit_ppm, out_path):
    data_blue = "#2F6C9E"
    fit_blue  = "#8FC2E8"
    plt.figure(figsize=(7, 5))
    plt.plot(T_plot * 1000, fit_ppm, '-',  label="Total",
             color=fit_blue, linewidth=2.2, zorder=5)
    plt.plot(T_plot * 1000, tls_ppm, '--', label="TLS term",
             color="green", linewidth=2.2, zorder=2)
    plt.plot(T_plot * 1000, qp_ppm,  '--', label="QP term",
             color="orange", linewidth=2.2, zorder=3)
    plt.plot(T * 1000, frac_ppm, 'o', label="Data",
             color=data_blue, markeredgecolor=data_blue, zorder=4)
    plt.xlabel("Temperature (mK)")
    plt.ylabel("δf/f0 (ppm)")
    ax = plt.gca()
    _axes_common(ax)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Saved components plot with data points: {out_path}")

# ---- Run ----
def main():
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
    data_path  = os.path.join(REPO_ROOT, "examples", "csv_data", "freqshift", "Freq_shift_nbonly.csv")
    plots_dir  = os.path.join(REPO_ROOT, "examples", "plots", "freqshift")
    os.makedirs(plots_dir, exist_ok=True)
    json_path  = os.path.join(plots_dir, "freqshift.json")

    if not os.path.exists(data_path):
        print("Data file missing:", data_path)
        sys.exit(1)

    T, frac = load_data(data_path)

    # Bounds: Q0_TLS ∈ [1e3, 1e7], Delta_K ∈ [6, 20] K, alpha ∈ [1e-12, 1]
    bounds = ((1e3, 6, 1e-12), (1e7, 20, 1.0))
    best_p, best_loss = fit(T, frac, bounds, json_path)
    if best_p is None:
        print("No fit found.")
        return

    # Dense plotting grid 0..1700 mK; avoid T=0 singularity
    T_plot = np.linspace(1e-6, 1.7, 1200)

    # Curves and ppm conversion (for readable y-axis)
    tls_curve = tls_shift(T_plot, best_p[0])
    qp_curve  = qp_shift(T_plot, best_p[1], best_p[2])
    fit_curve = tls_curve + qp_curve
    frac_ppm, tls_ppm, qp_ppm, fit_ppm = (arr * 1e6 for arr in (frac, tls_curve, qp_curve, fit_curve))
    plot_total(T, frac_ppm, T_plot, fit_ppm,
                os.path.join(plots_dir, "freqshift.png"))
    plot_components(T, frac_ppm, T_plot, tls_ppm, qp_ppm, fit_ppm,
                        os.path.join(plots_dir, "freqshift_loss_channels.png"))

# Move the entry point guard to module level
if __name__ == "__main__":
    main()
