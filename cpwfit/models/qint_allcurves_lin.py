import os
import json
import numpy as np
import pandas as pd
from scipy.special import kv as K_v
from scipy.optimize import least_squares
import matplotlib.pyplot as plt

# -------- Plot style --------
plt.rcParams.update({
    "font.size": 24,
    "axes.labelsize": 28,
    "axes.titlesize": 28,
    "xtick.labelsize": 22,
    "ytick.labelsize": 22,
    "legend.fontsize": 22,
    "legend.title_fontsize": 24,
})

# -------- Config --------
POWERS = ["80", "100", "120", "140", "160"]
# Base paths resolved relative to the repository root
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
CSV_PATTERN = "{p}dBm.csv"

# Directories per user request
CSV_DIR = os.path.join(REPO_ROOT, "examples", "csv_data", "qint")
OUTPUT_DIR = os.path.join(REPO_ROOT, "examples", "plots", "qint")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Weight 160 dBm 5×
DATASET_WEIGHTS = {"160": 5.0}

# Per-power nbar
NBAR_BY_POWER = {
    "80": 1.45466e8, "100": 1.43609183602e6, "120": 1.413910552e4,
    "140": 138.34054, "160": 1.35032,
}

# Bounds and initial guesses (shared across powers)
BOUNDS = {
    "Q_TLS_0": (1e5, 5e8), "beta1": (1e-3, 25), "beta2": (1e-5, 2e3),
    "D": (1.0, 5e4), "Q_other": (1e4, 5e9), "A_QP": (1e-8, 1e9), "Tc": (0.5, 12.0),
}
INIT = {
    "Q_TLS_0": 6.368547e5, "beta1": 0.25, "beta2": 0.2, "D": 800.0,
    "Q_other": 5e7, "A_QP": 1e3, "Tc": 7.0,
}

# Iteration controls
MAX_ITERS = 100
OUTER_MAX_ROUNDS = 50
IMPROVEMENT_TOL = 0.0
EXP_CLIP = 700

# -------- Physics --------
kB = 1.380649e-23
hbar = 6.62607015e-34 / (2 * np.pi)
f_res = 5.61989e9

# -------- Models --------
def Q_TLS(T, Q_TLS_0, beta1, beta2, D, nbar):
    T = np.maximum(T, 1e-6)
    arg = hbar * 2*np.pi*f_res / (2 * kB * T)
    tanh_arg = np.tanh(arg)
    core = 1.0 + (nbar**beta2 / (D * T**beta1)) * tanh_arg
    return Q_TLS_0 * (np.sqrt(np.maximum(core, 1e-30)) / tanh_arg)

def Q_QP(T, A_QP, Tc):
    T = np.maximum(T, 1e-6)
    arg = hbar * 2*np.pi*f_res / (2 * kB * T)
    denom = np.sinh(arg) * K_v(0, arg)
    denom = np.where(denom == 0, 1e-300, denom)
    expo = np.clip(1.76 * kB * Tc / (kB * T), -EXP_CLIP, EXP_CLIP)
    return A_QP * np.exp(expo) / denom

def Q_total_shared(T, nbar, params):
    Q_TLS_0, beta1, beta2, D, Q_other, A_QP, Tc = params
    tls = Q_TLS(T, Q_TLS_0, beta1, beta2, D, nbar)
    qp = Q_QP(T, A_QP, Tc)
    return 1.0 / (1.0 / tls + 1.0 / qp + 1.0 / Q_other)

# -------- Data I/O --------
def load_dataset(power):
    csv_path = os.path.join(CSV_DIR, CSV_PATTERN.format(p=power))
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"Missing CSV for power {power}: {csv_path}")
    df = pd.read_csv(csv_path, delimiter=';')
    df.columns = [c.strip() for c in df.columns]
    tcol = next((c for c in df.columns if 'temp' in c.lower()), None)
    qcol = next((c for c in df.columns if 'qint' in c.lower()), None)
    if tcol is None or qcol is None:
        raise RuntimeError(f"Missing temperature/Qint columns in {csv_path}")
    T = df[tcol].astype(str).str.replace(',', '.').astype(float).to_numpy() / 1000.0
    Q = df[qcol].astype(str).str.replace(',', '.').astype(float).to_numpy()
    mask = np.isfinite(T) & np.isfinite(Q) & (T > 0)
    return T[mask], Q[mask]

# -------- Helpers --------
PARAM_ORDER = ["Q_TLS_0","beta1","beta2","D","Q_other","A_QP","Tc"]
LOWER = np.array([BOUNDS[k][0] for k in PARAM_ORDER])
UPPER = np.array([BOUNDS[k][1] for k in PARAM_ORDER])
INIT_VEC = np.array([INIT[k] for k in PARAM_ORDER], dtype=float)

def total_loss(all_data, params):
    # Weighted mean of mean relative absolute errors
    accum = 0.0
    sum_w = 0.0
    for d in all_data:
        pred = Q_total_shared(d["T"], d["nbar"], params)
        err = np.mean(np.abs(pred - d["Q"]) / np.maximum(d["Q"], 1e-12))
        w = DATASET_WEIGHTS.get(d["power"], 1.0)
        accum += w * err
        sum_w += w
    return accum / max(sum_w, 1e-12)

def residual_vector(all_data, params):
    res = []
    for d in all_data:
        r = (Q_total_shared(d["T"], d["nbar"], params) - d["Q"]) / np.maximum(d["Q"], 1e-12)
        w = DATASET_WEIGHTS.get(d["power"], 1.0)
        res.append(np.sqrt(w) * r if w != 1.0 else r)
    return np.concatenate(res)

def refine_params(all_data, p_start):
    res = least_squares(lambda p: residual_vector(all_data, p),
                        p_start, bounds=(LOWER, UPPER), max_nfev=8000, verbose=0)
    return np.clip(res.x, LOWER, UPPER)

# -------- Optimization --------
def optimize_shared(start_params=None, max_iters=MAX_ITERS):
    # Load datasets
    all_data = []
    for p in POWERS:
        if p not in NBAR_BY_POWER:
            continue
        try:
            T, Q = load_dataset(p)
        except FileNotFoundError as e:
            print(e); continue
        all_data.append({"power": p, "T": T, "Q": Q, "nbar": NBAR_BY_POWER[p]})
    if not all_data:
        raise RuntimeError("No datasets loaded.")

    params = INIT_VEC.copy() if start_params is None else np.array(start_params, dtype=float)
    best_params, best_loss = params.copy(), np.inf

    JITTER0, MIN_JITTER = 0.35, 0.02
    for it in range(1, max_iters + 1):
        jitter_scale = max(MIN_JITTER, JITTER0 * (1 - it / max_iters))

        params = refine_params(all_data, params)
        loss = total_loss(all_data, params)

        # Explore with jitter then refine
        trial = params * (1 + jitter_scale * np.random.randn(params.size))
        trial = np.clip(trial, LOWER, UPPER)
        trial = refine_params(all_data, trial)
        loss_trial = total_loss(all_data, trial)
        if loss_trial < loss:
            params, loss = trial, loss_trial

        if loss < best_loss:
            best_params, best_loss = params.copy(), loss

        print(f"[Iter {it:03d}/{max_iters}] loss={loss:.6e} best={best_loss:.6e} jit={jitter_scale:.3f}")

    print(f"Finished. Best loss {best_loss:.6e}")
    return best_params, best_loss, all_data

# -------- Results I/O --------
def save_results(json_path, params, best_loss):
    data = {
        "best_loss_mean_rel_abs": float(best_loss),
        "weights": DATASET_WEIGHTS,
        "params": {k: float(v) for k, v in zip(PARAM_ORDER, params)},
        "param_order": PARAM_ORDER,
    }
    with open(json_path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"Saved params JSON: {json_path}")

def load_previous(json_path):
    if not os.path.exists(json_path):
        return None
    with open(json_path, "r") as f:
        d = json.load(f)
    p = d.get("params", {})
    if "Q_TLS_0" not in p:
        print("Old JSON without Q_TLS_0 detected; ignoring previous results.")
        return None
    vec = np.array([p[k] for k in PARAM_ORDER], dtype=float)
    return vec, d.get("best_loss_mean_rel_abs")

# -------- Plot --------
def plot_shared(all_data, params, best_loss):
    fig, ax = plt.subplots(figsize=(10, 8))
    cmap = plt.get_cmap("Blues", len(all_data) + 2)
    colors = [cmap(i + 2) for i in range(len(all_data))][::-1]

    Y_SCALE = 1e7
    ax.set_xlabel("Temperature (mK)")
    ax.set_ylabel(r"$Q_{\mathrm{int}} (\times 10^{7})$")


    all_scaled = []
    for i, d in enumerate(all_data):
        Tmin, Tmax = d["T"].min(), d["T"].max()
        if Tmax <= Tmin:
            continue
        T_fit = np.linspace(Tmin, Tmax, 600)
        y_dat = d["Q"] / Y_SCALE
        y_fit = Q_total_shared(T_fit, d["nbar"], params) / Y_SCALE
        all_scaled += [y_dat, y_fit]

        # Constant visual style (no visual weighting)
        ms = 6                              # marker size
        lw = 2.0                            # line width
        alpha = 0.65

        ax.plot(d["T"] * 1000, y_dat, 'o', color=colors[i], alpha=alpha,
                markersize=ms, label="_nolegend_")
        ax.plot(T_fit * 1000, y_fit, '-', color=colors[i], linewidth=lw,
                label=f"-{d['power']} dBm")

    # Auto y-limits (slightly padded), print them
    y = np.concatenate(all_scaled)
    y = y[np.isfinite(y) & (y > 0)]
    ymin, ymax = np.min(y)/1.05, np.max(y)*1.05
    ax.set_ylim(ymin, ymax)
    print(f"Y-limits (×1e7 units): {ymin:.3g} — {ymax:.3g}")

    # Linear ticks: use a reasonable number of major ticks
    try:
        from matplotlib.ticker import MaxNLocator
        ax.yaxis.set_major_locator(MaxNLocator(nbins=6))
    except Exception:
        pass

    leg = ax.legend(loc="upper right", frameon=True, handlelength=1.1, handletextpad=0.6, borderpad=0.4, labelspacing=0.5)
    leg.get_frame().set_alpha(0.9)
    for line in leg.get_lines():
        line.set_linewidth(2.0)

    ax.grid(False)
    plt.tight_layout()
    out_plot = os.path.join(OUTPUT_DIR, "qint_lin.png")
    plt.savefig(out_plot, dpi=300)
    plt.close()
    print(f"Saved plot: {out_plot}")

# -------- Main --------
def main():
    json_path = os.path.join(OUTPUT_DIR, "qint_lin.json")
    prev = load_previous(json_path)
    prev_params, prev_loss = (prev if prev is not None else (None, None))

    final_pack = None
    for round_idx in range(1, OUTER_MAX_ROUNDS + 1):
        print(f"\n=== Outer Round {round_idx} ===")
        params, loss, all_data = optimize_shared(start_params=prev_params, max_iters=MAX_ITERS)
        if prev_loss is None or loss < prev_loss - IMPROVEMENT_TOL:
            print(f"Improved loss {loss:.6e} (prev {prev_loss})")
            save_results(json_path, params, loss)
            prev_params, prev_loss = params.copy(), loss
            final_pack = (params, loss, all_data)
        else:
            print(f"No further improvement (loss {loss:.6e} >= prev {prev_loss:.6e}). Stop.")
            break

    if final_pack is None:
        loaded = load_previous(json_path)
        if loaded is None:
            raise RuntimeError("No results to plot.")
        params, best_loss = loaded
        all_data = []
        for pwr in POWERS:
            if pwr in NBAR_BY_POWER:
                try:
                    T, Q = load_dataset(pwr)
                except Exception:
                    continue
                all_data.append({"power": pwr, "T": T, "Q": Q, "nbar": NBAR_BY_POWER[pwr]})
    else:
        params, best_loss, all_data = final_pack

    print(f"Plotting with best loss {best_loss:.6e}")
    plot_shared(all_data, params, best_loss)

if __name__ == "__main__":
    main()