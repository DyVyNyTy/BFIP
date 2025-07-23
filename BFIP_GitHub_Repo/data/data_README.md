# 📊 BFIP Dataset Overview

This directory contains all key data outputs used in the BFIP research project. These CSV files are simulation results derived from threshold analysis, modulation sweeps, and multi-ion dynamics.

---

## Files

### 🧪 `bfip_sweep.csv`
- Contains full parameter sweep results across ions (Fe²⁺, Ca²⁺, H⁺)
- Columns:
  - `ion`, `K_d_scaled`, `n_H`, `L*`, `theta*`, `ΔG`, `MI_approx`, `bfip`

### 🌊 `bfip_pH_sweep.csv`
- Variation of pH conditions for H⁺ to assess stability of BFIP emergence
- Columns:
  - `pH`, `L*`, `theta*`, `ΔG`, `MI_approx`, `bfip`

### 🔁 `scaffold_breakpoint.csv`
- Sweep of scaffold-driven K₍d₎ fluctuation
- Measures phase coherence loss at high modulation amplitudes
- Columns:
  - `mod_amp`, `theta_mean`, `theta_std`, `dynamic_range`, `bfip_break`

---

## ⚛ Multi-Ion Dynamic Simulations

Each file includes time-series data for θ(t) of two ions under shared noisy ligand conditions.

| File | Description |
|------|-------------|
| `Fe2+_H+_BFIP_BFIP_dynamics.csv` | Both ions in BFIP state |
| `Fe2+_BFIP_H+_nonBFIP_dynamics.csv` | One BFIP, one degraded |
| `Fe2+_nonBFIP_H+_BFIP_dynamics.csv` | Reversed roles |
| `Ca2+_BFIP_H+_BFIP_dynamics.csv` | Dual BFIP |
| `Ca2+_nonBFIP_Fe2+_BFIP_dynamics.csv` | Ca²⁺ degraded, Fe²⁺ stable |

Columns in each:
- `time`, `θ_1`, `θ_2`

---

Each dataset supports analysis, plotting, and scientific validation of the BFIP phase behavior under multiple conditions.

