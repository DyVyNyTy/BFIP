
#!/usr/bin/env python3
"""
explore_bfip_dynamic_ds_contour.py

Sweep ΔH and ΔS, simulate dynamic θ(t), compute dynamic MI, then plot the BFIP region.
"""
import argparse
import yaml
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# Ensure paths to local modules
toplevel = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, toplevel)
sys.path.insert(0, os.path.join(toplevel, 'simulation'))

from models.simulation import run_dynamic
from models.information import compute_mutual_information
from models.thermodynamics import gibbs_free_energy

# Thresholds
THR_THETA   = 0.10   # θ* mean threshold
THR_MI      = 2.4   # Mutual Information threshold
MI_PERTURB  = 0.50   # 50% perturbation applied to k_off and amplitude
R_gas       = 8.314  # J/mol·K

def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)

def run_contour(ion, cfg, out_png):
    ion_cfg = cfg['ions'][ion]

    # Load simulation parameters
    dH0       = ion_cfg['Delta_H']['mean']
    dS0       = ion_cfg['Delta_S']['mean']
    lig_range = np.linspace(*ion_cfg['ligand_range'])
    T0        = np.mean(cfg.get('T_range', [298.0, 310.0]))
    t_span    = np.linspace(*cfg.get('t_span', [0.0, 120.0, 300]))

    # ΔH / ΔS sweep space (20×20 for speed and resolution)
    spanH = ion_cfg['Delta_H']['std'] * 2
    spanS = ion_cfg['Delta_S']['std'] * 2
    H_vals = np.linspace(dH0 - spanH, dH0 + spanH, 20)
    S_vals = np.linspace(dS0 - spanS, dS0 + spanS, 20)
    H, S = np.meshgrid(H_vals, S_vals, indexing='ij')

    Z = np.zeros_like(H, dtype=bool)

    for i, dH in enumerate(H_vals):
        for j, dS in enumerate(S_vals):
            print(f"→ Simulating ΔH = {dH:.1f}, ΔS = {dS:.2f}...", flush=True)

            # Baseline kinetic parameters
            kin_base = {
                'k_on': ion_cfg['k_on']['mean'],
                'k_off': ion_cfg['k_off']['mean'],
                'amplitude': ion_cfg.get('amplitude', 1.0)
            }
            t, theta_base = run_dynamic(kin_base, lig_range, T0, t_span)

            # Dual perturbation: increase k_off and decrease amplitude
            kin_pert = {
                'k_on': kin_base['k_on'],
                'k_off': kin_base['k_off'] * (1 + MI_PERTURB),
                'amplitude': kin_base['amplitude'] * (1 - MI_PERTURB)
            }
            _, theta_pert = run_dynamic(kin_pert, lig_range, T0, t_span)

            # Mutual Information
            hist2d, _, _ = np.histogram2d(theta_base, theta_pert, bins=30)
            pXY = hist2d / np.sum(hist2d)
            pX = np.sum(pXY, axis=1)
            pY = np.sum(pXY, axis=0)
            mask = pXY > 0
            denom = np.outer(pX, pY)
            mi_dyn = np.sum(pXY[mask] * np.log(pXY[mask] / denom[mask]))

            # Thermodynamic favorability
            theta_mean = np.mean(theta_base)
            G = gibbs_free_energy(theta_mean, dH, dS, T0)

            # Log diagnostics
            print(f"    θ̄ = {theta_mean:.4f}, MI = {mi_dyn:.4f}, ΔG = {G:.2f}")

            # Apply BFIP condition
            Z[i, j] = (
                (theta_mean > THR_THETA) and
                (mi_dyn > THR_MI) and
                (G < - (R_gas * T0 / 1000) * np.log(2))
            )

    # Plot BFIP region
    plt.figure(figsize=(6, 5))
    cs = plt.contourf(H, S, Z.T, levels=[-0.5, 0.5, 1.5], cmap='coolwarm')
    plt.colorbar(cs, ticks=[0, 1], label='Dynamic BFIP region')
    plt.xlabel('ΔH (kJ/mol)')
    plt.ylabel('ΔS (kJ/(mol·K))')
    plt.title(f'{ion}: Dynamic BFIP in ΔH–ΔS')
    plt.tight_layout()
    plt.savefig(out_png)
    print(f'\n✅ Saved dynamic BFIP contour to {out_png}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot dynamic BFIP contour in ΔH–ΔS')
    parser.add_argument('--ion',    required=True, help='Ion name, e.g. Fe2+')
    parser.add_argument('--config', default='pipeline/config.yaml', help='Path to config.yaml')
    parser.add_argument('--out',    default='dyn_dH_dS_contour.png', help='Output PNG filename')
    args = parser.parse_args()
    cfg = load_config(args.config)
    run_contour(args.ion, cfg, args.out)
