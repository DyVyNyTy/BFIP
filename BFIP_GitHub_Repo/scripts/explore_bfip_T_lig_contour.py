#!/usr/bin/env python3
"""
explore_bfip_T_lig_contour.py

Sweep temperature and ligand concentration for a given ion at fixed ΔH, ΔS, and K_d,
then plot the BFIP region based on thresholds for θ*, MI, and ΔG.

Usage:
  python explore_bfip_T_lig_contour.py \
    --ion Fe2+ \
    --config pipeline/config.yaml \
    --out Fe2_T_Lig_contour.png
"""
import argparse
import yaml
import numpy as np
import matplotlib.pyplot as plt
import os, sys

# ensure models/ is on the import path
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

from models.kinetics import hill_equation
from models.thermodynamics import gibbs_free_energy
from models.information import compute_mutual_information

# BFIP thresholds and MI perturbation
THR_THETA  = 0.05    # binding saturation threshold
THR_MI     = 0.01    # lowered MI threshold
MI_PERTURB = 0.10    # 10% perturbation
R_gas      = 8.314   # J/(mol·K)

def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)

def run_contour(ion, cfg, out_png):
    ion_cfg = cfg['ions'][ion]

    # fixed thermodynamics and affinity
    dH_mean = ion_cfg['Delta_H']['mean']
    dS_mean = ion_cfg['Delta_S']['mean']
    Kd_mean = ion_cfg['K_d']['mean']
    n_H     = ion_cfg['n_H']['mean']

    # ranges
    T_min, T_max = cfg['T_range']
    L_min, L_max, _ = ion_cfg['ligand_range']

    T_vals = np.linspace(T_min, T_max, 200)
    L_vals = np.linspace(L_min, L_max, 200)
    T_grid, L_grid = np.meshgrid(T_vals, L_vals, indexing='ij')

    RTln2_grid = (R_gas * T_grid / 1000) * np.log(2)

    Z = np.zeros_like(T_grid, dtype=bool)

    for i, T in enumerate(T_vals):
        for j, L in enumerate(L_vals):
            # θ* at this L
            theta = hill_equation(np.array([L]), n_H, Kd_mean)[0]
            # ΔG at θ*
            G = gibbs_free_energy(theta, dH_mean, dS_mean, T)
            # MI around θ*
            P = np.array([1-theta, theta])
            Q = np.array([1-theta*(1+MI_PERTURB), theta*(1+MI_PERTURB)])
            pJ = np.outer(P, Q)
            MI = compute_mutual_information(pJ, P, Q)

            Z[i,j] = (theta > THR_THETA) and (MI > THR_MI) and (G < -RTln2_grid[i,j])

    # Plotting
    plt.figure(figsize=(6,5))
    cs = plt.contourf(L_vals, T_vals, Z, levels=[-0.5, 0.5, 1.5], cmap='viridis')
    cbar = plt.colorbar(cs, ticks=[0,1], label='BFIP region')
    plt.xlabel('Ligand concentration [L]')
    plt.ylabel('Temperature (K)')
    plt.title(f'{ion}: BFIP in [L]–T at ΔH={dH_mean:.0f}, ΔS={dS_mean:.1f}, Kd={Kd_mean:.1f}')
    plt.tight_layout()
    plt.savefig(out_png)
    print(f'Saved T–[L] BFIP contour to {out_png}')

if __name__ == '__main__':
    p = argparse.ArgumentParser(description='Plot BFIP contour in T–ligand space')
    p.add_argument('--ion',    required=True, help='Ion name, e.g. Fe2+')
    p.add_argument('--config', default='pipeline/config.yaml', help='Path to config.yaml')
    p.add_argument('--out',    default='bfip_T_Lig_contour.png', help='Output PNG filename')
    args = p.parse_args()
    cfg = load_config(args.config)
    run_contour(args.ion, cfg, args.out)
