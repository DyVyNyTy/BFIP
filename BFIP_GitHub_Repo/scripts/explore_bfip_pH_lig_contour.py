#!/usr/bin/env python3
"""
explore_bfip_pH_lig_contour.py

Sweep pH and ligand concentration for a given ion at fixed ΔH, ΔS, and T,
then plot the BFIP region based on thresholds for θ*, MI, and ΔG.
Usage:
  python explore_bfip_pH_lig_contour.py \
    --ion Fe2+ \
    --config pipeline/config.yaml \
    --out Fe2_pH_Lig_contour.png
"""
import argparse
import yaml
import numpy as np
import matplotlib.pyplot as plt
import os, sys
# ensure models/ is on path
toplevel = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, toplevel)

from models.kinetics import hill_equation
from models.thermodynamics import gibbs_free_energy
from models.information import compute_mutual_information

# BFIP thresholds and MI perturbation
THR_THETA  = 0.05   # binding saturation
THR_MI     = 0.01   # lowered mutual-information threshold
MI_PERTURB = 0.10   # 10% perturbation around θ for MI
RT = 8.314 / 1000  # kJ/(mol·K)


def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)


def run_contour(ion, cfg, out_png):
    ion_cfg = cfg['ions'][ion]
    # fixed thermodynamics
    dH = ion_cfg.get('Delta_H',{}).get('mean', ion_cfg.get('dH',{}).get('mean'))
    dS = ion_cfg.get('Delta_S',{}).get('mean', ion_cfg.get('dS',{}).get('mean'))
    T_vals = cfg.get('temperature_range', cfg.get('T_range', [298.0, 310.0]))
    T0 = np.mean(T_vals)

    # nominal binding parameters
    Kd = ion_cfg.get('K_d', {}).get('mean', ion_cfg.get('Kd',{}).get('mean'))
    n_H = ion_cfg.get('n_H', {}).get('mean', ion_cfg.get('nH',{}).get('mean', 1.0))

    # detect pH range key
    if 'pH_range' in ion_cfg:
        pH_min, pH_max = ion_cfg['pH_range'][:2]
    elif 'pH_bounds' in ion_cfg:
        pH_min, pH_max = ion_cfg['pH_bounds'][:2]
    else:
        raise KeyError("pH range not found in config for ion")

    # detect ligand range key
    if 'ligand_range' in ion_cfg:
        L_min, L_max = ion_cfg['ligand_range'][:2]
    elif 'ligand_bounds' in ion_cfg:
        L_min, L_max = ion_cfg['ligand_bounds'][:2]
    else:
        raise KeyError("Ligand range not found in config for ion")

    pH_vals = np.linspace(pH_min, pH_max, 200)
    L_vals  = np.linspace(L_min, L_max, 200)

    Z = np.zeros((len(pH_vals), len(L_vals)), dtype=bool)
    for i, pH in enumerate(pH_vals):
        for j, Lg in enumerate(L_vals):
         for j, Lg in enumerate(L_vals):
            # binding saturation at Lg
            theta = hill_equation(np.array([Lg]), n_H, Kd)[0]
            # Gibbs free energy at theta
            G = gibbs_free_energy(theta, dH, dS, T0)
            # MI: discrete two-point around theta
            P_dist = np.array([1 - theta, theta])
            Q_dist = np.array([1 - theta * (1 + MI_PERTURB), theta * (1 + MI_PERTURB)])
            pJ = np.outer(P_dist, Q_dist)
            MI = compute_mutual_information(pJ, P_dist, Q_dist)
            # BFIP condition
            Z[i, j] = (theta > THR_THETA) and (MI > THR_MI) and (G < -RT * np.log(2))

    # plot
    plt.figure(figsize=(6,5))
    cs = plt.contourf(L_vals, pH_vals, Z, levels=[-0.5, 0.5, 1.5], cmap='plasma')
    plt.colorbar(cs, ticks=[0,1], label='BFIP region')
    plt.xlabel('Ligand concentration [L]')
    plt.ylabel('pH')
    plt.title(f'{ion}: BFIP in pH–[L] at ΔH={dH:.1f}, ΔS={dS:.2f}, T={T0:.0f}K')
    plt.tight_layout()
    plt.savefig(out_png)
    print(f'Saved BFIP pH–[L] contour to {out_png}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot BFIP contour in pH–ligand space')
    parser.add_argument('--ion',    required=True, help='Ion name, e.g. Fe2+')
    parser.add_argument('--config', default='pipeline/config.yaml', help='Path to config.yaml')
    parser.add_argument('--out',    default='bfip_pH_Lig_contour.png', help='Output PNG')
    args = parser.parse_args()
    cfg = load_config(args.config)
    run_contour(args.ion, cfg, args.out)
