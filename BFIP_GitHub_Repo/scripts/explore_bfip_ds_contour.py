#!/usr/bin/env python3
"""
explore_bfip_kd_lig_contour.py

Sweep Kd and ligand concentration for a given ion at fixed ΔH, ΔS, pH, and T,
then plot the BFIP region based on thresholds for θ*, MI, and ΔG.
Usage:
  python explore_bfip_kd_lig_contour.py \
    --ion Fe2+ \
    --config pipeline/config.yaml \
    --out Fe2_Kd_Lig_contour.png
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

# BFIP thresholds (hardcoded, adjusted for MI sensitivity)
THR_THETA   = 0.05    # binding saturation threshold
THR_MI      = 0.01    # lowered mutual‐information threshold
MI_PERTURB  = 0.10    # 10% perturbation around θ* for MI estimation
THR_n       = 1.0     # factor for RTln2



def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)


def run_contour(ion, cfg, out_png):
    ion_cfg = cfg['ions'][ion]
    # fixed thermodynamics
    dH0 = ion_cfg['Delta_H']['mean']
    dS0 = ion_cfg['Delta_S']['mean']
    T_vals = cfg.get('T_range', cfg.get('temperature_range', [298.0,310.0]))
    T0 = np.mean(T_vals)

    # sweep ligand and Kd
    lig_vals = ion_cfg['ligand_range']
    Lmin, Lmax = lig_vals[0], lig_vals[1]
    L = np.linspace(Lmin, Lmax, 200)
    Kd0 = ion_cfg['K_d']['mean']
    sigma_Kd = ion_cfg['K_d']['std']
    Kd = np.linspace(max(1e-6, Kd0 - 3*sigma_Kd), Kd0 + 3*sigma_Kd, 200)
    K, Lg = np.meshgrid(Kd, L, indexing='ij')

    # preparatory
    n_H = ion_cfg['n_H']['mean']
    RTln2 = 8.314 * T0/1000 * np.log(2)

    Z = np.zeros(K.shape, dtype=bool)
    for i, kd in enumerate(Kd):
        # theta*_nominal at kd
        theta_s = hill_equation(np.array([kd]), n_H, kd)[0]
        G_s = gibbs_free_energy(theta_s, dH0, dS0, T0)
        for j, lig in enumerate(L):
            # MI: compare theta curves at nominal vs slight Kd shift
            theta_ref = hill_equation(np.array([lig]), n_H, kd)[0]
            # small perturb of kd
            kp = kd * 1.01
            theta_samp = hill_equation(np.array([lig]), n_H, kp)[0]
            # discrete MI
            P = np.array([1-theta_s, theta_s])
            Q = np.array([1-theta_s*(1+MI_PERTURB), theta_s*(1+MI_PERTURB)])
            pJ = np.outer(P, Q)
            MI = compute_mutual_information(pJ, P, Q)

            Z[i,j] = (theta_s > THR_THETA) and (MI > THR_MI) and (G_s < -THR_n * RTln2)

    plt.figure(figsize=(6,5))
    cs = plt.contourf(K, Lg, Z.T, levels=[-0.5, 0.5, 1.5], cmap='plasma')
    plt.colorbar(cs, ticks=[0,1], label='BFIP region')
    plt.xlabel('K_d')
    plt.ylabel('Ligand concentration')
    plt.title(f'{ion}: BFIP region in K_d–[L] at ΔH={dH0:.1f}, ΔS={dS0:.2f}, T={T0:.0f}K')
    plt.tight_layout()
    plt.savefig(out_png)
    print(f'Saved BFIP contour map to {out_png}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot BFIP contour in Kd–Lig')
    parser.add_argument('--ion',   required=True, help='Ion name, e.g. Fe2+')
    parser.add_argument('--config',default='pipeline/config.yaml', help='Path to config.yaml')
    parser.add_argument('--out',   default='bfip_Kd_Lig_contour.png', help='Output PNG')
    args = parser.parse_args()
    cfg = load_config(args.config)
    run_contour(args.ion, cfg, args.out)
