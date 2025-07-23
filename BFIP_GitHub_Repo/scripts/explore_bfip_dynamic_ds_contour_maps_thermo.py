#!/usr/bin/env python3
"""
explore_bfip_dynamic_ds_contour_maps_thermo.py

Sweep ΔH and ΔS, simulate dynamic θ(t), compute MI and θ̄ maps, and plot full region truth surfaces.
"""
import argparse
import yaml
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

toplevel = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, toplevel)

from thermo_dynamic_model import run_dynamic_thermo
from models.information import compute_mutual_information
from models.thermodynamics import gibbs_free_energy

# Thresholds
THR_THETA   = 0.10
THR_MI      = 2.2
MI_PERTURB  = 0.50
R_gas       = 8.314

def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)

def run_contour_maps(ion, cfg, base_out):
    ion_cfg = cfg['ions'][ion]

    dH0       = ion_cfg['Delta_H']['mean']
    dS0       = ion_cfg['Delta_S']['mean']
    T0        = np.mean(cfg.get('T_range', [298.0, 310.0]))
    t_span    = np.linspace(*cfg.get('t_span', [0.0, 120.0, 300]))

    spanH = ion_cfg['Delta_H']['std'] * 2
    spanS = ion_cfg['Delta_S']['std'] * 2
    H_vals = np.linspace(dH0 - spanH, dH0 + spanH, 20)
    S_vals = np.linspace(dS0 - spanS, dS0 + spanS, 20)
    H, S = np.meshgrid(H_vals, S_vals, indexing='ij')

    Z = np.zeros_like(H, dtype=bool)
    MI_map = np.zeros_like(H)
    THETA_map = np.zeros_like(H)

    for i, dH in enumerate(H_vals):
        for j, dS in enumerate(S_vals):
            print(f"→ ΔH = {dH:.1f}, ΔS = {dS:.2f}", flush=True)

            kin_base = {
                'k_on': ion_cfg['k_on']['mean'],
                'amplitude': ion_cfg.get('amplitude', 1.0)
            }
            n_H = ion_cfg.get('n_H', {}).get('mean', 1.0)
            t, theta_base = run_dynamic_thermo(kin_base, None, T0, t_span, dH, dS, n_H)

            kin_pert = {
                'k_on': kin_base['k_on'],
                'amplitude': kin_base['amplitude'] * (1 - MI_PERTURB)
            }
            _, theta_pert = run_dynamic_thermo(kin_pert, None, T0, t_span, dH, dS, n_H)

            hist2d, _, _ = np.histogram2d(theta_base, theta_pert, bins=30)
            pXY = hist2d / np.sum(hist2d)
            pX = np.sum(pXY, axis=1)
            pY = np.sum(pXY, axis=0)
            mask = pXY > 0
            denom = np.outer(pX, pY)
            mi_dyn = np.sum(pXY[mask] * np.log(pXY[mask] / denom[mask]))

            theta_mean = np.mean(theta_base)
            G = gibbs_free_energy(theta_mean, dH, dS, T0)

            MI_map[i, j] = mi_dyn
            THETA_map[i, j] = theta_mean

            Z[i, j] = (
                (theta_mean > THR_THETA) and
                (mi_dyn > THR_MI) and
                (G < - (R_gas * T0 / 1000) * np.log(2))
            )

    # Plot results
    plt.figure(figsize=(6, 5))
    cs = plt.contourf(H, S, Z.T, levels=[-0.5, 0.5, 1.5], cmap='coolwarm')
    plt.colorbar(cs, ticks=[0, 1], label='Dynamic BFIP region')
    plt.xlabel('ΔH (kJ/mol)')
    plt.ylabel('ΔS (kJ/(mol·K))')
    plt.title(f'{ion}: BFIP Region (Z)')
    plt.tight_layout()
    plt.savefig(base_out + "_Z.png")

    plt.figure(figsize=(6, 5))
    cs = plt.contourf(H, S, MI_map.T, cmap='viridis')
    plt.colorbar(cs, label='Mutual Information (bits)')
    plt.xlabel('ΔH (kJ/mol)')
    plt.ylabel('ΔS (kJ/(mol·K))')
    plt.title(f'{ion}: Mutual Information Map')
    plt.tight_layout()
    plt.savefig(base_out + "_MI.png")

    plt.figure(figsize=(6, 5))
    cs = plt.contourf(H, S, THETA_map.T, cmap='plasma')
    plt.colorbar(cs, label='Mean θ')
    plt.xlabel('ΔH (kJ/mol)')
    plt.ylabel('ΔS (kJ/(mol·K))')
    plt.title(f'{ion}: Mean θ Map')
    plt.tight_layout()
    plt.savefig(base_out + "_THETA.png")

    # ✅ Save raw maps for overlap atlas reconstruction
    np.save(base_out + "_Z.npy", Z.astype(int))
    np.save(base_out + "_MI.npy", MI_map)
    np.save(base_out + "_THETA.npy", THETA_map)

    print(f"✅ Saved Z, MI, θ̄ plots and raw data arrays to '{base_out}_*.png' and '.npy'")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot BFIP + MI + θ maps using thermodynamic model')
    parser.add_argument('--ion',    required=True, help='Ion name, e.g. H+')
    parser.add_argument('--config', default='pipeline/config.yaml', help='Path to config.yaml')
    parser.add_argument('--out',    default='maps_thermo', help='Base output name (no extension)')
    args = parser.parse_args()
    cfg = load_config(args.config)
    run_contour_maps(args.ion, cfg, args.out)
