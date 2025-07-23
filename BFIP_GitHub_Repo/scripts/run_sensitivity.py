#!/usr/bin/env python3
"""
run_sensitivity.py

Compute Sobol sensitivity indices for BFIP metrics (MI, θ*, ΔG) using a KDE-based MI estimator.
Usage:
  python run_sensitivity.py --ion Fe2+ --metric mi --config pipeline/config.yaml --output sobol_Fe2.csv
"""
import os
import sys
# Ensure project root is on path for local imports
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

import argparse
import yaml
import numpy as np
from SALib.sample import sobol as salib_sobol_sample
from SALib.analyze import sobol as salib_sobol_analyze
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.feature_selection import mutual_info_regression

# Import BFIP model functions
from models.kinetics import hill_equation
from models.thermodynamics import gibbs_free_energy


def load_config(path):
    """Load and return the YAML config dict"""
    with open(path) as f:
        return yaml.safe_load(f)


def get_bounds(config, ion):
    """
    Extract parameter bounds [mean - std, mean + std] for BFIP parameters.
    Expects config['ions'][ion][param] to be dicts with 'mean' and 'std'.
    """
    ion_cfg = config['ions'].get(ion)
    if ion_cfg is None:
        raise KeyError(f"Ion '{ion}' not found in config.")
    names = ['Delta_H', 'Delta_S', 'K_d', 'k_on', 'k_off']
    bounds = []
    for key in names:
        stats = ion_cfg.get(key)
        if not isinstance(stats, dict) or 'mean' not in stats or 'std' not in stats:
            raise KeyError(f"Config for ion '{ion}' missing mean/std for '{key}'")
        mean = float(stats['mean'])
        std = float(stats['std'])
        bounds.append((mean - std, mean + std))
    return names, bounds


def mutual_information_histogram(x, y, bins=20):
    """
    Estimate mutual information between two continuous variables via 2D histogram.
    """
    hist_2d, xedges, yedges = np.histogram2d(x, y, bins=bins)
    p_xy = hist_2d / np.sum(hist_2d)
    p_x = np.sum(p_xy, axis=1)
    p_y = np.sum(p_xy, axis=0)
    mask = p_xy > 0
    den = np.outer(p_x, p_y)
    mi = np.sum(p_xy[mask] * np.log(p_xy[mask] / den[mask]))
    return mi


def compute_static_metrics(params, ion_cfg, temp_range):
    """
    Compute θ*, ΔG, and MI for a single parameter vector.
    Uses KDE-based MI between θ curves for nominal vs. sampled Kd.
    """
    lig_vals = ion_cfg['ligand_range']
    L_star = (lig_vals[0] + lig_vals[1]) / 2.0
    n_H = ion_cfg['n_H']['mean'] if isinstance(ion_cfg['n_H'], dict) else ion_cfg['n_H']
    Kd_nominal = ion_cfg['K_d']['mean'] if isinstance(ion_cfg['K_d'], dict) else ion_cfg['K_d']
    theta_nom = hill_equation(np.array([L_star]), n_H, Kd_nominal)[0]
    dH, dS, Kd_sample, k_on, k_off = params
    T0 = np.mean(temp_range)
    G_s = gibbs_free_energy(theta_nom, dH, dS, T0)
    lig_grid = np.linspace(lig_vals[0], lig_vals[1], 100)
    theta_ref  = hill_equation(lig_grid, n_H, Kd_nominal)
    theta_samp = hill_equation(lig_grid, n_H, Kd_sample)
    theta_ref_reshaped = theta_ref.reshape(-1, 1)
    MI_kde = mutual_info_regression(theta_ref_reshaped, theta_samp, discrete_features=False, random_state=0)
    MI = MI_kde[0]
    return {'mi': MI, 'theta': theta_nom, 'deltaG': G_s}


def main():
    parser = argparse.ArgumentParser(description="Sobol sensitivity for BFIP metrics")
    parser.add_argument('--ion',    required=True, help='Ion name (e.g. Fe2+, Ca2+, H+)')
    parser.add_argument('--metric', default='mi', choices=['mi','theta','deltaG'], help='Metric to analyze')
    parser.add_argument('--config', default='pipeline/config.yaml', help='Path to config.yaml')
    parser.add_argument('--output', default='sobol_results.csv', help='Base filename for CSV and PNG')
    args = parser.parse_args()

    cfg = load_config(args.config)
    names, bounds = get_bounds(cfg, args.ion)
    problem = {'num_vars': len(names), 'names': names, 'bounds': bounds}
    samples = salib_sobol_sample.sample(problem, N=256, calc_second_order=False)
    ion_cfg = cfg['ions'][args.ion]
    temp_range = cfg.get('temperature_range', [298.0, 310.0])
    Y = np.array([compute_static_metrics(s, ion_cfg, temp_range)[args.metric] for s in samples])
    Si = salib_sobol_analyze.analyze(problem, Y, calc_second_order=False, print_to_console=False)
    df = pd.DataFrame({'S1': Si['S1'], 'ST': Si['ST']}, index=names)
    out_csv = args.output
    df.to_csv(out_csv)
    print(f"Sobol indices saved to {out_csv}")
    title_map = {'mi':'MI','theta':'θ*','deltaG':'ΔG'}
    ax = df['S1'].plot.bar(title=f"{args.ion} {title_map[args.metric]} First-order Sobol")
    ax.set_ylabel('S1')
    plt.tight_layout()
    out_png = out_csv.replace('.csv', f"_{args.metric}.png")
    plt.savefig(out_png)
    print(f"Bar chart saved to {out_png}")

if __name__ == '__main__':
    main()
