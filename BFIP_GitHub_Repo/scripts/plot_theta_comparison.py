
#!/usr/bin/env python3
"""
plot_theta_comparison.py

Plot θ(t) for baseline and perturbed kinetics at a given ΔH and ΔS.
Usage:
  python plot_theta_comparison.py --ion Fe2+ --config pipeline/config.yaml
"""
import argparse
import yaml
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

toplevel = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, toplevel)
sys.path.insert(0, os.path.join(toplevel, 'simulation'))

from models.simulation import run_dynamic

MI_PERTURB = 0.50

def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)

def run_theta_plot(ion, cfg, out_png):
    ion_cfg = cfg['ions'][ion]

    # Use central thermodynamic values
    dH = ion_cfg['Delta_H']['mean']
    dS = ion_cfg['Delta_S']['mean']
    T0 = np.mean(cfg.get('T_range', [298.0, 310.0]))
    lig_range = np.linspace(*ion_cfg['ligand_range'])
    t_span = np.linspace(*cfg.get('t_span', [0.0, 120.0, 300]))

    # Baseline kinetics
    kin_base = {
        'k_on': ion_cfg['k_on']['mean'],
        'k_off': ion_cfg['k_off']['mean'],
        'amplitude': ion_cfg.get('amplitude', 1.0)
    }
    t, theta_base = run_dynamic(kin_base, lig_range, T0, t_span)

    # Perturbed kinetics
    kin_pert = {
        'k_on': kin_base['k_on'],
        'k_off': kin_base['k_off'] * (1 + MI_PERTURB),
        'amplitude': kin_base['amplitude'] * (1 - MI_PERTURB)
    }
    _, theta_pert = run_dynamic(kin_pert, lig_range, T0, t_span)

    # Plot comparison
    plt.figure(figsize=(7, 4))
    plt.plot(t, theta_base, label='Baseline θ(t)', lw=2)
    plt.plot(t, theta_pert, label='Perturbed θ(t)', lw=2, linestyle='--')
    plt.xlabel('Time (s)')
    plt.ylabel('θ(t)')
    plt.title(f'{ion} θ(t) Comparison at ΔH={dH}, ΔS={dS}')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(out_png)
    print(f'✅ Saved θ(t) comparison plot to {out_png}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compare θ(t) for baseline and perturbed kinetics')
    parser.add_argument('--ion',    required=True, help='Ion name, e.g. Fe2+')
    parser.add_argument('--config', default='pipeline/config.yaml', help='Path to config.yaml')
    parser.add_argument('--out',    default='theta_comparison.png', help='Output PNG filename')
    args = parser.parse_args()
    cfg = load_config(args.config)
    run_theta_plot(args.ion, cfg, args.out)
