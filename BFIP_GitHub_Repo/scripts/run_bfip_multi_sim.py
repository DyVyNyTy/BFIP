
import sys
import os
import yaml
import numpy as np
import matplotlib.pyplot as plt

# Path adjustment for local module access
sys.path.append(os.path.join(os.path.dirname(__file__), 'models', 'models'))

from kinetics import hill_equation
from thermodynamics import gibbs_free_energy
from simulation import run_dynamic

def load_config(path):
    with open(path, 'r') as f:
        return yaml.safe_load(f)

def run_multi_static(cfg):
    ligand_vals = np.linspace(cfg['ligand_range'][0], cfg['ligand_range'][1], 100)
    theta_curves = {}
    g_curves = {}

    for ion_cfg in cfg['ions']:
        ion = ion_cfg['ion'].replace('+', 'p')
        θ_vals = []
        G_vals = []
        for L in ligand_vals:
            θ = hill_equation(L, ion_cfg['parameters']['n_H'], ion_cfg['parameters']['K_d'])
            ΔG = gibbs_free_energy(θ, ion_cfg['parameters']['ΔH'], ion_cfg['parameters']['ΔS'], cfg['temperature'])
            θ_vals.append(θ)
            G_vals.append(ΔG)
        theta_curves[ion] = θ_vals
        g_curves[ion] = G_vals

    return ligand_vals, theta_curves, g_curves

def run_multi_dynamic(cfg):
    t_curves = {}
    θ_curves = {}

    for ion_cfg in cfg['ions']:
        ion = ion_cfg['ion'].replace('+', 'p')
        ion_cfg['parameters']['amplitude'] = cfg.get('amplitude', 1.0)
        t, θ_t = run_dynamic(ion_cfg['parameters'], cfg['ligand_range'], cfg['temperature'], cfg['t_span'])
        t_curves[ion] = t
        θ_curves[ion] = θ_t

    return t_curves, θ_curves

def run_and_plot(cfg_path):
    cfg = load_config(cfg_path)
    combo_name = "_".join([ion['ion'].replace('+', 'p') for ion in cfg['ions']])

    ligand_vals, theta_curves, g_curves = run_multi_static(cfg)

    # Static Plot
    fig, ax1 = plt.subplots()
    ax1.set_title(f"Multi-Ion Static Phase – {combo_name}")
    ax1.set_xlabel("Ligand Concentration")
    ax2 = ax1.twinx()

    for ion, θ_vals in theta_curves.items():
        ax1.plot(ligand_vals, θ_vals, label=f'{ion} θ')

    for ion, G_vals in g_curves.items():
        ax2.plot(ligand_vals, G_vals, label=f'{ion} ΔG', linestyle='dashed')

    ax1.set_ylabel("θ (saturation)")
    ax2.set_ylabel("ΔG")
    fig.legend(loc="lower center", bbox_to_anchor=(0.5, -0.2), ncol=2)
    fig.tight_layout()
    plt.savefig(f"{combo_name}_multi_static.png")
    plt.close()

    # Dynamic Plot
    t_curves, θ_curves = run_multi_dynamic(cfg)

    plt.figure()
    plt.title(f"Multi-Ion Dynamic θ(t) – {combo_name}")
    for ion, θ_t in θ_curves.items():
        plt.plot(t_curves[ion], θ_t, label=f'{ion} θ(t)')
    plt.xlabel("Time")
    plt.ylabel("θ(t)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{combo_name}_multi_dynamic.png")
    plt.close()

if __name__ == '__main__':
    run_and_plot(sys.argv[1])
