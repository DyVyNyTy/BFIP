
import sys
import os
import yaml
import numpy as np
import matplotlib.pyplot as plt

# Dynamically add the models directory to the path
sys.path.append(os.path.join(os.path.dirname(__file__), 'models', 'models'))

from kinetics import hill_equation
from thermodynamics import gibbs_free_energy
from simulation import run_dynamic

def load_config(path):
    with open(path, 'r') as f:
        return yaml.safe_load(f)

def run_static_simulation(cfg):
    ligand_vals = np.linspace(cfg['ligand_range'][0], cfg['ligand_range'][1], 100)
    θ_vals = []
    G_vals = []
    for L in ligand_vals:
        θ = hill_equation(L, cfg['parameters']['n_H'], cfg['parameters']['K_d'])
        ΔG = gibbs_free_energy(θ, cfg['parameters']['ΔH'], cfg['parameters']['ΔS'], cfg['temperature'])
        θ_vals.append(θ)
        G_vals.append(ΔG)
    return ligand_vals, np.array(θ_vals), np.array(G_vals)

def run_and_plot(cfg_path):
    cfg = load_config(cfg_path)
    ion_name = cfg['ion'].replace('+', 'p')

    ligand_vals, θ_vals, G_vals = run_static_simulation(cfg)

    # Plotting θ and ΔG
    fig, ax1 = plt.subplots()
    ax1.set_title(f"Static BFIP Phase – {ion_name}")
    ax1.plot(ligand_vals, θ_vals, label='θ (saturation)', color='blue')
    ax1.set_xlabel('Ligand Concentration')
    ax1.set_ylabel('θ', color='blue')
    ax2 = ax1.twinx()
    ax2.plot(ligand_vals, G_vals, label='ΔG', color='red')
    ax2.set_ylabel('ΔG', color='red')
    fig.tight_layout()
    plt.savefig(f"{ion_name}_static_phase.png")
    plt.close()

    # Run dynamic
    cfg['parameters']['amplitude'] = cfg.get('amplitude', 1.0)
    t, θ_t = run_dynamic(cfg['parameters'], cfg['ligand_range'], cfg['temperature'], cfg['t_span'])

    # Plot θ over time
    plt.figure()
    plt.title(f"Dynamic θ(t) – {ion_name}")
    plt.plot(t, θ_t, label='θ(t)', color='green')
    plt.xlabel('Time')
    plt.ylabel('θ(t)')
    plt.tight_layout()
    plt.savefig(f"{ion_name}_dynamic_theta.png")
    plt.close()

if __name__ == '__main__':
    run_and_plot(sys.argv[1])
