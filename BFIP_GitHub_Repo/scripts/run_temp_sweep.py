
import os
import sys
import yaml
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), 'models', 'models'))
from kinetics import hill_equation
from thermodynamics import gibbs_free_energy

def load_config(path):
    with open(path, 'r') as f:
        return yaml.safe_load(f)

def temperature_sweep(cfg_path, output_dir):
    cfg = load_config(cfg_path)
    ion = cfg['ion'].replace('+', 'p')
    T_vals = np.linspace(298, 315, 6)
    ligand_vals = np.linspace(cfg['ligand_range'][0], cfg['ligand_range'][1], 100)

    plt.figure()
    for T in T_vals:
        θ_vals = []
        G_vals = []
        for L in ligand_vals:
            θ = hill_equation(L, cfg['parameters']['n_H'], cfg['parameters']['K_d'])
            ΔG = gibbs_free_energy(θ, cfg['parameters']['ΔH'], cfg['parameters']['ΔS'], T)
            θ_vals.append(θ)
            G_vals.append(ΔG)
        plt.plot(ligand_vals, θ_vals, label=f"θ @ T={T:.1f}K")

    plt.title(f"Temperature Sweep – θ Curve for {ion}")
    plt.xlabel("Ligand Concentration")
    plt.ylabel("θ (Saturation)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    output_path = os.path.join(output_dir, f"{ion}_temp_sweep_theta.png")
    plt.savefig(output_path)
    plt.close()

    plt.figure()
    for T in T_vals:
        G_vals = []
        for L in ligand_vals:
            θ = hill_equation(L, cfg['parameters']['n_H'], cfg['parameters']['K_d'])
            ΔG = gibbs_free_energy(θ, cfg['parameters']['ΔH'], cfg['parameters']['ΔS'], T)
            G_vals.append(ΔG)
        plt.plot(ligand_vals, G_vals, label=f"ΔG @ T={T:.1f}K")

    plt.title(f"Temperature Sweep – ΔG Curve for {ion}")
    plt.xlabel("Ligand Concentration")
    plt.ylabel("ΔG (kcal/mol)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    output_path = os.path.join(output_dir, f"{ion}_temp_sweep_dG.png")
    plt.savefig(output_path)
    plt.close()

if __name__ == "__main__":
    temperature_sweep(sys.argv[1], os.path.dirname(__file__))
