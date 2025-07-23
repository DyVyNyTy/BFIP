
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

def entropy_sweep(cfg_path, output_dir):
    cfg = load_config(cfg_path)
    ion = cfg['ion'].replace('+', 'p')
    S_vals = np.linspace(0.005, 0.02, 6)
    ligand_vals = np.linspace(cfg['ligand_range'][0], cfg['ligand_range'][1], 100)
    T = cfg['temperature']

    plt.figure()
    for ΔS in S_vals:
        θ_vals = []
        G_vals = []
        for L in ligand_vals:
            θ = hill_equation(L, cfg['parameters']['n_H'], cfg['parameters']['K_d'])
            ΔG = gibbs_free_energy(θ, cfg['parameters']['ΔH'], ΔS, T)
            θ_vals.append(θ)
            G_vals.append(ΔG)
        plt.plot(ligand_vals, θ_vals, label=f"θ @ ΔS={ΔS:.3f}")

    plt.title(f"Entropy Sweep – θ Curve for {ion}")
    plt.xlabel("Ligand Concentration")
    plt.ylabel("θ (Saturation)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    output_path = os.path.join(output_dir, f"{ion}_entropy_sweep_theta.png")
    plt.savefig(output_path)
    plt.close()

    plt.figure()
    for ΔS in S_vals:
        G_vals = []
        for L in ligand_vals:
            θ = hill_equation(L, cfg['parameters']['n_H'], cfg['parameters']['K_d'])
            ΔG = gibbs_free_energy(θ, cfg['parameters']['ΔH'], ΔS, T)
            G_vals.append(ΔG)
        plt.plot(ligand_vals, G_vals, label=f"ΔG @ ΔS={ΔS:.3f}")

    plt.title(f"Entropy Sweep – ΔG Curve for {ion}")
    plt.xlabel("Ligand Concentration")
    plt.ylabel("ΔG (kcal/mol)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    output_path = os.path.join(output_dir, f"{ion}_entropy_sweep_dG.png")
    plt.savefig(output_path)
    plt.close()

if __name__ == "__main__":
    entropy_sweep(sys.argv[1], os.path.dirname(__file__))
