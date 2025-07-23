# BFIP Parameter Sweep Simulator
# Explores where BFIP activation occurs across ΔG₀ and ligand space

import numpy as np
import matplotlib.pyplot as plt
import os

def calculate_free_energy(ligand, deltaG0, RT):
    return deltaG0 + RT * np.log(ligand + 1e-8)

def calculate_theta(deltaG, RT):
    return 1 / (1 + np.exp(deltaG / RT))

def calculate_mutual_information(theta):
    eps = 1e-8
    theta = np.clip(theta, eps, 1 - eps)
    return -(theta * np.log2(theta) + (1 - theta) * np.log2(1 - theta))

def run_parameter_sweep(ion_name="NewIon", output_dir="sweep_results"):
    ligand_range = np.linspace(0.01, 10, 100)
    deltaG0_range = np.linspace(-30, 0, 100)
    temp = 310.15  # Kelvin
    RT = 8.314 * temp / 1000

    BFIP_map = np.zeros((len(deltaG0_range), len(ligand_range)))

    for i, deltaG0 in enumerate(deltaG0_range):
        for j, ligand in enumerate(ligand_range):
            dG = calculate_free_energy(ligand, deltaG0, RT)
            theta = calculate_theta(dG, RT)
            mi = calculate_mutual_information(theta)

            if (theta > 0.1) and (mi > 0.5) and (dG < -RT * np.log(2)):
                BFIP_map[i, j] = 1

    # Save as .npy and plot
    os.makedirs(output_dir, exist_ok=True)
    np.save(os.path.join(output_dir, f"{ion_name}_BFIP_sweep.npy"), BFIP_map)

    plt.figure(figsize=(8, 6))
    plt.imshow(BFIP_map, extent=[0.01, 10, -30, 0], aspect='auto', origin='lower', cmap='viridis')
    plt.colorbar(label="BFIP Activated (1=True, 0=False)")
    plt.xlabel("Ligand Concentration")
    plt.ylabel("ΔG₀ (kcal/mol)")
    plt.title(f"BFIP Phase Activation Sweep: {ion_name}")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{ion_name}_BFIP_sweep.png"))
    plt.show()

    print(f"✅ Sweep completed for {ion_name}. Data saved in '{output_dir}'")

if __name__ == "__main__":
    run_parameter_sweep()
