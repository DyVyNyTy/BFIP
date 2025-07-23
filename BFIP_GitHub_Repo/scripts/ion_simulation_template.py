# Ion BFIP Simulation Template
# Copy this script into any ion folder and modify appropriately

import numpy as np
import os

def calculate_free_energy(ligand_conc, params):
    # TODO: Insert ion-specific ΔG model
    # Example: ΔG = ΔG0 + RT * ln([ligand])
    return params["deltaG0"] + params["RT"] * np.log(ligand_conc + 1e-8)

def calculate_theta(deltaG, params):
    # Langmuir isotherm example: θ = 1 / (1 + exp(ΔG / RT))
    return 1 / (1 + np.exp(deltaG / params["RT"]))

def calculate_mutual_information(theta):
    # Simple proxy: MI = -θ log2(θ) - (1-θ) log2(1-θ)
    eps = 1e-8
    theta = np.clip(theta, eps, 1 - eps)
    return -(theta * np.log2(theta) + (1 - theta) * np.log2(1 - theta))

def run_simulation(output_dir, ion_name):
    ligand_range = np.linspace(0, 10, 20)  # Customize range
    temp = 310.15  # 37C
    RT = 8.314 * temp / 1000

    params = {
        "deltaG0": -20.0,  # Modify per ion
        "RT": RT
    }

    Z_map = np.zeros((20, 20))
    THETA_map = np.zeros((20, 20))
    MI_map = np.zeros((20, 20))

    for i, lig1 in enumerate(ligand_range):
        for j, lig2 in enumerate(ligand_range):
            ligand = (lig1 + lig2) / 2  # Mix or define separate axes
            dG = calculate_free_energy(ligand, params)
            theta = calculate_theta(dG, params)
            mi = calculate_mutual_information(theta)

            Z_map[j, i] = int((theta > 0.1) and (mi > 0.5) and (dG < -RT * np.log(2)))
            THETA_map[j, i] = theta
            MI_map[j, i] = mi

    os.makedirs(output_dir, exist_ok=True)
    np.save(os.path.join(output_dir, f"{ion_name}_thermo_Z.npy"), Z_map)
    np.save(os.path.join(output_dir, f"{ion_name}_thermo_THETA.npy"), THETA_map)
    np.save(os.path.join(output_dir, f"{ion_name}_thermo_MI.npy"), MI_map)

    print(f"✅ Maps saved for {ion_name} in {output_dir}")

if __name__ == "__main__":
    run_simulation(output_dir="maps", ion_name="NEWION")
