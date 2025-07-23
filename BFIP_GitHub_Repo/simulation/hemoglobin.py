#!/usr/bin/env python
import numpy as np
import logging
from bfip.kinetics import hill_equation

# Set up logging
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

try:
    logging.info("Starting simulate_hemoglobin.py")
    
    # Configuration & wildcards injected by Snakemake
    cfg = snakemake.config
    ion = snakemake.wildcards.ion
    params = cfg["ions"][ion]
    logging.info(f"Processing ion: {ion}")

    # Load the Monte Carlo samples
    logging.info(f"Loading params from {snakemake.input.params}")
    sampled = np.load(snakemake.input.params)
    n_samples = sampled["K_d"].shape[0]
    logging.info(f"Number of samples: {n_samples}")

    # Ligand axis and temperature
    lig_vals = np.linspace(*params["ligand_range"])
    T = cfg["T_range"][1]
    logging.info(f"Ligand values: {lig_vals[:5]}... (first 5), Temperature: {T}")

    # Allocate storage
    theta_grid = np.zeros((n_samples, len(lig_vals)))
    G_grid = np.zeros_like(theta_grid)
    logging.info("Allocated theta_grid and G_grid")

    # Static binding + free-energy calculation per sample
    for i in range(n_samples):
        Kd = float(sampled["K_d"][i])
        nH = params["n_H"]["mean"]
        DH = float(sampled["Delta_H"][i])
        DS = float(sampled["Delta_S"][i])
        logging.info(f"Sample {i}: Kd={Kd}, nH={nH}, Delta_H={DH}, Delta_S={DS}")

        theta = hill_equation(lig_vals, nH, Kd)
        G = DH - T * DS * 0.001

        theta_grid[i] = theta
        G_grid[i] = G

    # Save output (fix: use first index of snakemake.output)
    output_path = snakemake.output[0]
    logging.info(f"Saving output to {output_path}")
    np.savez(output_path, theta=theta_grid, G=G_grid)
    logging.info("Simulation completed successfully")

except Exception as e:
    logging.error(f"Error in simulate_hemoglobin.py: {str(e)}")
    raise


def main():
    pass

if __name__ == '__main__':
    main()
