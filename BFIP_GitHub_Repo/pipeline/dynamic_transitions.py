#!/usr/bin/env python
import numpy as np
from bfip.simulation import run_dynamic

# Snakemake context
ion = snakemake.wildcards.ion
cfg = snakemake.config
t_span = cfg["t_span"]
T = cfg["T_range"][1]

# Get fixed values from config
n_H = cfg["ions"][ion]["n_H"]["mean"]
amplitude = cfg["ions"][ion].get("amplitude", 1.0)

# Load sampled params
archive = np.load(snakemake.input.params)
kinetic_names = ["K_d", "k_on", "k_off"]
thermo_names = ["Delta_H", "Delta_S"]

n_samples = archive[kinetic_names[0]].shape[0]
samples = []
for i in range(n_samples):
    kp = {k: float(archive[k][i]) for k in kinetic_names}
    kp["n_H"] = n_H
    kp["amplitude"] = amplitude
    tp = {k: float(archive[k][i]) for k in thermo_names}
    samples.append((kp, tp))

# Pull ligand range spec directly from config
ligand_range_spec = cfg["ions"][ion]["ligand_range"]

# Run dynamic transitions
theta_dynamics = []
time_grid = None
for kp, tp in samples:
    t, theta_t = run_dynamic(kp, ligand_range_spec, T, t_span)
    if time_grid is None:
        time_grid = t
    theta_dynamics.append(theta_t)

# Save output
theta_dynamics = np.array(theta_dynamics)
np.savez(snakemake.output[0], time=time_grid, theta=theta_dynamics)


def main():
    pass

if __name__ == '__main__':
    main()
