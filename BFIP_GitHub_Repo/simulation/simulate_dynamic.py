#!/usr/bin/env python
import yaml, argparse
import numpy as np
from bfip.simulation import run_dynamic

parser = argparse.ArgumentParser()
parser.add_argument("--config", required=True)
parser.add_argument("--ion",    required=True)
parser.add_argument("--out",    required=True)
args = parser.parse_args()

# Load config + params
cfg = yaml.safe_load(open(args.config))
ion_cfg = cfg["ions"][args.ion]
samples = np.load(f"sampling/{args.ion}/params.npz")
# For simplicity, take the *mean* kinetics from your sampled arrays:
kin = {
    "k_on":  float(samples["k_on"].mean()),
    "k_off": float(samples["k_off"].mean()),
    "n_H":   float(ion_cfg["n_H"]["mean"]),
    "amplitude": 0.1  # or pull from config if you add it
}
T      = cfg["T_range"][1]
t_span = cfg["t_span"]

# Run dynamic sim
t, theta_t = run_dynamic(kin, ion_cfg["conc_range"][0], T, t_span)

# Save
np.savez(args.out, t=t, theta=theta_t)


def main():
    pass

if __name__ == '__main__':
    main()
