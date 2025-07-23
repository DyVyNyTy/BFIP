#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import os

# Define ions and file paths
ions = {
    "Fe2+": "dynamic/Fe2+/theta_dynamic.npz",
    "Ca2+": "dynamic/Ca2+/theta_dynamic.npz",
    "H+": "dynamic/H+/theta_dynamic.npz"
}

for ion, path in ions.items():
    if not os.path.exists(path):
        print(f"[WARN] File not found: {path}")
        continue

    # Load NPZ
    data = np.load(path)
    time = data["time"]
    theta = data["theta"]  # shape: (n_samples, len(time))

    # Compute stats
    theta_mean = np.mean(theta, axis=0)
    theta_std = np.std(theta, axis=0)

    # Plot
    plt.figure(figsize=(8, 5))
    plt.plot(time, theta_mean, label="Mean θ(t)", color="blue")
    plt.fill_between(time, theta_mean - theta_std, theta_mean + theta_std,
                     color="blue", alpha=0.3, label="±1 SD")

    plt.title(f"Dynamic θ(t) Over Time: {ion}")
    plt.xlabel("Time (s)")
    plt.ylabel("θ Occupancy")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"theta_dynamics_{ion.replace('+','plus')}.png")
    plt.close()

print("✅ Dynamic θ(t) plots saved.")


def main():
    pass

if __name__ == '__main__':
    main()
