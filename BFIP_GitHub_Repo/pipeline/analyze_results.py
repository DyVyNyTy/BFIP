#!/usr/bin/env python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logging
from bfip.analysis import mask_bfip
from scipy.stats import entropy

# ——— Setup logging to the Snakemake log file ———
logging.basicConfig(
    filename=snakemake.log[0],
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s'
)

try:
    logging.info("Starting analyze_results.py")

    # ——— Load config & wildcards ———
    cfg = snakemake.config
    ion = snakemake.wildcards.ion
    T = cfg["T_range"][1]
    logging.info(f"Ion = {ion}, Temperature = {T} K")

    # ——— Load the combined .npz (theta & G) ———
    archive = np.load(snakemake.input.data)
    theta_grid = archive["theta"]  # shape: (n_samples, n_ligand)
    G_grid     = archive["G"]      # same shape
    logging.info(f"Loaded θ grid shape {theta_grid.shape}, G grid shape {G_grid.shape}")

    # ——— Build ligand axis ———
    lig_vals = np.linspace(*cfg["ions"][ion]["ligand_range"])
    logging.info(f"Ligand range: {lig_vals[0]} → {lig_vals[-1]} ({len(lig_vals)} points)")

    # ——— Compute per-ligand means ———
    theta_mean = theta_grid.mean(axis=0)
    G_mean     = G_grid.mean(axis=0)
    logging.info(f"θ_mean min/max = {theta_mean.min():.3f}/{theta_mean.max():.3f}")
    logging.info(f"G_mean min/max = {G_mean.min():.3f}/{G_mean.max():.3f}")

    # ——— Mutual information per ligand ———
    col_sums = theta_grid.sum(axis=0, keepdims=True)
    P = np.where(col_sums > 0, theta_grid / col_sums, 1e-10)
    Q = np.full(len(lig_vals), 1 / len(lig_vals))
    I = entropy(P + 1e-10, axis=0) - entropy(Q)
    logging.info(f"I (mutual info) min/max = {I.min():.3f}/{I.max():.3f}")

    # ——— BFIP mask ———
    bfip_mask = mask_bfip(theta_mean, I, G_mean, T)
    logging.info(f"BFIP mask (True count) = {bfip_mask.sum()} / {len(bfip_mask)}")

    # ——— Write summary CSV ———
    df = pd.DataFrame({
        "ligand":     lig_vals,
        "theta_mean": theta_mean,
        "G_mean":     G_mean,
        "I":          I,
        "bfip":       bfip_mask
    })
    df.to_csv(snakemake.output.csv, index=False)
    logging.info(f"Wrote summary to {snakemake.output.csv}")

    # ——— Plot ligand vs. θ_mean ———
    plt.figure(figsize=(8,5))
    if bfip_mask.any():
        plt.scatter(
            lig_vals[bfip_mask],
            theta_mean[bfip_mask],
            c="red", s=50, marker="o", edgecolors="none",
            label=f"{ion} BFIP"
        )
    else:
        plt.scatter(
            lig_vals,
            theta_mean,
            c="red", s=50, marker="o", edgecolors="none",
            label=ion
        )

    plt.xlabel("Ligand Concentration")
    plt.ylabel("Mean θ")
    plt.title(f"BFIP Regions: {ion}")
    plt.ticklabel_format(useOffset=False, style="plain", axis="y")
    plt.legend()
    plt.tight_layout()
    plt.savefig(snakemake.output.png)
    logging.info(f"Saved plot to {snakemake.output.png}")
    plt.close()

    logging.info("Finished analyze_results.py successfully")

except Exception as e:
    logging.error(f"analyze_results.py FAILED: {e}")
    raise


def main():
    pass

if __name__ == '__main__':
    main()
