
import numpy as np
import matplotlib.pyplot as plt
from thermo_dynamic_model import run_dynamic_thermo
from models.thermodynamics import gibbs_free_energy

R_gas = 8.314
THR_THETA = 0.10
THR_MI = 2.2
MI_PERTURB = 0.50

def compute_mi(theta_base, theta_pert):
    hist2d, _, _ = np.histogram2d(theta_base, theta_pert, bins=30)
    pXY = hist2d / np.sum(hist2d)
    pX = np.sum(pXY, axis=1)
    pY = np.sum(pXY, axis=0)
    mask = pXY > 0
    denom = np.outer(pX, pY)
    return np.sum(pXY[mask] * np.log(pXY[mask] / denom[mask]))

def sweep_amplitude_threshold(dH, dS, T0=300.0):
    t_span = np.linspace(0, 120.0, 300)
    n_H = 1.0
    k_on = 1e5

    amp_vals = np.linspace(0.1, 0.5, 20)
    results = []

    for amp in amp_vals:
        kin_base = {'k_on': k_on, 'amplitude': amp}
        t, theta_base = run_dynamic_thermo(kin_base, None, T0, t_span, dH, dS, n_H)

        kin_pert = {'k_on': k_on, 'amplitude': amp * (1 - MI_PERTURB)}
        _, theta_pert = run_dynamic_thermo(kin_pert, None, T0, t_span, dH, dS, n_H)

        theta_mean = np.mean(theta_base)
        mi = compute_mi(theta_base, theta_pert)
        dG = gibbs_free_energy(theta_mean, dH, dS, T0)
        bfip = (theta_mean > THR_THETA) and (mi > THR_MI) and (dG < -(R_gas * T0 / 1000) * np.log(2))

        results.append({
            'amplitude': amp,
            'theta_mean': theta_mean,
            'MI': mi,
            'dG': dG,
            'BFIP': bfip
        })

    return results

def plot_flicker_threshold(results):
    amps = [r['amplitude'] for r in results]
    thetas = [r['theta_mean'] for r in results]
    mis = [r['MI'] for r in results]
    dgs = [r['dG'] for r in results]
    bfip_status = [r['BFIP'] for r in results]

    plt.figure(figsize=(10, 6))
    plt.plot(amps, thetas, label="θ̄", marker='o')
    plt.plot(amps, mis, label="MI", marker='o')
    plt.plot(amps, dgs, label="ΔG", marker='o')
    plt.axhline(THR_THETA, color='gray', linestyle='--', label='θ* threshold')
    plt.axhline(THR_MI, color='green', linestyle='--', label='MI threshold')
    plt.xlabel("Amplitude")
    plt.ylabel("Value")
    plt.title("BFIP Metrics vs. Amplitude")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("BFIP_flicker_threshold.png", dpi=600)
    plt.show()

    plt.figure(figsize=(8, 4))
    plt.plot(amps, bfip_status, drawstyle='steps-post', marker='o')
    plt.xlabel("Amplitude")
    plt.ylabel("BFIP Active (1=Yes, 0=No)")
    plt.title("BFIP Phase Flicker Threshold")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("BFIP_flicker_binary.png", dpi=600)
    plt.show()

if __name__ == '__main__':
    dH = -46.8421
    dS = -0.1021
    results = sweep_amplitude_threshold(dH, dS)
    plot_flicker_threshold(results)
