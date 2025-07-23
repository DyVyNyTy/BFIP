
import numpy as np
import matplotlib.pyplot as plt
from thermo_dynamic_model import run_dynamic_thermo
from models.thermodynamics import gibbs_free_energy
import os

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

def simulate_bfip_flicker(dH, dS, amp_list, T0=300.0):
    t_span = np.linspace(0, 120.0, 300)
    n_H = 1.0
    k_on = 1e5  # generic value

    results = []

    for amp in amp_list:
        kin_base = {'k_on': k_on, 'amplitude': amp}
        t, theta_base = run_dynamic_thermo(kin_base, None, T0, t_span, dH, dS, n_H)

        kin_pert = {'k_on': k_on, 'amplitude': amp * (1 - MI_PERTURB)}
        _, theta_pert = run_dynamic_thermo(kin_pert, None, T0, t_span, dH, dS, n_H)

        mi = compute_mi(theta_base, theta_pert)
        theta_mean = np.mean(theta_base)
        dG = gibbs_free_energy(theta_mean, dH, dS, T0)

        bfip_active = (theta_mean > THR_THETA) and (mi > THR_MI) and (dG < -(R_gas * T0 / 1000) * np.log(2))

        results.append({
            'amplitude': amp,
            'theta_mean': theta_mean,
            'MI': mi,
            'dG': dG,
            'BFIP': bfip_active,
            't': t,
            'theta': theta_base
        })

    return results

def plot_results(results):
    plt.figure(figsize=(10, 6))
    for res in results:
        label = f"Amp={res['amplitude']:.2f} | BFIP={'Yes' if res['BFIP'] else 'No'}"
        plt.plot(res['t'], res['theta'], label=label)
    plt.xlabel("Time (s)")
    plt.ylabel("θ(t)")
    plt.title("θ(t) under varying amplitude at ΔH/ΔS")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("BFIP_theta_dynamics.png", dpi=600)
    plt.show()

if __name__ == '__main__':
    dH = -46.8421
    dS = -0.1021
    amp_list = [0.5, 0.75, 1.0, 1.25, 1.5]
    results = simulate_bfip_flicker(dH, dS, amp_list)
    plot_results(results)
