
import numpy as np
import csv
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

def test_bfip(dH, dS, amp, kon, nH, T0=300.0):
    t_span = np.linspace(0, 120.0, 300)
    kin_base = {'k_on': kon, 'amplitude': amp}
    _, theta_base = run_dynamic_thermo(kin_base, None, T0, t_span[:2], dH, dS, nH)

    kin_pert = {'k_on': kon, 'amplitude': amp * (1 - MI_PERTURB)}
    _, theta_pert = run_dynamic_thermo(kin_pert, None, T0, t_span[:2], dH, dS, nH)

    theta_mean = np.mean(theta_base)
    mi = compute_mi(theta_base, theta_pert)
    dG = gibbs_free_energy(theta_mean, dH, dS, T0)

    is_bfip = (theta_mean > THR_THETA) and (mi > THR_MI) and (dG < -(R_gas * T0 / 1000) * np.log(2))
    return is_bfip, theta_mean, mi, dG

def scan_parameters():
    dH = -46.0
    dS = -0.11

    amplitudes = [0.75, 0.85, 0.95, 1.0]
    kon_vals = [1e4, 1e5, 1e6]
    nH_vals = [1.0, 1.5, 2.0]

    results = []

    for amp in amplitudes:
        for kon in kon_vals:
            for nH in nH_vals:
                active, theta, mi, dG = test_bfip(dH, dS, amp, kon, nH)
                logic = int(active)
                results.append([
                    amp, kon, nH, logic, theta, mi, dG
                ])

    with open("BFIP_param_scan_Ca.csv", "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Amplitude', 'k_on', 'n_H', 'BFIP', 'Theta_mean', 'MI', 'dG'])
        writer.writerows(results)

    print("✅ Parameter scan complete. Results saved to BFIP_param_scan_Ca.csv")

if __name__ == '__main__':
    scan_parameters()
