
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

def check_bfip(dH, dS, amp=0.30, T0=300.0):
    t_span = np.linspace(0, 120.0, 300)
    kin_base = {'k_on': 1e5, 'amplitude': amp}
    n = 1.0
    _, theta_base = run_dynamic_thermo(kin_base, None, T0, t_span[:2], dH, dS, n)

    kin_pert = {'k_on': 1e5, 'amplitude': amp * (1 - MI_PERTURB)}
    _, theta_pert = run_dynamic_thermo(kin_pert, None, T0, t_span[:2], dH, dS, n)

    theta_mean = np.mean(theta_base)
    mi = compute_mi(theta_base, theta_pert)
    dG = gibbs_free_energy(theta_mean, dH, dS, T0)

    return (theta_mean > THR_THETA) and (mi > THR_MI) and (dG < -(R_gas * T0 / 1000) * np.log(2))

def evaluate_quadrants():
    # Define candidate ΔH/ΔS coordinates from earlier scans
    logic_points = {
        '00 (OFF, OFF)':  ((-20.0, -0.05), (-20.0, -0.05)),
        '10 (ON, OFF)':   ((-46.84, -0.102), (-20.0, -0.05)),
        '01 (OFF, ON)':   ((-20.0, -0.05), (-46.0, -0.110)),
        '11 (ON, ON)':    ((-46.84, -0.102), (-46.0, -0.110))
    }

    results = []
    for label, (h_coords, ca_coords) in logic_points.items():
        H_on = check_bfip(*h_coords)
        Ca_on = check_bfip(*ca_coords)
        pattern = f"{int(H_on)}{int(Ca_on)}"
        results.append((label, pattern, h_coords, ca_coords))

    return results

def print_results(results):
    print("\n🧪 BFIP Logic Quadrant Test Results:")
    for label, pattern, h_coords, ca_coords in results:
        print(f"• Target: {label}")
        print(f"  ➤ Result: {pattern} | H⁺: {h_coords} | Ca²⁺: {ca_coords}\n")

if __name__ == '__main__':
    results = evaluate_quadrants()
    print_results(results)
