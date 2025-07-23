
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

def test_bfip(dH, dS, amp=0.60, T0=300.0):
    t_span = np.linspace(0, 120.0, 300)
    kin_base = {'k_on': 1e5, 'amplitude': amp}
    _, theta_base = run_dynamic_thermo(kin_base, None, T0, t_span[:2], dH, dS, 1.0)

    kin_pert = {'k_on': 1e5, 'amplitude': amp * (1 - MI_PERTURB)}
    _, theta_pert = run_dynamic_thermo(kin_pert, None, T0, t_span[:2], dH, dS, 1.0)

    theta_mean = np.mean(theta_base)
    mi = compute_mi(theta_base, theta_pert)
    dG = gibbs_free_energy(theta_mean, dH, dS, T0)

    is_bfip = (theta_mean > THR_THETA) and (mi > THR_MI) and (dG < -(R_gas * T0 / 1000) * np.log(2))
    return is_bfip, theta_mean, mi, dG

def scan_and_log():
    sweep = np.linspace(-20, 20, 9)  # -20 to +20 in steps of 5
    base_points = {
        '10': ((-46.84, -0.102), (-20.0, -0.05)),
        '01': ((-20.0, -0.05), (-46.0, -0.110)),
        '11': ((-46.84, -0.102), (-46.0, -0.110))
    }

    results = []
    for label, (h_base, ca_base) in base_points.items():
        for dx in sweep:
            for dy in sweep:
                dH_H, dS_H = h_base[0] + dx, h_base[1] + dy
                dH_Ca, dS_Ca = ca_base[0] + dx, ca_base[1] + dy

                H_on, theta_H, mi_H, dg_H = test_bfip(dH_H, dS_H)
                Ca_on, theta_Ca, mi_Ca, dg_Ca = test_bfip(dH_Ca, dS_Ca)

                logic = f"{int(H_on)}{int(Ca_on)}"
                results.append([
                    label, logic,
                    dH_H, dS_H, H_on, theta_H, mi_H, dg_H,
                    dH_Ca, dS_Ca, Ca_on, theta_Ca, mi_Ca, dg_Ca
                ])

    header = [
        'TargetLogic', 'ResultLogic',
        'dH_H', 'dS_H', 'H_on', 'theta_H', 'MI_H', 'dG_H',
        'dH_Ca', 'dS_Ca', 'Ca_on', 'theta_Ca', 'MI_Ca', 'dG_Ca'
    ]

    with open('BFIP_logic_scan_results.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(results)

    print(f"✅ Scan complete. Results saved to BFIP_logic_scan_results.csv")

if __name__ == '__main__':
    scan_and_log()
