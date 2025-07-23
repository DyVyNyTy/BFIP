
import numpy as np
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

def test_bfip(dH, dS, amp=0.50, T0=300.0):
    t_span = np.linspace(0, 120.0, 300)
    kin_base = {'k_on': 1e5, 'amplitude': amp}
    _, theta_base = run_dynamic_thermo(kin_base, None, T0, t_span[:2], dH, dS, 1.0)

    kin_pert = {'k_on': 1e5, 'amplitude': amp * (1 - MI_PERTURB)}
    _, theta_pert = run_dynamic_thermo(kin_pert, None, T0, t_span[:2], dH, dS, 1.0)

    theta_mean = np.mean(theta_base)
    mi = compute_mi(theta_base, theta_pert)
    dG = gibbs_free_energy(theta_mean, dH, dS, T0)

    return (theta_mean > THR_THETA) and (mi > THR_MI) and (dG < -(R_gas * T0 / 1000) * np.log(2))

def scan_for_bfip_targets():
    # Each tuple: (label, initial H⁺, initial Ca²⁺)
    targets = {
        "10 (H⁺ ON, Ca²⁺ OFF)":  ((-46.84, -0.102), (-20.0, -0.05)),
        "01 (H⁺ OFF, Ca²⁺ ON)":  ((-20.0, -0.05), (-46.0, -0.11)),
        "11 (H⁺ ON, Ca²⁺ ON)":   ((-46.84, -0.102), (-46.0, -0.11))
    }

    sweep = np.linspace(-10, 10, 5)  # ±10 kJ/mol or kJ/mol·K
    found = []

    for label, (h_guess, ca_guess) in targets.items():
        for dx in sweep:
            for dy in sweep:
                dH_H, dS_H = h_guess[0] + dx, h_guess[1] + dy
                dH_Ca, dS_Ca = ca_guess[0] + dx, ca_guess[1] + dy

                H_on = test_bfip(dH_H, dS_H)
                Ca_on = test_bfip(dH_Ca, dS_Ca)

                logic = f"{int(H_on)}{int(Ca_on)}"
                if logic == label.split()[0]:
                    found.append((label, logic, (dH_H, dS_H), (dH_Ca, dS_Ca)))
                    print(f"✅ Found match for {label}: logic {logic}")
                    break
            if found and found[-1][0] == label:
                break

    return found

if __name__ == '__main__':
    results = scan_for_bfip_targets()
    print("\n🎯 BFIP Active Pair Matches:")
    for label, logic, h_coords, ca_coords in results:
        print(f"• {label} → Logic: {logic}")
        print(f"  ➤ H⁺: {h_coords}")
        print(f"  ➤ Ca²⁺: {ca_coords}\n")
