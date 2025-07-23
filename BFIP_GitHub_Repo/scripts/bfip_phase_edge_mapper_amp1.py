
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

def test_bfip(dH, dS, amp=1.0, T0=300.0):
    t_span = np.linspace(0, 120.0, 300)
    kin_base = {'k_on': 1e5, 'amplitude': amp}
    _, theta_base = run_dynamic_thermo(kin_base, None, T0, t_span[:2], dH, dS, 1.0)

    kin_pert = {'k_on': 1e5, 'amplitude': amp * (1 - MI_PERTURB)}
    _, theta_pert = run_dynamic_thermo(kin_pert, None, T0, t_span[:2], dH, dS, 1.0)

    theta_mean = np.mean(theta_base)
    mi = compute_mi(theta_base, theta_pert)
    dG = gibbs_free_energy(theta_mean, dH, dS, T0)

    is_bfip = (theta_mean > THR_THETA) and (mi > THR_MI) and (dG < -(R_gas * T0 / 1000) * np.log(2))
    return int(is_bfip)

def map_phase_space(dH_center=-46.0, dS_center=-0.11, span=5.0, step=0.5, amp=1.0):
    dH_vals = np.arange(dH_center - span, dH_center + span + step, step)
    dS_vals = np.arange(dS_center - span, dS_center + span + step, step)
    Z = np.zeros((len(dH_vals), len(dS_vals)), dtype=int)

    for i, dH in enumerate(dH_vals):
        for j, dS in enumerate(dS_vals):
            Z[i, j] = test_bfip(dH, dS, amp)

    return dH_vals, dS_vals, Z

def plot_heatmap(dH_vals, dS_vals, Z):
    plt.figure(figsize=(8, 6))
    plt.imshow(Z.T, origin='lower', aspect='auto',
               extent=[dH_vals[0], dH_vals[-1], dS_vals[0], dS_vals[-1]],
               cmap='coolwarm', interpolation='nearest')
    plt.colorbar(label='BFIP State (1=ON, 0=OFF)')
    plt.xlabel('ΔH (kJ/mol)')
    plt.ylabel('ΔS (kJ/mol·K)')
    plt.title('Refined BFIP Phase Edge Mapping (Amp = 1.0)')
    plt.tight_layout()
    plt.savefig('BFIP_phase_edge_heatmap_amp1.0.png', dpi=600)
    plt.show()

if __name__ == '__main__':
    dH_vals, dS_vals, Z = map_phase_space()
    plot_heatmap(dH_vals, dS_vals, Z)
