
import numpy as np
import matplotlib.pyplot as plt
from thermo_dynamic_model import run_dynamic_thermo
from models.thermodynamics import gibbs_free_energy

R = 8.314
THR_THETA = 0.10
THR_MI = 2.2
MI_PERTURB = 0.50

def compute_mi(base, pert):
    hist2d, _, _ = np.histogram2d(base, pert, bins=30)
    pXY = hist2d / np.sum(hist2d)
    pX = np.sum(pXY, axis=1)
    pY = np.sum(pXY, axis=0)
    mask = pXY > 0
    denom = np.outer(pX, pY)
    return np.sum(pXY[mask] * np.log(pXY[mask] / denom[mask]))

def simulate_gate_chain(n_gates=3, amp_base=0.50, amp_high=0.80, T0=300.0):
    t_span = np.linspace(0, 120, 300)
    dH = -46.84
    dS = -0.102

    theta_history = []
    bfip_states = []

    for gate in range(n_gates):
        amp = amp_base + gate * (amp_high - amp_base) / max(n_gates - 1, 1)
        kin_base = {'k_on': 1e5, 'amplitude': amp}
        _, theta_base = run_dynamic_thermo(kin_base, None, T0, t_span, dH, dS, 1.0)

        kin_pert = {'k_on': 1e5, 'amplitude': amp * (1 - MI_PERTURB)}
        _, theta_pert = run_dynamic_thermo(kin_pert, None, T0, t_span, dH, dS, 1.0)

        theta_mean = np.mean(theta_base)
        mi = compute_mi(theta_base, theta_pert)
        dG = gibbs_free_energy(theta_mean, dH, dS, T0)
        is_bfip = int((theta_mean > THR_THETA) and (mi > THR_MI) and (dG < -(R * T0 / 1000) * np.log(2)))

        theta_history.append(theta_base)
        bfip_states.append(is_bfip)

    return t_span, theta_history, bfip_states

def plot_logic_chain(t, theta_data, states):
    plt.figure(figsize=(10, 6))
    for i, theta in enumerate(theta_data):
        plt.plot(t, theta, label=f'Gate {i+1} (BFIP={states[i]})')
    plt.xlabel('Time (s)')
    plt.ylabel('θ(t)')
    plt.title('Chained BFIP Logic Gates (H⁺)')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('BFIP_chained_logic_gates.png', dpi=600)
    plt.show()

if __name__ == '__main__':
    t, thetas, logic = simulate_gate_chain()
    plot_logic_chain(t, thetas, logic)
