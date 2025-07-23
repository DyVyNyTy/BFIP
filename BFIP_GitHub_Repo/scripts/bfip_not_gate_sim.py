
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

def control_amp(t, write=0.72, base=0.45, period=60):
    return base + (write - base) * ((np.sin(2 * np.pi * t / period) > 0).astype(float))

def simulate_not_gate(T0=300.0):
    t_span = np.linspace(0, 240, 480)
    dH = -46.84
    dS = -0.102

    a_states = []
    b_states = []
    a_amps = []

    for t in t_span:
        amp_A = control_amp(t)
        amp_B = 1.0 - (amp_A - 0.4)  # Inverted amplitude logic

        a_amps.append(amp_A)

        def gate_response(amp):
            kin_base = {'k_on': 1e5, 'amplitude': amp}
            _, theta_base = run_dynamic_thermo(kin_base, None, T0, [t, t + 1], dH, dS, 1.0)
            kin_pert = {'k_on': 1e5, 'amplitude': amp * (1 - MI_PERTURB)}
            _, theta_pert = run_dynamic_thermo(kin_pert, None, T0, [t, t + 1], dH, dS, 1.0)

            theta_mean = np.mean(theta_base)
            mi = compute_mi(theta_base, theta_pert)
            dG = gibbs_free_energy(theta_mean, dH, dS, T0)
            return int((theta_mean > THR_THETA) and (mi > THR_MI) and (dG < -(R * T0 / 1000) * np.log(2)))

        a_states.append(gate_response(amp_A))
        b_states.append(gate_response(amp_B))

    return t_span, np.array(a_amps), np.array(a_states), np.array(b_states)

def plot_not_gate(t, amps, a_logic, b_logic):
    plt.figure(figsize=(12, 6))
    plt.step(t, a_logic, label='Gate A (Input)', where='mid')
    plt.step(t, b_logic + 1.1, label='Gate B (NOT Output)', where='mid')
    plt.plot(t, amps, 'k--', label='Input Amplitude A')
    plt.xlabel('Time (s)')
    plt.ylabel('Logic State (0/1)')
    plt.title('BFIP NOT Gate Simulation (H⁺)')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('BFIP_not_gate_sim.png', dpi=600)
    plt.show()

if __name__ == '__main__':
    t, amps, a, b = simulate_not_gate()
    plot_not_gate(t, amps, a, b)
