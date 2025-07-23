
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

def logic_pulse(t, period=60):
    phase = (t % period) / period
    return 0.75 if phase < 0.5 else 0.45  # High half of period

def simulate_and_gate(T0=300.0):
    t_span = np.linspace(0, 240, 480)
    dH = -46.84
    dS = -0.102

    a_logic, b_logic, c_logic = [], [], []

    def gate_response(amp, t):
        kin_base = {'k_on': 1e5, 'amplitude': amp}
        _, theta_base = run_dynamic_thermo(kin_base, None, T0, [t, t + 1], dH, dS, 1.0)
        kin_pert = {'k_on': 1e5, 'amplitude': amp * (1 - MI_PERTURB)}
        _, theta_pert = run_dynamic_thermo(kin_pert, None, T0, [t, t + 1], dH, dS, 1.0)
        theta_mean = np.mean(theta_base)
        mi = compute_mi(theta_base, theta_pert)
        dG = gibbs_free_energy(theta_mean, dH, dS, T0)
        return int((theta_mean > THR_THETA) and (mi > THR_MI) and (dG < -(R * T0 / 1000) * np.log(2)))

    for t in t_span:
        amp_A = logic_pulse(t, period=80)       # A high 0–40s, low 40–80s
        amp_B = logic_pulse(t + 20, period=80)  # B delayed → creates various combinations

        A_on = gate_response(amp_A, t)
        B_on = gate_response(amp_B, t)

        amp_C = 0.85 if (A_on and B_on) else 0.45
        C_on = gate_response(amp_C, t)

        a_logic.append(A_on)
        b_logic.append(B_on)
        c_logic.append(C_on)

    return t_span, np.array(a_logic), np.array(b_logic), np.array(c_logic)

def plot_and_gate(t, A, B, C):
    plt.figure(figsize=(12, 6))
    plt.step(t, A, label='Gate A', where='mid')
    plt.step(t, B + 1.1, label='Gate B', where='mid')
    plt.step(t, C + 2.2, label='AND Output (C)', where='mid')
    plt.xlabel('Time (s)')
    plt.ylabel('Logic State')
    plt.title('BFIP AND Gate Simulation (H⁺)')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('BFIP_and_gate_sim.png', dpi=600)
    plt.show()

if __name__ == '__main__':
    t, A, B, C = simulate_and_gate()
    plot_and_gate(t, A, B, C)
