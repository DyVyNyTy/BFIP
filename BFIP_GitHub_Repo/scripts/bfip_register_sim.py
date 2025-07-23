
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

def control_pulse(t, base=0.45, write=0.70, lock=0.85, period=80):
    phase = (t % period) / period
    if phase < 0.25:
        return base  # Watch
    elif phase < 0.5:
        return write  # Write
    else:
        return lock  # Lock

def simulate_register(n_gates=3, T0=300.0):
    t_span = np.linspace(0, 240, 480)
    dH = -46.84
    dS = -0.102

    states = []
    amplitudes = []

    for t in t_span:
        amp = control_pulse(t)
        amplitudes.append(amp)
        gate_outputs = []
        for g in range(n_gates):
            amp_g = amp + g * 0.03
            kin_base = {'k_on': 1e5, 'amplitude': amp_g}
            _, theta_base = run_dynamic_thermo(kin_base, None, T0, [t, t + 1], dH, dS, 1.0)

            kin_pert = {'k_on': 1e5, 'amplitude': amp_g * (1 - MI_PERTURB)}
            _, theta_pert = run_dynamic_thermo(kin_pert, None, T0, [t, t + 1], dH, dS, 1.0)

            theta_mean = np.mean(theta_base)
            mi = compute_mi(theta_base, theta_pert)
            dG = gibbs_free_energy(theta_mean, dH, dS, T0)
            bfip = int((theta_mean > THR_THETA) and (mi > THR_MI) and (dG < -(R * T0 / 1000) * np.log(2)))

            gate_outputs.append(bfip)
        states.append(gate_outputs)

    return t_span, np.array(amplitudes), np.array(states)

def plot_register(t, amps, data):
    n_gates = data.shape[1]
    plt.figure(figsize=(12, 7))

    for g in range(n_gates):
        plt.step(t, data[:, g] + g, label=f'Gate {g+1} Bit', where='mid')

    plt.plot(t, amps, 'k--', label='Input Amplitude (Control)')
    plt.xlabel('Time (s)')
    plt.ylabel('Logic Bit (+ Index)')
    plt.title('BFIP Phase-Based Register Simulation (3-Gate)')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('BFIP_phase_register.png', dpi=600)
    plt.show()

if __name__ == '__main__':
    t, amps, state_matrix = simulate_register()
    plot_register(t, amps, state_matrix)
