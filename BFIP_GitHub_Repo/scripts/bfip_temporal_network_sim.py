
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

def pulse_amplitude(t, base=0.4, peak=0.75, period=60):
    return base + (peak - base) * (np.sin(2 * np.pi * t / period) > 0).astype(float)

def simulate_temporal_chain(n_gates=3, T0=300.0):
    t_span = np.linspace(0, 300, 600)
    dH = -46.84
    dS = -0.102

    results = []
    amplitudes = []
    for t in t_span:
        amp = pulse_amplitude(t)
        amplitudes.append(amp)
        gate_results = []
        for g in range(n_gates):
            gate_amp = amp + g * 0.05
            kin_base = {'k_on': 1e5, 'amplitude': gate_amp}
            _, theta_base = run_dynamic_thermo(kin_base, None, T0, [t, t + 1], dH, dS, 1.0)

            kin_pert = {'k_on': 1e5, 'amplitude': gate_amp * (1 - MI_PERTURB)}
            _, theta_pert = run_dynamic_thermo(kin_pert, None, T0, [t, t + 1], dH, dS, 1.0)

            theta_mean = np.mean(theta_base)
            mi = compute_mi(theta_base, theta_pert)
            dG = gibbs_free_energy(theta_mean, dH, dS, T0)
            bfip = int((theta_mean > THR_THETA) and (mi > THR_MI) and (dG < -(R * T0 / 1000) * np.log(2)))

            gate_results.append((theta_mean, mi, dG, bfip))
        results.append(gate_results)

    return t_span, np.array(amplitudes), np.array(results)

def plot_temporal_results(t, amps, data):
    n_gates = data.shape[1]
    bfip_matrix = data[:, :, 3]

    plt.figure(figsize=(12, 8))

    for g in range(n_gates):
        plt.plot(t, bfip_matrix[:, g] + g, label=f'Gate {g+1} BFIP')

    plt.plot(t, amps, 'k--', label='Input Amplitude (shared)')
    plt.xlabel('Time (s)')
    plt.ylabel('BFIP State (+ Gate Index)')
    plt.title('Temporal BFIP Logic Network Response')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('BFIP_temporal_logic_network.png', dpi=600)
    plt.show()

if __name__ == '__main__':
    t, amps, data = simulate_temporal_chain()
    plot_temporal_results(t, amps, data)
