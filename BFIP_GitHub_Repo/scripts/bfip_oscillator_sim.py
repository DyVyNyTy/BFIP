
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

def sinusoidal_amp(t, base=0.45, amp=0.3, freq=1/60):
    return base + amp * np.sin(2 * np.pi * freq * t)

def simulate_oscillator(T0=300.0):
    t_span = np.linspace(0, 240, 480)
    dH = -46.84
    dS = -0.102

    amplitudes = []
    bfip_states = []

    for t in t_span:
        amp = sinusoidal_amp(t)
        amplitudes.append(amp)

        kin_base = {'k_on': 1e5, 'amplitude': amp}
        _, theta_base = run_dynamic_thermo(kin_base, None, T0, [t, t + 1], dH, dS, 1.0)
        kin_pert = {'k_on': 1e5, 'amplitude': amp * (1 - MI_PERTURB)}
        _, theta_pert = run_dynamic_thermo(kin_pert, None, T0, [t, t + 1], dH, dS, 1.0)

        theta_mean = np.mean(theta_base)
        mi = compute_mi(theta_base, theta_pert)
        dG = gibbs_free_energy(theta_mean, dH, dS, T0)
        bfip = int((theta_mean > THR_THETA) and (mi > THR_MI) and (dG < -(R * T0 / 1000) * np.log(2)))

        bfip_states.append(bfip)

    return t_span, np.array(amplitudes), np.array(bfip_states)

def plot_oscillator(t, amps, bfip):
    plt.figure(figsize=(12, 6))
    plt.plot(t, amps, 'k--', label='Input Amplitude (sinusoidal)')
    plt.step(t, bfip + 1.2, label='BFIP Logic State (0=OFF, 1=ON)', where='mid')
    plt.xlabel('Time (s)')
    plt.ylabel('Amplitude / Logic State')
    plt.title('BFIP Oscillator Simulation (H⁺ Clock)')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('BFIP_oscillator_sim.png', dpi=600)
    plt.show()

if __name__ == '__main__':
    t, amps, states = simulate_oscillator()
    plot_oscillator(t, amps, states)
