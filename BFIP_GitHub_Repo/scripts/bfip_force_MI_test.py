
import numpy as np
import matplotlib.pyplot as plt
from thermo_dynamic_model import run_dynamic_thermo
from models.thermodynamics import gibbs_free_energy

R = 8.314
THR_THETA = 0.10
MI_FORCED = 3.0  # Force MI artificially above threshold
MI_REQUIRED = 2.2

def sinusoidal_amp(t, base=0.20, amp=0.45, freq=1/60):
    return base + amp * np.sin(2 * np.pi * freq * t)

def simulate_force_mi(T0=300.0):
    t_span = np.linspace(0, 240, 480)
    dH = -40
    dS = -0.05

    amps = []
    thetas = []
    deltaGs = []
    bfips = []

    for t in t_span:
        amp = sinusoidal_amp(t)
        amps.append(amp)

        kin = {'k_on': 1e5, 'amplitude': amp}
        _, theta_base = run_dynamic_thermo(kin, None, T0, [t, t + 1], dH, dS, 1.0)
        theta_mean = np.mean(theta_base)
        dG = gibbs_free_energy(theta_mean, dH, dS, T0)

        bfip = int((theta_mean > THR_THETA) and (MI_FORCED > MI_REQUIRED) and (dG < -(R * T0 / 1000) * np.log(2)))

        thetas.append(theta_mean)
        deltaGs.append(dG)
        bfips.append(bfip)

    return t_span, np.array(amps), np.array(thetas), np.array(deltaGs), np.array(bfips)

def plot_force_mi(t, amps, thetas, dGs, bfips):
    plt.figure(figsize=(12, 10))

    plt.subplot(4, 1, 1)
    plt.plot(t, amps, 'k--')
    plt.ylabel('Amplitude')
    plt.title('Forced MI > 2.2 — BFIP Flicker Test (H⁺)')

    plt.subplot(4, 1, 2)
    plt.plot(t, thetas, 'b')
    plt.axhline(0.10, color='gray', linestyle='--', label='θ* Threshold')
    plt.ylabel('Mean θ(t)')
    plt.legend()

    plt.subplot(4, 1, 3)
    plt.plot(t, dGs, 'c')
    plt.axhline(-(R * 300 / 1000) * np.log(2), color='gray', linestyle='--', label='ΔG < RTln2')
    plt.ylabel('ΔG')
    plt.legend()

    plt.subplot(4, 1, 4)
    plt.step(t, bfips, where='mid', color='purple')
    plt.ylabel('BFIP')
    plt.xlabel('Time (s)')

    plt.tight_layout()
    plt.savefig('BFIP_force_MI_test.png', dpi=600)
    plt.show()

if __name__ == '__main__':
    t, amps, thetas, dGs, bfips = simulate_force_mi()
    plot_force_mi(t, amps, thetas, dGs, bfips)
