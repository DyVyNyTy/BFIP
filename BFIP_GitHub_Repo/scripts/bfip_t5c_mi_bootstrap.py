
import numpy as np
import matplotlib.pyplot as plt
from thermo_dynamic_model import run_dynamic_thermo
from models.thermodynamics import gibbs_free_energy

R = 8.314
THR_THETA = 0.10
THR_MI = 2.2
T0 = 300.0

def bootstrap_amp(t, mi_trace, base=0.20, amp=0.45, freq=1/60):
    """Feedback amplitude adjusted by smoothed MI"""
    mi_feedback = np.clip(np.mean(mi_trace[-5:]), 0.5, 2.5)
    gain = (mi_feedback - 0.5) / 2.0  # Normalize to [0,1]
    mod = gain * 0.2  # Max +0.2 adjustment
    return base + amp * np.sin(2 * np.pi * freq * t) + mod

def compute_mi(base, pert):
    hist2d, _, _ = np.histogram2d(base, pert, bins=30)
    pXY = hist2d / np.sum(hist2d)
    pX = np.sum(pXY, axis=1)
    pY = np.sum(pXY, axis=0)
    mask = pXY > 0
    denom = np.outer(pX, pY)
    return np.sum(pXY[mask] * np.log(pXY[mask] / denom[mask]))

def simulate_mi_bootstrap():
    t_span = np.linspace(0, 240, 480)
    dH = -46.84
    dS = -0.102

    amps = []
    thetas = []
    mis = []
    dGs = []
    bfips = []

    for i, t in enumerate(t_span):
        amp = bootstrap_amp(t, mis) if i >= 5 else 0.20 + 0.45 * np.sin(2 * np.pi * (1/60) * t)
        amps.append(amp)

        kin_base = {'k_on': 1e5, 'amplitude': amp}
        _, theta_base = run_dynamic_thermo(kin_base, None, T0, [t, t + 1], dH, dS, 1.0)

        kin_pert = {'k_on': 1e5, 'amplitude': amp * 0.9}
        _, theta_pert = run_dynamic_thermo(kin_pert, None, T0, [t, t + 1], dH, dS, 1.0)

        theta_mean = np.mean(theta_base)
        mi = compute_mi(theta_base, theta_pert)
        dG = gibbs_free_energy(theta_mean, dH, dS, T0)

        bfip = int((theta_mean > THR_THETA) and (mi > THR_MI) and (dG < -(R * T0 / 1000) * np.log(2)))

        thetas.append(theta_mean)
        mis.append(mi)
        dGs.append(dG)
        bfips.append(bfip)

    return t_span, np.array(amps), np.array(thetas), np.array(mis), np.array(dGs), np.array(bfips)

def plot_mi_bootstrap(t, amps, thetas, mis, dGs, bfips):
    plt.figure(figsize=(12, 10))

    plt.subplot(4, 1, 1)
    plt.plot(t, amps, 'k--')
    plt.ylabel('Amplitude')
    plt.title('BFIP T5.C: MI Feedback Bootstrap')

    plt.subplot(4, 1, 2)
    plt.plot(t, thetas, 'b')
    plt.axhline(0.10, color='gray', linestyle='--', label='θ* Threshold')
    plt.ylabel('Mean θ(t)')
    plt.legend()

    plt.subplot(4, 1, 3)
    plt.plot(t, mis, 'g')
    plt.axhline(2.2, color='gray', linestyle='--', label='MI Threshold')
    plt.ylabel('Mutual Info (bits)')
    plt.legend()

    plt.subplot(4, 1, 4)
    plt.step(t, bfips, where='mid', color='purple')
    plt.ylabel('BFIP State')
    plt.xlabel('Time (s)')

    plt.tight_layout()
    plt.savefig('BFIP_t5c_mi_bootstrap.png', dpi=600)
    plt.show()

if __name__ == '__main__':
    t, amps, thetas, mis, dGs, bfips = simulate_mi_bootstrap()
    plot_mi_bootstrap(t, amps, thetas, mis, dGs, bfips)
