
import numpy as np
import matplotlib.pyplot as plt
from thermo_dynamic_model import run_dynamic_thermo
from models.thermodynamics import gibbs_free_energy

R = 8.314
THR_THETA = 0.10
THR_MI = 2.2
MI_PERTURB = 0.50
T0 = 300.0

def pulsed_amp(t, base=0.20, amp=0.45, pulse_width=60, total_time=240):
    """Return a pulse that starts ON and falls to 0 after pulse_width seconds."""
    return base + amp if t < pulse_width else 0.0

def simulate_phase_lock():
    t_span = np.linspace(0, 240, 480)
    dH = -46.84
    dS = -0.102

    amps = []
    thetas = []
    mis = []
    dGs = []
    bfips = []
    latch_state = 0

    for t in t_span:
        amp = pulsed_amp(t)
        amps.append(amp)

        kin_base = {'k_on': 1e5, 'amplitude': amp}
        _, theta_base = run_dynamic_thermo(kin_base, None, T0, [t, t + 1], dH, dS, 1.0)

        kin_pert = {'k_on': 1e5, 'amplitude': amp * (1 - MI_PERTURB)}
        _, theta_pert = run_dynamic_thermo(kin_pert, None, T0, [t, t + 1], dH, dS, 1.0)

        theta_mean = np.mean(theta_base)
        mi = compute_mi(theta_base, theta_pert)
        dG = gibbs_free_energy(theta_mean, dH, dS, T0)

        # Latch logic: if BFIP was ever ON and ΔG & θ̄ remain favorable, stay ON
        if (theta_mean > THR_THETA and mi > THR_MI and dG < -(R * T0 / 1000) * np.log(2)):
            latch_state = 1
        elif latch_state and (theta_mean > THR_THETA and dG < -(R * T0 / 1000) * np.log(2)):
            latch_state = 1  # sustain without MI
        else:
            latch_state = 0

        thetas.append(theta_mean)
        mis.append(mi)
        dGs.append(dG)
        bfips.append(latch_state)

    return t_span, np.array(amps), np.array(thetas), np.array(mis), np.array(dGs), np.array(bfips)

def compute_mi(base, pert):
    hist2d, _, _ = np.histogram2d(base, pert, bins=30)
    pXY = hist2d / np.sum(hist2d)
    pX = np.sum(pXY, axis=1)
    pY = np.sum(pXY, axis=0)
    mask = pXY > 0
    denom = np.outer(pX, pY)
    return np.sum(pXY[mask] * np.log(pXY[mask] / denom[mask]))

def plot_phase_lock(t, amps, thetas, mis, dGs, bfips):
    plt.figure(figsize=(12, 10))

    plt.subplot(4, 1, 1)
    plt.plot(t, amps, 'k--')
    plt.ylabel('Amplitude')
    plt.title('BFIP Phase Lock Register (Set Once, Retain ON)')

    plt.subplot(4, 1, 2)
    plt.plot(t, thetas, 'b')
    plt.axhline(0.10, color='gray', linestyle='--', label='θ* Threshold')
    plt.ylabel('Mean θ(t)')
    plt.legend()

    plt.subplot(4, 1, 3)
    plt.plot(t, dGs, 'c')
    plt.axhline(-(R * T0 / 1000) * np.log(2), color='gray', linestyle='--', label='ΔG < RTln2')
    plt.ylabel('ΔG')
    plt.legend()

    plt.subplot(4, 1, 4)
    plt.step(t, bfips, where='mid', color='purple')
    plt.ylabel('BFIP State')
    plt.xlabel('Time (s)')

    plt.tight_layout()
    plt.savefig('BFIP_phase_lock_register.png', dpi=600)
    plt.show()

if __name__ == '__main__':
    t, amps, thetas, mis, dGs, bfips = simulate_phase_lock()
    plot_phase_lock(t, amps, thetas, mis, dGs, bfips)
