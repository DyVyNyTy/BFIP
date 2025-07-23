
import numpy as np
import matplotlib.pyplot as plt
from thermo_dynamic_model import run_dynamic_thermo
from models.thermodynamics import gibbs_free_energy

R = 8.314
THR_THETA = 0.10
THR_MI = 2.2
T0 = 300.0

def coupling_amp(t, base=0.22, amp=0.40, freq=1/60):
    return base + amp * np.sin(2 * np.pi * freq * t)

def compute_mi(base, pert):
    hist2d, _, _ = np.histogram2d(base, pert, bins=30)
    pXY = hist2d / np.sum(hist2d)
    pX = np.sum(pXY, axis=1)
    pY = np.sum(pXY, axis=0)
    mask = pXY > 0
    denom = np.outer(pX, pY)
    return np.sum(pXY[mask] * np.log(pXY[mask] / denom[mask]))

def simulate_dual_ion_coupling():
    t_span = np.linspace(0, 240, 480)
    dH_H = -46.84
    dS_H = -0.102
    dH_Ca = -46.00
    dS_Ca = -0.11

    amps = []
    theta_H = []
    theta_Ca = []
    mis = []
    bfip_logic = []

    for t in t_span:
        amp = coupling_amp(t)
        amps.append(amp)

        kin = {'k_on': 1e5, 'amplitude': amp}
        _, th_H = run_dynamic_thermo(kin, None, T0, [t, t + 1], dH_H, dS_H, 1.0)
        _, th_Ca = run_dynamic_thermo(kin, None, T0, [t, t + 1], dH_Ca, dS_Ca, 1.0)

        theta_H_mean = np.mean(th_H)
        theta_Ca_mean = np.mean(th_Ca)
        theta_H.append(theta_H_mean)
        theta_Ca.append(theta_Ca_mean)

        mi = compute_mi(th_H, th_Ca)
        mis.append(mi)

        dG_H = gibbs_free_energy(theta_H_mean, dH_H, dS_H, T0)
        bfip = int((theta_H_mean > THR_THETA) and (mi > THR_MI) and (dG_H < -(R * T0 / 1000) * np.log(2)))
        bfip_logic.append(bfip)

    return t_span, np.array(amps), np.array(theta_H), np.array(theta_Ca), np.array(mis), np.array(bfip_logic)

def plot_dual_ion_coupling(t, amps, theta_H, theta_Ca, mis, bfip_logic):
    plt.figure(figsize=(12, 10))

    plt.subplot(4, 1, 1)
    plt.plot(t, amps, 'k--')
    plt.ylabel('Amplitude')
    plt.title('BFIP T5.B: Dual-Ion Coherence (H⁺ + Ca²⁺)')

    plt.subplot(4, 1, 2)
    plt.plot(t, theta_H, label='H⁺', color='b')
    plt.plot(t, theta_Ca, label='Ca²⁺', color='orange')
    plt.axhline(0.10, color='gray', linestyle='--', label='θ* Threshold')
    plt.ylabel('θ̄ (binding)')
    plt.legend()

    plt.subplot(4, 1, 3)
    plt.plot(t, mis, 'g')
    plt.axhline(2.2, color='gray', linestyle='--', label='MI Threshold')
    plt.ylabel('Mutual Info (bits)')
    plt.legend()

    plt.subplot(4, 1, 4)
    plt.step(t, bfip_logic, where='mid', color='purple')
    plt.ylabel('BFIP State')
    plt.xlabel('Time (s)')

    plt.tight_layout()
    plt.savefig('BFIP_t5b_dual_ion_coupling.png', dpi=600)
    plt.show()

if __name__ == '__main__':
    t, amps, thH, thCa, mis, bfip = simulate_dual_ion_coupling()
    plot_dual_ion_coupling(t, amps, thH, thCa, mis, bfip)
