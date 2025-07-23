
import numpy as np
import matplotlib.pyplot as plt
from thermo_dynamic_model import run_dynamic_thermo
from models.thermodynamics import gibbs_free_energy

R = 8.314
THR_THETA = 0.10
THR_MI = 2.2
T0 = 300.0

def pulsed_signal(t, phase_shift=0, high=0.65, low=0.25, repeat=80, width=40):
    return high if int(t + phase_shift) % repeat < width else low

def compute_mi(base, pert):
    hist2d, _, _ = np.histogram2d(base, pert, bins=30)
    pXY = hist2d / np.sum(hist2d)
    pX = np.sum(pXY, axis=1)
    pY = np.sum(pXY, axis=0)
    mask = pXY > 0
    denom = np.outer(pX, pY)
    return np.sum(pXY[mask] * np.log(pXY[mask] / denom[mask]))

def simulate_comparator():
    t_span = np.linspace(0, 240, 480)
    dH = -46.84
    dS = -0.102

    theta_A = []
    theta_B = []
    delta_theta = []
    mis = []
    bfips = []

    for t in t_span:
        aA = pulsed_signal(t, phase_shift=0)
        aB = pulsed_signal(t, phase_shift=40)

        kinA = {'k_on': 1e5, 'amplitude': aA}
        _, thA = run_dynamic_thermo(kinA, None, T0, [t, t + 1], dH, dS, 1.0)

        kinB = {'k_on': 1e5, 'amplitude': aB}
        _, thB = run_dynamic_thermo(kinB, None, T0, [t, t + 1], dH, dS, 1.0)

        θA = np.mean(thA)
        θB = np.mean(thB)
        θΔ = abs(θA - θB)
        mi = compute_mi(thA, thB)

        dG = gibbs_free_energy((θA + θB)/2, dH, dS, T0)
        bfip = int((θΔ > 0.07) and (mi > THR_MI) and (dG < -(R * T0 / 1000) * np.log(2)))

        theta_A.append(θA)
        theta_B.append(θB)
        delta_theta.append(θΔ)
        mis.append(mi)
        bfips.append(bfip)

    return t_span, np.array(theta_A), np.array(theta_B), np.array(delta_theta), np.array(mis), np.array(bfips)

def plot_comparator(t, θA, θB, θΔ, mis, bfips):
    plt.figure(figsize=(12, 10))

    plt.subplot(4, 1, 1)
    plt.plot(t, θA, label='θ̄ A', color='blue')
    plt.plot(t, θB, label='θ̄ B', color='orange')
    plt.axhline(0.10, color='gray', linestyle='--', label='θ*')
    plt.ylabel('θ̄(t)')
    plt.title('BFIP T6.C: Comparator (θ̄ Difference + MI)')
    plt.legend()

    plt.subplot(4, 1, 2)
    plt.plot(t, θΔ, label='|Δθ̄|', color='black')
    plt.axhline(0.07, color='gray', linestyle='--', label='Δθ̄ Threshold')
    plt.ylabel('Δθ̄')
    plt.legend()

    plt.subplot(4, 1, 3)
    plt.plot(t, mis, 'g')
    plt.axhline(2.2, color='gray', linestyle='--', label='MI Threshold')
    plt.ylabel('MI (bits)')
    plt.legend()

    plt.subplot(4, 1, 4)
    plt.step(t, bfips, where='mid', color='purple')
    plt.ylabel('BFIP')
    plt.xlabel('Time (s)')

    plt.tight_layout()
    plt.savefig('BFIP_t6c_comparator.png', dpi=600)
    plt.show()

if __name__ == '__main__':
    t, θA, θB, θΔ, mi, bfip = simulate_comparator()
    plot_comparator(t, θA, θB, θΔ, mi, bfip)
