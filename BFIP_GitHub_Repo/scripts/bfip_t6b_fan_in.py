
import numpy as np
import matplotlib.pyplot as plt
from thermo_dynamic_model import run_dynamic_thermo
from models.thermodynamics import gibbs_free_energy

R = 8.314
THR_THETA = 0.10
THR_MI = 2.2
T0 = 300.0

def weak_pulse(t, phase=0, width=40, repeat=80, base=0.18, high=0.26):
    """Generate weak staggered pulses from multiple sources"""
    return high if int(t + phase) % repeat < width else base

def compute_mi(base, pert):
    hist2d, _, _ = np.histogram2d(base, pert, bins=30)
    pXY = hist2d / np.sum(hist2d)
    pX = np.sum(pXY, axis=1)
    pY = np.sum(pXY, axis=0)
    mask = pXY > 0
    denom = np.outer(pX, pY)
    return np.sum(pXY[mask] * np.log(pXY[mask] / denom[mask]))

def simulate_fan_in():
    t_span = np.linspace(0, 240, 480)
    dH = -46.84
    dS = -0.102

    amp_inputs = [[], [], []]
    amp_sum = []
    thetas = []
    mis = []
    dGs = []
    bfips = []

    for t in t_span:
        a1 = weak_pulse(t, phase=0)
        a2 = weak_pulse(t, phase=20)
        a3 = weak_pulse(t, phase=40)
        a_total = (a1 + a2 + a3) / 3

        amp_inputs[0].append(a1)
        amp_inputs[1].append(a2)
        amp_inputs[2].append(a3)
        amp_sum.append(a_total)

        kin = {'k_on': 1e5, 'amplitude': a_total}
        _, th_base = run_dynamic_thermo(kin, None, T0, [t, t + 1], dH, dS, 1.0)
        kin_pert = {'k_on': 1e5, 'amplitude': a_total * 0.9}
        _, th_pert = run_dynamic_thermo(kin_pert, None, T0, [t, t + 1], dH, dS, 1.0)

        th_mean = np.mean(th_base)
        mi = compute_mi(th_base, th_pert)
        dG = gibbs_free_energy(th_mean, dH, dS, T0)

        bfip = int((th_mean > THR_THETA) and (mi > THR_MI) and (dG < -(R * T0 / 1000) * np.log(2)))
        thetas.append(th_mean)
        mis.append(mi)
        dGs.append(dG)
        bfips.append(bfip)

    return t_span, np.array(amp_inputs), np.array(amp_sum), np.array(thetas), np.array(mis), np.array(bfips)

def plot_fan_in(t, amp_inputs, amp_sum, thetas, mis, bfips):
    plt.figure(figsize=(12, 10))

    plt.subplot(4, 1, 1)
    for i in range(3):
        plt.plot(t, amp_inputs[i], label=f'Source {i+1}')
    plt.plot(t, amp_sum, 'k--', label='Summed Input')
    plt.ylabel('Amplitude')
    plt.title('BFIP T6.B: Convergent Fan-In (Weak Inputs → Strong Signal)')
    plt.legend()

    plt.subplot(4, 1, 2)
    plt.plot(t, thetas, 'b')
    plt.axhline(0.10, color='gray', linestyle='--', label='θ* Threshold')
    plt.ylabel('θ̄(t)')
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
    plt.savefig('BFIP_t6b_fan_in.png', dpi=600)
    plt.show()

if __name__ == '__main__':
    t, amps_in, amp_sum, theta, mi, bfip = simulate_fan_in()
    plot_fan_in(t, amps_in, amp_sum, theta, mi, bfip)
