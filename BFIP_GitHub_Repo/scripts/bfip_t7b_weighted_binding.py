
import numpy as np
import matplotlib.pyplot as plt
from thermo_dynamic_model import run_dynamic_thermo
from models.thermodynamics import gibbs_free_energy

R = 8.314
T0 = 300.0
THR_THETA = 0.10
THR_MI = 2.2

def input_sequence(t, on=0.65, off=0.25, repeat=50, width=25):
    return on if int(t) % repeat < width else off

def compute_mi(x, y):
    hist2d, _, _ = np.histogram2d(x, y, bins=30)
    pXY = hist2d / np.sum(hist2d)
    pX = np.sum(pXY, axis=1)
    pY = np.sum(pXY, axis=0)
    mask = pXY > 0
    denom = np.outer(pX, pY)
    return np.sum(pXY[mask] * np.log(pXY[mask] / denom[mask]))

def simulate_weighted_binding():
    t_span = np.linspace(0, 240, 480)
    dH = -46.84
    dS = -0.102

    θs = []
    MIs = []
    Gs = []
    BFIP = []

    repetition_count = 0
    memory_bias = 1.0

    for t in t_span:
        amp = input_sequence(t)
        if amp > 0.6:
            repetition_count += 1
            memory_bias = 1.0 + 0.05 * repetition_count
        else:
            repetition_count = max(0, repetition_count - 1)
            memory_bias = 1.0 + 0.05 * repetition_count

        kin = {'k_on': 1e5, 'amplitude': amp}
        _, θ = run_dynamic_thermo(kin, None, T0, [t, t + 1], dH, dS, memory_bias)
        kin_p = {'k_on': 1e5, 'amplitude': amp * 0.95}
        _, θp = run_dynamic_thermo(kin_p, None, T0, [t, t + 1], dH, dS, memory_bias)

        th_mean = np.mean(θ)
        mi = compute_mi(θ, θp)
        dG = gibbs_free_energy(th_mean, dH, dS, T0)
        active = th_mean > THR_THETA and mi > THR_MI and dG < -(R*T0/1000)*np.log(2)

        θs.append(th_mean)
        MIs.append(mi)
        Gs.append(dG)
        BFIP.append(1 if active else 0)

    return t_span, np.array(θs), np.array(MIs), np.array(Gs), np.array(BFIP)

def plot_binding(t, θ, MI, G, BFIP):
    plt.figure(figsize=(12, 10))

    plt.subplot(4, 1, 1)
    plt.plot(t, θ, label='θ̄(t)', color='blue')
    plt.axhline(THR_THETA, linestyle='--', color='gray', label='θ*')
    plt.ylabel('θ̄(t)')
    plt.title('BFIP T7.B: Repetition-Weighted Binding Memory')
    plt.legend()

    plt.subplot(4, 1, 2)
    plt.plot(t, MI, label='MI', color='green')
    plt.axhline(2.2, linestyle='--', color='gray', label='MI Threshold')
    plt.ylabel('Mutual Info')
    plt.legend()

    plt.subplot(4, 1, 3)
    plt.plot(t, G, label='ΔG', color='cyan')
    plt.axhline(-(R*T0/1000)*np.log(2), linestyle='--', color='gray', label='ΔG < RTln2')
    plt.ylabel('ΔG')
    plt.legend()

    plt.subplot(4, 1, 4)
    plt.step(t, BFIP, label='BFIP', color='purple')
    plt.ylabel('BFIP')
    plt.xlabel('Time (s)')
    plt.legend()

    plt.tight_layout()
    plt.savefig('BFIP_t7b_weighted_binding.png', dpi=600)
    plt.show()

if __name__ == '__main__':
    t, θ, mi, dG, bfip = simulate_weighted_binding()
    plot_binding(t, θ, mi, dG, bfip)
