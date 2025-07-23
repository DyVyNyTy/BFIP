
import numpy as np
import matplotlib.pyplot as plt
from thermo_dynamic_model import run_dynamic_thermo
from models.thermodynamics import gibbs_free_energy

R = 8.314
THR_THETA = 0.10
THR_MI = 2.2
T0 = 300.0

def chain_input(t, pulse_width=40, repeat=80):
    """Returns pulsed input at defined intervals"""
    return 0.65 if int(t) % repeat < pulse_width else 0.20

def compute_mi(base, pert):
    hist2d, _, _ = np.histogram2d(base, pert, bins=30)
    pXY = hist2d / np.sum(hist2d)
    pX = np.sum(pXY, axis=1)
    pY = np.sum(pXY, axis=0)
    mask = pXY > 0
    denom = np.outer(pX, pY)
    return np.sum(pXY[mask] * np.log(pXY[mask] / denom[mask]))

def simulate_bfip_chain():
    t_span = np.linspace(0, 240, 480)
    dH = -46.84
    dS = -0.102

    gates = 3
    theta_chain = [[] for _ in range(gates)]
    bfip_chain = [[] for _ in range(gates)]
    inputs = []

    states = [0] * gates

    for t in t_span:
        amp = chain_input(t)
        inputs.append(amp)

        prev_amp = amp
        for i in range(gates):
            kin = {'k_on': 1e5, 'amplitude': prev_amp}
            _, theta_base = run_dynamic_thermo(kin, None, T0, [t, t + 1], dH, dS, 1.0)

            kin_pert = {'k_on': 1e5, 'amplitude': prev_amp * 0.9}
            _, theta_pert = run_dynamic_thermo(kin_pert, None, T0, [t, t + 1], dH, dS, 1.0)

            theta_mean = np.mean(theta_base)
            mi = compute_mi(theta_base, theta_pert)
            dG = gibbs_free_energy(theta_mean, dH, dS, T0)

            if (theta_mean > THR_THETA and mi > THR_MI and dG < -(R * T0 / 1000) * np.log(2)):
                states[i] = 1
            elif states[i] == 1 and theta_mean > THR_THETA and dG < -(R * T0 / 1000) * np.log(2):
                states[i] = 1
            else:
                states[i] = 0

            theta_chain[i].append(theta_mean)
            bfip_chain[i].append(states[i])
            prev_amp = theta_mean  # Pass to next gate

    return t_span, np.array(inputs), np.array(theta_chain), np.array(bfip_chain)

def plot_bfip_chain(t, inputs, thetas, bfips):
    plt.figure(figsize=(12, 10))

    plt.subplot(3, 1, 1)
    plt.plot(t, inputs, 'k--')
    plt.ylabel('Input Amplitude')
    plt.title('BFIP T6.A: Memory Chain (Sequential Phase Lock)')

    plt.subplot(3, 1, 2)
    for i in range(thetas.shape[0]):
        plt.plot(t, thetas[i], label=f'Gate {i+1}')
    plt.axhline(0.10, color='gray', linestyle='--', label='θ* Threshold')
    plt.ylabel('θ̄(t)')
    plt.legend()

    plt.subplot(3, 1, 3)
    for i in range(bfips.shape[0]):
        plt.step(t, bfips[i] + i, where='mid', label=f'Gate {i+1}')
    plt.ylabel('BFIP States')
    plt.xlabel('Time (s)')
    plt.legend()

    plt.tight_layout()
    plt.savefig('BFIP_t6a_memory_chain.png', dpi=600)
    plt.show()

if __name__ == '__main__':
    t, amps, thetas, bfips = simulate_bfip_chain()
    plot_bfip_chain(t, amps, thetas, bfips)
