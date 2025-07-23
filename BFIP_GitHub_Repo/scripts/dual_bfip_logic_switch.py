
import numpy as np
import matplotlib.pyplot as plt
from thermo_dynamic_model import run_dynamic_thermo
from models.thermodynamics import gibbs_free_energy

R_gas = 8.314
THR_THETA = 0.10
THR_MI = 2.2
MI_PERTURB = 0.50

def compute_mi(theta_base, theta_pert):
    hist2d, _, _ = np.histogram2d(theta_base, theta_pert, bins=30)
    pXY = hist2d / np.sum(hist2d)
    pX = np.sum(pXY, axis=1)
    pY = np.sum(pXY, axis=0)
    mask = pXY > 0
    denom = np.outer(pX, pY)
    return np.sum(pXY[mask] * np.log(pXY[mask] / denom[mask]))

def generate_pulse_wave(t_span, low=0.10, high=0.30, period=60):
    pulse_wave = []
    for t in t_span:
        if (t % period) < (period / 2):
            pulse_wave.append(high)
        else:
            pulse_wave.append(low)
    return np.array(pulse_wave)

def simulate_dual_bfip(dH_H, dS_H, dH_Ca, dS_Ca, T0=300.0):
    t_span = np.linspace(0, 120.0, 300)
    amp_wave = generate_pulse_wave(t_span)

    logic_states = []
    logic_labels = []

    for idx, amp in enumerate(amp_wave):
        def run_logic(dH, dS):
            kin_base = {'k_on': 1e5, 'amplitude': amp}
            n = 1.0
            _, theta_base = run_dynamic_thermo(kin_base, None, T0, t_span[:2], dH, dS, n)
            kin_pert = {'k_on': 1e5, 'amplitude': amp * (1 - MI_PERTURB)}
            _, theta_pert = run_dynamic_thermo(kin_pert, None, T0, t_span[:2], dH, dS, n)
            theta_mean = np.mean(theta_base)
            mi = compute_mi(theta_base, theta_pert)
            dG = gibbs_free_energy(theta_mean, dH, dS, T0)
            return (theta_mean > THR_THETA) and (mi > THR_MI) and (dG < -(R_gas * T0 / 1000) * np.log(2))

        H_on = run_logic(dH_H, dS_H)
        Ca_on = run_logic(dH_Ca, dS_Ca)

        logic = (int(H_on) << 1) | int(Ca_on)
        logic_states.append(logic)
        logic_labels.append(f"{int(H_on)}{int(Ca_on)}")

    return t_span, amp_wave, logic_states, logic_labels

def plot_logic_trace(t, amps, logic_states):
    plt.figure(figsize=(12, 5))
    plt.plot(t, amps, label='Amplitude', color='orange')
    plt.step(t, logic_states, where='mid', label='Logic State (H⁺, Ca²⁺)', color='purple')
    plt.yticks([0, 1, 2, 3], ['00', '01', '10', '11'])
    plt.xlabel("Time (s)")
    plt.title("Dual Ion BFIP Logic Switching (H⁺ + Ca²⁺)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("BFIP_dual_logic_switch.png", dpi=600)
    plt.show()

if __name__ == '__main__':
    # Pick known BFIP-stable points from previous simulations
    dH_H = -46.8421
    dS_H = -0.1021
    dH_Ca = -46.0
    dS_Ca = -0.110

    t, amps, logic_states, logic_labels = simulate_dual_bfip(dH_H, dS_H, dH_Ca, dS_Ca)
    plot_logic_trace(t, amps, logic_states)
