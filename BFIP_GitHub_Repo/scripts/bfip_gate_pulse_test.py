
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

def simulate_pulse_amplitude(dH, dS, T0=300.0):
    t_span = np.linspace(0, 120.0, 300)
    amp_wave = generate_pulse_wave(t_span)

    theta_all = []
    bfip_status = []
    mi_all = []
    amps = []

    for idx, amp in enumerate(amp_wave):
        kin_base = {'k_on': 1e5, 'amplitude': amp}
        n_H = 1.0
        t_local, theta_base = run_dynamic_thermo(kin_base, None, T0, t_span[:2], dH, dS, n_H)

        kin_pert = {'k_on': kin_base['k_on'], 'amplitude': amp * (1 - MI_PERTURB)}
        _, theta_pert = run_dynamic_thermo(kin_pert, None, T0, t_span[:2], dH, dS, n_H)

        theta_mean = np.mean(theta_base)
        mi = compute_mi(theta_base, theta_pert)
        dG = gibbs_free_energy(theta_mean, dH, dS, T0)

        bfip = (theta_mean > THR_THETA) and (mi > THR_MI) and (dG < -(R_gas * T0 / 1000) * np.log(2))

        theta_all.append(theta_mean)
        mi_all.append(mi)
        bfip_status.append(1 if bfip else 0)
        amps.append(amp)

    return t_span, amps, theta_all, mi_all, bfip_status

def plot_pulse_response(t, amps, thetas, mis, bfip_status):
    fig, ax1 = plt.subplots(figsize=(12, 6))

    ax1.plot(t, amps, label="Amplitude", color='orange')
    ax1.set_ylabel("Amplitude", color='orange')
    ax1.tick_params(axis='y', labelcolor='orange')

    ax2 = ax1.twinx()
    ax2.plot(t, thetas, label="θ̄", color='blue')
    ax2.plot(t, mis, label="MI", color='green')
    ax2.plot(t, bfip_status, label="BFIP (binary)", color='black', linestyle='--')
    ax2.set_ylabel("θ̄ / MI / BFIP")
    ax2.legend(loc='upper right')

    plt.title("BFIP Gate Response to Pulsed Amplitude")
    fig.tight_layout()
    plt.savefig("BFIP_pulse_response.png", dpi=600)
    plt.show()

if __name__ == '__main__':
    dH = -46.8421
    dS = -0.1021
    t, amps, thetas, mis, bfip_status = simulate_pulse_amplitude(dH, dS)
    plot_pulse_response(t, amps, thetas, mis, bfip_status)
