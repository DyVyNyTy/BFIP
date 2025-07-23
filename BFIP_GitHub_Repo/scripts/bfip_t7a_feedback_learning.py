
import numpy as np
import matplotlib.pyplot as plt
from thermo_dynamic_model import run_dynamic_thermo
from models.thermodynamics import gibbs_free_energy

R = 8.314
THR_THETA = 0.10
THR_MI = 2.2
T0 = 300.0

def training_pulse(t, cycle=40, on_width=20):
    return 0.65 if int(t) % cycle < on_width else 0.25

def compute_mi(base, pert):
    hist2d, _, _ = np.histogram2d(base, pert, bins=30)
    pXY = hist2d / np.sum(hist2d)
    pX = np.sum(pXY, axis=1)
    pY = np.sum(pXY, axis=0)
    mask = pXY > 0
    denom = np.outer(pX, pY)
    return np.sum(pXY[mask] * np.log(pXY[mask] / denom[mask]))

def simulate_feedback_learning():
    t_span = np.linspace(0, 240, 480)
    dH = -46.84
    dS = -0.102

    theta_vals = []
    mi_vals = []
    dG_vals = []
    bfip_vals = []

    retained = False
    boost_factor = 1.0

    for t in t_span:
        amp = training_pulse(t)

        kin = {'k_on': 1e5, 'amplitude': amp}
        _, theta = run_dynamic_thermo(kin, None, T0, [t, t + 1], dH, dS, boost_factor)

        kin_pert = {'k_on': 1e5, 'amplitude': amp * 0.9}
        _, theta_pert = run_dynamic_thermo(kin_pert, None, T0, [t, t + 1], dH, dS, boost_factor)

        th_mean = np.mean(theta)
        mi = compute_mi(theta, theta_pert)
        dG = gibbs_free_energy(th_mean, dH, dS, T0)

        passed = th_mean > THR_THETA and mi > THR_MI and dG < -(R * T0 / 1000) * np.log(2)

        if passed:
            retained = True
            boost_factor = 1.25  # simulate learning retention
        elif retained:
            boost_factor *= 0.999  # gradual decay

        theta_vals.append(th_mean)
        mi_vals.append(mi)
        dG_vals.append(dG)
        bfip_vals.append(1 if passed else 0)

    return t_span, np.array(theta_vals), np.array(mi_vals), np.array(dG_vals), np.array(bfip_vals)

def plot_feedback_learning(t, theta, mi, dG, bfip):
    plt.figure(figsize=(12, 10))

    plt.subplot(4, 1, 1)
    plt.plot(t, theta, label='θ̄(t)', color='blue')
    plt.axhline(0.10, linestyle='--', color='gray', label='θ*')
    plt.ylabel('θ̄(t)')
    plt.title('BFIP T7.A: Feedback Learning (Retained Phase after Training)')
    plt.legend()

    plt.subplot(4, 1, 2)
    plt.plot(t, mi, label='MI', color='green')
    plt.axhline(2.2, linestyle='--', color='gray', label='MI Threshold')
    plt.ylabel('Mutual Info (bits)')
    plt.legend()

    plt.subplot(4, 1, 3)
    plt.plot(t, dG, label='ΔG', color='cyan')
    plt.axhline(-(R*T0/1000)*np.log(2), linestyle='--', color='gray', label='ΔG < RTln2')
    plt.ylabel('ΔG (kJ/mol)')
    plt.legend()

    plt.subplot(4, 1, 4)
    plt.step(t, bfip, where='mid', label='BFIP', color='purple')
    plt.ylabel('BFIP')
    plt.xlabel('Time (s)')
    plt.legend()

    plt.tight_layout()
    plt.savefig('BFIP_t7a_feedback_learning.png', dpi=600)
    plt.show()

if __name__ == '__main__':
    t, θ, mi, dG, bfip = simulate_feedback_learning()
    plot_feedback_learning(t, θ, mi, dG, bfip)
