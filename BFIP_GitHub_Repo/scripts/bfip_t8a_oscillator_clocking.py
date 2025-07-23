
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import square
import pandas as pd

t = np.linspace(0, 10, 1000)
I_t = np.sin(2 * np.pi * 1.5 * t) ** 2

O_sine = (np.sin(2 * np.pi * 1.5 * t) + 1) / 2
O_square = (square(2 * np.pi * 1.5 * t) + 1) / 2
O_phase_aligned = (np.sin(2 * np.pi * 1.5 * t + np.pi/4) + 1) / 2

G_sine = I_t * O_sine
G_square = I_t * O_square
G_phase = I_t * O_phase_aligned

def theta_bar(input_signal):
    return np.cumsum(input_signal) / len(input_signal)

def delta_G(input_signal):
    return -np.gradient(input_signal)

def mutual_info_proxy(I, O):
    return np.correlate(I - np.mean(I), O - np.mean(O), mode='same') / len(I)

theta_sine = theta_bar(G_sine)
theta_square = theta_bar(G_square)
theta_phase = theta_bar(G_phase)

dG_sine = delta_G(theta_sine)
dG_square = delta_G(theta_square)
dG_phase = delta_G(theta_phase)

MI_sine = mutual_info_proxy(I_t, O_sine)
MI_square = mutual_info_proxy(I_t, O_square)
MI_phase = mutual_info_proxy(I_t, O_phase_aligned)

BFIP_sine = (theta_sine > 0.15) & (MI_sine > 0.005)
BFIP_square = (theta_square > 0.15) & (MI_square > 0.005)
BFIP_phase = (theta_phase > 0.15) & (MI_phase > 0.005)

fig, axs = plt.subplots(3, 1, figsize=(12, 12), sharex=True)

axs[0].plot(t, I_t, label='Ion Input I(t)', linestyle='--', alpha=0.5)
axs[0].plot(t, G_sine, label='Gated Input (Sine)')
axs[0].plot(t, theta_sine, label='θ̄(t)', linewidth=2)
axs[0].plot(t, dG_sine, label='ΔG(t)')
axs[0].plot(t, BFIP_sine * 0.3, label='BFIP Active', color='green', linestyle='dotted')
axs[0].set_title('Sine-Gated Clock Logic')
axs[0].legend()

axs[1].plot(t, I_t, label='Ion Input I(t)', linestyle='--', alpha=0.5)
axs[1].plot(t, G_square, label='Gated Input (Square)')
axs[1].plot(t, theta_square, label='θ̄(t)', linewidth=2)
axs[1].plot(t, dG_square, label='ΔG(t)')
axs[1].plot(t, BFIP_square * 0.3, label='BFIP Active', color='green', linestyle='dotted')
axs[1].set_title('Square-Gated Clock Logic')
axs[1].legend()

axs[2].plot(t, I_t, label='Ion Input I(t)', linestyle='--', alpha=0.5)
axs[2].plot(t, G_phase, label='Gated Input (Phase-Aligned)')
axs[2].plot(t, theta_phase, label='θ̄(t)', linewidth=2)
axs[2].plot(t, dG_phase, label='ΔG(t)')
axs[2].plot(t, BFIP_phase * 0.3, label='BFIP Active', color='green', linestyle='dotted')
axs[2].set_title('Phase-Aligned Gated Clock Logic')
axs[2].legend()

plt.xlabel('Time (s)')
plt.tight_layout()
plt.show()

summary = pd.DataFrame({
    'Gate Type': ['Sine', 'Square', 'Phase-Aligned'],
    'Final θ̄': [theta_sine[-1], theta_square[-1], theta_phase[-1]],
    'Max ΔG': [np.max(dG_sine), np.max(dG_square), np.max(dG_phase)],
    'Max MI Proxy': [np.max(MI_sine), np.max(MI_square), np.max(MI_phase)],
    'BFIP Triggered?': [BFIP_sine.any(), BFIP_square.any(), BFIP_phase.any()]
})

print("\n=== Tier 8A Summary ===\n", summary)


def main():
    pass

if __name__ == '__main__':
    main()
