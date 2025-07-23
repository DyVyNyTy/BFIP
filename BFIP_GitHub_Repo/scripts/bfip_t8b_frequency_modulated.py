
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import square
import pandas as pd

# Time base
t = np.linspace(0, 10, 1000)

# Ion pulse input (same Tier 7 carryover)
I_t = np.sin(2 * np.pi * 1.5 * t) ** 2

# Adaptive oscillator function: starts mismatched, adjusts every X seconds based on MI
def adaptive_gate(t, input_signal, window=100):
    freq = 1.0  # Start at 1 Hz, lower than pulse freq
    gate = np.zeros_like(t)
    current_phase = 0
    for i in range(0, len(t), window):
        segment = input_signal[i:i+window]
        t_segment = t[i:i+window]
        if len(segment) == 0: break
        phase_shift = np.pi / 4  # small bias
        osc = (np.sin(2 * np.pi * freq * t_segment + phase_shift) + 1) / 2
        gate[i:i+window] = osc

        # Update frequency using MI heuristic (cross-correlation proxy)
        corr = np.correlate(segment - np.mean(segment), osc - np.mean(osc), mode='valid')
        if corr.size > 0 and corr[0] > 0.005:
            freq += 0.05  # tune up slightly
        else:
            freq -= 0.025  # back off slightly
        freq = np.clip(freq, 0.5, 3.0)  # avoid runaway
    return gate

# Generate adaptive gate
O_adaptive = adaptive_gate(t, I_t)

# Apply gate
G_adaptive = I_t * O_adaptive

# θ̄(t), ΔG(t), MI(t)
def theta_bar(input_signal):
    return np.cumsum(input_signal) / len(input_signal)

def delta_G(input_signal):
    return -np.gradient(input_signal)

def mutual_info_proxy(I, O):
    return np.correlate(I - np.mean(I), O - np.mean(O), mode='same') / len(I)

theta_adaptive = theta_bar(G_adaptive)
dG_adaptive = delta_G(theta_adaptive)
MI_adaptive = mutual_info_proxy(I_t, O_adaptive)

# BFIP Activation Rule
BFIP_adaptive = (theta_adaptive > 0.15) & (MI_adaptive > 0.005)

# Plot
plt.figure(figsize=(12, 6))
plt.plot(t, I_t, label='Ion Input I(t)', linestyle='--', alpha=0.5)
plt.plot(t, G_adaptive, label='Gated Input (Adaptive)')
plt.plot(t, theta_adaptive, label='θ̄(t)', linewidth=2)
plt.plot(t, dG_adaptive, label='ΔG(t)')
plt.plot(t, BFIP_adaptive * 0.3, label='BFIP Active', color='green', linestyle='dotted')
plt.title('Adaptive Frequency-Modulated Gate Logic')
plt.xlabel('Time (s)')
plt.legend()
plt.tight_layout()
plt.savefig("BFIP_t8b_frequency_modulated.png", dpi=300)
plt.show()

# Summary
summary = pd.DataFrame({
    'Final θ̄': [theta_adaptive[-1]],
    'Max ΔG': [np.max(dG_adaptive)],
    'Max MI Proxy': [np.max(MI_adaptive)],
    'BFIP Triggered?': [BFIP_adaptive.any()]
})

print("\n=== Tier 8B Summary ===\n", summary)


def main():
    pass

if __name__ == '__main__':
    main()
