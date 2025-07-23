
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import pandas as pd

# Time vector
t = np.linspace(0, 10, 1000)

# Ion pulse input (Tier 7 model)
I_t = np.sin(2 * np.pi * 1.5 * t) ** 2

# Phase Lock System:
# Simulates a gate that adjusts in real-time to stabilize coherence after BFIP activation
def phase_lock_gate(t, input_signal, lock_threshold=0.005, memory_window=100):
    gate = np.zeros_like(t)
    phase_state = 0.0
    freq = 1.5
    locked = False

    for i in range(memory_window, len(t)):
        # Predict using memory
        memory = input_signal[i-memory_window:i]
        mean_mem = np.mean(memory)
        signal = (np.sin(2 * np.pi * freq * t[i] + phase_state) + 1) / 2

        # Gate opens based on memory resonance
        if mean_mem > 0.2:
            gate[i] = signal

        # Locking logic
        if not locked and mean_mem > 0.25:
            locked = True
            phase_state += 0.1  # Shift to stabilize
        elif locked:
            phase_state += 0.02  # Sustain lock
        else:
            phase_state -= 0.01  # Drift if not stable

        freq = np.clip(freq, 1.0, 2.5)
        phase_state = np.clip(phase_state, 0, 2 * np.pi)

    return gaussian_filter1d(gate, sigma=2)

# Apply phase lock
O_lock = phase_lock_gate(t, I_t)
G_lock = I_t * O_lock

# θ̄, ΔG, MI
def theta_bar(input_signal):
    return np.cumsum(input_signal) / len(input_signal)

def delta_G(input_signal):
    return -np.gradient(input_signal)

def mutual_info_proxy(I, O):
    return np.correlate(I - np.mean(I), O - np.mean(O), mode='same') / len(I)

theta = theta_bar(G_lock)
dG = delta_G(theta)
MI_proxy = mutual_info_proxy(I_t, O_lock)
BFIP = (theta > 0.15) & (MI_proxy > 0.005)

# Plot
plt.figure(figsize=(12, 6))
plt.plot(t, I_t, label='Ion Input I(t)', linestyle='--', alpha=0.5)
plt.plot(t, G_lock, label='Phase-Locked Gated Input')
plt.plot(t, theta, label='θ̄(t)', linewidth=2)
plt.plot(t, dG, label='ΔG(t)')
plt.plot(t, BFIP * 0.3, label='BFIP Active', color='green', linestyle='dotted')
plt.title('Tier 9 — Phase Lock & Functional Integration')
plt.xlabel('Time (s)')
plt.legend()
plt.tight_layout()
plt.savefig("BFIP_t9_phase_lock.png", dpi=300)
plt.show()

# Summary
summary = pd.DataFrame({
    'Final θ̄': [theta[-1]],
    'Max ΔG': [np.max(dG)],
    'Max MI Proxy': [np.max(MI_proxy)],
    'BFIP Triggered?': [BFIP.any()]
})

print("\n=== Tier 9 Summary ===\n", summary)


def main():
    pass

if __name__ == '__main__':
    main()
