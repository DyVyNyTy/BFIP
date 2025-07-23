
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import pandas as pd

t = np.linspace(0, 10, 1000)

# Ion input signal
I_t = np.sin(2 * np.pi * 1.5 * t) ** 2

# Gate 1: Standard θ̄₁ logic
G1 = (np.sin(2 * np.pi * 1.5 * t) + 1) / 2 * I_t
θ1 = np.cumsum(G1) / len(G1)

# Gate 2: Encoded with memory logic (Tier 12)
def encoder_gate(θ_input, t, window=100, encode_threshold=0.18):
    gate = np.zeros_like(t)
    for i in range(window, len(t)):
        avg = np.mean(θ_input[i - window:i])
        if avg > encode_threshold:
            gate[i] = (np.sin(2 * np.pi * 1.2 * t[i]) + 1) / 2 * avg
    return gaussian_filter1d(gate, sigma=2)

G2 = encoder_gate(θ1, t)
θ2 = np.cumsum(G2) / len(G2)

# Tier 13: G2 enters a feedback rehearsal loop before triggering G3
def rehearsal_gate(G_encoded, t, memory_feedback_gain=0.5, window=50):
    rehearsal = np.zeros_like(t)
    for i in range(window, len(t)):
        avg_mem = np.mean(G_encoded[i - window:i])
        if avg_mem > 0.01:
            loop_signal = (np.sin(2 * np.pi * 1.1 * t[i]) + 1) / 2
            rehearsal[i] = loop_signal * avg_mem * memory_feedback_gain
    return gaussian_filter1d(rehearsal, sigma=3)

G2_rehearsed = rehearsal_gate(G2, t)
θ2_rehearsed = np.cumsum(G2_rehearsed) / len(G2_rehearsed)

# Gate 3: Triggers only if rehearsal crosses coherence
def final_gate(G_looped, t, threshold=0.01):
    gate = np.zeros_like(t)
    for i in range(len(t)):
        if G_looped[i] > threshold:
            gate[i] = G_looped[i] * ((np.sin(2 * np.pi * 1.0 * t[i]) + 1) / 2)
    return gaussian_filter1d(gate, sigma=2)

G3 = final_gate(G2_rehearsed, t)
θ3 = np.cumsum(G3) / len(G3)

# ΔG and MI
ΔG3 = -np.gradient(θ3)
MI3 = np.correlate(G2_rehearsed - np.mean(G2_rehearsed), G3 - np.mean(G3), mode='same') / len(G2_rehearsed)
BFIP3 = (θ3 > 0.15) & (MI3 > 0.005)

# Plot
plt.figure(figsize=(14, 7))
plt.plot(t, I_t, '--', alpha=0.3, label='Ion Input I(t)')
plt.plot(t, θ1, label='θ̄₁(t)', linewidth=2)
plt.plot(t, θ2, label='θ̄₂(t)', linewidth=2)
plt.plot(t, θ2_rehearsed, label='θ̄₂ Rehearsed(t)', linewidth=2)
plt.plot(t, θ3, label='θ̄₃(t)', linewidth=2)
plt.plot(t, ΔG3, label='ΔG₃(t)', linestyle='--')
plt.plot(t, BFIP3 * 0.3, 'g:', label='BFIP₃ Active')
plt.title('Tier 13 — Memory Rehearsal & Resonant Recall')
plt.xlabel('Time (s)')
plt.legend()
plt.tight_layout()
plt.savefig("BFIP_t13_rehearsal_recall.png", dpi=300)
plt.show()

# Summary
summary = pd.DataFrame({
    'Final θ̄₁': [θ1[-1]],
    'Final θ̄₂': [θ2[-1]],
    'Final θ̄₂ Rehearsed': [θ2_rehearsed[-1]],
    'Final θ̄₃': [θ3[-1]],
    'Max ΔG₃': [np.max(ΔG3)],
    'Max MI₃ Proxy': [np.max(MI3)],
    'BFIP₃ Triggered?': [BFIP3.any()]
})

print("\n=== Tier 13 Summary ===\n", summary)


def main():
    pass

if __name__ == '__main__':
    main()
