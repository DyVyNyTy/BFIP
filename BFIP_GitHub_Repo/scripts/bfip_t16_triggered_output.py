
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import pandas as pd

t = np.linspace(0, 10, 1000)
pulse_freq = 1.5

# Input with recall region reintroduced
I_t = np.sin(2 * np.pi * pulse_freq * t) ** 2
I_t[t > 5] += np.sin(2 * np.pi * (pulse_freq + 0.1) * t[t > 5]) ** 2
I_t = np.clip(I_t, 0, 1)

# Gate 1: Standard ionic integration
G1 = (np.sin(2 * np.pi * pulse_freq * t) + 1) / 2 * I_t
θ1 = np.cumsum(G1) / len(G1)

# Gate 2: Encoded pattern detection
def encode_gate(θ_input, t, threshold=0.2, window=100):
    output = np.zeros_like(t)
    for i in range(window, len(t)):
        mem_avg = np.mean(θ_input[i - window:i])
        if mem_avg > threshold:
            signal = (np.sin(2 * np.pi * 1.2 * t[i]) + 1) / 2
            output[i] = signal * mem_avg
    return gaussian_filter1d(output, sigma=2)

G2 = encode_gate(θ1, t)
θ2 = np.cumsum(G2) / len(G2)

# Gate 3: Reentrant pattern recognition gate
def reentrant_recognition_gate(G_encoded, mem_threshold=0.015):
    return gaussian_filter1d(np.where(G_encoded > mem_threshold, G_encoded, 0), sigma=2)

G3 = reentrant_recognition_gate(G2)
θ3 = np.cumsum(G3) / len(G3)

# BFIP₄ Determination
ΔG3 = -np.gradient(θ3)
MI3 = np.correlate(G2 - np.mean(G2), G3 - np.mean(G3), mode='same') / len(G2)
BFIP4 = (θ3 > 0.15) & (MI3 > 0.005)

# Gate 4: Functional output triggered only when BFIP₄ holds
G4 = np.zeros_like(t)
for i in range(len(t)):
    if BFIP4[i]:
        G4[i] = (np.sin(2 * np.pi * 0.8 * t[i]) + 1) / 2 * θ3[i]

θ4 = np.cumsum(G4) / len(G4)

# Plot
plt.figure(figsize=(14, 7))
plt.plot(t, I_t, '--', alpha=0.3, label='Ion Input I(t)')
plt.plot(t, θ1, label='θ̄₁(t)', linewidth=2)
plt.plot(t, θ2, label='θ̄₂(t)', linewidth=2)
plt.plot(t, θ3, label='θ̄₃(t)', linewidth=2)
plt.plot(t, θ4, label='θ̄₄ (Triggered Output)', linewidth=2)
plt.plot(t, BFIP4 * 0.3, 'g:', label='BFIP₄ Active')
plt.title('Tier 16 — Functional Recall with Triggered Output')
plt.xlabel('Time (s)')
plt.legend()
plt.tight_layout()
plt.savefig("BFIP_t16_triggered_output.png", dpi=300)
plt.show()

# Summary
summary = pd.DataFrame({
    'Final θ̄₁': [θ1[-1]],
    'Final θ̄₂': [θ2[-1]],
    'Final θ̄₃': [θ3[-1]],
    'Final θ̄₄': [θ4[-1]],
    'Max ΔG₃': [np.max(ΔG3)],
    'Max MI₃ Proxy': [np.max(MI3)],
    'BFIP₄ Triggered?': [BFIP4.any()]
})

print("\n=== Tier 16 Summary ===\n", summary)


def main():
    pass

if __name__ == '__main__':
    main()
