
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import pandas as pd

t = np.linspace(0, 10, 1000)
pulse_freq = 1.5

# Ion input with repeating pattern (simulate memory reappearance)
I_t = np.sin(2 * np.pi * pulse_freq * t) ** 2
I_t[t > 5] += np.sin(2 * np.pi * (pulse_freq + 0.1) * t[t > 5]) ** 2  # introduce a phase recall ripple
I_t = np.clip(I_t, 0, 1)

# Gate 1: Base encoding
G1 = (np.sin(2 * np.pi * pulse_freq * t) + 1) / 2 * I_t
θ1 = np.cumsum(G1) / len(G1)

# Gate 2: Pattern detector + memory storage
def pattern_memory_gate(θ_input, t, threshold=0.2, window=80):
    memory = []
    out = np.zeros_like(t)
    for i in range(window, len(t)):
        avg = np.mean(θ_input[i - window:i])
        if avg > threshold:
            pattern = np.sin(2 * np.pi * 1.1 * t[i])
            out[i] = (pattern + 1) / 2 * avg
            memory.append(avg)
    return gaussian_filter1d(out, sigma=2), memory

G2, mem_store = pattern_memory_gate(θ1, t)
θ2 = np.cumsum(G2) / len(G2)

# Gate 3: Reentrant Oscillation — tries to match previous memory signatures
def reentrant_gate(G_input, memory_trace, t, match_thresh=0.18):
    gate = np.zeros_like(t)
    for i in range(len(t)):
        echo_signal = (np.sin(2 * np.pi * 1.2 * t[i]) + 1) / 2
        if any([abs(m - echo_signal) < match_thresh for m in memory_trace]):
            gate[i] = echo_signal * 0.5
    return gaussian_filter1d(gate, sigma=3)

G3 = reentrant_gate(G2, mem_store, t)
θ3 = np.cumsum(G3) / len(G3)

# BFIP₄ Activation Metrics
ΔG3 = -np.gradient(θ3)
MI3 = np.correlate(G2 - np.mean(G2), G3 - np.mean(G3), mode='same') / len(G2)
BFIP4 = (θ3 > 0.15) & (MI3 > 0.005)

# Plot
plt.figure(figsize=(14, 7))
plt.plot(t, I_t, '--', alpha=0.3, label='Ion Input I(t)')
plt.plot(t, θ1, label='θ̄₁(t)', linewidth=2)
plt.plot(t, θ2, label='θ̄₂(t)', linewidth=2)
plt.plot(t, θ3, label='θ̄₃(t)', linewidth=2)
plt.plot(t, ΔG3, label='ΔG₃(t)', linestyle='--')
plt.plot(t, BFIP4 * 0.3, 'g:', label='BFIP₄ Active')
plt.title('Tier 15 — Ionic Pattern Recognition via Reentrant Oscillation')
plt.xlabel('Time (s)')
plt.legend()
plt.tight_layout()
plt.savefig("BFIP_t15_pattern_recognition.png", dpi=300)
plt.show()

# Summary
summary = pd.DataFrame({
    'Final θ̄₁': [θ1[-1]],
    'Final θ̄₂': [θ2[-1]],
    'Final θ̄₃': [θ3[-1]],
    'Max ΔG₃': [np.max(ΔG3)],
    'Max MI₃ Proxy': [np.max(MI3)],
    'BFIP₄ Triggered?': [BFIP4.any()]
})

print("\n=== Tier 15 Summary ===\n", summary)


def main():
    pass

if __name__ == '__main__':
    main()
