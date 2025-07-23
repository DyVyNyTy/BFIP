
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import pandas as pd

t = np.linspace(0, 10, 1000)

# Ion input
I_t = np.sin(2 * np.pi * 1.5 * t) ** 2

# Gate 1: Normal ionic logic
G1 = (np.sin(2 * np.pi * 1.5 * t) + 1) / 2 * I_t
θ1 = np.cumsum(G1) / len(G1)

# Gate 2: Encodes function from θ̄₁ pattern
def encoder_gate(θ_input, t, window=100, encode_threshold=0.18):
    pattern_memory = []
    gate = np.zeros_like(t)
    for i in range(window, len(t)):
        pattern = θ_input[i - window:i]
        avg = np.mean(pattern)
        pattern_memory.append(avg)

        # Fire if a rising encoded structure is recognized
        if avg > encode_threshold:
            gate[i] = (np.sin(2 * np.pi * 1.2 * t[i]) + 1) / 2 * avg
    return gaussian_filter1d(gate, sigma=2)

G2 = encoder_gate(θ1, t)
θ2 = np.cumsum(G2) / len(G2)

# Gate 3: Receives from encoder and applies logic filter
def logic_filter_gate(G_encoded, t, logic_threshold=0.015):
    logic_gate = np.zeros_like(t)
    for i in range(len(t)):
        signal = G_encoded[i]
        if signal > logic_threshold:
            logic_gate[i] = signal * ((np.sin(2 * np.pi * 1.0 * t[i]) + 1) / 2)
    return gaussian_filter1d(logic_gate, sigma=2)

G3 = logic_filter_gate(G2, t)
θ3 = np.cumsum(G3) / len(G3)

# ΔG and MI
ΔG3 = -np.gradient(θ3)
MI3 = np.correlate(G2 - np.mean(G2), G3 - np.mean(G3), mode='same') / len(G2)
BFIP3 = (θ3 > 0.15) & (MI3 > 0.005)

# Plot
plt.figure(figsize=(14, 7))
plt.plot(t, I_t, '--', alpha=0.3, label='Ion Input I(t)')
plt.plot(t, θ1, label='θ̄₁(t)', linewidth=2)
plt.plot(t, θ2, label='θ̄₂(t)', linewidth=2)
plt.plot(t, θ3, label='θ̄₃(t)', linewidth=2)
plt.plot(t, ΔG3, label='ΔG₃(t)', linestyle='--')
plt.plot(t, BFIP3 * 0.3, 'g:', label='BFIP₃ Active')
plt.title('Tier 12 — Functional Encoding Integration')
plt.xlabel('Time (s)')
plt.legend()
plt.tight_layout()
plt.savefig("BFIP_t12_functional_encoding.png", dpi=300)
plt.show()

# Summary
summary = pd.DataFrame({
    'Final θ̄₁': [θ1[-1]],
    'Final θ̄₂': [θ2[-1]],
    'Final θ̄₃': [θ3[-1]],
    'Max ΔG₃': [np.max(ΔG3)],
    'Max MI₃ Proxy': [np.max(MI3)],
    'BFIP₃ Triggered?': [BFIP3.any()]
})

print("\n=== Tier 12 Summary ===\n", summary)


def main():
    pass

if __name__ == '__main__':
    main()
