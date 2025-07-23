
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import pandas as pd

t = np.linspace(0, 10, 1000)
pulse_freq = 1.5

# Ion input
I_t = np.sin(2 * np.pi * pulse_freq * t) ** 2
I_t[t > 5] += np.sin(2 * np.pi * (pulse_freq + 0.1) * t[t > 5]) ** 2
I_t = np.clip(I_t, 0, 1)

# Gate 1: Base logic
G1 = (np.sin(2 * np.pi * pulse_freq * t) + 1) / 2 * I_t
θ1 = np.cumsum(G1) / len(G1)

# Gate 2A and 2B: Two parallel encoding pathways
def encode_pathway(signal, t, shift, gain):
    return gaussian_filter1d(((np.sin(2 * np.pi * (1.2 + shift) * t) + 1) / 2) * signal * gain, sigma=2)

G2a = encode_pathway(θ1, t, shift=0.0, gain=0.7)
G2b = encode_pathway(θ1, t, shift=0.2, gain=0.5)
θ2a = np.cumsum(G2a) / len(G2a)
θ2b = np.cumsum(G2b) / len(G2b)

# Gate 3: Scaffold merge of parallel memory pathways
G3 = gaussian_filter1d((G2a + G2b) / 2, sigma=2)
θ3 = np.cumsum(G3) / len(G3)

# ΔG and MI
ΔG3 = -np.gradient(θ3)
MI3 = np.correlate((G2a + G2b)/2 - np.mean((G2a + G2b)/2), G3 - np.mean(G3), mode='same') / len(G3)
BFIP4 = (θ3 > 0.15) & (MI3 > 0.005)

# Gate 4: Functional response from scaffolded memory state
G4 = np.zeros_like(t)
for i in range(len(t)):
    if BFIP4[i]:
        G4[i] = (np.sin(2 * np.pi * 0.8 * t[i]) + 1) / 2 * θ3[i]

θ4 = np.cumsum(G4) / len(G4)

# Plot
plt.figure(figsize=(14, 7))
plt.plot(t, I_t, '--', alpha=0.3, label='Ion Input I(t)')
plt.plot(t, θ1, label='θ̄₁ (Base)', linewidth=2)
plt.plot(t, θ2a, label='θ̄₂A (Channel A)', linewidth=2)
plt.plot(t, θ2b, label='θ̄₂B (Channel B)', linewidth=2)
plt.plot(t, θ3, label='θ̄₃ (Scaffold Merge)', linewidth=2)
plt.plot(t, θ4, label='θ̄₄ (Triggered Output)', linewidth=2)
plt.plot(t, BFIP4 * 0.3, 'g:', label='BFIP₄ Active')
plt.title('Tier 18 — Multi-Channel Scaffolded Memory (BFIP-Networked Architecture)')
plt.xlabel('Time (s)')
plt.legend()
plt.tight_layout()
plt.savefig("BFIP_t18_scaffolded_memory.png", dpi=300)
plt.show()

# Summary
summary = pd.DataFrame({
    'Final θ̄₁': [θ1[-1]],
    'Final θ̄₂A': [θ2a[-1]],
    'Final θ̄₂B': [θ2b[-1]],
    'Final θ̄₃': [θ3[-1]],
    'Final θ̄₄': [θ4[-1]],
    'Max ΔG₃': [np.max(ΔG3)],
    'Max MI₃ Proxy': [np.max(MI3)],
    'BFIP₄ Triggered?': [BFIP4.any()]
})

print("\n=== Tier 18 Summary ===\n", summary)


def main():
    pass

if __name__ == '__main__':
    main()
