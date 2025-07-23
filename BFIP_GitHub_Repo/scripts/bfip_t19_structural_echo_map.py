
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import pandas as pd

t = np.linspace(0, 10, 1000)
pulse_freq = 1.5

# Ion input (same echo structure)
I_t = np.sin(2 * np.pi * pulse_freq * t) ** 2
I_t[t > 5] += np.sin(2 * np.pi * (pulse_freq + 0.1) * t[t > 5]) ** 2
I_t = np.clip(I_t, 0, 1)

# Oscillatory Subnets (3 branches + merge)
def subnet_phase(t, offset, mod, gain):
    return ((np.sin(2 * np.pi * (1.2 + mod) * t + offset) + 1) / 2) * gain

P1 = subnet_phase(t, offset=0, mod=0.0, gain=0.4)
P2 = subnet_phase(t, offset=1, mod=0.1, gain=0.3)
P3 = subnet_phase(t, offset=2, mod=0.15, gain=0.3)

# Combine with base signal
G1 = I_t * ((P1 + P2 + P3) / 3)
θ1 = np.cumsum(G1) / len(G1)

# Scaffold Echo Memory
Echo = gaussian_filter1d(θ1 * ((np.sin(2 * np.pi * 0.85 * t) + 1) / 2), sigma=3)
θ_echo = np.cumsum(Echo) / len(Echo)

# Final Resonant Phase Merge
G2 = gaussian_filter1d((θ1 + θ_echo) / 2, sigma=2)
θ2 = np.cumsum(G2) / len(G2)

# ΔG and MI
ΔG2 = -np.gradient(θ2)
MI2 = np.correlate((Echo - np.mean(Echo)), (G2 - np.mean(G2)), mode='same') / len(Echo)
BFIP5 = (θ2 > 0.15) & (MI2 > 0.005)

# Output Trigger
G3 = np.zeros_like(t)
for i in range(len(t)):
    if BFIP5[i]:
        G3[i] = θ2[i] * ((np.sin(2 * np.pi * 0.65 * t[i]) + 1) / 2)
θ3 = np.cumsum(G3) / len(G3)

# Plot
plt.figure(figsize=(14, 7))
plt.plot(t, I_t, '--', alpha=0.3, label='Ion Input I(t)')
plt.plot(t, θ1, label='θ̄₁ (Subnet Merge)', linewidth=2)
plt.plot(t, θ_echo, label='θ̄ Echo Memory', linewidth=2)
plt.plot(t, θ2, label='θ̄₂ (Resonant Merge)', linewidth=2)
plt.plot(t, θ3, label='θ̄₃ (Final Output)', linewidth=2)
plt.plot(t, BFIP5 * 0.3, 'g:', label='BFIP₅ Active')
plt.title('Tier 19 — Oscillatory Network Integration (Final Structural Echo Map)')
plt.xlabel('Time (s)')
plt.legend()
plt.tight_layout()
plt.savefig("BFIP_t19_structural_echo_map.png", dpi=300)
plt.show()

# Summary
summary = pd.DataFrame({
    'Final θ̄₁': [θ1[-1]],
    'Final θ̄ Echo': [θ_echo[-1]],
    'Final θ̄₂': [θ2[-1]],
    'Final θ̄₃': [θ3[-1]],
    'Max ΔG₂': [np.max(ΔG2)],
    'Max MI₂ Proxy': [np.max(MI2)],
    'BFIP₅ Triggered?': [BFIP5.any()]
})

print("\n=== Tier 19 Summary ===\n", summary)


def main():
    pass

if __name__ == '__main__':
    main()
