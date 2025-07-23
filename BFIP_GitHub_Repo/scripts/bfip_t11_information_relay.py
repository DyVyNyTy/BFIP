
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import pandas as pd

t = np.linspace(0, 10, 1000)

# Ion input
I_t = np.sin(2 * np.pi * 1.5 * t) ** 2

# Gate 1: Initial processing
G1 = (np.sin(2 * np.pi * 1.5 * t) + 1) / 2 * I_t
θ1 = np.cumsum(G1) / len(G1)

# Gate 2: Cascaded from θ̄₁
def relay_encode_gate(θ_input, t, base_freq=1.6, threshold=0.1):
    gate = np.zeros_like(t)
    freq_mod = base_freq
    for i in range(100, len(t)):
        avg_theta = np.mean(θ_input[i-100:i])
        # Modulate based on θ̄₁ signal strength
        if avg_theta > threshold:
            freq_mod = base_freq + (avg_theta - threshold) * 5  # encode higher frequency
            signal = (np.sin(2 * np.pi * freq_mod * t[i]) + 1) / 2
            gate[i] = signal * avg_theta
    return gaussian_filter1d(gate, sigma=2)

G2 = relay_encode_gate(θ1, t)
θ2 = np.cumsum(G2) / len(G2)

# Gate 3: Receives from G2 output
G3 = gaussian_filter1d(G2, sigma=4) * ((np.sin(2 * np.pi * 1.3 * t) + 1) / 2)
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
plt.title('Tier 11 — Information Relay Encoding')
plt.xlabel('Time (s)')
plt.legend()
plt.tight_layout()
plt.savefig("BFIP_t11_information_relay.png", dpi=300)
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

print("\n=== Tier 11 Summary ===\n", summary)


def main():
    pass

if __name__ == '__main__':
    main()
