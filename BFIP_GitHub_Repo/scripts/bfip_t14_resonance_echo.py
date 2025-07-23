
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import pandas as pd

t = np.linspace(0, 10, 1000)

# Ion input
I_t = np.sin(2 * np.pi * 1.5 * t) ** 2
G1 = (np.sin(2 * np.pi * 1.5 * t) + 1) / 2 * I_t
θ1 = np.cumsum(G1) / len(G1)

# Tier 12 logic preserved
def encoder_gate(θ_input, t, window=100, encode_threshold=0.18):
    gate = np.zeros_like(t)
    for i in range(window, len(t)):
        avg = np.mean(θ_input[i - window:i])
        if avg > encode_threshold:
            gate[i] = (np.sin(2 * np.pi * 1.2 * t[i]) + 1) / 2 * avg
    return gaussian_filter1d(gate, sigma=2)

G2 = encoder_gate(θ1, t)
θ2 = np.cumsum(G2) / len(G2)

# Echo loop to reinforce weak θ̄₂ memory
def echo_loop(G_input, t, echo_gain=0.8, echo_window=50, echo_trigger=0.01):
    echo = np.zeros_like(t)
    for i in range(echo_window, len(t)):
        recent = G_input[i - echo_window:i]
        if np.mean(recent) > echo_trigger:
            phase = np.sin(2 * np.pi * 1.0 * t[i])
            echo[i] = echo_gain * np.mean(recent) * ((phase + 1) / 2)
    return gaussian_filter1d(echo, sigma=3)

G_echo = echo_loop(G2, t)
θ_echo = np.cumsum(G_echo) / len(G_echo)

# G3 receives both direct G2 and echo G2
G3 = gaussian_filter1d(G2 + G_echo, sigma=2)
θ3 = np.cumsum(G3) / len(G3)

# ΔG and MI
ΔG3 = -np.gradient(θ3)
MI3 = np.correlate((G2 + G_echo) - np.mean(G2 + G_echo), G3 - np.mean(G3), mode='same') / len(G3)
BFIP3 = (θ3 > 0.15) & (MI3 > 0.005)

# Plot
plt.figure(figsize=(14, 7))
plt.plot(t, I_t, '--', alpha=0.3, label='Ion Input I(t)')
plt.plot(t, θ1, label='θ̄₁(t)', linewidth=2)
plt.plot(t, θ2, label='θ̄₂(t)', linewidth=2)
plt.plot(t, θ_echo, label='θ̄ Echo(t)', linewidth=2)
plt.plot(t, θ3, label='θ̄₃(t)', linewidth=2)
plt.plot(t, ΔG3, label='ΔG₃(t)', linestyle='--')
plt.plot(t, BFIP3 * 0.3, 'g:', label='BFIP₃ Active')
plt.title('Tier 14 — Resonance Lock & Echo Memory')
plt.xlabel('Time (s)')
plt.legend()
plt.tight_layout()
plt.savefig("BFIP_t14_resonance_echo.png", dpi=300)
plt.show()

# Summary
summary = pd.DataFrame({
    'Final θ̄₁': [θ1[-1]],
    'Final θ̄₂': [θ2[-1]],
    'Final θ̄ Echo': [θ_echo[-1]],
    'Final θ̄₃': [θ3[-1]],
    'Max ΔG₃': [np.max(ΔG3)],
    'Max MI₃ Proxy': [np.max(MI3)],
    'BFIP₃ Triggered?': [BFIP3.any()]
})

print("\n=== Tier 14 Summary ===\n", summary)


def main():
    pass

if __name__ == '__main__':
    main()
