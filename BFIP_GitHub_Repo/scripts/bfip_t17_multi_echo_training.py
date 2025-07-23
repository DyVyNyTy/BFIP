
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import pandas as pd

t = np.linspace(0, 10, 1000)
pulse_freq = 1.5

# Ion input: repeating signature with phase shift
I_t = np.sin(2 * np.pi * pulse_freq * t) ** 2
I_t[t > 5] += np.sin(2 * np.pi * (pulse_freq + 0.1) * t[t > 5]) ** 2
I_t = np.clip(I_t, 0, 1)

# Gate 1: Base response
G1 = (np.sin(2 * np.pi * pulse_freq * t) + 1) / 2 * I_t
θ1 = np.cumsum(G1) / len(G1)

# Gate 2: Encode and store pattern
def encode_gate(θ_input, t, threshold=0.2, window=100):
    encoded = np.zeros_like(t)
    for i in range(window, len(t)):
        avg = np.mean(θ_input[i - window:i])
        if avg > threshold:
            encoded[i] = (np.sin(2 * np.pi * 1.2 * t[i]) + 1) / 2 * avg
    return gaussian_filter1d(encoded, sigma=2)

G2 = encode_gate(θ1, t)
θ2 = np.cumsum(G2) / len(G2)

# Gate 3: Echo training from both θ̄₁ and θ̄₂
def echo_train_gate(t, θ1, θ2, training_gain=0.4, window=80):
    echo = np.zeros_like(t)
    for i in range(window, len(t)):
        echo_avg = (np.mean(θ1[i - window:i]) + np.mean(θ2[i - window:i])) / 2
        if echo_avg > 0.01:
            echo[i] = (np.sin(2 * np.pi * 1.1 * t[i]) + 1) / 2 * echo_avg * training_gain
    return gaussian_filter1d(echo, sigma=2)

G3 = echo_train_gate(t, θ1, θ2)
θ3 = np.cumsum(G3) / len(G3)

# Gate 4: Output fires if BFIP₄ thresholds are achieved after training
ΔG3 = -np.gradient(θ3)
MI3 = np.correlate(G2 - np.mean(G2), G3 - np.mean(G3), mode='same') / len(G3)
BFIP4 = (θ3 > 0.15) & (MI3 > 0.005)

G4 = np.zeros_like(t)
for i in range(len(t)):
    if BFIP4[i]:
        G4[i] = ((np.sin(2 * np.pi * 0.8 * t[i]) + 1) / 2) * θ3[i]

θ4 = np.cumsum(G4) / len(G4)

# Plot
plt.figure(figsize=(14, 7))
plt.plot(t, I_t, '--', alpha=0.3, label='Ion Input I(t)')
plt.plot(t, θ1, label='θ̄₁(t)', linewidth=2)
plt.plot(t, θ2, label='θ̄₂(t)', linewidth=2)
plt.plot(t, θ3, label='θ̄₃ (Echo Trained)', linewidth=2)
plt.plot(t, θ4, label='θ̄₄ (Triggered Output)', linewidth=2)
plt.plot(t, BFIP4 * 0.3, 'g:', label='BFIP₄ Active')
plt.title('Tier 17 — Memory Amplification via Multi-Echo Phase Training')
plt.xlabel('Time (s)')
plt.legend()
plt.tight_layout()
plt.savefig("BFIP_t17_multi_echo_training.png", dpi=300)
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

print("\n=== Tier 17 Summary ===\n", summary)


def main():
    pass

if __name__ == '__main__':
    main()
