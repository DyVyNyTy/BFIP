
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import pandas as pd

# Time base
t = np.linspace(0, 10, 1000)

# Ion input
I_t = np.sin(2 * np.pi * 1.5 * t) ** 2

# Gating Layer
def base_gate(input_signal, freq, phase_shift):
    return (np.sin(2 * np.pi * freq * t + phase_shift) + 1) / 2 * input_signal

# Cascade Logic Layer
def cascade_gate(previous_theta, freq=1.5, phase_shift=0.0, memory_window=50):
    gate = np.zeros_like(previous_theta)
    for i in range(memory_window, len(previous_theta)):
        mem_avg = np.mean(previous_theta[i - memory_window:i])
        signal = (np.sin(2 * np.pi * freq * t[i] + phase_shift) + 1) / 2
        gate[i] = signal * mem_avg
    return gaussian_filter1d(gate, sigma=2)

# First gate: base ion gating
G1 = base_gate(I_t, freq=1.5, phase_shift=0)
θ1 = np.cumsum(G1) / len(G1)
ΔG1 = -np.gradient(θ1)
MI1 = np.correlate(I_t - np.mean(I_t), G1 - np.mean(G1), mode='same') / len(I_t)
BFIP1 = (θ1 > 0.15) & (MI1 > 0.005)

# Second gate: cascaded information from θ1
G2 = cascade_gate(θ1, freq=1.6, phase_shift=np.pi/4)
θ2 = np.cumsum(G2) / len(G2)
ΔG2 = -np.gradient(θ2)
MI2 = np.correlate(G1 - np.mean(G1), G2 - np.mean(G2), mode='same') / len(G1)
BFIP2 = (θ2 > 0.15) & (MI2 > 0.005)

# Plot
plt.figure(figsize=(14, 7))

plt.plot(t, I_t, '--', alpha=0.4, label='Ion Input I(t)')
plt.plot(t, G1, label='Gate 1 Output')
plt.plot(t, G2, label='Gate 2 Output (Cascaded)')
plt.plot(t, θ1, label='θ̄₁(t)', linewidth=2)
plt.plot(t, θ2, label='θ̄₂(t)', linewidth=2)
plt.plot(t, ΔG2, label='ΔG₂(t)', color='brown', linestyle='--')
plt.plot(t, BFIP2 * 0.3, 'g:', label='BFIP₂ Active')

plt.title('Tier 10 — Information Cascade Logic')
plt.xlabel('Time (s)')
plt.legend()
plt.tight_layout()
plt.savefig("BFIP_t10_information_cascade.png", dpi=300)
plt.show()

# Summary
summary = pd.DataFrame({
    'Final θ̄₁': [θ1[-1]],
    'Final θ̄₂': [θ2[-1]],
    'Max ΔG₂': [np.max(ΔG2)],
    'Max MI₂ Proxy': [np.max(MI2)],
    'BFIP₂ Triggered?': [BFIP2.any()]
})

print("\n=== Tier 10 Summary ===\n", summary)


def main():
    pass

if __name__ == '__main__':
    main()
