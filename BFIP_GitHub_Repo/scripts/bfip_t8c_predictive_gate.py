
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import square
import pandas as pd

# Time base
t = np.linspace(0, 10, 1000)
dt = t[1] - t[0]

# Ion input: same repeating pulse structure
I_t = np.sin(2 * np.pi * 1.5 * t) ** 2

# Memory-based predictor gate: looks back at pulse history to pre-activate
def predictive_gate(t, input_signal, memory_window=50, prediction_horizon=25):
    gate = np.zeros_like(t)
    pred_pattern = np.zeros(memory_window)

    for i in range(memory_window, len(t) - prediction_horizon):
        # Capture past window
        past_window = input_signal[i - memory_window:i]
        # Predict a future spike if past has strong structure
        if np.sum(past_window) > memory_window * 0.25:  # threshold memory pattern
            gate[i + prediction_horizon] = 1.0  # fire preemptively

    # Smooth the gate
    from scipy.ndimage import gaussian_filter1d
    return gaussian_filter1d(gate, sigma=5)

# Build gate
O_predict = predictive_gate(t, I_t)

# Apply to input
G_predict = I_t * O_predict

# θ̄, ΔG, MI
def theta_bar(input_signal):
    return np.cumsum(input_signal) / len(input_signal)

def delta_G(input_signal):
    return -np.gradient(input_signal)

def mutual_info_proxy(I, O):
    return np.correlate(I - np.mean(I), O - np.mean(O), mode='same') / len(I)

theta_predict = theta_bar(G_predict)
dG_predict = delta_G(theta_predict)
MI_predict = mutual_info_proxy(I_t, O_predict)

# BFIP thresholding
BFIP_predict = (theta_predict > 0.15) & (MI_predict > 0.005)

# Plot
plt.figure(figsize=(12, 6))
plt.plot(t, I_t, label='Ion Input I(t)', linestyle='--', alpha=0.5)
plt.plot(t, G_predict, label='Gated Input (Predictive)')
plt.plot(t, theta_predict, label='θ̄(t)', linewidth=2)
plt.plot(t, dG_predict, label='ΔG(t)')
plt.plot(t, BFIP_predict * 0.3, label='BFIP Active', color='green', linestyle='dotted')
plt.title('Tier 8.C — Predictive Gate Anticipation')
plt.xlabel('Time (s)')
plt.legend()
plt.tight_layout()
plt.savefig("BFIP_t8c_predictive_gate.png", dpi=300)
plt.show()

# Summary
summary = pd.DataFrame({
    'Final θ̄': [theta_predict[-1]],
    'Max ΔG': [np.max(dG_predict)],
    'Max MI Proxy': [np.max(MI_predict)],
    'BFIP Triggered?': [BFIP_predict.any()]
})

print("\n=== Tier 8C Summary ===\n", summary)


def main():
    pass

if __name__ == '__main__':
    main()
