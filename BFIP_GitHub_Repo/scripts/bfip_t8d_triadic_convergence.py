
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import square
from scipy.ndimage import gaussian_filter1d
import pandas as pd

t = np.linspace(0, 10, 1000)
dt = t[1] - t[0]

# Ion input pulses
I_t = np.sin(2 * np.pi * 1.5 * t) ** 2

# ⏱️ Tier 8.D Triadic Oscillator Function
def triadic_gate(t, I_t, memory_window=50, prediction_horizon=25, freq_start=1.0):
    gate = np.zeros_like(t)
    freq = freq_start
    current_phase = 0
    MI_trace = []

    for i in range(memory_window, len(t) - prediction_horizon, 10):
        # Phase gate
        t_seg = t[i:i+prediction_horizon]
        if len(t_seg) == 0:
            continue
        osc = (np.sin(2 * np.pi * freq * t_seg + current_phase) + 1) / 2

        # Memory predictor
        memory = I_t[i-memory_window:i]
        if np.sum(memory) > memory_window * 0.25:
            end = min(i + prediction_horizon + len(osc), len(gate))
            segment_len = end - (i + prediction_horizon)
            gate[i + prediction_horizon:end] += osc[:segment_len]

        # MI proxy (past window to oscillator)
        past = I_t[i-memory_window:i]
        mi_val = np.correlate(past - np.mean(past), osc - np.mean(osc), mode='valid') / len(past)
        MI_trace.append(mi_val[0] if len(mi_val) else 0)

        # Frequency adjust
        if mi_val[0] > 0.005:
            freq += 0.05
        else:
            freq -= 0.025
        freq = np.clip(freq, 0.5, 3.0)

        # Phase adjust
        current_phase += 0.05

    gate = gaussian_filter1d(gate, sigma=4)
    return np.clip(gate, 0, 1), np.array(MI_trace)

# Run triadic control gate
O_triadic, MI_series = triadic_gate(t, I_t)

# Gated signal
G_triadic = I_t * O_triadic

# θ̄, ΔG, MI
def theta_bar(input_signal):
    return np.cumsum(input_signal) / len(input_signal)

def delta_G(input_signal):
    return -np.gradient(input_signal)

def mutual_info_proxy(I, O):
    return np.correlate(I - np.mean(I), O - np.mean(O), mode='same') / len(I)

theta = theta_bar(G_triadic)
dG = delta_G(theta)
MI_proxy = mutual_info_proxy(I_t, O_triadic)
BFIP = (theta > 0.15) & (MI_proxy > 0.005)

# Plot
plt.figure(figsize=(12, 6))
plt.plot(t, I_t, label='Ion Input I(t)', linestyle='--', alpha=0.5)
plt.plot(t, G_triadic, label='Triadic Gated Input')
plt.plot(t, theta, label='θ̄(t)', linewidth=2)
plt.plot(t, dG, label='ΔG(t)')
plt.plot(t, BFIP * 0.3, label='BFIP Active', color='green', linestyle='dotted')
plt.title('Tier 8.D — Triadic Convergence Trigger')
plt.xlabel('Time (s)')
plt.legend()
plt.tight_layout()
plt.savefig("BFIP_t8d_triadic_convergence.png", dpi=300)
plt.show()

# Summary
summary = pd.DataFrame({
    'Final θ̄': [theta[-1]],
    'Max ΔG': [np.max(dG)],
    'Max MI Proxy': [np.max(MI_proxy)],
    'BFIP Triggered?': [BFIP.any()]
})

print("\n=== Tier 8D Summary ===\n", summary)


def main():
    pass

if __name__ == '__main__':
    main()
