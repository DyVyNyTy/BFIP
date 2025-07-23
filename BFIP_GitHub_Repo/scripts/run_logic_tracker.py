
import os
import sys
import yaml
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), 'models', 'models'))
from simulation import run_dynamic

def load_config(path):
    with open(path, 'r') as f:
        return yaml.safe_load(f)

def run_logic_tracker(cfg_path, output_dir):
    cfg = load_config(cfg_path)
    ions = [ion_cfg['ion'].replace('+', 'p') for ion_cfg in cfg['ions']]
    logic_threshold = 0.05
    time_steps = None
    logic_codes = []

    for idx, ion_cfg in enumerate(cfg['ions']):
        ion_cfg['parameters']['amplitude'] = cfg.get('amplitude', 1.0)
        t, θ_t = run_dynamic(ion_cfg['parameters'], cfg['ligand_range'], cfg['temperature'], cfg['t_span'])
        if time_steps is None:
            time_steps = t
        if idx == 0:
            θ_matrix = np.zeros((len(cfg['ions']), len(θ_t)))
        θ_matrix[idx, :] = θ_t

    for i in range(len(time_steps)):
        logic_code = ''.join(['1' if θ_matrix[j, i] > logic_threshold else '0' for j in range(len(cfg['ions']))])
        logic_codes.append(logic_code)

    # Plot logic code transitions over time
    numeric_codes = [int(code, 2) for code in logic_codes]
    plt.figure(figsize=(10, 4))
    plt.plot(time_steps, numeric_codes, drawstyle='steps-post', color='purple')
    plt.title(f"Real-Time Logic Code Transitions – {'_'.join(ions)}")
    plt.xlabel("Time")
    plt.ylabel("Logic Code (Binary → Int)")
    plt.grid(True)
    plt.tight_layout()
    output_path = os.path.join(output_dir, f"{'_'.join(ions)}_logic_flips.png")
    plt.savefig(output_path)
    plt.close()

if __name__ == "__main__":
    run_logic_tracker(sys.argv[1], os.path.dirname(__file__))
