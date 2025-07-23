# bfip/sampling.py
import numpy as np

def sample_parameters(param_dict, n_samples):
    """
    Sample kinetic and thermodynamic parameters from normal distributions.
    Expects param_dict to have nested dicts for keys “K_d”, “k_on”, “k_off”, “Delta_H”, “Delta_S”,
    each with “mean” and “std” entries.
    Returns a dict mapping those keys to length-n_samples numpy arrays.
    """
    samples = {}
    param_keys = ["K_d", "k_on", "k_off", "Delta_H", "Delta_S"]
    for key in param_keys:
        stats = param_dict.get(key)
        if stats is None or "mean" not in stats or "std" not in stats:
            raise ValueError(f"Missing or malformed stats for '{key}' in config.yaml")
        mean = float(stats["mean"])
        std  = float(stats["std"])
        samples[key] = np.random.normal(mean, std, size=n_samples)
    return samples


def main():
    pass

if __name__ == '__main__':
    main()
