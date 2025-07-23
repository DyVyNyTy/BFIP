from .kinetics import hill_equation
from .thermodynamics import gibbs_free_energy
from scipy.integrate import odeint
from scipy.interpolate import interp1d
import numpy as np

def run_single(ligand_vals, kinetic_params, thermo_params, T, t_span):
    θ_grid = np.zeros_like(ligand_vals, dtype=float)
    G_grid = np.zeros_like(ligand_vals, dtype=float)
    for i, L in enumerate(ligand_vals):
        θ = hill_equation(L, kinetic_params["n_H"], kinetic_params["K_d"])
        θ_grid[i] = θ
        G_grid[i] = gibbs_free_energy(θ, thermo_params["ΔH"], thermo_params["ΔS"], T)
    return θ_grid, G_grid

def run_dynamic(kinetic_params, ligand_range_vals, T, t_span):
    t = t_span  # Already a precomputed array
    amplitude = kinetic_params["amplitude"]

    # Interpret ligand range spec: [start, end, num_points]
    L_start = ligand_range_vals[0]
    L_end = ligand_range_vals[1]

    # Create ligand profile over time
    L_vals = np.linspace(L_start, L_end, num=len(t)) * amplitude
    ligand_interp = interp1d(t, L_vals, kind="linear", fill_value="extrapolate")

    def wrapped_binding_dynamics(θ, t, k_on, k_off):
        L = float(ligand_interp(t))
        return k_on * L * (1 - θ) - k_off * θ

    θ_t = odeint(wrapped_binding_dynamics, 0.0, t,
                 args=(kinetic_params["k_on"], kinetic_params["k_off"])).flatten()
    return t, θ_t


def main():
    pass

if __name__ == '__main__':
    main()
