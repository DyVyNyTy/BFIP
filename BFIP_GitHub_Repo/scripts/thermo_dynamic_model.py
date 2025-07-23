
import numpy as np
from scipy.integrate import odeint

R = 8.314  # J/(mol·K)

def hill_equation(theta, t, k_on, amplitude, K_d, n_H):
    """ODE for Hill-like cooperative binding"""
    L_t = amplitude * np.sin(0.1 * t) + amplitude  # ligand varies with time
    dtheta_dt = k_on * L_t**n_H * (1 - theta) - (k_on * K_d) * theta
    return dtheta_dt

def run_dynamic_thermo(kinetic_params, ligand_range, T, t_span, dH, dS, n_H=1.0):
    """
    Run dynamic theta simulation with thermodynamic-dependent kinetics.

    Inputs:
    - kinetic_params: dict with 'k_on', 'amplitude'
    - ligand_range: not used here (only for compatibility)
    - T: temperature (K)
    - t_span: np.array of time points
    - dH: enthalpy (J/mol)
    - dS: entropy (J/mol·K)
    - n_H: Hill coefficient
    """
    k_on = kinetic_params['k_on']
    amplitude = kinetic_params.get('amplitude', 1.0)

    # Compute ΔG and K_d
    dG = dH - T * dS
    K_d = np.exp(dG / (R * T))  # unitless

    # Integrate θ(t)
    theta0 = 0.0
    theta_t = odeint(hill_equation, theta0, t_span, args=(k_on, amplitude, K_d, n_H))
    return t_span, theta_t.flatten()


def main():
    pass

if __name__ == '__main__':
    main()
