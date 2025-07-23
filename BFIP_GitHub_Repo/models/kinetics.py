import numpy as np
from scipy.integrate import odeint

def hill_equation(ligand, n_H, K_d):
    return ligand**n_H / (K_d**n_H + ligand**n_H)

def binding_dynamics(theta, t, ligand_func, k_on, k_off):
    conc = ligand_func(t)
    return k_on * conc * (1 - theta) - k_off * theta

def dynamic_ligand_base(base, amp, period=100):
    def fn(t, θ):
        feedback = θ*(1-θ)
        return base*(1 + amp*np.sin(2*np.pi*t/period)*feedback)
    return fn


def main():
    pass

if __name__ == '__main__':
    main()
