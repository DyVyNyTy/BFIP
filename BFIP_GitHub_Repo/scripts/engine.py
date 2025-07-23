import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.stats import norm
import csv, time, os

from ion_phase_lab.models.kinetics import hill_equation, binding_dynamics, dynamic_ligand_base
from ion_phase_lab.models.thermodynamics import gibbs_free_energy
from ion_phase_lab.models.information import compute_mutual_information

import os

def run_bfip_engine(output_path=None):
    if output_path is None:
        base_dir = os.path.dirname(__file__)
        output_path = os.path.join(base_dir, 'results', 'bfip_points.csv')


    start = time.time()
    print(f"Start: {time.time()-start:.2f}s")

    time_series = np.linspace(0, 200, 200)  # seconds
    temperatures = np.array([305.0, 310.0, 315.0])  # Kelvin

    targets = {
        'Fe2+': {'MI': 0.05, 'theta': 0.05, 'n_RTln2': 1},
        'Ca2+': {'MI': 0.05, 'theta': 0.05, 'n_RTln2': 1},
        'H+':  {'MI': 0.0, 'theta': 0.6, 'n_RTln2': 2},
    }

    def get_ion_params():
        return {
            'Fe2+': {
                'conc_range': (0.01, 5.0, 20), 'pH_range': (6.5, 8.0, 10), 'lig_range': (0.0, 50.0, 20),
                'n_H': 2.0, 'K_d': 26.0, 'k_on': 2e6, 'k_off': 5e3, 'dH': -65.0, 'dS': -0.2
            },
            'Ca2+': {
                'conc_range': (1e-5, 0.1, 20), 'pH_range': (6.5, 8.0, 10), 'lig_range': (4.0, 8.0, 20),
                'n_H': 2.0, 'K_d': 3e-6, 'k_on': 2e7, 'k_off': 5e2, 'dH': -30.0, 'dS': -0.1
            },
            'H+': {
                'conc_range': (1e-5, 1e-3, 10), 'pH_range': (5.0, 8.0, 10), 'lig_range': (5.0, 8.0, 20),
                'n_H': 1.5, 'K_d': 1e-7, 'k_on': 1e8, 'k_off': 1e4, 'dH': -20.0, 'dS': -0.05
            }
        }

    def bind_ode(theta, t, ion, lig_fn, k_on, k_off, pH):
        L = lig_fn(t, theta)
        conc = (L / 760.0) if ion == 'Fe2+' else 10**(-L)
        return k_on * conc * (1 - theta) - k_off * theta

    ion_params = get_ion_params()
    csv_rows = []

    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    for ion, params in ion_params.items():
        thr = targets[ion]
        concs = np.linspace(*params['conc_range'])
        pHs = np.linspace(*params['pH_range'])
        Ls = np.linspace(*params['lig_range'])

        for c in concs:
            for pH in pHs:
                for T in temperatures:
                    k_on  = params['k_on'] * (1 - 0.1*(pH-7.4)**2)
                    k_off = params['k_off'] * (1 + 0.1*(pH-7.4)**2)
                    baseL = Ls[len(Ls)//2]
                    lig_fn = dynamic_ligand_base(baseL, 1.0)


                    theta_arr = np.array([hill_equation(L, params['n_H'], params['K_d']) for L in Ls])
                    G_arr = np.array([gibbs_free_energy(theta, params['dH'], params['dS'], T) for theta in theta_arr])

                    dtheta_dL = np.gradient(theta_arr, Ls)
                    idx = np.argmax(dtheta_dL)
                    L_star = Ls[idx]
                    theta_s = hill_equation(L_star, params['n_H'], params['K_d'])
                    G_s = gibbs_free_energy(theta_s, params['dH'], params['dS'], T)

                    grid = np.linspace(0, 1, 100)
                    pP = norm.pdf(grid, np.mean(theta_arr), 0.05); pP /= pP.sum()
                    dyn = odeint(bind_ode, 0.5, time_series, args=(ion, lig_fn, k_on, k_off, pH)).flatten()
                    hist, _ = np.histogram(dyn, bins=50, range=(0,1), density=True)
                    pF = np.clip(hist, 1e-12, None)
                    pJ = np.outer(pP, pF)

                    MI = compute_mutual_information(pJ, pP, pF)

                    RTln2 = 8.314 * T/1000 * np.log(2)
                    if MI > thr['MI'] and theta_s > thr['theta'] and G_s < -thr['n_RTln2']*RTln2:
                        csv_rows.append([ion, c, pH, T, L_star, theta_s, MI, G_s])

    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['ion', 'conc', 'pH', 'T', 'L*', 'theta*', 'MI', 'G_S'])
        writer.writerows(csv_rows)

    print(f"Exported {len(csv_rows)} BFIP points to {output_path}")

if __name__ == "__main__":
    run_bfip_engine()
