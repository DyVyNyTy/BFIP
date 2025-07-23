
#!/usr/bin/env python3
import os
import numpy as np
import time
from scipy.integrate import solve_ivp
import plotly.graph_objects as go

print("Working directory:", os.getcwd())

ions = {
    'hbs_cushion': {
        'conc_range': (0.01, 1.0, 20), 'pH_range': (6.8, 7.4, 5),
        'n_H': 1.0, 'K_d': 5e-4, 'k_on': 8e4, 'k_off': 3e3,
        'Delta_H': -6000.0, 'Delta_S': 10.0
    },
    'Ca2+': {
        'conc_range': (0.0001, 0.01, 3), 'pH_range': (6.5, 8.0, 3),
        'n_H': 2.0, 'K_d': 1e-6, 'k_on': 1e4, 'k_off': 30,
        'Delta_H': -3000.0, 'Delta_S': 7.0
    },
    'Fe2+': {
        'conc_range': (0.01, 1.0, 8), 'pH_range': (6.5, 8.0, 4),
        'n_H': 1.0, 'K_d': 1e-5, 'k_on': 5e4, 'k_off': 500,
        'Delta_H': -5000.0, 'Delta_S': 8.0
    }
}

def o2_dynamics(t):
    if t < 100:
        return 0.95
    elif t < 300:
        return 0.05 + 0.03 * np.sin(0.04 * (t - 100))
    else:
        return 0.2 + 0.1 * np.cos(0.01 * t)  # simulate stress

def gibbs_free_energy(theta, ΔH, ΔS, T):
    return ΔH - T * ΔS * theta

def binding_dynamics(t, y, conc, pH, PCO2, T, ion_type, perturbation=0.0):
    theta, MI = y[0], y[1]
    theta = np.clip(theta, 0.0, 1.0)
    ion = ions[ion_type]
    O2_level = o2_dynamics(t)
    effective_ligand = conc * (1 - 0.7 * (1 - O2_level))
    k_on = ion['k_on']
    k_off = ion['k_off'] * (1 + 0.3 * PCO2 / 50 + 0.1 * (7.4 - pH)**2 + 0.02 * (T - 310.15))

    dtheta_dt = k_on * effective_ligand * (1 - theta) - k_off * theta

    dMI_dt = (
        -0.1 * (1 - O2_level) * theta +  # stress effect
        0.1 * np.sin(0.05 * t) +         # rhythmic fluctuations
        perturbation                     # direct perturbation
    )
    MI = max(min(MI + dMI_dt, 1000), -200)
    return [dtheta_dt, MI]

def simulate():
    T = 310.15
    t = np.linspace(0, 800, 400)
    results = {}

    for ion_type in ions.keys():
        concs = np.linspace(*ions[ion_type]['conc_range'])
        pHs = np.linspace(*ions[ion_type]['pH_range'])
        PCO2 = 50
        rupture_mask = np.zeros((len(concs), len(pHs)), dtype=bool)
        bfip_mask = np.zeros((len(concs), len(pHs)), dtype=bool)

        for i, c in enumerate(concs):
            for j, pH in enumerate(pHs):
                y0 = [0.4, 1000.0]
                sol = solve_ivp(binding_dynamics, [0, 800], y0,
                    args=(c, pH, PCO2, T, ion_type, -0.3), t_eval=t,
                    method='Radau', rtol=1e-9, atol=1e-12)
                theta = sol.y[0][-1]
                MI = sol.y[1][-1]
                dG = gibbs_free_energy(theta, ions[ion_type]['Delta_H'], ions[ion_type]['Delta_S'], T)

                bfip = (theta > 0.3 and 10 < MI < 900 and dG < -3000)
                rupture = (theta > 0.85 and MI < 20)

                bfip_mask[i, j] = bfip
                rupture_mask[i, j] = rupture
                print(f"[{ion_type}] θ={theta:.3f}, MI={MI:.2f}, ΔG={dG:.2f}, BFIP={bfip}, RUPTURE={rupture}")

        results[ion_type] = {
            'conc': concs,
            'pH': pHs,
            'bfip': bfip_mask,
            'rupture': rupture_mask
        }

    return results

def plot(results):
    for ion_type, data in results.items():
        Ion, PH = np.meshgrid(data['conc'], data['pH'], indexing='ij')
        fig = go.Figure()
        fig.add_trace(go.Scatter3d(
            x=Ion.flatten(), y=PH.flatten(), z=data['bfip'].flatten().astype(int),
            mode='markers', marker=dict(size=4, color='blue'), name='BFIP'
        ))
        fig.add_trace(go.Scatter3d(
            x=Ion.flatten(), y=PH.flatten(), z=data['rupture'].flatten().astype(int),
            mode='markers', marker=dict(size=4, color='red'), name='Rupture'
        ))
        fig.update_layout(
            title=f'Perturbation Sweep Map: {ion_type}',
            scene=dict(xaxis_title='[Ion] (mM)', yaxis_title='pH', zaxis_title='State Flag'),
        )
        fig.write_html(f'perturbation_sweep_map_{ion_type}.html')

if __name__ == "__main__":
    start = time.time()
    result = simulate()
    plot(result)
    print(f"Simulation complete in {time.time() - start:.2f} seconds.")
