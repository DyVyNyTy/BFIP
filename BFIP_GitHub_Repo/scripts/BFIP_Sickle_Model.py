#!/usr/bin/env python3
import os
import time
import numpy as np
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
        return 0.3

def gibbs_free_energy(theta, ΔH, ΔS, T):
    return ΔH - T * ΔS * theta

def binding_dynamics(t, y, conc, pH, PCO2, T, ion_type):
    theta, MI = y[0], y[1]
    theta = np.clip(theta, 0.0, 1.0)
    ion = ions[ion_type]
    O2_level = o2_dynamics(t)
    effective_ligand = conc * (1 - 0.7 * (1 - O2_level))

    k_on = ion['k_on']
    k_off = ion['k_off'] * (
        1 + 0.3 * PCO2 / 50 + 0.1 * (7.4 - pH)**2 + 0.02 * (T - 310.15)
    )

    dtheta_dt = k_on * effective_ligand * (1 - theta) - k_off * theta

    if O2_level < 0.3 and theta < 0.15:
        dMI_dt = 0.05 * (1 - O2_level) * (1 - theta)
    else:
        decay_strength = 0.06 if MI > 50 else 0.12
        dMI_dt = -decay_strength * MI - 0.005 * (t / 800)

    MI = max(min(MI + dMI_dt, 1000), -100)
    return [dtheta_dt, MI]

def simulate():
    T = 310.15
    t = np.linspace(0, 800, 400)
    PCO2 = 50
    results = {}

    for ion_type in ions.keys():
        concs = np.linspace(*ions[ion_type]['conc_range'])
        pHs = np.linspace(*ions[ion_type]['pH_range'])

        θ_grid = np.zeros((len(concs), len(pHs)))
        MI_grid = np.zeros((len(concs), len(pHs)))
        G_grid = np.zeros((len(concs), len(pHs)))
        bfip_mask = np.zeros((len(concs), len(pHs)), dtype=bool)
        bfip_rupture_mask = np.zeros((len(concs), len(pHs)), dtype=bool)

        for i, c in enumerate(concs):
            for j, pH in enumerate(pHs):
                y0 = [0.4, 10.0]
                sol = solve_ivp(binding_dynamics, [0, 800], y0,
                                args=(c, pH, PCO2, T, ion_type),
                                t_eval=t, method='Radau', rtol=1e-9, atol=1e-12)
                θ = sol.y[0][-1]
                MI = max(min(sol.y[1][-1], 1000), -100)
                ΔG = gibbs_free_energy(θ, ions[ion_type]['Delta_H'], ions[ion_type]['Delta_S'], T)

                current_BFIP = (θ > 0.3 and 10 < MI <= 90 and ΔG < -3000)
                bfip_mask[i, j] = current_BFIP

                if current_BFIP and MI < 10:
                    bfip_rupture_mask[i, j] = True

                θ_grid[i, j] = θ
                MI_grid[i, j] = MI
                G_grid[i, j] = ΔG

                print(f"[{ion_type}] θ={θ:.3f}, MI={MI:.3f}, ΔG={ΔG:.2f}, BFIP={current_BFIP}, RUPTURE={bfip_rupture_mask[i,j]}")

        results[ion_type] = {
            'conc': concs,
            'pH': pHs,
            'theta_grid': θ_grid,
            'MI_grid': MI_grid,
            'G_grid': G_grid,
            'bfip_mask': bfip_mask,
            'bfip_rupture_mask': bfip_rupture_mask
        }

    print(f"Simulation done in {time.time() - start:.2f}s")
    return results

def plot(results):
    for ion_type, data in results.items():
        conc = data['conc']
        pH = data['pH']
        θ = data['theta_grid']
        MI = data['MI_grid']
        ΔG = data['G_grid']
        BFIP = data['bfip_mask']
        RUPTURE = data['bfip_rupture_mask']
        Ion, PH = np.meshgrid(conc, pH, indexing='ij')

        fig = go.Figure(data=[
            go.Scatter3d(x=Ion.flatten(), y=PH.flatten(), z=θ.flatten(), mode='markers',
                         marker=dict(size=4, color=θ.flatten(), colorscale='Viridis', colorbar=dict(title='θ'))),
            go.Scatter3d(x=Ion.flatten(), y=PH.flatten(), z=MI.flatten(), mode='markers',
                         marker=dict(size=4, color=MI.flatten(), colorscale='Reds', colorbar=dict(title='MI'))),
            go.Scatter3d(x=Ion.flatten(), y=PH.flatten(), z=ΔG.flatten(), mode='markers',
                         marker=dict(size=4, color=ΔG.flatten(), colorscale='Cividis', colorbar=dict(title='ΔG'))),
            go.Scatter3d(x=Ion[BFIP], y=PH[BFIP], z=ΔG[BFIP], mode='markers',
                         marker=dict(size=6, color='magenta', symbol='diamond'), name='BFIP Zone'),
            go.Scatter3d(x=Ion[RUPTURE], y=PH[RUPTURE], z=ΔG[RUPTURE], mode='markers',
                         marker=dict(size=6, color='red', symbol='x'), name='BFIP Rupture')
        ])
        fig.update_layout(title=f'BFIP Mapping: {ion_type}',
                          scene=dict(xaxis_title='[Ion] (mM)', yaxis_title='pH', zaxis_title='Value'))
        fig.write_html(f"bfip_phase_map_{ion_type}.html", include_plotlyjs='cdn', auto_open=False)

        # ➕ New 2D ΔG vs θ trace with MI coloring
        trace_fig = go.Figure()
        trace_fig.add_trace(go.Scatter(
            x=data['theta_grid'].flatten(),
            y=data['G_grid'].flatten(),
            mode='markers',
            marker=dict(
                size=6,
                color=data['MI_grid'].flatten(),
                colorscale='Turbo',
                colorbar=dict(title='MI'),
                showscale=True
            ),
            name='ΔG vs θ (MI-colored)'
        ))
        trace_fig.update_layout(
            title=f"ΔG vs θ — {ion_type}",
            xaxis_title="θ",
            yaxis_title="ΔG (kJ/mol)",
        )
        trace_fig.write_html(f"deltaG_vs_theta_trace_{ion_type}.html", include_plotlyjs='cdn', auto_open=False)
        print(f"Saved: bfip_phase_map_{ion_type}.html and deltaG_vs_theta_trace_{ion_type}.html")

if __name__ == "__main__":
    start = time.time()
    result = simulate()
    plot(result)
