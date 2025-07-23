
# Experimental Validation Protocol for BFIP Model
## Objective
Validate BFIP regions for Fe²⁺, Ca²⁺, and H⁺ by measuring phase boundaries, functional coherence, and thermodynamic potentials.

## Methods
1. **Fe²⁺ in Hemoglobin**:
   - **Setup**: Purified hemoglobin (Sigma-Aldrich H7379) in a gas-controlled chamber (OxyGraph-2k, Oroboros).
   - **Measurements**:
     - Oxygen saturation (θ) via oximetry (OxyLite, Oxford Optronix) at pO₂ = 0–150 mmHg, pH = 6.5–8.0, [Fe²⁺] = 0.1–2 mM.
     - Redox state via UV-Vis spectroscopy (Agilent Cary 60, 400–600 nm) to confirm >90% Fe²⁺.
     - Dynamic pO₂ oscillations (10% amplitude, 100 s period) to measure θ(t).
     - Enthalpy (ΔH) and entropy (ΔS) via isothermal titration calorimetry (ITC, Malvern MicroCal PEAQ-ITC).
   - **Analysis**: Fit θ to Hill equation, compute I(P; F) using Fe²⁺ states and θ, estimate ΔG = ΔH - T ΔS.
2. **Ca²⁺ in Calmodulin**:
   - **Setup**: Recombinant calmodulin (Thermo Fisher P0DP29) in buffer (20 mM HEPES, pH 6.5–8.0).
   - **Measurements**:
     - Calcium binding via fluorescence (Fura-2, Invitrogen F1221) at pCa = 4–8, [Ca²⁺] = 0.1–10 μM.
     - Confirm >90% Ca²⁺ via mass spectrometry (Q Exactive, Thermo Fisher).
     - Dynamic pCa changes (0.1 unit amplitude, 100 s period) to measure θ(t).
     - ΔH and ΔS via ITC.
   - **Analysis**: Fit binding curves, compute I(P; F), estimate ΔG.
3. **H⁺ in Proton Gradients**:
   - **Setup**: Mitochondrial membranes (isolated from HEK293 cells) or liposomes with pH gradients.
   - **Measurements**:
     - Proton flux via pH-sensitive dyes (BCECF, Thermo Fisher B1151) at pH = 5–8, [H⁺] = 1 μM–1 mM.
     - H⁺ uniformity via NMR (Bruker Avance 600 MHz).
     - Dynamic pH oscillations (0.1 unit amplitude, 100 s period).
     - ΔH and ΔS via ITC or differential scanning calorimetry (DSC, TA Instruments Nano DSC).
   - **Analysis**: Fit flux curves, compute I(P; F), estimate ΔG.

## Expected Outcomes
- **Fe²⁺**: BFIP at pO₂ ~20–40 mmHg, pH ~7.4, [Fe²⁺] ~0.5 mM, I(P; F) > 0.5, ΔG < 0.
- **Ca²⁺**: BFIP at pCa ~6, pH ~7.0, [Ca²⁺] ~1 μM, I(P; F) > 0.5, ΔG < 0.
- **H⁺**: BFIP at pH ~6.5, [H⁺] ~0.1 mM, I(P; F) > 0.5, ΔG < 0.

## Data Analysis
- Compare simulated vs. experimental θ(t), I(P; F), and ΔG.
- Validate phase boundaries (ligand*) and uniformity constraints.
- Compute Ψ = G - T I to identify stable BFIP regions.
