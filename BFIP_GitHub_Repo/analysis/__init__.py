# ion_phase_lab/__init__.py
# Marks the root of the lab as a package

# ion_phase_lab/models/__init__.py
# Imports for shared use
from .kinetics import hill_equation, binding_dynamics, dynamic_ligand_base
from .thermodynamics import gibbs_free_energy, psi_function
from .sampling import sample_parameters
from .simulation import run_single, run_dynamic

# ion_phase_lab/analysis/__init__.py
# Future shared analysis imports go here

# ion_phase_lab/engine.py
# Formerly bfip_model.py — high-level model orchestrator
# This is where simulation parameters are defined and orchestrated
# (Kept modular for use in CLI or Snakemake)

# Placeholder for upcoming content migration

def main():
    pass

if __name__ == '__main__':
    main()
