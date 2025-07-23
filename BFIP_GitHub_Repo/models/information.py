# ion_phase_lab/__init__.py
# Root lab package

# ion_phase_lab/models/__init__.py
from .kinetics import hill_equation, binding_dynamics, dynamic_ligand_base
from .thermodynamics import gibbs_free_energy, psi_function
from .sampling import sample_parameters
from .simulation import run_single, run_dynamic


# ion_phase_lab/models/information.py
def compute_mutual_information(pJ, pP, pF, eps=1e-12):
    """
    Compute mutual information (I(P;F)) using outer product formulation.
    
    Parameters:
    - pJ: joint distribution matrix (P × F)
    - pP: marginal over parameter space
    - pF: marginal over functional space
    - eps: small value to prevent log(0)
    
    Returns:
    - MI value in bits
    """
    import numpy as np

    pJ = np.clip(pJ, eps, None)
    pP = np.clip(pP, eps, None)
    pF = np.clip(pF, eps, None)

    ratio = pJ / (pP[:, None] * pF[None, :])
    mi = max(np.nansum(pJ * np.log2(ratio)), 0.0)
    return mi


def main():
    pass

if __name__ == '__main__':
    main()
