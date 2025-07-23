import numpy as np
from scipy.stats import entropy

def mutual_information(P, Q):
    """Compute I(P;Q) robustly via discrete distributions + entropies."""
    P = P/np.sum(P); Q = Q/np.sum(Q)
    P_joint = np.outer(P, Q)
    return np.sum(P_joint * np.log2(P_joint / (P[:,None]*Q[None,:])), where=P_joint>0)

def mask_bfip(theta_star, info, G_star, T):
    RT = 8.314*T/1000
    return (info > 0.5) & (theta_star > 0.1) & (G_star < -RT*np.log(2))


def main():
    pass

if __name__ == '__main__':
    main()
