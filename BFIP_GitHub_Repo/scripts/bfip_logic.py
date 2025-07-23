
import numpy as np

def load_bfip_maps(h_path, ca_path, fe_path):
    Z_H = np.load(h_path)
    Z_Ca = np.load(ca_path)
    Z_Fe = np.load(fe_path)
    return Z_H, Z_Ca, Z_Fe

def classify_bfip_region(h, ca, fe):
    pattern = (h << 2) | (ca << 1) | fe
    interpretations = {
        0b000: "No BFIP activity — phase silent",
        0b100: "Proton-driven logic: pH gating or charge relay",
        0b010: "Scaffold memory / Ca²⁺ binding logic",
        0b001: "Redox-encoded gate: Fe²⁺ electron logic",
        0b110: "Dual signal: pH + structure memory",
        0b101: "Acid-redox switch: dynamic electron flow control",
        0b011: "Redox-structure pair: signal buffering or integrity check",
        0b111: "Consensus phase: Bio-functional control node (high-order logic)"
    }
    return pattern, interpretations.get(pattern, "Unknown pattern")

def analyze_overlap(Z_H, Z_Ca, Z_Fe):
    results = []
    shape = Z_H.shape
    for i in range(shape[0]):
        for j in range(shape[1]):
            h, ca, fe = Z_H[i, j], Z_Ca[i, j], Z_Fe[i, j]
            pattern, meaning = classify_bfip_region(h, ca, fe)
            results.append({
                'i': i, 'j': j,
                'H+': h, 'Ca2+': ca, 'Fe2+': fe,
                'Pattern': format(pattern, '03b'),
                'Interpretation': meaning
            })
    return results


def main():
    pass

if __name__ == '__main__':
    main()
