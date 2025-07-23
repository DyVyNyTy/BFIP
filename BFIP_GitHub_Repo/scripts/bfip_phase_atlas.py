
#!/usr/bin/env python3
"""
bfip_phase_atlas.py (final version)

Generate a BFIP Phase Atlas using safe filenames without '+' characters.
"""
import os
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# Mapping: ion display name → filename prefix (safe)
ION_MAP = {
    'H+': 'H',
    'Ca2+': 'Ca2',
    'Fe2+': 'Fe2'
}
LABELS = ['Z', 'MI', 'THETA']
OUTPUT_DIR = '.'

def plot_atlas():
    fig, axes = plt.subplots(nrows=len(ION_MAP), ncols=3, figsize=(12, 10))

    for i, (ion, prefix) in enumerate(ION_MAP.items()):
        for j, label in enumerate(LABELS):
            filename = f"{prefix}_maps_thermo_{label}.png"
            path = os.path.join(OUTPUT_DIR, filename)
            ax = axes[i, j]
            if os.path.exists(path):
                img = mpimg.imread(path)
                ax.imshow(img)
                ax.set_title(f"{ion} - {label}")
                ax.axis('off')
            else:
                ax.text(0.5, 0.5, f"Missing: {filename}", ha='center', va='center')
                ax.axis('off')

    plt.tight_layout()
    out_path = os.path.join(OUTPUT_DIR, "BFIP_Phase_Atlas.png")
    plt.savefig(out_path)
    print(f"✅ Saved atlas to: {out_path}")

if __name__ == "__main__":
    plot_atlas()
