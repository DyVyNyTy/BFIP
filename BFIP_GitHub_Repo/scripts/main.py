# ion_phase_lab/main.py

from ion_phase_lab.engine import run_bfip_engine

if __name__ == "__main__":
    print("🔬 Ion Phase Lab CLI")
    print("1. Run BFIP Engine")
    print("2. Export results")
    choice = input("Choose an option [1–2]: ").strip()

    if choice == "1":
        run_bfip_engine()
    elif choice == "2":
        out_path = input("Enter output path (default: results/bfip_points.csv): ").strip() or "results/bfip_points.csv"
        run_bfip_engine(output_path=out_path)
    else:
        print("Invalid selection. Exiting.")
