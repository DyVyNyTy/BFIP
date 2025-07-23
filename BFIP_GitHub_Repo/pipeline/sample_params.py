#!/usr/bin/env python
import numpy as np
import logging
from bfip.sampling import sample_parameters

# Set up logging
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

try:
    logging.info("Starting sample_params.py")
    
    # Snakemake injects this object automatically
    ion = snakemake.wildcards.ion
    config = snakemake.config
    logging.info(f"Processing ion: {ion}")

    # Read the per-ion parameter dictionary
    ion_params = config["ions"][ion]
    logging.info(f"Ion parameters: {ion_params}")

    # Number of Monte Carlo samples (fallback to 100 if not set in config)
    n_samples = config.get("n_samples", 100)
    logging.info(f"Number of samples: {n_samples}")

    # Generate the samples
    logging.info("Generating parameter samples")
    samples = sample_parameters(ion_params, n_samples=n_samples)
    logging.info(f"Sampled parameters: {list(samples.keys())}")

    # Save as a single .npz
    output_path = snakemake.output[0]
    logging.info(f"Saving samples to {output_path}")
    np.savez(output_path, **samples)
    logging.info("Parameter sampling completed successfully")

except Exception as e:
    logging.error(f"Error in sample_params.py: {str(e)}")
    raise

def main():
    pass

if __name__ == '__main__':
    main()
