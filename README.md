# Methane FEP Hydration Free Energy

## Difficulty
Computing solvation free energies via alchemical decoupling is a standard but critical workflow in computational drug discovery, used daily by researchers to estimate binding affinities and partition coefficients. The data (methane, TIP3P water, Amber14) are real experimental and standard force-field parameters. The level of challenge is moderate-to-high: the protocol requires correctly setting up soft-core potentials, statistical decorrelation, and MBAR estimation.

The key insight is that one must decouple electrostatics first, then sterics, and evaluate the reduced potential of every sampled frame at all 12 states to ensure the MBAR overlap matrix remains well-conditioned.

While the system is a classic benchmark, the task is novel in its rigorous verification constraints. It relies on mathematical recomputation rather than trusting text declarations to verify the physics results.

## Reference solution
`solution/solve.py` executes a full 12-state OpenMM simulation. It checkpoints per state, subsamples by statistical inefficiency, and applies `pymbar` to compute the final free energy. The script uses a fixed random seed (42) to improve reproducibility and writes the final solvated system to `final_structure.pdb`.

## Verification
`tests/test_fep.py` parses the output PDB to ensure the methane solute, solvent box, and periodic boundaries are physically present. It re-runs MBAR on the submitted `u_kn` matrix to prevent hardcoding. 

### Calibration Evidence
The verifier gates were calibrated by running the full 25,000-step reference solution (matching the default `solve.py` arguments) in the same containerized environment. The measured results were:
- **Delta G:** 2.560 kcal/mol
- **Uncertainty:** 0.251 kcal/mol
- **Minimum neighbor overlap:** 0.079
- **Minimum samples kept per state:** 13 (out of 250 frames)

Based on these measured results, the verification thresholds are set as follows:
- **Overlap > 0.03**: The test run produced a minimum adjacent overlap of 0.079. The 0.03 threshold allows a safe margin for statistical noise.
- **Uncertainty range (0.01 - 0.5)**: The test run yielded an uncertainty of 0.251. The 0.5 upper limit allows for the higher noise expected from the agent's implementation.
- **MBAR agreement (0.02)**: Typo-level tolerance allowing for floating-point differences between the agent's MBAR run and the verifier's MBAR run.
- **Samples per state (>= 8)**: The test run kept a minimum of 13 samples. The 8-sample floor is safely below this.
- **Free energy tolerance (1.0)**: The measured Delta G is 2.560 kcal/mol, which is 0.46 away from the experimental reference (2.10). A tolerance of 1.0 accommodates this force-field deviation.

### Ground Truth Provenance
The reference value (2.10 kcal/mol) is the experimental value from the FreeSolv database (Mobley and Guthrie, 2014). 
