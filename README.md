# Methane FEP Hydration Free Energy

## Difficulty
Computing solvation free energies via alchemical decoupling is a standard but critical workflow in computational drug discovery, used daily by researchers to estimate binding affinities and partition coefficients. The data (methane, TIP3P water, Amber14) are real experimental and standard force-field parameters. The level of challenge is moderate-to-high: the protocol requires correctly setting up soft-core potentials, statistical decorrelation, and MBAR estimation.

The key insight is that one must decouple electrostatics first, then sterics, and evaluate the reduced potential of every sampled frame at all 12 states to ensure the MBAR overlap matrix remains well-conditioned.

While the system is a classic benchmark, the task is novel in its rigorous verification constraints. It relies on mathematical recomputation rather than trusting text declarations, forcing the agent to genuinely execute the full multi-state FEP pipeline.

## Reference solution
`solution/solve.py` executes a full 12-state OpenMM simulation. It checkpoints per state, subsamples by statistical inefficiency, and applies `pymbar` to compute the final free energy. The script uses a fixed random seed (42) for both the integrator and the barostat. This ensures deterministic reproducibility: repeated runs on identical hardware will produce the exact same trajectory and thermodynamic outputs. The script also writes the final solvated system to `final_structure.pdb`.

## Verification
`tests/test_fep.py` parses the output PDB to ensure the methane solute and solvent box are physically present (anti-cheat). It re-runs MBAR on the submitted `u_kn` matrix to prevent hardcoding.

Calibration evidence for the verification gates:
| Gate | Threshold | Reason / Evidence |
| :--- | :--- | :--- |
| Adjacent overlap | > 0.03 | Standard threshold in the pymbar community (Klimovich et al., 2015) to ensure lambda states are statistically connected. |
| MBAR agreement | 0.02 | Typo-level tolerance allowing for floating-point differences between the agent's MBAR run and the verifier's MBAR run. |
| Uncertainty range | 0.01 - 0.5 | 0.01 is the noise floor for this simulation length; 0.5 is the maximum acceptable error for a physical result. |
| Samples per state | >= 8 | Ensures enough independent data points to construct a valid overlap matrix. |
| Free energy tolerance | 1.0 | TIP3P/GAFF simulations typically deviate from experiment by 0.3-0.8 kcal/mol. 1.0 accounts for this. |

Ground Truth Provenance:
The reference value (2.10 kcal/mol) is the experimental value from the FreeSolv database (Mobley and Guthrie, 2014, J. Comput. Aided Mol. Des.), Entry 451.
