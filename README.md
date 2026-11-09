# Methane-fep-hydration
Calculate the hydration free energy of methane with OpenMM, `openmmtools.alchemy` and MBAR.
The agent-facing description is in `instruction.md`.

## Layout

| Path | Purpose |
|---|---|
| `task.toml` | Task config and metadata |
| `instruction.md` | Prompt given to the agent (full protocol, output schema) |
| `environment/Dockerfile` | Agent image (conda-forge: OpenMM 8.1.1, openmmtools 0.24.0, pymbar 4.0.3) |
| `environment/inputs/` | `methane.pdb`, `methane.xml` (GAFF c3/hc, AM1-BCC charges) copied to `/app/inputs` |
| `solution/solve.py`, `solve.sh` | Oracle. Runs the real simulation (multi-hour, resumable) |
| `tests/Dockerfile` | Verifier image (pytest, numpy, pymbar) |
| `tests/test.sh`, `tests/test_fep.py` | Verifier |
| `tests/ground_truth.json` | Reference value and tolerance (never in the agent image) |
| `dev/make_synthetic_submission.py` | Builds a synthetic submission to test the verifier without OpenMM |

## What the verifier checks

1. `fep.json` and `mbar_input.npz` exist; all required keys present; values finite.
2. Declared protocol: 12 states, 500,000 steps, seed 42, 298.15 K, 1 atm, 2 fs, 0.9 nm cutoff.
3. Lambda schedule: 12 distinct states, endpoints (1,1) and (0,0), non-increasing, electrostatics off before sterics drop below 1.
4. `u_kn` has shape `(12, N)`, finite, depends on lambda; every state has at least 50 decorrelated samples.
5. MBAR is re-run on `u_kn`. The reported free energy (within 0.02 kcal/mol), uncertainty (within 25%) and overlap matrix (within 0.02) must match.
6. Overlap matrix is 12x12, non-negative, rows sum to 1, and adjacent-state overlap is above 0.03.
7. Uncertainty is in [0.01, 0.5] kcal/mol.
8. Free energy is within `tolerance_kcal_mol` of the reference.

Limit: a submission that fabricates a self-consistent `u_kn` would pass checks 4-7. The declared-protocol fields are declarations, not measurements. The task relies on the agent environment having no network and on the cost of faking a believable `u_kn` that also lands on the reference.

## Difficulty

Alchemical decoupling has to be set up correctly: soft-core sterics, a lambda schedule that keeps neighbouring states overlapping (electrostatics off before sterics), reduced potentials evaluated at every state for every sample, decorrelation, and MBAR with the right sign convention.

## Reference solution

`solution/solve.sh` runs the full multi-hour simulation in `solve.py` to produce the reference outputs.

## Verification

`tests/test.sh` runs `test_fep.py`, which checks the schema, protocol, lambda schedule, and re-runs MBAR on the submitted raw data to verify the reported free energy and overlap matrix.
