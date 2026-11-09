Calculate the hydration free energy of methane (CH4) with OpenMM, openmmtools and MBAR, using alchemical decoupling of the solute from TIP3P water.

**Inputs**

`/app/inputs/` contains the solute:

- `methane.pdb`: methane, residue name `CH4`.
- `methane.xml`: an OpenMM force-field file for methane (GAFF `c3`/`hc` parameters, AM1-BCC charges). Load it together with `amber14-all.xml` and `amber14/tip3p.xml`.

**Protocol**

Use these settings exactly. The reference value was produced with them.

| Item | Setting |
|---|---|
| Solvent | TIP3P (`amber14/tip3p.xml`), rigid, cubic box, 1.0 nm padding, no added ions |
| Nonbonded | PME, 0.9 nm cutoff, dispersion correction on, Ewald tolerance 5e-4 |
| Constraints | `HBonds`, 2 fs time step |
| Ensemble | NPT, 298.15 K, 1 atm; `LangevinMiddleIntegrator` (1/ps friction), `MonteCarloBarostat` every 25 steps |
| Alchemical method | `openmmtools.alchemy` (`AbsoluteAlchemicalFactory`, `AlchemicalRegion` on the 5 methane atoms), soft-core potential, openmmtools defaults for the soft-core parameters |
| States | exactly 12 lambda states, state 0 fully interacting, state 11 fully decoupled |
| Schedule rule | `lambda_electrostatics` reaches 0 before `lambda_sterics` drops below 1; both are non-increasing along the state index. The rest of the schedule is up to you |
| Equilibration | at least 50,000 steps per state, not used in the analysis |
| Production | exactly 500,000 steps per state |
| Random seed | 42 for integrator and barostat |
| Estimator | MBAR (`pymbar` 4), after subsampling each state by its statistical inefficiency |

Sign and standard state: the hydration free energy is G(solute interacting with water) minus G(solute decoupled), so it is `f_0 - f_11` in units of kT, converted to kcal/mol. No standard-state correction is needed.

Budget: 12 states of 500,000 steps is a multi-hour CPU job (4 cores). Time a short run first, run states so that completed ones are saved, and keep going if one state needs a rerun. You have 28,800 seconds in total.

**Output**

Write two files to `/app/output/`.

**`mbar_input.npz`** with arrays:

- `u_kn`: shape `(12, N)`, reduced potentials `u = beta * (E + P*V)` of all retained samples evaluated at each of the 12 states.
- `N_k`: shape `(12,)`, integer number of retained samples that came from each state. Samples in `u_kn` are ordered by the state they came from, and `N_k.sum() == N`.

**`fep.json`** with exactly these keys:

```json
{
  "delta_G_hyd_kcal_mol": 0.0,
  "uncertainty_kcal_mol": 0.0,
  "overlap_matrix": [[0.0]],
  "n_states": 12,
  "steps_per_state": 500000,
  "equilibration_steps": 50000,
  "seed": 42,
  "temperature_K": 298.15,
  "pressure_atm": 1.0,
  "timestep_fs": 2.0,
  "nonbonded_cutoff_nm": 0.9,
  "lambda_electrostatics": [],
  "lambda_sterics": []
}
