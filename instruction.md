Calculate the hydration free energy of methane (CH4) with OpenMM, openmmtools and MBAR, using alchemical decoupling of the solute from TIP3P water.

Inputs:
The directory `/app/inputs/` contains the solute:
- `/app/inputs/methane.pdb`: methane, residue name `CH4`.
- `/app/inputs/methane.xml`: an OpenMM force-field file for methane (GAFF `c3`/`hc` parameters, AM1-BCC charges). Load it together with `amber14-all.xml` and `amber14/tip3p.xml`.

Protocol:
Use these settings: Solvent is TIP3P, rigid, cubic box, 1.0 nm padding, no added ions. Nonbonded interactions use PME, 0.9 nm cutoff, dispersion correction on, Ewald tolerance 5e-4. Constraints are HBonds with a 2 fs time step. Ensemble is NPT at 298.15 K and 1 atm with LangevinMiddleIntegrator (1/ps friction) and MonteCarloBarostat every 25 steps. Use openmmtools.alchemy with a soft-core potential on the 5 methane atoms. Use exactly 12 lambda states, state 0 fully interacting, state 11 fully decoupled. The rule is lambda_electrostatics reaches 0 before lambda_sterics drops below 1; both are non-increasing. Equilibration is at least 5,000 steps per state. Production is exactly 25,000 steps per state. Random seed is 42 for integrator and barostat. Use MBAR (pymbar 4) after subsampling each state by its statistical inefficiency.

The hydration free energy is G(solute interacting) minus G(solute decoupled), so it is f_0 - f_11 in units of kT, converted to kcal/mol.

Acceptance criteria:
To ensure statistical convergence and pass verification, your final u_kn must yield at least 8 decorrelated samples per state, a minimum adjacent-state overlap greater than 0.03, and an uncertainty between 0.01 and 0.5 kcal/mol. The reported free energy, uncertainty, and overlap matrix must match an independent MBAR recomputation within 0.02 kcal/mol (or 25% for uncertainty). The final_structure.pdb must contain the solvated system.

Output:
Write three files to `/app/output/`:
1. `mbar_input.npz` containing exactly `u_kn` and `N_k`. `u_kn` must be shape `(12, N)`, the reduced potentials u = beta * (E + P*V) of all retained samples evaluated at each of the 12 states. `N_k` must be shape `(12,)`, the integer number of retained samples that came from each state. Samples in `u_kn` are ordered by the state they came from, and `N_k.sum() == N`.
2. `fep.json` containing exactly these keys: `delta_G_hyd_kcal_mol`, `uncertainty_kcal_mol`, `overlap_matrix`, `n_states`, `steps_per_state`, `equilibration_steps`, `seed`, `temperature_K`, `pressure_atm`, `timestep_fs`, `nonbonded_cutoff_nm`, `lambda_electrostatics`, `lambda_sterics`.
3. `final_structure.pdb` containing the final coordinates of the solvated system.

The result must be within an undisclosed tolerance of an undisclosed reference value. Do not use online solutions or hints specific to this task.

You have 28800 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
