#!/usr/bin/env python3
"""Reference solution: hydration free energy of methane by alchemical decoupling."""
import argparse
import json
import os
import time

import numpy as np
import openmm as mm
from openmm import app, unit
from openmmtools import alchemy, states
from pymbar import MBAR, timeseries

SEED = 42
TEMPERATURE = 298.15 * unit.kelvin
PRESSURE = 1.0 * unit.atmosphere
TIMESTEP = 2.0 * unit.femtoseconds
FRICTION = 1.0 / unit.picosecond
CUTOFF = 0.9 * unit.nanometer
PADDING = 1.0 * unit.nanometer
BAROSTAT_FREQ = 25

LAMBDA_ELEC = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
LAMBDA_STERICS = [1.0, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]
K = len(LAMBDA_ELEC)

KJ = unit.kilojoule_per_mole
R_KCAL = 0.0019872041

def log(msg):
    print(time.strftime("%H:%M:%S"), msg, flush=True)

def build_system(inputs):
    pdb = app.PDBFile(os.path.join(inputs, "methane.pdb"))
    ff = app.ForceField("amber14-all.xml", "amber14/tip3p.xml", os.path.join(inputs, "methane.xml"))
    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.addSolvent(ff, model="tip3p", padding=PADDING, neutralize=False)
    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=CUTOFF,
        constraints=app.HBonds,
        rigidWater=True,
        ewaldErrorTolerance=5e-4,
    )
    for force in system.getForces():
        if isinstance(force, mm.NonbondedForce):
            force.setUseDispersionCorrection(True)
    barostat = mm.MonteCarloBarostat(PRESSURE, TEMPERATURE, BAROSTAT_FREQ)
    barostat.setRandomNumberSeed(SEED)
    system.addForce(barostat)
    solute = [a.index for a in modeller.topology.atoms() if a.residue.name == "CH4"]
    if len(solute) != 5:
        raise RuntimeError(f"expected 5 methane atoms, found {len(solute)}")
    return modeller, system, solute

def build_alchemical_state(system, solute):
    factory = alchemy.AbsoluteAlchemicalFactory()
    region = alchemy.AlchemicalRegion(alchemical_atoms=solute)
    alch_system = factory.create_alchemical_system(system, region)
    alch_state = alchemy.AlchemicalState.from_system(alch_system)
    thermo = states.ThermodynamicState(alch_system, temperature=TEMPERATURE, pressure=PRESSURE)
    return states.CompoundThermodynamicState(thermo, composable_states=[alch_state])

def set_box(context, box):
    context.setPeriodicBoxVectors(mm.Vec3(*box[0]), mm.Vec3(*box[1]), mm.Vec3(*box[2]))

def set_lambda(compound, context, k):
    compound.lambda_electrostatics = LAMBDA_ELEC[k]
    compound.lambda_sterics = LAMBDA_STERICS[k]
    compound.apply_to_context(context)

def run_state(k, context, integrator, compound, positions, box, n_eq, n_prod, sample_every, first):
    set_lambda(compound, context, k)
    set_box(context, box)
    context.setPositions(positions)
    if first:
        log("minimising")
        mm.LocalEnergyMinimizer.minimize(context, maxIterations=1000)
    context.setVelocitiesToTemperature(TEMPERATURE, SEED)
    log(f"state {k}: equilibrating {n_eq} steps")
    integrator.step(n_eq)
    n_samples = n_prod // sample_every
    pos_frames = np.empty((n_samples, context.getSystem().getNumParticles(), 3))
    box_frames = np.empty((n_samples, 3, 3))
    log(f"state {k}: production {n_prod} steps, {n_samples} frames")
    for i in range(n_samples):
        integrator.step(sample_every)
        st = context.getState(getPositions=True)
        pos_frames[i] = st.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        box_frames[i] = np.array(st.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer))
    st = context.getState(getPositions=True)
    return (
        pos_frames,
        box_frames,
        st.getPositions(asNumpy=True).value_in_unit(unit.nanometer),
        np.array(st.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer)),
    )

def reduced_potentials(compound, context, pos_all, box_all):
    kT = (unit.MOLAR_GAS_CONSTANT_R * TEMPERATURE).value_in_unit(KJ)
    n = len(pos_all)
    u = np.empty((K, n))
    for l in range(K):
        set_lambda(compound, context, l)
        for i in range(n):
            set_box(context, box_all[i])
            context.setPositions(pos_all[i])
            st = context.getState(getEnergy=True)
            e = st.getPotentialEnergy().value_in_unit(KJ)
            vol = st.getPeriodicBoxVolume().value_in_unit(unit.nanometer**3)
            pv = (PRESSURE * vol * unit.nanometer**3 * unit.AVOGADRO_CONSTANT_NA).value_in_unit(KJ)
            u[l, i] = (e + pv) / kT
        log(f"reduced potentials at state {l} done")
    return u

def analyse(u_kn, N_k):
    kT_kcal = R_KCAL * 298.15
    mbar = MBAR(u_kn, N_k)
    res = mbar.compute_free_energy_differences()
    dG = -float(res["Delta_f"][0, -1]) * kT_kcal
    err = float(res["dDelta_f"][0, -1]) * kT_kcal
    overlap = np.asarray(mbar.compute_overlap()["matrix"])
    return dG, err, overlap

def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", default="/app/inputs")
    ap.add_argument("--outdir", default="/app/output")
    ap.add_argument("--workdir", default="/app/work")
    ap.add_argument("--steps", type=int, default=500_000)
    ap.add_argument("--eq-steps", type=int, default=50_000)
    ap.add_argument("--sample-every", type=int, default=500)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.workdir, exist_ok=True)

    modeller, system, solute = build_system(args.inputs)
    log(f"system: {system.getNumParticles()} particles, solute atoms {solute}")
    compound = build_alchemical_state(system, solute)
    integrator = mm.LangevinMiddleIntegrator(TEMPERATURE, FRICTION, TIMESTEP)
    integrator.setRandomNumberSeed(SEED)
    context = compound.create_context(integrator)

    positions = modeller.positions.value_in_unit(unit.nanometer)
    positions = np.array(positions)
    box = np.array(modeller.topology.getPeriodicBoxVectors().value_in_unit(unit.nanometer))

    for k in range(K):
        ckpt = os.path.join(args.workdir, f"state_{k:02d}.npz")
        if os.path.exists(ckpt):
            z = np.load(ckpt)
            positions, box = z["final_pos"], z["final_box"]
            log(f"state {k}: loaded checkpoint")
            continue
        pf, bf, positions, box = run_state(
            k, context, integrator, compound, positions, box,
            args.eq_steps, args.steps, args.sample_every, first=(k == 0),
        )
        np.savez(ckpt, pos=pf, box=bf, final_pos=positions, final_box=box)

    pos_all, box_all, owner = [], [], []
    for k in range(K):
        z = np.load(os.path.join(args.workdir, f"state_{k:02d}.npz"))
        pos_all.append(z["pos"])
        box_all.append(z["box"])
        owner += [k] * len(z["pos"])
    pos_all = np.concatenate(pos_all)
    box_all = np.concatenate(box_all)
    owner = np.array(owner)

    u_path = os.path.join(args.workdir, "u_all.npy")
    if os.path.exists(u_path):
        u_all = np.load(u_path)
    else:
        u_all = reduced_potentials(compound, context, pos_all, box_all)
        np.save(u_path, u_all)

    blocks, N_k = [], []
    for k in range(K):
        cols = np.where(owner == k)[0]
        series = u_all[k, cols]
        g = timeseries.statistical_inefficiency(series)
        idx = timeseries.subsample_correlated_data(series, g=g)
        log(f"state {k}: g = {g:.1f}, kept {len(idx)} of {len(cols)} samples")
        blocks.append(u_all[:, cols[idx]])
        N_k.append(len(idx))
    u_kn = np.concatenate(blocks, axis=1)
    N_k = np.array(N_k, dtype=int)

    dG, err, overlap = analyse(u_kn, N_k)
    log(f"delta G hyd = {dG:.3f} +/- {err:.3f} kcal/mol")
    nb = min(min(overlap[i, i + 1], overlap[i + 1, i]) for i in range(K - 1))
    log(f"minimum neighbour overlap = {nb:.3f}")

    np.savez(os.path.join(args.outdir, "mbar_input.npz"), u_kn=u_kn, N_k=N_k)
    out = {
        "delta_G_hyd_kcal_mol": dG,
        "uncertainty_kcal_mol": err,
        "overlap_matrix": overlap.tolist(),
        "n_states": K,
        "steps_per_state": args.steps,
        "equilibration_steps": args.eq_steps,
        "seed": SEED,
        "temperature_K": 298.15,
        "pressure_atm": 1.0,
        "timestep_fs": 2.0,
        "nonbonded_cutoff_nm": 0.9,
        "lambda_electrostatics": LAMBDA_ELEC,
        "lambda_sterics": LAMBDA_STERICS,
    }
    with open(os.path.join(args.outdir, "fep.json"), "w") as f:
        json.dump(out, f, indent=2)
    log("wrote fep.json and mbar_input.npz")

if __name__ == "__main__":
    main()
