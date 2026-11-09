"""Build a synthetic 12-state harmonic submission for testing tests/test_fep.py."""
import json, sys, numpy as np
rng = np.random.default_rng(1)
K = 12; n = 200
kappa = np.geomspace(4.0, 0.5, K)
kT = 0.0019872041 * 298.15
true = 0.5*np.log(kappa[0]/kappa[-1]) * kT
x = np.concatenate([rng.normal(0, 1/np.sqrt(kappa[k]), n) for k in range(K)])
u_kn = 0.5 * kappa[:, None] * x[None, :]**2
N_k = np.full(K, n)
import pymbar
m = pymbar.MBAR(u_kn, N_k)
r = m.compute_free_energy_differences()
dG = -float(r["Delta_f"][0, -1]) * kT; err = float(r["dDelta_f"][0, -1]) * kT
O = m.compute_overlap()["matrix"]
print("analytic", true, "mbar", dG, "+/-", err, "min nbr overlap", min(min(O[i,i+1],O[i+1,i]) for i in range(K-1)))
out_dir = sys.argv[1]
le = [1.0] + [0.0]*11
ls = [1.0, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]
fep = dict(delta_G_hyd_kcal_mol=dG, uncertainty_kcal_mol=err, overlap_matrix=O.tolist(), n_states=12,
           steps_per_state=500000, equilibration_steps=50000, seed=42, temperature_K=298.15, pressure_atm=1.0,
           timestep_fs=2.0, nonbonded_cutoff_nm=0.9, lambda_electrostatics=le, lambda_sterics=ls)
json.dump(fep, open(out_dir + "/fep.json", "w"))
np.savez(out_dir + "/mbar_input.npz", u_kn=u_kn, N_k=N_k)
json.dump({"delta_G_hyd_kcal_mol": dG, "tolerance_kcal_mol": 0.5}, open(out_dir + "/truth.json", "w"))
