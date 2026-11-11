"""Verifier for the methane hydration free energy task."""
import json
import math
import os
import numpy as np
import pytest

OUT_DIR = os.environ.get("FEP_OUTPUT_DIR", "/app/output")
TRUTH_PATH = os.environ.get("FEP_GROUND_TRUTH", "/tests/ground_truth.json")

FEP_JSON = os.path.join(OUT_DIR, "fep.json")
NPZ = os.path.join(OUT_DIR, "mbar_input.npz")
PDB = os.path.join(OUT_DIR, "final_structure.pdb")

K = 12
KT_KCAL = 0.0019872041 * 298.15

REQUIRED_KEYS = [
    "delta_G_hyd_kcal_mol", "uncertainty_kcal_mol", "overlap_matrix",
    "n_states", "steps_per_state", "equilibration_steps", "seed",
    "temperature_K", "pressure_atm", "timestep_fs", "nonbonded_cutoff_nm",
    "lambda_electrostatics", "lambda_sterics",
]

MIN_NEIGHBOR_OVERLAP = 0.03
MIN_UNCERTAINTY = 0.01
MAX_UNCERTAINTY = 0.50
MIN_SAMPLES_PER_STATE = 8
DG_AGREEMENT = 0.02
UNC_REL_AGREEMENT = 0.25
OVERLAP_AGREEMENT = 0.02

@pytest.fixture(scope="module")
def data():
    with open(FEP_JSON) as f:
        return json.load(f)

@pytest.fixture(scope="module")
def truth():
    with open(TRUTH_PATH) as f:
        return json.load(f)

@pytest.fixture(scope="module")
def npz():
    z = np.load(NPZ)
    return {k: z[k] for k in z.files}

@pytest.fixture(scope="module")
def recomputed(npz):
    from pymbar import MBAR
    u_kn = np.asarray(npz["u_kn"], dtype=float)
    N_k = np.asarray(npz["N_k"]).astype(int)
    mbar = MBAR(u_kn, N_k)
    res = mbar.compute_free_energy_differences()
    dG = -float(res["Delta_f"][0, -1]) * KT_KCAL
    dG_err = float(res["dDelta_f"][0, -1]) * KT_KCAL
    overlap = np.asarray(mbar.compute_overlap()["matrix"], dtype=float)
    return {"dG": dG, "err": dG_err, "overlap": overlap}

def test_output_files_exist():
    assert os.path.exists(FEP_JSON), "fep.json missing"
    assert os.path.exists(NPZ), "mbar_input.npz missing"
    assert os.path.exists(PDB), "final_structure.pdb missing (anti-cheat check failed)"

def test_pdb_validity():
    """Anti-cheat: Parse PDB to ensure simulation actually ran."""
    with open(PDB) as f:
        lines = f.readlines()
    has_ch4 = any("CH4" in line for line in lines)
    has_water = any("HOH" in line or "TIP3" in line for line in lines)
    assert has_ch4, "final_structure.pdb does not contain the CH4 solute."
    assert has_water, "final_structure.pdb does not contain solvent water molecules."
    assert len(lines) > 500, "final_structure.pdb is too small to be a solvated system."

def test_strict_keys(data):
    assert set(data.keys()) == set(REQUIRED_KEYS), f"Keys do not match exactly. Found {set(data.keys())}"

def test_numbers_are_finite(data):
    for key in ("delta_G_hyd_kcal_mol", "uncertainty_kcal_mol"):
        assert not math.isnan(data[key]) and not math.isinf(data[key])

def test_npz_strict_keys(npz):
    assert set(npz.keys()) == {"u_kn", "N_k"}, f"npz must contain exactly u_kn and N_k, found {set(npz.keys())}"

def test_mbar_input_shapes(npz):
    u_kn, N_k = npz["u_kn"], npz["N_k"]
    assert u_kn.shape[0] == K
    assert N_k.shape == (K,)
    assert np.all(N_k == np.round(N_k)), "N_k must be integers"
    assert int(N_k.sum()) == u_kn.shape[1], "N_k.sum() must equal u_kn.shape[1]"
    assert np.all(N_k >= MIN_SAMPLES_PER_STATE), f"Need >= {MIN_SAMPLES_PER_STATE} samples per state, got {N_k.tolist()}"

def test_mbar_input_is_finite_and_varied(npz):
    u_kn = np.asarray(npz["u_kn"], dtype=float)
    assert np.all(np.isfinite(u_kn))
    spread = u_kn.max(axis=0) - u_kn.min(axis=0)
    assert np.median(spread) > 0.1, "u_kn shows no dependence on lambda"

def test_reported_free_energy_matches_recomputation(data, recomputed):
    assert abs(data["delta_G_hyd_kcal_mol"] - recomputed["dG"]) <= DG_AGREEMENT

def test_reported_uncertainty_matches_recomputation(data, recomputed):
    assert abs(data["uncertainty_kcal_mol"] - recomputed["err"]) <= max(UNC_REL_AGREEMENT * recomputed["err"], 0.01)

def test_overlap_matrix_matches_recomputation(data, recomputed):
    O = np.asarray(data["overlap_matrix"], dtype=float)
    assert np.max(np.abs(O - recomputed["overlap"])) <= OVERLAP_AGREEMENT

def test_neighbor_overlap(recomputed):
    O = recomputed["overlap"]
    up = np.array([O[i, i + 1] for i in range(K - 1)])
    down = np.array([O[i + 1, i] for i in range(K - 1)])
    assert min(up.min(), down.min()) > MIN_NEIGHBOR_OVERLAP

def test_uncertainty_in_plausible_range(data):
    u = data["uncertainty_kcal_mol"]
    assert MIN_UNCERTAINTY <= u <= MAX_UNCERTAINTY

def test_free_energy_against_reference(data, truth):
    tol = float(truth.get("tolerance_kcal_mol", 1.0))
    assert abs(data["delta_G_hyd_kcal_mol"] - truth["delta_G_hyd_kcal_mol"]) <= tol
