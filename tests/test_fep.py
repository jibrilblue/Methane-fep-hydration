"""Verifier for the methane hydration free energy task.

Checks, in order: output files and schema, declared protocol, lambda schedule,
reduced-potential file, an independent MBAR recomputation that must agree with
what was reported, overlap-matrix quality, a plausible uncertainty, and
finally the free energy against the reference.
"""
import json
import math
import os

import numpy as np
import pytest

OUT_DIR = os.environ.get("FEP_OUTPUT_DIR", "/app/output")
TRUTH_PATH = os.environ.get("FEP_GROUND_TRUTH", "/tests/ground_truth.json")

FEP_JSON = os.path.join(OUT_DIR, "fep.json")
NPZ = os.path.join(OUT_DIR, "mbar_input.npz")

K = 12
KT_KCAL = 0.0019872041 * 298.15  # kcal/mol at 298.15 K (gas constant in kcal/mol/K)

REQUIRED_KEYS = [
    "delta_G_hyd_kcal_mol",
    "uncertainty_kcal_mol",
    "overlap_matrix",
    "n_states",
    "steps_per_state",
    "equilibration_steps",
    "seed",
    "temperature_K",
    "pressure_atm",
    "timestep_fs",
    "nonbonded_cutoff_nm",
    "lambda_electrostatics",
    "lambda_sterics",
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
    assert os.path.exists(FEP_JSON), f"Output file missing: {FEP_JSON}"
    with open(FEP_JSON) as f:
        try:
            return json.load(f)
        except json.JSONDecodeError as e:
            pytest.fail(f"fep.json is not valid JSON: {e}")


@pytest.fixture(scope="module")
def truth():
    with open(TRUTH_PATH) as f:
        return json.load(f)


@pytest.fixture(scope="module")
def npz():
    assert os.path.exists(NPZ), f"Output file missing: {NPZ}"
    try:
        z = np.load(NPZ)
        return {k: z[k] for k in z.files}
    except Exception as e:
        pytest.fail(f"mbar_input.npz could not be read: {e}")


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
    assert os.path.exists(FEP_JSON), f"Output file missing: {FEP_JSON}"
    assert os.path.exists(NPZ), f"Output file missing: {NPZ}"


def test_required_keys(data):
    missing = [k for k in REQUIRED_KEYS if k not in data]
    assert not missing, f"fep.json is missing keys: {missing}"


def test_numbers_are_finite(data):
    for key in ("delta_G_hyd_kcal_mol", "uncertainty_kcal_mol"):
        v = data[key]
        assert isinstance(v, (int, float)) and not isinstance(v, bool), f"{key} must be a number"
        assert not math.isnan(v), f"{key} is NaN"
        assert not math.isinf(v), f"{key} is Inf"


def test_declared_protocol(data):
    assert data["n_states"] == K, f"n_states must be {K}"
    assert data["steps_per_state"] == 500000, "steps_per_state must be 500000"
    assert data["equilibration_steps"] >= 50000, "equilibration_steps must be >= 50000"
    assert data["seed"] == 42, "seed must be 42"
    assert abs(data["temperature_K"] - 298.15) < 1e-6, "temperature_K must be 298.15"
    assert abs(data["pressure_atm"] - 1.0) < 1e-6, "pressure_atm must be 1.0"
    assert abs(data["timestep_fs"] - 2.0) < 1e-6, "timestep_fs must be 2.0"
    assert abs(data["nonbonded_cutoff_nm"] - 0.9) < 1e-6, "nonbonded_cutoff_nm must be 0.9"


def test_lambda_schedule(data):
    le = np.asarray(data["lambda_electrostatics"], dtype=float)
    ls = np.asarray(data["lambda_sterics"], dtype=float)
    assert le.shape == (K,) and ls.shape == (K,), f"lambda lists must have {K} entries"
    assert np.all((le >= 0) & (le <= 1)) and np.all((ls >= 0) & (ls <= 1)), "lambdas must lie in [0, 1]"
    assert le[0] == 1.0 and ls[0] == 1.0, "state 0 must be fully interacting"
    assert le[-1] == 0.0 and ls[-1] == 0.0, "last state must be fully decoupled"
    assert np.all(np.diff(le) <= 1e-12), "lambda_electrostatics must be non-increasing"
    assert np.all(np.diff(ls) <= 1e-12), "lambda_sterics must be non-increasing"
    below = ls < 1.0
    assert np.all(le[below] == 0.0), "electrostatics must be fully off before sterics drop below 1"
    pairs = set(zip(le.round(8), ls.round(8)))
    assert len(pairs) == K, "all 12 lambda states must be distinct"


def test_mbar_input_shapes(npz):
    assert "u_kn" in npz and "N_k" in npz, "mbar_input.npz must contain u_kn and N_k"
    u_kn, N_k = npz["u_kn"], npz["N_k"]
    assert u_kn.ndim == 2 and u_kn.shape[0] == K, f"u_kn must have shape ({K}, N)"
    assert N_k.shape == (K,), f"N_k must have shape ({K},)"
    assert np.all(N_k == np.round(N_k)), "N_k must be integers"
    assert int(N_k.sum()) == u_kn.shape[1], "N_k.sum() must equal the number of columns of u_kn"
    assert np.all(N_k >= MIN_SAMPLES_PER_STATE), (
        f"every state needs at least {MIN_SAMPLES_PER_STATE} decorrelated samples, got {N_k.tolist()}"
    )


def test_mbar_input_is_finite_and_varied(npz):
    u_kn = np.asarray(npz["u_kn"], dtype=float)
    assert np.all(np.isfinite(u_kn)), "u_kn contains NaN or Inf"
    spread = u_kn.max(axis=0) - u_kn.min(axis=0)
    assert np.median(spread) > 0.1, "u_kn shows no dependence on lambda"
    assert np.unique(np.round(u_kn[0], 6)).size > 0.5 * u_kn.shape[1], "u_kn rows look duplicated or constant"


def test_reported_free_energy_matches_recomputation(data, recomputed):
    diff = abs(data["delta_G_hyd_kcal_mol"] - recomputed["dG"])
    assert diff <= DG_AGREEMENT, (
        f"Reported dG {data['delta_G_hyd_kcal_mol']:.4f} differs from MBAR on mbar_input.npz "
        f"({recomputed['dG']:.4f}) by {diff:.4f} kcal/mol."
    )


def test_reported_uncertainty_matches_recomputation(data, recomputed):
    rep, rec = data["uncertainty_kcal_mol"], recomputed["err"]
    assert abs(rep - rec) <= max(UNC_REL_AGREEMENT * rec, 0.01), (
        f"Reported uncertainty {rep:.4f} differs from MBAR recomputation {rec:.4f}"
    )


def test_overlap_matrix_structure(data):
    O = np.asarray(data["overlap_matrix"], dtype=float)
    assert O.shape == (K, K), f"overlap_matrix must be {K}x{K}, got {O.shape}"
    assert np.all(np.isfinite(O)), "overlap_matrix contains NaN or Inf"
    assert np.all(O >= -1e-8), "overlap_matrix has negative entries"
    assert np.allclose(O.sum(axis=1), 1.0, atol=0.01), "overlap_matrix rows must sum to 1"


def test_overlap_matrix_matches_recomputation(data, recomputed):
    O = np.asarray(data["overlap_matrix"], dtype=float)
    assert O.shape == (K, K), f"overlap_matrix must be {K}x{K}"
    diff = np.max(np.abs(O - recomputed["overlap"]))
    assert diff <= OVERLAP_AGREEMENT, f"Reported overlap matrix differs from recomputation by {diff:.4f}"


def test_neighbor_overlap(recomputed):
    O = recomputed["overlap"]
    up = np.array([O[i, i + 1] for i in range(K - 1)])
    down = np.array([O[i + 1, i] for i in range(K - 1)])
    worst = float(min(up.min(), down.min()))
    assert worst > MIN_NEIGHBOR_OVERLAP, (
        f"Adjacent-state overlap {worst:.4f} is below {MIN_NEIGHBOR_OVERLAP}"
    )


def test_uncertainty_in_plausible_range(data):
    u = data["uncertainty_kcal_mol"]
    assert u > 0.0, "Uncertainty must be non-zero"
    assert MIN_UNCERTAINTY <= u <= MAX_UNCERTAINTY, (
        f"Uncertainty {u:.4f} kcal/mol outside plausible range [{MIN_UNCERTAINTY}, {MAX_UNCERTAINTY}]"
    )


def test_free_energy_against_reference(data, truth):
    tol = float(truth.get("tolerance_kcal_mol", 0.5))
    diff = abs(data["delta_G_hyd_kcal_mol"] - truth["delta_G_hyd_kcal_mol"])
    assert diff <= tol, f"Energy out of tolerance: |{data['delta_G_hyd_kcal_mol']:.3f} - ref| = {diff:.3f} > {tol}"
