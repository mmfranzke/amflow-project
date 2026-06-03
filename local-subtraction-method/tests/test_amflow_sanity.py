import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "two_loop_method3_mc_compare.py"
sys.path.insert(0, str(SCRIPT_PATH.parent))
spec = importlib.util.spec_from_file_location("two_loop_method3_mc_compare", SCRIPT_PATH)
compare = importlib.util.module_from_spec(spec)
spec.loader.exec_module(compare)


FIXTURE = """
BEGIN_AMFLOW_METADATA
point = equal_mass_offshell_positive
pPlus = 1
pMinus = 5
pPerp2 = 0
p2 = 5
x = 1/4
y = 1/4
ml2 = 1
Ml2 = 1
mk2 = 1
Mk2 = 1
AMFLOW_EPS_ORDER = 4
AMFLOW_PRECISION_GOAL = 10
END_AMFLOW_METADATA
AMFLOW_COEFF_POWER=-2 RAW=1 NORMALIZED=0.6666666666666666
AMFLOW_COEFF_POWER=-1 RAW=1 NORMALIZED=3.005686100222
AMFLOW_COEFF_POWER=0 RAW=1 NORMALIZED=2.226118637208
"""


def args(**overrides):
    base = {
        "amflow_eps_order": "4",
        "precision_goal": "10",
        "amflow_sanity_check": True,
        "allow_high_order_amflow": False,
    }
    base.update(overrides)
    return SimpleNamespace(**base)


def test_amflow_metadata_parses_and_matches_point():
    metadata = compare.parse_amflow_metadata_block(FIXTURE)
    assert metadata["point"] == "equal_mass_offshell_positive"
    assert compare.amflow_metadata_mismatches(metadata, "equal_mass_offshell_positive", args()) == []


def test_amflow_metadata_mismatch_is_detected():
    metadata = compare.parse_amflow_metadata_block(FIXTURE.replace("x = 1/4", "x = 1/5"))
    mismatches = compare.amflow_metadata_mismatches(metadata, "equal_mass_offshell_positive", args())
    assert any(key == "x" for key, _actual, _expected in mismatches)


def test_known_equal_mass_offshell_amflow_sanity_passes(monkeypatch, tmp_path):
    _raw, normalized = compare.parse_imported_amflow_coefficients(FIXTURE)
    monkeypatch.setattr(compare, "sanity_marker_path", lambda point_name: tmp_path / f"{point_name}.json")
    status, _messages = compare.amflow_sanity_check("equal_mass_offshell_positive", normalized, args())
    assert status == "PASSED"


def test_absurd_amflow_coefficients_fail_sanity():
    coeffs = {-2: 1e6 + 1e5j, -1: 0, 0: 0}
    status, messages = compare.amflow_sanity_check("equal_mass_offshell_positive", coeffs, args())
    assert status == "FAILED"
    assert messages


def test_high_order_amflow_requires_explicit_allow_or_marker():
    try:
        compare.check_high_order_amflow_allowed(
            "four_mass_offshell_branch_safe_B",
            args(amflow_eps_order="6", amflow="fresh", allow_high_order_amflow=False),
        )
    except SystemExit as exc:
        assert "Refusing high-order AMFlow" in str(exc)
    else:
        raise AssertionError("high-order AMFlow was not rejected")

    compare.check_high_order_amflow_allowed(
        "four_mass_offshell_branch_safe_B",
        args(amflow_eps_order="6", amflow="fresh", allow_high_order_amflow=True),
    )
