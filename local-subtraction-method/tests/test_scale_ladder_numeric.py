import importlib.util
import sys
from fractions import Fraction
from pathlib import Path


SCRIPT_PATH = Path(__file__).resolve().parents[1] / "scripts" / "two_loop_method3_mc_compare.py"
sys.path.insert(0, str(SCRIPT_PATH.parent))
spec = importlib.util.spec_from_file_location("two_loop_method3_mc_compare", SCRIPT_PATH)
compare = importlib.util.module_from_spec(spec)
spec.loader.exec_module(compare)


def test_ksubloop_numeric_matches_short_scale_for_branch_safe_q2():
    point = compare.get_point("equal_mass_offshell_positive")
    eps = Fraction(1, 10)
    eta = Fraction(1, 10**8)
    q2 = Fraction(0)

    numeric, _raw, _pref = compare.ksubloop_numeric_mpmath(point, eps, eta, q2, dps=40)
    short_value = compare.ksubloop_analytic(point, eps, eta, q2, "short")
    full_value = compare.ksubloop_analytic(point, eps, eta, q2, "full")

    assert abs(numeric - short_value) / abs(short_value) < 1e-10
    assert abs(numeric - full_value) / abs(full_value) > 1e-2
