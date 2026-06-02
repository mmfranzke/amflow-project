import mpmath as mp

from lsmethod.kernels import (
    l_residue_A1,
    l_residue_A1_short_without_extra_term,
    l_residue_B1,
    l_residue_delta0,
    l_residue_extra_offshell_term,
    l_residue_z,
    method1_l_residue_closed,
)
from lsmethod.two_loop_points import get_point


def test_l_residue_ingredients_equal_mass_offshell_positive():
    mp.mp.dps = 80
    point = get_point("equal_mass_offshell_positive")

    assert abs(l_residue_delta0(point) - mp.mpf(1) / 16) < mp.mpf("1e-70")
    assert abs(l_residue_A1_short_without_extra_term(point) - mp.mpf(5) / 6) < mp.mpf("1e-70")
    assert abs(l_residue_extra_offshell_term(point) + mp.mpf(55) / 36) < mp.mpf("1e-70")
    assert abs(l_residue_A1(point) + mp.mpf(25) / 36) < mp.mpf("1e-70")
    assert abs(l_residue_B1(point) - mp.mpf(8) / 9) < mp.mpf("1e-70")
    assert abs(l_residue_z(point) - mp.mpf(27) / 25) < mp.mpf("1e-70")
    assert abs(l_residue_z(point, short_C_diagnostic=True) - mp.mpf(14) / 15) < mp.mpf("1e-70")


def test_l_residue_full_minus_short_is_extra_offshell_term():
    offshell = get_point("equal_mass_offshell_positive")
    onshell = get_point("equal_mass_onshell_branch")

    assert abs(
        l_residue_A1(offshell)
        - l_residue_A1_short_without_extra_term(offshell)
        - l_residue_extra_offshell_term(offshell)
    ) < mp.mpf("1e-70")
    assert abs(l_residue_A1(offshell) - l_residue_A1_short_without_extra_term(offshell)) > mp.mpf("1")
    assert abs(l_residue_A1(onshell) - l_residue_A1_short_without_extra_term(onshell)) < mp.mpf("1e-70")


def test_l_residue_closed_leading_pole_is_half_method12_leading_pole():
    point = get_point("equal_mass_offshell_positive")
    x = mp.mpf(point.x.numerator) / point.x.denominator
    expected_lres = 1 / (2 * (1 - x))
    expected_method12 = 1 / (1 - x)

    eps_values = [mp.mpf(1) / n for n in (400, 500, 650, 800, 1000)]
    estimates = [eps**2 * method1_l_residue_closed(point, eps) for eps in eps_values]

    assert abs(estimates[-1] - expected_lres) < mp.mpf("1e-2")
    assert expected_method12 == 2 * expected_lres


def test_l_residue_closed_eps_point_is_close_to_direct_numeric_reference():
    mp.mp.dps = 50
    point = get_point("equal_mass_offshell_positive")
    eps = mp.mpf("0.1")
    eta = mp.mpf("1e-6")
    x = mp.mpf(point.x.numerator) / point.x.denominator
    delta0 = l_residue_delta0(point)
    a1 = l_residue_A1(point)
    b1 = l_residue_B1(point)
    d = 4 - 2 * eps
    omega = 2 * mp.pi ** ((2 - 2 * eps) / 2) / mp.gamma((2 - 2 * eps) / 2)
    pref = mp.gamma(eps) / (1 - x) * (1 / (1j * mp.pi ** (d / 2))) * mp.mpf("0.5") * (2 * mp.pi * 1j) * omega / 2

    def integrand_t(t):
        return t ** (-eps) * (a1 + b1 * t - 1j * eta) ** (-eps) / (t + delta0 - 1j * eta)

    root = -a1 / b1
    splits = [mp.mpf("0"), root * (1 - mp.mpf("1e-8")), root, root * (1 + mp.mpf("1e-8")), 2, 10, 100, mp.inf]
    numeric = pref * mp.quad(integrand_t, splits)
    closed = method1_l_residue_closed(point, eps, eta=eta)

    assert abs(closed - numeric) / abs(closed) < mp.mpf("5e-3")
