import mpmath as mp

from lsmethod.closed_forms import (
    product_closed_form,
    product_laurent_coefficients,
    residue_closed_form,
    residue_laurent_coefficients,
    residue_ingredients,
)
from lsmethod.two_loop_points import get_point


def test_closed_form_api_leading_poles_equal_mass_offshell_positive():
    point = get_point("equal_mass_offshell_positive")
    x = mp.mpf(point.x.numerator) / point.x.denominator
    X = 1 - x

    product = product_laurent_coefficients(point, max_order=0)
    residue = residue_laurent_coefficients(point, max_order=0)

    assert abs(product[-2] - complex(1 / X)) < 1e-12
    assert abs(residue[-2] - complex(1 / (2 * X))) < 1e-12


def test_residue_ingredients_expose_full_and_short_C():
    point = get_point("equal_mass_offshell_positive")
    ingredients = residue_ingredients(point)
    full_by_name = residue_ingredients(point, scale="full")
    short_by_name = residue_ingredients(point, scale="short")
    short = residue_ingredients(point, short_C_diagnostic=True)

    assert ingredients["C_full_includes_extra_offshell_term"] is True
    assert full_by_name["C_full_includes_extra_offshell_term"] is True
    assert short_by_name["C_full_includes_extra_offshell_term"] is False
    assert short["C_full_includes_extra_offshell_term"] is False
    assert abs(ingredients["A1"] + mp.mpf(25) / 36) < mp.mpf("1e-70")
    assert abs(full_by_name["A1"] + mp.mpf(25) / 36) < mp.mpf("1e-70")
    assert abs(short_by_name["A1"] - mp.mpf(5) / 6) < mp.mpf("1e-70")
    assert abs(short["A1"] - mp.mpf(5) / 6) < mp.mpf("1e-70")


def test_closed_form_values_are_callable():
    point = get_point("equal_mass_offshell_positive")
    eps = mp.mpf("0.1")

    assert abs(product_closed_form(point, eps)) > 0
    assert abs(residue_closed_form(point, eps)) > 0
    assert abs(residue_closed_form(point, eps, scale="short")) > 0
    assert abs(residue_closed_form(point, eps, scale="short") - residue_closed_form(point, eps, scale="full")) > 1


def test_four_mass_offshell_branch_safe_point_is_branch_safe_for_residue_scales():
    for point_name in [
        "four_mass_offshell_branch_safe_A",
        "four_mass_offshell_branch_safe_B",
        "four_mass_offshell_branch_safe_C",
    ]:
        point = get_point(point_name)
        short = residue_ingredients(point, scale="short")
        full = residue_ingredients(point, scale="full")

        assert point.in_support()
        assert point.x != point.y
        assert point.p2 != 0
        assert len({point.ml2, point.Ml2, point.mk2, point.Mk2}) > 1
        assert short["Delta0"] > 0
        assert short["A1"] > 0
        assert full["A1"] > 0
        assert abs(short["A1"] - full["A1"]) > mp.mpf("1e-30")
        assert abs(short["z"] - full["z"]) > mp.mpf("1e-30")


def test_new_four_mass_branch_safe_points_have_clean_z_values():
    for point_name in [
        "four_mass_offshell_branch_safe_B",
        "four_mass_offshell_branch_safe_C",
    ]:
        point = get_point(point_name)
        short = residue_ingredients(point, scale="short")
        full = residue_ingredients(point, scale="full")

        assert 0 < short["z"] < 1
        assert 0 < full["z"] < 1
