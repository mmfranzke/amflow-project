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
    short = residue_ingredients(point, short_C_diagnostic=True)

    assert ingredients["C_full_includes_extra_offshell_term"] is True
    assert short["C_full_includes_extra_offshell_term"] is False
    assert abs(ingredients["A1"] + mp.mpf(25) / 36) < mp.mpf("1e-70")
    assert abs(short["A1"] - mp.mpf(5) / 6) < mp.mpf("1e-70")


def test_closed_form_values_are_callable():
    point = get_point("equal_mass_offshell_positive")
    eps = mp.mpf("0.1")

    assert abs(product_closed_form(point, eps)) > 0
    assert abs(residue_closed_form(point, eps)) > 0

