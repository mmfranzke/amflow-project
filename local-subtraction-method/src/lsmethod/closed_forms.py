import cmath
import math

import mpmath as mp

from .kernels import (
    l_residue_closed_ingredients,
    method1_l_residue_closed,
)

EULER_GAMMA = 0.577215664901532860606512090082402431
ZETA3 = 1.20205690315959428539973816151144999


def _as_mpf(value):
    if hasattr(value, "numerator") and hasattr(value, "denominator"):
        return mp.mpf(value.numerator) / value.denominator
    return mp.mpf(value)


def _as_float(value):
    return float(_as_mpf(value))


def _delta_log(delta, branch="-i0"):
    if delta < 0:
        phase = -math.pi if branch == "-i0" else math.pi
        return math.log(abs(float(delta))) + 1j * phase
    return math.log(float(delta))


def product_delta0(point):
    return (1 - point.x) * point.ml2 + point.x * point.Ml2 - point.x * (1 - point.x) * point.p2


def product_delta1(point):
    lam = point.lam
    mu_k_lam = (1 - lam) * point.mk2 + lam * point.Mk2
    return -mu_k_lam + lam * (1 - lam) * point.Ml2 + lam * (1 - lam * point.x) * point.p2


def product_closed_form(point, eps, eta=mp.mpf("1e-30"), branch="-i0"):
    eps_f = _as_float(eps)
    if branch == "-i0":
        d0 = complex(float(product_delta0(point)), -float(eta))
        d1 = complex(float(product_delta1(point)), -float(eta))
        return math.gamma(eps_f) ** 2 / float(point.p_plus - point.x) * d0 ** (-eps_f) * d1 ** (-eps_f)
    d0_log = _delta_log(product_delta0(point), branch)
    d1_log = _delta_log(product_delta1(point), branch)
    return math.gamma(eps_f) ** 2 / float(point.p_plus - point.x) * cmath.exp(-eps_f * (d0_log + d1_log))


method12_closed = product_closed_form


def residue_closed_form(point, eps, eta=mp.mpf("1e-30"), short_C_diagnostic=False):
    return complex(
        method1_l_residue_closed(
            point,
            _as_mpf(eps),
            eta=_as_mpf(eta),
            short_C_diagnostic=short_C_diagnostic,
        )
    )


def residue_ingredients(point, short_C_diagnostic=False):
    return l_residue_closed_ingredients(point, short_C_diagnostic=short_C_diagnostic)


def product_laurent_coefficients(point, max_order=2, eta=mp.mpf("1e-30"), min_power=-2, branch="-i0"):
    if min_power < -2:
        raise ValueError("product coefficients are implemented from c[-2] upward.")
    if max_order > 2:
        raise ValueError("product coefficients are implemented through c[2].")

    s_log = _delta_log(product_delta0(point), branch) + _delta_log(product_delta1(point), branch)
    q1 = -(2.0 * EULER_GAMMA + s_log)
    q2 = math.pi**2 / 6.0
    q3 = -2.0 * ZETA3 / 3.0
    q4 = math.pi**4 / 180.0
    exp_coeffs = {
        0: 1.0 + 0j,
        1: q1,
        2: q2 + q1**2 / 2.0,
        3: q3 + q1 * q2 + q1**3 / 6.0,
        4: q4 + q1 * q3 + q2**2 / 2.0 + q1**2 * q2 / 2.0 + q1**4 / 24.0,
    }
    X = float(point.p_plus - point.x)
    return {power: exp_coeffs[power + 2] / X for power in range(min_power, max_order + 1)}


def residue_laurent_coefficients(point, max_order=2, eta=mp.mpf("1e-40"), min_power=-2, short_C_diagnostic=False):
    if min_power < -2:
        raise ValueError("residue coefficients are implemented from c[-2] upward.")
    mp.mp.dps = 100
    degree = max_order + 10
    poly_powers = list(range(0, degree + 1))
    eps_nodes = [
        mp.mpf(1) / n
        for n in (30, 35, 42, 50, 60, 72, 86, 103, 124, 150, 180, 220, 270, 330, 410, 520, 660, 850, 1100, 1450, 1900, 2500)
    ][: degree + 6]
    matrix = mp.matrix([[eps ** power for power in poly_powers] for eps in eps_nodes])
    vector = mp.matrix([
        eps**2
        * method1_l_residue_closed(
            point,
            eps,
            eta=_as_mpf(eta),
            short_C_diagnostic=short_C_diagnostic,
        )
        for eps in eps_nodes
    ])
    solved = mp.lu_solve(matrix.T * matrix, matrix.T * vector)
    return {
        power: complex(float(mp.re(solved[power + 2])), float(mp.im(solved[power + 2])))
        for power in range(min_power, max_order + 1)
    }

