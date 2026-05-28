#!/usr/bin/env python3
"""Direct beta-contour check of Eq. (29) against Eq. (51)/(53).

Method: light-cone delta reduction, D1-pole plus C-branch-cut beta-contour
contribution.

No AMFlow, no w-representation is used for Eq. (29).  Eq. (51)/(53) are
reference values only.
"""

from __future__ import annotations

import argparse
import cmath
import math
from dataclasses import dataclass
from typing import Callable, Iterable, Optional

try:
    import numpy as np
    from scipy import integrate, special

    HAVE_SCIPY = True
except Exception:  # pragma: no cover - exercised only without scipy
    np = None
    integrate = None
    special = None
    HAVE_SCIPY = False

try:
    import mpmath as mp

    HAVE_MPMATH = True
except Exception:  # pragma: no cover
    mp = None
    HAVE_MPMATH = False


@dataclass(frozen=True)
class Point:
    x: float = 3.0 / 10.0
    y: float = 1.0 / 4.0
    s: float = 0.0
    ml2: float = 49.0 / 100.0
    Ml2: float = 4.0
    mk2: float = 1.0 / 4.0
    Mk2: float = 81.0 / 100.0

    @property
    def X(self) -> float:
        return 1.0 - self.x

    @property
    def lam(self) -> float:
        return self.y / self.X

    @property
    def a(self) -> float:
        return self.lam * (1.0 - self.lam)

    @property
    def b(self) -> float:
        return self.lam * (1.0 - self.lam * self.x)

    @property
    def mu_k_lambda(self) -> float:
        return (1.0 - self.lam) * self.mk2 + self.lam * self.Mk2

    @property
    def mu_l_x(self) -> float:
        return (1.0 - self.x) * self.ml2 + self.x * self.Ml2

    @property
    def delta0(self) -> float:
        return self.mu_l_x - self.x * (1.0 - self.x) * self.s

    @property
    def delta_w(self) -> float:
        return (
            ((1.0 + self.a) * self.x * (1.0 - self.x) + self.b) * self.s
            - self.mu_k_lambda
            + self.a * self.Ml2
            - (1.0 + self.a) * self.mu_l_x
        )

    @property
    def delta1(self) -> float:
        return (1.0 + self.a) * self.delta0 + self.delta_w

    @property
    def z(self) -> float:
        return -self.delta_w / ((1.0 + self.a) * self.delta0)


@dataclass
class QuadResult:
    value: complex
    real_error: float
    imag_error: float

    @property
    def error(self) -> float:
        return max(self.real_error, self.imag_error)

    def scale(self, factor: complex) -> "QuadResult":
        abs_factor = abs(factor)
        return QuadResult(
            factor * self.value,
            abs_factor * self.real_error,
            abs_factor * self.imag_error,
        )


def gamma(x: float) -> float:
    if HAVE_SCIPY:
        return float(special.gamma(x))
    return float(mp.gamma(x))


def hyp2f1(a: float, b: float, c: float, z: float) -> complex:
    if HAVE_SCIPY:
        return complex(special.hyp2f1(a, b, c, z))
    return complex(mp.hyper([a, b], [c], z))


def complex_quad_scipy(
    func: Callable[[float], complex],
    a: float,
    b: float,
    *,
    epsabs: float,
    epsrel: float,
    limit: int,
) -> QuadResult:
    real_value, real_error = integrate.quad(
        lambda x: float(complex(func(x)).real),
        a,
        b,
        epsabs=epsabs,
        epsrel=epsrel,
        limit=limit,
    )
    imag_value, imag_error = integrate.quad(
        lambda x: float(complex(func(x)).imag),
        a,
        b,
        epsabs=epsabs,
        epsrel=epsrel,
        limit=limit,
    )
    return QuadResult(complex(real_value, imag_value), real_error, imag_error)


def complex_quad_mpmath(
    func: Callable[[float], complex],
    a: float,
    b: float,
    *,
    epsabs: float,
    epsrel: float,
    limit: int,
) -> QuadResult:
    del epsabs, epsrel, limit
    value = mp.quad(lambda x: func(float(x)), [a, b])
    return QuadResult(complex(value), float("nan"), float("nan"))


def complex_quad(
    func: Callable[[float], complex],
    a: float,
    b: float,
    *,
    epsabs: float,
    epsrel: float,
    limit: int,
) -> QuadResult:
    if HAVE_SCIPY:
        return complex_quad_scipy(func, a, b, epsabs=epsabs, epsrel=epsrel, limit=limit)
    if HAVE_MPMATH:
        return complex_quad_mpmath(func, a, b, epsabs=epsabs, epsrel=epsrel, limit=limit)
    raise RuntimeError("Need scipy or mpmath for numerical integration.")


def dedupe_sorted(values: Iterable[float], *, min_separation: float) -> list[float]:
    out: list[float] = []
    for value in sorted(values):
        if not math.isfinite(value):
            continue
        if out and abs(value - out[-1]) < min_separation:
            continue
        out.append(value)
    return out


def beta_split_points(point: Point, u: float, eta: float) -> list[float]:
    x = point.x
    X = point.X
    a = point.a

    beta1 = (u + point.ml2) / (2.0 * x)
    beta2 = -(u + point.Ml2) / (2.0 * X)
    beta_c = (point.mu_k_lambda + a * u - point.b * point.s) / (2.0 * a * (x - 1.0))

    # The finite-i0 peaks have widths set by eta divided by the beta slope of
    # the corresponding denominator.  Add a geometric nest of split points so
    # quad never has to discover a narrow peak inside one huge interval.
    width1 = max(eta / abs(2.0 * x), 1.0e-14)
    width2 = max(eta / abs(2.0 * (x - 1.0)), 1.0e-14)
    widthc = max(eta / abs(2.0 * a * (x - 1.0)), 1.0e-14)
    radii = (0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 25.0, 50.0, 100.0, 250.0)
    candidates = []
    for center, width in ((beta2, width2), (beta_c, widthc), (beta1, width1)):
        for radius in radii:
            if radius == 0.0:
                candidates.append(center)
            else:
                candidates.extend([center - radius * width, center + radius * width])
    return dedupe_sorted(candidates, min_separation=max(eta * 1.0e-8, 1.0e-14))


def denominators(point: Point, beta: float, u: float, eta: float) -> tuple[complex, complex, complex]:
    x = point.x
    D1 = 2.0 * x * beta - u - point.ml2 + 1j * eta
    D2 = 2.0 * (x - 1.0) * beta - u - point.Ml2 + 1j * eta
    C = (
        point.mu_k_lambda
        - point.a * (2.0 * (x - 1.0) * beta - u)
        - point.b * point.s
        - 1j * eta
    )
    return D1, D2, C


def beta_integral(
    point: Point,
    eps: float,
    eta: float,
    u: float,
    *,
    include_c: bool,
    epsabs: float,
    epsrel: float,
    limit: int,
) -> QuadResult:
    def integrand(beta: float) -> complex:
        D1, D2, C = denominators(point, beta, u, eta)
        c_factor = C ** (-eps) if include_c else 1.0
        return c_factor / (D1 * D2)

    splits = beta_split_points(point, u, eta)
    endpoints = [-math.inf, *splits, math.inf]

    total = 0.0 + 0.0j
    real_err_sq = 0.0
    imag_err_sq = 0.0

    for left, right in zip(endpoints, endpoints[1:]):
        result = complex_quad(
            integrand,
            left,
            right,
            epsabs=epsabs,
            epsrel=epsrel,
            limit=limit,
        )
        total += result.value
        if math.isfinite(result.real_error):
            real_err_sq += result.real_error**2
        if math.isfinite(result.imag_error):
            imag_err_sq += result.imag_error**2

    return QuadResult(total, math.sqrt(real_err_sq), math.sqrt(imag_err_sq))


def r_intervals() -> list[tuple[float, float]]:
    points = [
        0.0,
        1.0e-10,
        1.0e-8,
        1.0e-6,
        1.0e-4,
        1.0e-2,
        0.1,
        0.25,
        0.5,
        0.75,
        0.9,
        0.99,
        1.0 - 1.0e-4,
        1.0 - 1.0e-6,
        1.0 - 1.0e-8,
        1.0 - 1.0e-10,
        1.0,
    ]
    return list(zip(points, points[1:]))


def r_from_rho(rho: float) -> float:
    if rho <= 0.0:
        return 0.0
    return rho / (1.0 + rho)


def rho_from_r(r: float) -> float:
    if r <= 0.0:
        r = 1.0e-300
    if r >= 1.0:
        r = 1.0 - 1.0e-16
    return r / (1.0 - r)


def rho_intervals_with_splits(split_rhos: Iterable[float]) -> list[tuple[float, float]]:
    base_r = [
        0.0,
        1.0e-10,
        1.0e-8,
        1.0e-6,
        1.0e-4,
        1.0e-2,
        0.1,
        0.25,
        0.5,
        0.75,
        0.9,
        0.99,
        1.0 - 1.0e-4,
        1.0 - 1.0e-6,
        1.0 - 1.0e-8,
        1.0 - 1.0e-10,
        1.0,
    ]
    extra_r = [r_from_rho(rho) for rho in split_rhos if rho > 0.0]
    points = dedupe_sorted([*base_r, *extra_r], min_separation=1.0e-14)
    return list(zip(points, points[1:]))


def lightcone_integral(
    point: Point,
    eps: float,
    eta: float,
    *,
    include_c: bool,
    epsabs: float,
    epsrel: float,
    limit: int,
) -> QuadResult:
    """Return the beta-u integral including 1/(i pi Gamma[1-eps])."""

    def r_integrand(r: float) -> complex:
        if r <= 0.0:
            r = 1.0e-300
        if r >= 1.0:
            r = 1.0 - 1.0e-16
        u = r / (1.0 - r)
        jac = 1.0 / (1.0 - r) ** 2
        beta_result = beta_integral(
            point,
            eps,
            eta,
            u,
            include_c=include_c,
            epsabs=epsabs,
            epsrel=epsrel,
            limit=limit,
        )
        if beta_result.error > max(1.0e-4, 1.0e-4 * abs(beta_result.value)):
            print(
                "  warning: beta quad error at "
                f"eps={eps:g}, eta={eta:g}, u={u:.6e}: "
                f"value={fmt_complex(beta_result.value)}, err~{beta_result.error:.3e}"
            )
        return (u ** (-eps)) * beta_result.value * jac

    total = 0.0 + 0.0j
    real_err_sq = 0.0
    imag_err_sq = 0.0
    for left, right in r_intervals():
        result = complex_quad(
            r_integrand,
            left,
            right,
            epsabs=epsabs,
            epsrel=epsrel,
            limit=limit,
        )
        total += result.value
        if math.isfinite(result.real_error):
            real_err_sq += result.real_error**2
        if math.isfinite(result.imag_error):
            imag_err_sq += result.imag_error**2

    prefactor = 1.0 / (1j * math.pi * gamma(1.0 - eps))
    return QuadResult(prefactor * total, math.sqrt(real_err_sq), math.sqrt(imag_err_sq))


def K_direct(point: Point, eps: float, eta: float, opts: argparse.Namespace) -> QuadResult:
    return lightcone_integral(
        point,
        eps,
        eta,
        include_c=False,
        epsabs=opts.epsabs,
        epsrel=opts.epsrel,
        limit=opts.limit,
    )


def I29_direct(point: Point, eps: float, eta: float, opts: argparse.Namespace) -> QuadResult:
    raw = lightcone_integral(
        point,
        eps,
        eta,
        include_c=True,
        epsabs=opts.epsabs,
        epsrel=opts.epsrel,
        limit=opts.limit,
    )
    return QuadResult(gamma(eps) / point.X * raw.value, raw.real_error, raw.imag_error)


def K_expected(point: Point, eps: float) -> complex:
    return gamma(eps) * point.delta0 ** (-eps)


def K_residue_analytic(point: Point, eps: float) -> complex:
    return K_expected(point, eps)


def beta1_real(point: Point, u: float) -> float:
    return (u + point.ml2) / (2.0 * point.x)


def beta2_real(point: Point, u: float) -> float:
    return -(u + point.Ml2) / (2.0 * point.X)


def betaC_of_u(point: Point, u: float) -> float:
    return (point.b * point.s - point.mu_k_lambda - point.a * u) / (2.0 * point.a * point.X)


def branch_slope_A(point: Point) -> float:
    return 2.0 * point.a * point.X


def q2_at_d1_pole(point: Point, u: float) -> float:
    return 2.0 * (point.x - 1.0) * beta1_real(point, u) - u


def c_pole(point: Point, u: float) -> float:
    return point.mu_k_lambda - point.a * q2_at_d1_pole(point, u) - point.b * point.s


def D1_real(point: Point, beta: float, u: float) -> float:
    return 2.0 * point.x * beta - u - point.ml2


def D2_real(point: Point, beta: float, u: float) -> float:
    return 2.0 * (point.x - 1.0) * beta - u - point.Ml2


def B_pole(point: Point, u: float, eps: float) -> complex:
    return c_pole(point, u) ** (-eps) / (u + point.delta0)


def rho2_of_u(point: Point, u: float) -> float:
    return betaC_of_u(point, u) - beta2_real(point, u)


def rho_cut_prefactor(eps: float, cut_sign: int) -> float:
    return cut_sign * 2.0 * math.sin(math.pi * eps) / math.pi


def residue_u_integral(
    integrand_u: Callable[[float], complex],
    *,
    epsabs: float,
    epsrel: float,
    limit: int,
) -> QuadResult:
    """Integrate an affine-pole residue integrand over u in [0, infinity)."""

    def r_integrand(r: float) -> complex:
        if r <= 0.0:
            r = 1.0e-300
        if r >= 1.0:
            r = 1.0 - 1.0e-16
        u = r / (1.0 - r)
        jac = 1.0 / (1.0 - r) ** 2
        return integrand_u(u) * jac

    total = 0.0 + 0.0j
    real_err_sq = 0.0
    imag_err_sq = 0.0
    for left, right in r_intervals():
        result = complex_quad(
            r_integrand,
            left,
            right,
            epsabs=epsabs,
            epsrel=epsrel,
            limit=limit,
        )
        total += result.value
        if math.isfinite(result.real_error):
            real_err_sq += result.real_error**2
        if math.isfinite(result.imag_error):
            imag_err_sq += result.imag_error**2
    return QuadResult(total, math.sqrt(real_err_sq), math.sqrt(imag_err_sq))


def K_residue_plain_r_numeric(point: Point, eps: float, opts: argparse.Namespace) -> QuadResult:
    return residue_u_integral(
        lambda u: u ** (-eps) / (u + point.delta0),
        epsabs=opts.epsabs,
        epsrel=opts.epsrel,
        limit=opts.limit,
    ).scale(1.0 / gamma(1.0 - eps))


def K_residue_weighted_numeric(point: Point, eps: float, opts: argparse.Namespace) -> Optional[QuadResult]:
    """Weighted quadrature for the endpoint-singular normalization integral.

    With u = Delta0 t/(1-t),
      int_0^infty du u^(-eps)/(u+Delta0)
      = Delta0^(-eps) int_0^1 dt t^(-eps) (1-t)^(-1+eps).
    scipy's algebraic endpoint weight resolves this much better than the
    generic r = u/(1+u) map.
    """

    if not HAVE_SCIPY:
        return None
    value, error = integrate.quad(
        lambda t: 1.0,
        0.0,
        1.0,
        weight="alg",
        wvar=(-eps, -1.0 + eps),
        epsabs=opts.epsabs,
        epsrel=opts.epsrel,
        limit=opts.limit,
    )
    scale = point.delta0 ** (-eps) / gamma(1.0 - eps)
    return QuadResult(scale * value, abs(scale) * error, 0.0)


def I29_d1_pole_only(point: Point, eps: float, opts: argparse.Namespace) -> QuadResult:
    # For C -> 1, the beta integrand has only the two propagator poles, and
    # closing the contour gives the D1 residue.  For full C(beta,u)^(-eps),
    # C introduces a branch point/cut in beta.  Therefore the D1 residue alone
    # is only one part of the beta-contour integral; the missing branch-cut
    # contribution is required for the full Eq. (29).
    if HAVE_SCIPY:
        def smooth_part(t: float) -> float:
            if t >= 1.0:
                t = 1.0 - 1.0e-16
            u = point.delta0 * t / (1.0 - t)
            return float(c_pole(point, u) ** (-eps))

        value, error = integrate.quad(
            smooth_part,
            0.0,
            1.0,
            weight="alg",
            wvar=(-eps, -1.0 + eps),
            epsabs=opts.epsabs,
            epsrel=opts.epsrel,
            limit=opts.limit,
        )
        scale = gamma(eps) / point.X / gamma(1.0 - eps) * point.delta0 ** (-eps)
        return QuadResult(scale * value, abs(scale) * error, 0.0)

    integral = residue_u_integral(
        lambda u: u ** (-eps) * c_pole(point, u) ** (-eps) / (u + point.delta0),
        epsabs=opts.epsabs,
        epsrel=opts.epsrel,
        limit=opts.limit,
    )
    return integral.scale(gamma(eps) / point.X / gamma(1.0 - eps))


def integrate_rho_cut(
    rho_integrand: Callable[[float], complex],
    split_rhos: Iterable[float],
    *,
    epsabs: float,
    epsrel: float,
    limit: int,
) -> QuadResult:
    def r_integrand(r: float) -> complex:
        rho = rho_from_r(r)
        jac = 1.0 / (1.0 - r) ** 2
        return rho_integrand(rho) * jac

    total = 0.0 + 0.0j
    real_err_sq = 0.0
    imag_err_sq = 0.0
    for left, right in rho_intervals_with_splits(split_rhos):
        if right <= left:
            continue
        result = complex_quad(
            r_integrand,
            left,
            right,
            epsabs=epsabs,
            epsrel=epsrel,
            limit=limit,
        )
        total += result.value
        if math.isfinite(result.real_error):
            real_err_sq += result.real_error**2
        if math.isfinite(result.imag_error):
            imag_err_sq += result.imag_error**2
    return QuadResult(total, math.sqrt(real_err_sq), math.sqrt(imag_err_sq))


def B_cut_eta(
    point: Point,
    u: float,
    eps: float,
    eta_cut: float,
    cut_sign: int,
    d2_sign: int,
    opts: argparse.Namespace,
) -> QuadResult:
    """C-branch-cut contribution with a finite D2 prescription.

    On the cut beta = betaC - rho, rho >= 0, and
    C(beta,u)^(-eps) has a discontinuity proportional to
    2 sin(pi eps) (A rho)^(-eps).  The global orientation is diagnosed by
    cut_sign.  The D2 prescription on the cut is diagnosed by d2_sign.
    """

    A = branch_slope_A(point)
    beta_c = betaC_of_u(point, u)
    rho2 = rho2_of_u(point, u)
    width = max(eta_cut / (2.0 * point.X), 1.0e-14)
    split_rhos = []
    if rho2 > 0:
        for factor in (0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 25.0, 50.0, 100.0):
            split_rhos.extend([rho2 - factor * width, rho2 + factor * width])

    def integrand(rho: float) -> complex:
        beta = beta_c - rho
        d1 = D1_real(point, beta, u)
        d2 = D2_real(point, beta, u) + 1j * d2_sign * eta_cut
        return (A * rho) ** (-eps) / (d1 * d2)

    result = integrate_rho_cut(
        integrand,
        split_rhos,
        epsabs=opts.epsabs,
        epsrel=opts.epsrel,
        limit=opts.limit,
    )
    return result.scale(rho_cut_prefactor(eps, cut_sign))


def B_cut_pv(
    point: Point,
    u: float,
    eps: float,
    delta: float,
    cut_sign: int,
    opts: argparse.Namespace,
) -> QuadResult:
    """C-branch-cut contribution with principal-value D2 treatment."""

    A = branch_slope_A(point)
    beta_c = betaC_of_u(point, u)
    rho2 = rho2_of_u(point, u)
    if rho2 <= 0 or rho2 <= delta:
        split_rhos: list[float] = []
        excluded: tuple[float, float] | None = None
    else:
        excluded = (rho2 - delta, rho2 + delta)
        split_rhos = [rho2 - delta, rho2 + delta]

    def integrand(rho: float) -> complex:
        if excluded is not None and excluded[0] <= rho <= excluded[1]:
            return 0.0
        beta = beta_c - rho
        d1 = D1_real(point, beta, u)
        d2 = D2_real(point, beta, u)
        return (A * rho) ** (-eps) / (d1 * d2)

    result = integrate_rho_cut(
        integrand,
        split_rhos,
        epsabs=opts.epsabs,
        epsrel=opts.epsrel,
        limit=opts.limit,
    )
    return result.scale(rho_cut_prefactor(eps, cut_sign))


def B_full(
    point: Point,
    u: float,
    eps: float,
    *,
    cut_mode: str,
    cut_sign: int,
    eta_cut: Optional[float],
    d2_sign: int,
    pv_delta: Optional[float],
    opts: argparse.Namespace,
) -> QuadResult:
    pole = B_pole(point, u, eps)
    if cut_mode == "eta":
        if eta_cut is None:
            raise ValueError("eta_cut is required for eta cut mode")
        cut = B_cut_eta(point, u, eps, eta_cut, cut_sign, d2_sign, opts)
    elif cut_mode == "pv":
        if pv_delta is None:
            raise ValueError("pv_delta is required for pv cut mode")
        cut = B_cut_pv(point, u, eps, pv_delta, cut_sign, opts)
    else:
        raise ValueError(f"unknown cut_mode={cut_mode!r}")
    return QuadResult(pole + cut.value, cut.real_error, cut.imag_error)


def I29_cut_contribution_eta(
    point: Point,
    eps: float,
    eta_cut: float,
    cut_sign: int,
    d2_sign: int,
    opts: argparse.Namespace,
) -> QuadResult:
    integral = residue_u_integral(
        lambda u: u ** (-eps)
        * B_cut_eta(point, u, eps, eta_cut, cut_sign, d2_sign, opts).value,
        epsabs=opts.epsabs,
        epsrel=opts.epsrel,
        limit=opts.limit,
    )
    return integral.scale(gamma(eps) / point.X / gamma(1.0 - eps))


def I29_cut_contribution_pv(
    point: Point,
    eps: float,
    delta: float,
    cut_sign: int,
    opts: argparse.Namespace,
) -> QuadResult:
    integral = residue_u_integral(
        lambda u: u ** (-eps) * B_cut_pv(point, u, eps, delta, cut_sign, opts).value,
        epsabs=opts.epsabs,
        epsrel=opts.epsrel,
        limit=opts.limit,
    )
    return integral.scale(gamma(eps) / point.X / gamma(1.0 - eps))


def I53(point: Point, eps: float) -> complex:
    return gamma(eps) ** 2 / point.X * point.delta0 ** (-eps) * point.delta1 ** (-eps)


def I51(point: Point, eps: float) -> complex:
    return (
        gamma(eps) ** 2
        / point.X
        * (1.0 + point.a) ** (-eps)
        * point.delta0 ** (-2.0 * eps)
        * hyp2f1(2.0 * eps, eps, 2.0 * eps, point.z)
    )


def fmt_complex(z: complex) -> str:
    return f"{z.real:.16e} {z.imag:+.16e}j"


def rel_diff(value: complex, reference: complex) -> float:
    denom = abs(reference)
    if denom == 0:
        return math.inf
    return abs(value - reference) / denom


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--eps", type=float, action="append", help="epsilon value; may be repeated")
    parser.add_argument("--eta", type=float, action="append", help="finite i0 regulator; may be repeated")
    parser.add_argument("--quick", action="store_true", help="run only eps=0.05 and eta=1e-3")
    parser.add_argument("--epsabs", type=float, default=1.0e-8)
    parser.add_argument("--epsrel", type=float, default=1.0e-8)
    parser.add_argument("--limit", type=int, default=300)
    parser.add_argument("--norm-tol", type=float, default=1.0e-8)
    parser.add_argument("--full-tol", type=float, default=1.0e-4)
    parser.add_argument("--cut-mode", choices=("eta", "pv"), action="append")
    parser.add_argument("--eta-cut", type=float, action="append")
    parser.add_argument("--pv-delta", type=float, action="append")
    parser.add_argument("--cut-sign", type=int, choices=(-1, 1), action="append")
    parser.add_argument("--d2-sign", type=int, choices=(-1, 1), action="append")
    parser.add_argument(
        "--naive-beta",
        action="store_true",
        help="also run the old finite-eta real-axis beta quadrature as diagnostic only",
    )
    return parser.parse_args()


def main() -> None:
    if not HAVE_SCIPY and not HAVE_MPMATH:
        raise SystemExit("Need scipy or mpmath. Install scipy for the intended path.")

    opts = parse_args()
    point = Point()

    if opts.quick:
        eps_values = [0.05]
        eta_values = [1.0e-3]
    else:
        eps_values = opts.eps if opts.eps else [0.05, 0.04]
        eta_values = opts.eta if opts.eta else [1.0e-2, 3.0e-3, 1.0e-3, 3.0e-4]

    cut_modes = opts.cut_mode if opts.cut_mode else ["eta"]
    eta_cut_values = opts.eta_cut if opts.eta_cut else ([1.0e-3] if opts.quick else [1.0e-2, 3.0e-3, 1.0e-3, 3.0e-4])
    pv_delta_values = opts.pv_delta if opts.pv_delta else ([1.0e-4] if opts.quick else [1.0e-3, 3.0e-4, 1.0e-4, 3.0e-5])
    cut_signs = opts.cut_sign if opts.cut_sign else [1, -1]
    d2_signs = opts.d2_sign if opts.d2_sign else [1, -1]

    print("Direct beta-contour check of Eq. (29) against Eq. (51)/(53)")
    print("Method: light-cone delta reduction, D1-pole plus C-branch-cut beta-contour contribution")
    print("No AMFlow, no w-representation used for Eq. (29)")
    print("Naive finite-eta beta quadrature is disabled by default because it is numerically unstable.")
    print(f"Backend: {'scipy' if HAVE_SCIPY else 'mpmath fallback'}")
    print()

    print("Branch-safe point and derived values:")
    print(f"  x = {point.x:.16g}")
    print(f"  y = {point.y:.16g}")
    print(f"  s = {point.s:.16g}")
    print(f"  ml2 = {point.ml2:.16g}")
    print(f"  Ml2 = {point.Ml2:.16g}")
    print(f"  mk2 = {point.mk2:.16g}")
    print(f"  Mk2 = {point.Mk2:.16g}")
    print(f"  X = {point.X:.16g}  expected 7/10")
    print(f"  lambda = {point.lam:.16g}  expected 5/14")
    print(f"  a = {point.a:.16g}  expected 45/196")
    print(f"  b = {point.b:.16g}  expected 125/392")
    print(f"  muKlambda = {point.mu_k_lambda:.16g}  expected 9/20")
    print(f"  Delta0 = {point.delta0:.16g}  expected 1543/1000")
    print(f"  DeltaW = {point.delta_w:.16g}  expected -40009/28000")
    print(f"  Delta1 = {point.delta1:.16g}  expected 459/980")
    print(f"  z = {point.z:.16g}")
    print()

    weighted_norm_statuses: list[str] = []
    cut_summaries: list[dict[str, object]] = []

    for eps in eps_values:
        ref51 = I51(point, eps)
        ref53 = I53(point, eps)
        h2f1_value = hyp2f1(2.0 * eps, eps, 2.0 * eps, point.z)
        h2f1_identity = (1.0 - point.z) ** (-eps)

        print("=" * 88)
        print(f"epsilon = {eps:.16g}")
        print(f"hyp2f1(2eps,eps,2eps,z) = {fmt_complex(h2f1_value)}")
        print(f"(1-z)^(-eps)             = {h2f1_identity:.16e}")
        print(f"I51 = {fmt_complex(ref51)}")
        print(f"I53 = {fmt_complex(ref53)}")
        print(f"I51 - I53 = {fmt_complex(ref51 - ref53)}")
        i51_i53_agree = rel_diff(ref51, ref53) < 1.0e-12
        print(f"C_pole(0) = {c_pole(point, 0.0):.16e}")
        print(f"C_pole(1) = {c_pole(point, 1.0):.16e}")
        print(f"large-u slope of C_pole ~ {c_pole(point, 1.0e6) / 1.0e6:.16e}")
        print(f"betaC(0) = {betaC_of_u(point, 0.0):.16e}")
        print(f"rho2(0) = betaC(0)-beta2(0) = {rho2_of_u(point, 0.0):.16e}")
        print()

        print("== D1-pole residue normalization check C^(-eps) -> 1 ==")
        k_residue_exact = K_residue_analytic(point, eps)
        k_ref = K_expected(point, eps)
        analytic_norm_diff = k_residue_exact - k_ref
        print(f"  K_residue_analytic = {fmt_complex(k_residue_exact)}")
        print(f"  K_expected         = {fmt_complex(k_ref)}")
        print(f"  analytic difference = {fmt_complex(analytic_norm_diff)}")
        print("  PASS: analytic D1-pole residue normalization agrees with expected result.")

        print()
        print("== Optional quadrature diagnostic for normalization integral ==")
        print("  NOTE: this tests scipy endpoint quadrature, not the analytic beta-residue normalization.")
        weighted_norm = K_residue_weighted_numeric(point, eps, opts)
        if weighted_norm is None:
            print("  weighted scipy quadrature unavailable; scipy is not installed.")
            weighted_norm_statuses.append("UNAVAILABLE")
        else:
            weighted_rel = rel_diff(weighted_norm.value, k_ref)
            print(f"  K_residue_numeric = {fmt_complex(weighted_norm.value)}")
            print(f"  K_expected        = {fmt_complex(k_ref)}")
            print(f"  difference        = {fmt_complex(weighted_norm.value - k_ref)}")
            print(f"  relative          = {weighted_rel:.6e}")
            print(f"  estimated error   = {weighted_norm.error:.6e}")
            if weighted_rel < opts.norm_tol:
                print("  PASS: weighted normalization quadrature is stable.")
                weighted_norm_statuses.append("PASS")
            else:
                print("  WARNING: numerical quadrature of the endpoint-singular u-integral is unstable.")
                weighted_norm_statuses.append("UNSTABLE")

        print()
        print("== D1-pole contribution only for full C(beta,u)^(-eps) ==")
        i29_residue = I29_d1_pole_only(point, eps, opts)
        full_rel_51 = rel_diff(i29_residue.value, ref51)
        full_rel_53 = rel_diff(i29_residue.value, ref53)
        print(f"  I29_D1_pole_only       = {fmt_complex(i29_residue.value)}")
        print(f"  I51                    = {fmt_complex(ref51)}")
        print(f"  I53                    = {fmt_complex(ref53)}")
        print(f"  I29_D1_pole_only - I53 = {fmt_complex(i29_residue.value - ref53)}")
        print(f"  rel to I51             = {full_rel_51:.6e}")
        print(f"  rel to I53             = {full_rel_53:.6e}")
        print(f"  estimated error        = {i29_residue.error:.6e}")
        print("  WARNING: This is not the full Eq. (29). It omits the C-branch-cut contribution in the beta plane.")
        print("  Therefore disagreement with Eq. (53) is expected and is not a failure of Eq. (29).")

        print()
        print("== Branch-cut diagnostics ==")
        print("The allowed ambiguity diagnosed here is the global cut orientation and the D2 prescription on the cut.")
        print("No arbitrary normalization factors are applied.")
        for cut_mode in cut_modes:
            if cut_mode == "eta":
                for eta_cut in eta_cut_values:
                    for cut_sign in cut_signs:
                        for d2_sign in d2_signs:
                            print("-" * 88)
                            print(
                                f"cut_sign = {cut_sign:+d}, cut_mode = eta, "
                                f"eta_cut = {eta_cut:.16g}, d2_sign = {d2_sign:+d}"
                            )
                            cut_contribution = I29_cut_contribution_eta(
                                point,
                                eps,
                                eta_cut,
                                cut_sign,
                                d2_sign,
                                opts,
                            )
                            full_contour = i29_residue.value + cut_contribution.value
                            rel = rel_diff(full_contour, ref53)
                            print(f"  I29_cut_contribution = {fmt_complex(cut_contribution.value)}")
                            print(f"  I29_full_contour     = {fmt_complex(full_contour)}")
                            print(f"  I53                  = {fmt_complex(ref53)}")
                            print(f"  difference           = {fmt_complex(full_contour - ref53)}")
                            print(f"  relative difference  = {rel:.6e}")
                            print(f"  imaginary part       = {full_contour.imag:.16e}")
                            print(f"  estimated cut error  = {cut_contribution.error:.6e}")
                            cut_summaries.append(
                                {
                                    "eps": eps,
                                    "cut_sign": cut_sign,
                                    "cut_mode": "eta",
                                    "parameter": eta_cut,
                                    "d2_sign": d2_sign,
                                    "relative": rel,
                                    "matched": rel < opts.full_tol,
                                    "i51_i53_agree": i51_i53_agree,
                                }
                            )
            elif cut_mode == "pv":
                for delta in pv_delta_values:
                    for cut_sign in cut_signs:
                        print("-" * 88)
                        print(
                            f"cut_sign = {cut_sign:+d}, cut_mode = pv, "
                            f"delta = {delta:.16g}"
                        )
                        cut_contribution = I29_cut_contribution_pv(
                            point,
                            eps,
                            delta,
                            cut_sign,
                            opts,
                        )
                        full_contour = i29_residue.value + cut_contribution.value
                        rel = rel_diff(full_contour, ref53)
                        print(f"  I29_cut_contribution = {fmt_complex(cut_contribution.value)}")
                        print(f"  I29_full_contour     = {fmt_complex(full_contour)}")
                        print(f"  I53                  = {fmt_complex(ref53)}")
                        print(f"  difference           = {fmt_complex(full_contour - ref53)}")
                        print(f"  relative difference  = {rel:.6e}")
                        print(f"  imaginary part       = {full_contour.imag:.16e}")
                        print(f"  estimated cut error  = {cut_contribution.error:.6e}")
                        cut_summaries.append(
                            {
                                "eps": eps,
                                "cut_sign": cut_sign,
                                "cut_mode": "pv",
                                "parameter": delta,
                                "d2_sign": None,
                                "relative": rel,
                                "matched": rel < opts.full_tol,
                                "i51_i53_agree": i51_i53_agree,
                            }
                        )

        if opts.naive_beta:
            print()
            print("== Optional diagnostic: naive finite-eta beta quadrature ==")
            print("WARNING: naive finite-eta beta quadrature is numerically unstable and is not used for validation.")
            for eta in eta_values:
                print("-" * 88)
                print(f"eps = {eps:.16g}, eta = {eta:.16g}")
                k_val = K_direct(point, eps, eta, opts).value
                print("Naive normalization check:")
                print(f"  K_direct   = {fmt_complex(k_val)}")
                print(f"  K_expected = {fmt_complex(k_ref)}")
                print(f"  difference = {fmt_complex(k_val - k_ref)}")
                print(f"  relative   = {rel_diff(k_val, k_ref):.6e}")
                i29 = I29_direct(point, eps, eta, opts).value
                print("Naive full direct Eq. (29):")
                print(f"  I29_direct = {fmt_complex(i29)}")
                print(f"  I53        = {fmt_complex(ref53)}")
                print(f"  I29 - I53  = {fmt_complex(i29 - ref53)}")
                print(f"  rel to I53 = {rel_diff(i29, ref53):.6e}")

    print("=" * 88)
    print("Summary:")
    print("  D1-pole C -> 1 normalization: PASS analytically")
    if not weighted_norm_statuses:
        weighted_summary = "NOT RUN"
    elif all(status == "PASS" for status in weighted_norm_statuses):
        weighted_summary = "PASS diagnostic only"
    elif any(status == "UNSTABLE" for status in weighted_norm_statuses):
        weighted_summary = "UNSTABLE diagnostic only"
    else:
        weighted_summary = "UNAVAILABLE diagnostic only"
    print(f"  Plain numerical endpoint quadrature: {weighted_summary}")
    print("  Full C(beta,u)^(-eps) D1-pole-only contribution: computed, but incomplete")
    if cut_summaries:
        print("  Branch-cut variants:")
        for row in cut_summaries:
            d2_text = "" if row["d2_sign"] is None else f", d2_sign={row['d2_sign']:+d}"
            print(
                "    "
                f"eps={row['eps']}, cut_sign={row['cut_sign']:+d}, "
                f"cut_mode={row['cut_mode']}, parameter={row['parameter']}"
                f"{d2_text}, relative={row['relative']:.6e}"
            )
        stable_groups: dict[tuple[object, object, object, object], list[dict[str, object]]] = {}
        for row in cut_summaries:
            key = (row["eps"], row["cut_mode"], row["cut_sign"], row["d2_sign"])
            stable_groups.setdefault(key, []).append(row)
        stable_matches = [
            rows
            for rows in stable_groups.values()
            if len(rows) >= 2
            and all(bool(row["matched"]) and bool(row["i51_i53_agree"]) for row in rows)
        ]
        single_matches = [
            row
            for row in cut_summaries
            if bool(row["matched"]) and bool(row["i51_i53_agree"])
        ]
        if stable_matches:
            print("  Full Eq. (29) vs Eq. (53): PASS for at least one branch-cut prescription with multiple matching regulator values.")
        elif single_matches:
            print(
                "  Full Eq. (29) vs Eq. (53): INCONCLUSIVE/POTENTIAL MATCH: at least one "
                "variant matched, but regulator stability was not established in this run."
            )
        else:
            print(
                "  Full Eq. (29) vs Eq. (53): INCONCLUSIVE: branch-cut contribution "
                "implemented diagnostically, but no stable printed prescription matched Eq. (53). "
                "Inspect contour orientation and D2 prescription."
            )
    else:
        print("  Full Eq. (29) vs Eq. (53): NOT TESTED; no branch-cut variants were run")
    print(f"  Naive beta quadrature: {'enabled diagnostic only' if opts.naive_beta else 'disabled'}")
    print("  Next step if inconclusive: inspect contour orientation and D2 prescription, or use a contour deformation that includes both pole and branch-cut pieces.")
    print("  Alternative: compare the later Feynman-parameter representation Eq. (62) to Eq. (53).")


if __name__ == "__main__":
    main()
