import math

import numpy as np


def omega(m):
    return 2 * math.pi ** (float(m) / 2) / math.gamma(float(m) / 2)


def method3_prefactor(point, eps):
    eps = float(eps)
    d = 4.0 - 2.0 * eps
    n = d - 2
    return (
        (2 * math.pi * 1j) ** 2
        * omega(n - 1)
        * omega(n - 2)
        / (1j * math.pi ** (d / 2)) ** 2
        / (float(point.x) * float(point.y))
    )


def A_equal_mass_pperp0(point, ell, k):
    x = float(point.x)
    y = float(point.y)
    p_plus = float(point.p_plus)
    p_minus = float(point.p_minus)
    M2 = float(point.M2)
    return (
        ((x - p_plus) / y) * ell**2
        + ((y - p_plus) / x) * k**2
        + M2 * (((x + y) * (x + y - p_plus) - x * y) / (x * y))
        - (x + y - p_plus) * p_minus
    )


def D2_equal_mass_pperp0(point, ell, eta):
    y = float(point.y)
    p_plus = float(point.p_plus)
    p_minus = float(point.p_minus)
    M2 = float(point.M2)
    return (
        (y - p_plus) * ((ell**2 + M2) / y - p_minus)
        - (ell**2 + M2)
        + 1j * eta * p_plus / y
    )


def gamma_xy(point):
    return (
        (point.x + point.y)
        * (point.p_plus - point.x - point.y)
        / (point.x * point.y)
    )


def _radical_inverse(indices, base):
    n = indices.copy()
    out = np.zeros_like(n, dtype=np.float64)
    factor = 1.0 / base
    while np.any(n):
        out += (n % base) * factor
        n //= base
        factor /= base
    return out


def _halton(n, dim, start=1):
    bases = [2, 3, 5, 7, 11, 13]
    if dim > len(bases):
        raise ValueError("Increase the Halton base list for higher dimensions.")
    indices = np.arange(start, start + n, dtype=np.uint64)
    return np.column_stack([_radical_inverse(indices, b) for b in bases[:dim]])


def _qmc_points(n, scramble, seed):
    try:
        from scipy.stats import qmc  # type: ignore

        m = int(round(math.log2(n)))
        if 2**m != n:
            raise ValueError("Sobol requires N to be a power of two.")
        sampler = qmc.Sobol(d=3, scramble=scramble, seed=seed)
        return sampler.random_base2(m)
    except ImportError:
        points = _halton(n, 3)
        if scramble:
            rng = np.random.default_rng(seed)
            points = (points + rng.random(3)) % 1.0
        return points


def method3_stripped(point, eps, eta, n_samples, rho=None, scrambles=4, seed=12345):
    if not point.in_support():
        return 0j, 0.0

    if point.p_perp2 != 0:
        raise ValueError("Method 3 currently implements p_perp2 = 0.")

    rho = math.sqrt(float(point.M2)) if rho is None else float(rho)
    eps_f = float(eps)
    eta_f = float(eta)
    gamma_f = float(gamma_xy(point))
    values = []

    for scramble_id in range(scrambles):
        pts = _qmc_points(n_samples, scramble=True, seed=seed + scramble_id)
        u = np.clip(pts[:, 0], np.finfo(float).tiny, 1.0 - np.finfo(float).eps)
        v = np.clip(pts[:, 1], np.finfo(float).tiny, 1.0 - np.finfo(float).eps)
        w = np.clip(pts[:, 2], np.finfo(float).tiny, 1.0 - np.finfo(float).eps)

        k = rho * u / (1.0 - u)
        ell = rho * v / (1.0 - v)
        theta = np.pi * w

        jac = rho / (1.0 - u) ** 2 * rho / (1.0 - v) ** 2 * np.pi
        measure = k ** (1.0 - 2.0 * eps_f) * ell ** (1.0 - 2.0 * eps_f)
        measure *= np.sin(theta) ** (-2.0 * eps_f)

        d2 = D2_equal_mass_pperp0(point, ell, eta_f)
        d4 = (
            A_equal_mass_pperp0(point, ell, k)
            - 2.0 * ell * k * np.cos(theta)
            + 1j * eta_f * gamma_f
        )

        values.append(np.mean(jac * measure / (d2 * d4)))

    mean = complex(np.mean(values))
    if len(values) < 2:
        return mean, float("nan")
    err = float(np.std(values, ddof=1) / math.sqrt(len(values)))
    return mean, err


def method3_normalized(point, eps, eta, n_samples, rho=None, scrambles=4, seed=12345):
    stripped, error = method3_stripped(
        point=point,
        eps=eps,
        eta=eta,
        n_samples=n_samples,
        rho=rho,
        scrambles=scrambles,
        seed=seed,
    )
    pref = complex(method3_prefactor(point, eps))
    return pref * stripped, abs(pref) * error, stripped
