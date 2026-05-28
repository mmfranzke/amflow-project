import math

import numpy as np


def omega(m):
    return 2 * math.pi ** (float(m) / 2) / math.gamma(float(m) / 2)


def method3_prefactor(point, eps):
    components = method3_prefactor_components(point, eps)
    return (
        components["contour_factor"]
        * components["omega_n_minus_1"]
        * components["omega_n_minus_2"]
        / components["loop_normalization"]
        * components["one_over_xy"]
    )


def method3_prefactor_components(point, eps):
    eps = float(eps)
    d = 4.0 - 2.0 * eps
    n = d - 2
    return {
        "d": d,
        "n": n,
        "omega_n_minus_1": omega(n - 1),
        "omega_n_minus_2": omega(n - 2),
        "contour_factor": (2 * math.pi * 1j) ** 2,
        "loop_normalization": (1j * math.pi ** (d / 2)) ** 2,
        "one_over_xy": 1.0 / (float(point.x) * float(point.y)),
    }


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
    Ml2 = float(point.Ml2)
    return (
        (y - p_plus) * ((ell**2 + Ml2) / y - p_minus)
        - (ell**2 + Ml2)
        + 1j * eta * p_plus / y
    )


def gamma_xy(point):
    return (
        (point.x + point.y)
        * (point.p_plus - point.x - point.y)
        / (point.x * point.y)
    )


def mapped_integrand_probe(point, eps, eta, u, v, w, rho=None, cutoff=None):
    eps_f = float(eps)
    eta_f = float(eta)
    u = float(u)
    v = float(v)
    w = float(w)

    if cutoff is None:
        rho = math.sqrt(float(point.M2)) if rho is None else float(rho)
        k = rho * u / (1.0 - u)
        ell = rho * v / (1.0 - v)
        jacobian = rho / (1.0 - u) ** 2 * rho / (1.0 - v) ** 2 * math.pi
        map_name = "infinite"
    else:
        cutoff = float(cutoff)
        k = cutoff * u
        ell = cutoff * v
        jacobian = cutoff * cutoff * math.pi
        map_name = "cutoff"

    theta = math.pi * w
    measure = k ** (1.0 - 2.0 * eps_f) * ell ** (1.0 - 2.0 * eps_f)
    measure *= math.sin(theta) ** (-2.0 * eps_f)
    d2 = D2_equal_mass_pperp0(point, ell, eta_f)
    a_value = A_equal_mass_pperp0(point, ell, k)
    gamma_value = gamma_xy(point)
    d4 = a_value - 2.0 * ell * k * math.cos(theta) + 1j * eta_f * float(gamma_value)
    integrand_without_jacobian = measure / (d2 * d4)

    return {
        "map": map_name,
        "u": u,
        "v": v,
        "w": w,
        "k": k,
        "ell": ell,
        "theta": theta,
        "jacobian": jacobian,
        "D2": d2,
        "A": a_value,
        "Gamma_xy": gamma_value,
        "D4": d4,
        "measure": measure,
        "integrand_without_jacobian": integrand_without_jacobian,
        "integrand_with_jacobian": jacobian * integrand_without_jacobian,
        "D2_real_sign": math.copysign(1.0, d2.real) if d2.real != 0 else 0.0,
        "D4_real_sign": math.copysign(1.0, d4.real) if d4.real != 0 else 0.0,
    }


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


def method3_stripped_cutoff(point, eps, eta, n_samples, cutoff, scrambles=4, seed=12345):
    if not point.in_support():
        return 0j, 0.0

    if point.p_perp2 != 0:
        raise ValueError("Method 3 currently implements p_perp2 = 0.")

    cutoff = float(cutoff)
    values = []

    for scramble_id in range(scrambles):
        pts = _qmc_points(n_samples, scramble=True, seed=seed + scramble_id)
        u = np.clip(pts[:, 0], np.finfo(float).tiny, 1.0 - np.finfo(float).eps)
        v = np.clip(pts[:, 1], np.finfo(float).tiny, 1.0 - np.finfo(float).eps)
        w = np.clip(pts[:, 2], np.finfo(float).tiny, 1.0 - np.finfo(float).eps)

        k = cutoff * u
        ell = cutoff * v
        theta = np.pi * w

        jac = cutoff * cutoff * np.pi
        measure = k ** (1.0 - 2.0 * float(eps)) * ell ** (1.0 - 2.0 * float(eps))
        measure *= np.sin(theta) ** (-2.0 * float(eps))

        d2 = D2_equal_mass_pperp0(point, ell, float(eta))
        d4 = (
            A_equal_mass_pperp0(point, ell, k)
            - 2.0 * ell * k * np.cos(theta)
            + 1j * float(eta) * float(gamma_xy(point))
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
