import mpmath as mp

from .numeric import quad_split


# Closed-form and parameter-integral representations of the kernel.


def delta0(kin):
    """LaTeX: eq:Delta0-final."""
    # First effective scale from the remaining l-collinear pair.
    return kin.mu_l_x() - kin.x * (1 - kin.x) * kin.p2


def delta1(kin):
    """LaTeX: eq:Delta1-final."""
    # Second compact scale in the final kernel convention.
    lam = kin.lam
    return (
        -kin.mu_k_lam()
        + lam * (1 - lam) * kin.Ml**2
        + lam * (1 - lam * kin.x) * kin.p2
    )


def closed_form(kin, eps, eta=mp.mpf("1e-30")):
    """LaTeX: eq:closed-form-kernel."""
    # Outside support the theta functions set the kernel to zero.
    if not kin.in_support():
        return mp.mpc(0)

    # Numerical implementation of the i0 prescription for complex powers.
    d0 = delta0(kin) - 1j * eta
    d1 = delta1(kin) - 1j * eta

    return mp.gamma(eps) ** 2 / kin.X * d0 ** (-eps) * d1 ** (-eps)


def delta_w(kin):
    """LaTeX: eq:method1-Delta0-Deltaw-def."""
    # Slope of Delta_x(w) in the Method 1 one-dimensional representation.
    a = kin.kappa
    b = kin.lam * (1 - kin.lam * kin.x)

    return (
        ((1 + a) * kin.x * (1 - kin.x) + b) * kin.p2
        - kin.mu_k_lam()
        + a * kin.Ml**2
        - (1 + a) * kin.mu_l_x()
    )


def delta_x_w(kin, w):
    """LaTeX: eq:method1-Delta-x-def."""
    # Linear scale appearing in the w-integral.
    return delta0(kin) + w * delta_w(kin)


def method1_w_integral(kin, eps, eta=mp.mpf("1e-30"), parts=8):
    """LaTeX: eq:method1-w-integral-before-euler."""
    # Numerical version of the w-parameter representation before Euler reduction.
    if not kin.in_support():
        return mp.mpc(0)

    a = kin.kappa
    upper = 1 / (1 + a)

    def integrand(w):
        # U vanishes at the upper endpoint and generates one endpoint pole.
        U = 1 - (1 + a) * w
        scale = delta_x_w(kin, w)
        return w ** (-1 + eps) * U ** (-1 + eps) * (scale - 1j * eta) ** (-2 * eps)

    return mp.gamma(2 * eps) / kin.X * quad_split(integrand, mp.mpf("0"), upper, parts=parts)


def C_u(kin, u):
    """LaTeX: eq:C-u-def."""
    X = kin.X
    lam = kin.lam
    kappa = kin.kappa
    au = u - X

    mass_part = -2 * X * kappa * kin.mu_l_u(u) - 2 * au * kin.mu_k_lam()

    p2_part = 2 * (
        X * kappa * (au**2 + u * (1 - u))
        + au * kappa * (X**2 - u**2)
        + au * lam * (1 - lam * kin.x)
    ) * kin.p2

    return mass_part + p2_part


def find_u_root_in_interval(kin, X, root_guess):
    try:
        root = mp.findroot(lambda u: C_u(kin, u), root_guess)
    except Exception:
        return None

    root = mp.mpf(root)

    if mp.mpf("0") < root < X:
        return root

    return None


def can_split_around_root(root, X, delta):
    return mp.mpf("0") < root - delta < root + delta < X


def principal_value_integral(integrand, a, b, root, delta):
    left = mp.quad(integrand, [a, root - delta])
    right = mp.quad(integrand, [root + delta, b])
    return left + right


def method2_u_integral(
    kin,
    eps,
    eta=mp.mpf("1e-30"),
    parts=8,
    method="split",
    delta=mp.mpf("1e-12"),
    root_guess=(mp.mpf("0.50"), mp.mpf("0.55")),
    scale_sign=-1,
):
    """LaTeX: eq:method1-u-representation."""

    if not kin.in_support():
        return mp.mpc(0)

    X = kin.X
    kappa = kin.kappa

    def integrand(u):
        return (
            u ** (-1 + eps)
            * (X - u) ** (-1 + eps)
            * (scale_sign * C_u(kin, u) - 1j * eta) ** (-2 * eps)
        )

    def std_integral():
        return quad_split(integrand, mp.mpf("0"), X, parts=parts)

    if method == "std":
        integral = std_integral()

    elif method in {"split", "exclude", "pv"}:
        u_root = find_u_root_in_interval(kin, X, root_guess)

        if u_root is None or not can_split_around_root(u_root, X, delta):
            integral = std_integral()

        elif method == "split":
            integral = mp.quad(
                integrand,
                [
                    mp.mpf("0"),
                    u_root - delta,
                    u_root + delta,
                    X,
                ],
            )

        elif method == "exclude":
            integral = (
                mp.quad(integrand, [mp.mpf("0"), u_root - delta])
                + mp.quad(integrand, [u_root + delta, X])
            )

        elif method == "pv":
            integral = principal_value_integral(
                integrand=integrand,
                a=mp.mpf("0"),
                b=X,
                root=u_root,
                delta=delta,
            )

    else:
        raise ValueError(
            f"Unknown method={method!r}. Use 'std', 'split', 'exclude', or 'pv'."
        )

    return mp.gamma(2 * eps) * (4 * kappa) ** eps * integral
