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


def _as_mpf(value):
    if hasattr(value, "numerator") and hasattr(value, "denominator"):
        return mp.mpf(value.numerator) / value.denominator
    return mp.mpf(value)


def _squared_mass(obj, square_name, mass_name):
    if hasattr(obj, square_name):
        return _as_mpf(getattr(obj, square_name))
    return _as_mpf(getattr(obj, mass_name)) ** 2


def _l_residue_supported(obj):
    p_plus = _as_mpf(getattr(obj, "p_plus", 1))
    p_perp2 = _as_mpf(getattr(obj, "p_perp2", 0))
    return p_plus == 1 and p_perp2 == 0


def l_residue_delta0(point):
    """Delta0 for the direct l^- residue diagnostic."""
    if not _l_residue_supported(point):
        raise ValueError("l-residue closed form currently implemented for pPlus=1, pPerp2=0.")
    x = _as_mpf(point.x)
    p2 = _as_mpf(point.p2)
    ml2 = _squared_mass(point, "ml2", "ml")
    Ml2 = _squared_mass(point, "Ml2", "Ml")
    return (1 - x) * ml2 + x * Ml2 - x * (1 - x) * p2


def l_residue_A1_short_without_extra_term(point):
    """Short diagnostic A1 omitting the current-PDF off-shell term."""
    if not _l_residue_supported(point):
        raise ValueError("l-residue closed form currently implemented for pPlus=1, pPerp2=0.")
    x = _as_mpf(point.x)
    y = _as_mpf(point.y)
    X = 1 - x
    lam = y / X
    kappa = lam * (1 - lam)
    p2 = _as_mpf(point.p2)
    ml2 = _squared_mass(point, "ml2", "ml")
    mk2 = _squared_mass(point, "mk2", "mk")
    Mk2 = _squared_mass(point, "Mk2", "Mk")
    mu_k = (1 - lam) * mk2 + lam * Mk2
    return mu_k - kappa * X * p2 + kappa * X / x * ml2


def l_residue_extra_offshell_term(point):
    """Extra current-PDF term -lambda*(1-lambda*x)*p2 in C(l)."""
    if not _l_residue_supported(point):
        raise ValueError("l-residue closed form currently implemented for pPlus=1, pPerp2=0.")
    x = _as_mpf(point.x)
    y = _as_mpf(point.y)
    lam = y / (1 - x)
    p2 = _as_mpf(point.p2)
    return -lam * (1 - lam * x) * p2


def l_residue_A1(point):
    """Full current-PDF A1 in C(l^-_1,t)=A1+B1 t."""
    return l_residue_A1_short_without_extra_term(point) + l_residue_extra_offshell_term(point)


def l_residue_B1(point):
    """B1 in C(l^-_1,t)=A1+B1 t for the direct l^- residue diagnostic."""
    if not _l_residue_supported(point):
        raise ValueError("l-residue closed form currently implemented for pPlus=1, pPerp2=0.")
    x = _as_mpf(point.x)
    y = _as_mpf(point.y)
    lam = y / (1 - x)
    return lam * (1 - lam) / x


def l_residue_z(point, short_C_diagnostic=False):
    """Hypergeometric argument z = 1 - B1 Delta0/A1."""
    a1 = l_residue_A1_short_without_extra_term(point) if short_C_diagnostic else l_residue_A1(point)
    return 1 - l_residue_B1(point) * l_residue_delta0(point) / a1


def l_residue_closed_ingredients(point, short_C_diagnostic=False):
    """Ingredients for the direct l^- residue closed-form diagnostic."""
    if not _l_residue_supported(point):
        raise ValueError("l-residue closed form currently implemented for pPlus=1, pPerp2=0.")
    x = _as_mpf(point.x)
    y = _as_mpf(point.y)
    X = 1 - x
    lam = y / X
    kappa = lam * (1 - lam)
    mk2 = _squared_mass(point, "mk2", "mk")
    Mk2 = _squared_mass(point, "Mk2", "Mk")
    a1_short = l_residue_A1_short_without_extra_term(point)
    extra = l_residue_extra_offshell_term(point)
    a1 = a1_short if short_C_diagnostic else a1_short + extra
    return {
        "X": X,
        "lambda": lam,
        "kappa": kappa,
        "mu_k_lambda": (1 - lam) * mk2 + lam * Mk2,
        "Delta0": l_residue_delta0(point),
        "A1": a1,
        "A1_full": a1_short + extra,
        "A1_short_without_extra_term": a1_short,
        "extra_offshell_term": extra,
        "C_full_includes_extra_offshell_term": not short_C_diagnostic,
        "B1": l_residue_B1(point),
        "z": 1 - l_residue_B1(point) * l_residue_delta0(point) / a1,
    }


def method1_l_residue_closed(point, eps, eta=mp.mpf("1e-30"), short_C_diagnostic=False):
    """Candidate closed form for the direct l^- residue representation.

    Diagnostic only: this is not the default Method-1/2 compact product.
    """
    if not _l_residue_supported(point):
        raise ValueError("l-residue closed form currently implemented for pPlus=1, pPerp2=0.")
    eps = _as_mpf(eps)
    eta = _as_mpf(eta)
    ingredients = l_residue_closed_ingredients(point, short_C_diagnostic=short_C_diagnostic)
    X = ingredients["X"]
    d0 = ingredients["Delta0"] - 1j * eta
    a1 = ingredients["A1"] - 1j * eta
    z = ingredients["z"]
    hyper_arg = z
    branch_conjugate = False
    if ingredients["A1"] < 0 and z > 1:
        # The full current-PDF off-shell scale can cross zero.  The direct
        # t-integral with C(t)-i0 corresponds to the opposite side of the
        # hypergeometric branch cut from mpmath's default continuation.
        hyper_arg = z + 1j * max(eta, mp.mpf("1e-30"))
        branch_conjugate = True
    value = (
        mp.gamma(2 * eps)
        / (X * eps)
        * (d0 * a1) ** (-eps)
        * mp.hyper([eps, 1 - eps], [1 + eps], hyper_arg)
    )
    if branch_conjugate:
        return mp.conj(value)
    return (
        value
    )


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
