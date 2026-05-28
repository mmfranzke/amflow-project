import mpmath as mp

import _path  # noqa: F401
from lsmethod.kinematics import Kinematics
from lsmethod.kernels import (
    C_u,
    closed_form,
    delta0,
    delta1,
    delta_w,
    find_u_root_in_interval,
    method1_w_integral,
    method2_u_integral,
)
from lsmethod.numeric import relative_error


def print_section(title):
    print()
    print(f"== {title} ==")


def print_kinematics(name, kin):
    """LaTeX: eq:support-theta, eq:lambda, eq:X-kappa-def."""
    print_section("Kinematics")
    print("point =", name)
    print("in_support =", kin.in_support())
    print("x =", kin.x)
    print("y =", kin.y)
    print("s = p2 =", kin.p2)
    print("ml2 =", kin.ml**2)
    print("Ml2 =", kin.Ml**2)
    print("mk2 =", kin.mk**2)
    print("Mk2 =", kin.Mk**2)
    print("X = 1 - x =", kin.X)
    print("lambda =", kin.lam)
    print("kappa =", kin.kappa)


def print_mass_scale_diagnostics(kin):
    """LaTeX: eq:mu-l-def, eq:method1-theta-z-def."""
    a = kin.kappa
    b = kin.lam * (1 - kin.lam * kin.x)
    d0 = delta0(kin)
    dw = delta_w(kin)
    d1 = delta1(kin)
    z = -dw / ((1 + a) * d0)

    print_section("Mass and Euler Scales")
    print("a = lambda (1 - lambda) =", a)
    print("b = lambda (1 - lambda x) =", b)
    print("muLx =", kin.mu_l_x())
    print("muKlambda =", kin.mu_k_lam())
    print("Delta0 =", d0)
    print("DeltaW =", dw)
    print("Delta1 =", d1)
    print("z = -DeltaW / ((1 + a) Delta0) =", z)


def print_w_scale_diagnostics(kin):
    """LaTeX: eq:method1-Delta-x-def, eq:method1-Delta0-Deltaw-def."""
    kappa = kin.kappa
    upper_w = 1 / (1 + kappa)

    d0 = delta0(kin)
    dw = delta_w(kin)
    endpoint = d0 + upper_w * dw

    print_section("Method 1 w Scale")
    print("Delta0 =", d0)
    print("Delta_w =", dw)
    print("w range = [0,", upper_w, "]")
    print("Delta_x(0) =", d0)
    print("Delta_x(wmax) =", endpoint)

    if dw == 0:
        print("Delta_x root = none, Delta_w is zero")
        return

    w_root = -d0 / dw
    print("Delta_x root =", w_root)
    print("root inside w range =", mp.mpf("0") < w_root < upper_w)


def print_C_u_endpoint_diagnostics(kin):
    """LaTeX: eq:C-endpoints."""
    X = kin.X
    kappa = kin.kappa
    u_root = find_u_root_in_interval(kin, X, (mp.mpf("0.40"), mp.mpf("0.45")))

    lhs_0 = C_u(kin, mp.mpf("0"))
    rhs_0 = -2 * X * delta1(kin)
    lhs_X = C_u(kin, X)
    rhs_X = -2 * X * kappa * delta0(kin)

    print_section("C_u Endpoints")
    print("C_u(0) =", lhs_0)
    print("-2 X Delta1 =", rhs_0)
    print("abs difference at u=0 =", abs(lhs_0 - rhs_0))
    print("C_u(X) =", lhs_X)
    print("-2 X kappa Delta0 =", rhs_X)
    print("abs difference at u=X =", abs(lhs_X - rhs_X))
    print("C_u root in (0, X) =", u_root)


def print_delta1_sign_diagnostics(kin):
    """LaTeX: eq:method1-Delta1-from-Deltaw, eq:Delta1-final."""
    S = (1 + kin.kappa) * delta0(kin) + delta_w(kin)
    d1 = delta1(kin)

    print_section("Delta1 Sign Convention")
    print("S = (1 + kappa) Delta0 + Delta_w =", S)
    print("Delta1 =", d1)
    print("S - Delta1 =", S - d1)
    print("S = Delta1 =", abs(S - d1) < mp.mpf("1e-40"))


def print_w_integral_vs_closed_form(kin, eps):
    """LaTeX: eq:method1-w-integral-before-euler, eq:closed-form-kernel."""
    w_integral = method1_w_integral(kin, eps)
    exact = closed_form(kin, eps)

    print_section("w Integral vs Closed Form")
    print("method1_w_integral =", w_integral)
    print("closed_form =", exact)
    print("rel_error =", relative_error(w_integral, exact))
    print("closed_form / 2 =", exact / 2)
    print("rel_error vs closed_form / 2 =", relative_error(w_integral, exact / 2))
    return w_integral


def print_w_integral_vs_u_integral(kin, eps, w_integral):
    """LaTeX: eq:method1-u-representation."""
    # This is a regression check between the two reduced representations.
    u_integral_plus = method2_u_integral(kin, eps, scale_sign=1)
    u_integral_minus = method2_u_integral(kin, eps, scale_sign=-1)

    print_section("w Integral vs u Integral")
    print("method1_w_integral =", w_integral)
    print("method2_u_integral, C_u scale =", u_integral_plus)
    print("rel_error, C_u scale =", relative_error(w_integral, u_integral_plus))
    print("method2_u_integral, -C_u scale =", u_integral_minus)
    print("rel_error, -C_u scale =", relative_error(w_integral, u_integral_minus))


def make_points():
    return {
        # No Delta_x root in the w interval.
        "safe": Kinematics(
            x=mp.mpf("0.30"),
            y=mp.mpf("0.25"),
            p2=mp.mpf("-1.00"),
            ml=mp.mpf("0.70"),
            Ml=mp.mpf("3.00"),
            mk=mp.mpf("0.50"),
            Mk=mp.mpf("0.90"),
        ),
        # Closer to the original point, but still mild.
        "mild": Kinematics(
            x=mp.mpf("0.30"),
            y=mp.mpf("0.25"),
            p2=mp.mpf("-0.10"),
            ml=mp.mpf("0.70"),
            Ml=mp.mpf("2.00"),
            mk=mp.mpf("0.50"),
            Mk=mp.mpf("0.90"),
        ),
        # Branch stress test: Delta_x crosses zero in the w interval.
        "original": Kinematics(
            x=mp.mpf("0.30"),
            y=mp.mpf("0.25"),
            p2=mp.mpf("-3.00"),
            ml=mp.mpf("0.70"),
            Ml=mp.mpf("1.10"),
            mk=mp.mpf("0.50"),
            Mk=mp.mpf("0.90"),
        ),
        # Rational branch-safe point used by the analytic Mathematica check.
        "branch_safe_rational": Kinematics(
            x=mp.mpf("0.30"),
            y=mp.mpf("0.25"),
            p2=mp.mpf("0.00"),
            ml=mp.mpf("0.70"),
            Ml=mp.mpf("2.00"),
            mk=mp.mpf("0.50"),
            Mk=mp.mpf("0.90"),
        ),
    }


def main():
    mp.mp.dps = 50

    eps = mp.mpf("0.08")
    points = make_points()

    # Try "branch_safe_rational" for the analytic-check point, or "original"
    # for the branch stress test.
    point_name = "branch_safe_rational"
    kin = points[point_name]

    print_kinematics(point_name, kin)
    print_mass_scale_diagnostics(kin)
    print_w_scale_diagnostics(kin)
    print_C_u_endpoint_diagnostics(kin)
    print_delta1_sign_diagnostics(kin)
    w_integral = print_w_integral_vs_closed_form(kin, eps)
    print_w_integral_vs_u_integral(kin, eps, w_integral)


if __name__ == "__main__":
    main()
