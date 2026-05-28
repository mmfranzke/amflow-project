import mpmath as mp

from lsmethod.kinematics import Kinematics
from lsmethod.kernels import C_u


def test_C_u_matches_D_u_definition_at_interior_point():
    """LaTeX: eq:D-u-def, eq:D-u-linear."""
    mp.mp.dps = 50

    kin = Kinematics(
        x=mp.mpf("0.30"),
        y=mp.mpf("0.25"),
        p2=mp.mpf("-1.00"),
        ml=mp.mpf("0.70"),
        Ml=mp.mpf("3.00"),
        mk=mp.mpf("0.50"),
        Mk=mp.mpf("0.90"),
    )

    u = mp.mpf("0.37")
    X = kin.X
    lam = kin.lam
    kappa = kin.kappa
    au = u - X
    beta_u = 2 * X * kappa
    Lperp2 = mp.mpf("0")

    Lambda_u = (
        au**2 * kin.p2
        + Lperp2
        - kin.mu_l_u(u)
        + u * (1 - u) * kin.p2
    )

    rho_u = (
        kin.mu_k_lam()
        - kappa * (((kin.x - 1 + u) * (kin.x - 1 - u)) * kin.p2 + Lperp2)
        - lam * (1 - lam * kin.x) * kin.p2
    )

    D_u = beta_u * Lambda_u - 2 * au * rho_u

    assert abs(C_u(kin, u) - D_u) < mp.mpf("1e-40")
