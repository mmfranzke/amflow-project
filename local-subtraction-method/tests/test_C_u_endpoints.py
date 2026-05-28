import mpmath as mp

from lsmethod.kinematics import Kinematics
from lsmethod.kernels import C_u, delta0, delta1


def test_C_u_endpoints():
    """LaTeX: eq:C-endpoints."""
    mp.mp.dps = 50

    kin = Kinematics(
        x=mp.mpf("0.30"),
        y=mp.mpf("0.25"),
        p2=mp.mpf("-3.00"),
        ml=mp.mpf("0.70"),
        Ml=mp.mpf("1.10"),
        mk=mp.mpf("0.50"),
        Mk=mp.mpf("0.90"),
    )

    X = kin.X
    kappa = kin.kappa

    assert abs(C_u(kin, mp.mpf("0")) + 2 * X * delta1(kin)) < mp.mpf("1e-40")
    assert abs(C_u(kin, X) + 2 * X * kappa * delta0(kin)) < mp.mpf("1e-40")
