import mpmath as mp

from lsmethod.kinematics import Kinematics
from lsmethod.kernels import closed_form, method1_w_integral


def test_closed_and_w_integral_are_evaluable():
    mp.mp.dps = 50

    kin = Kinematics(
        x=mp.mpf("0.30"),
        y=mp.mpf("0.25"),
        p2=mp.mpf("-3.0"),
        ml=mp.mpf("0.70"),
        Ml=mp.mpf("1.10"),
        mk=mp.mpf("0.50"),
        Mk=mp.mpf("0.90"),
    )
    eps = mp.mpf("0.08")

    exact = closed_form(kin, eps)
    numeric = method1_w_integral(kin, eps)

    assert mp.isfinite(abs(exact))
    assert mp.isfinite(abs(numeric))
