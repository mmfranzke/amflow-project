import mpmath as mp

import _path  # noqa: F401
from lsmethod.kinematics import Kinematics
from lsmethod.kernels import closed_form, method1_w_integral
from lsmethod.numeric import relative_error


# Compare the closed form with the Method 1 w-integral at one Euclidean point.
def main():
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

    print("closed_form =", exact)
    print("method1_w   =", numeric)
    print("rel_error   =", relative_error(numeric, exact))


if __name__ == "__main__":
    main()
