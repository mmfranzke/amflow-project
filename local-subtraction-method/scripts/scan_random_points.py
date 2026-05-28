import argparse
import csv
import random
from pathlib import Path

import mpmath as mp

import _path  # noqa: F401
from lsmethod.kinematics import Kinematics
from lsmethod.kernels import closed_form, method1_w_integral
from lsmethod.numeric import relative_error


# Generate one random point inside the support region.
def random_kinematics():
    x = mp.mpf(str(random.uniform(0.05, 0.80)))
    y = mp.mpf(str(random.uniform(0.02, float(0.95 * (1 - x)))))
    return Kinematics(
        x=x,
        y=y,
        p2=mp.mpf("-3.0"),
        ml=mp.mpf("0.70"),
        Ml=mp.mpf("1.10"),
        mk=mp.mpf("0.50"),
        Mk=mp.mpf("0.90"),
    )


# Run a small random scan and write the results to CSV.
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--npoints", type=int, default=20)
    parser.add_argument("--eps", type=str, default="0.08")
    parser.add_argument("--out", type=str, default="results/numeric_checks/scan.csv")
    args = parser.parse_args()

    mp.mp.dps = 50
    eps = mp.mpf(args.eps)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    with out.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["x", "y", "closed_real", "w_real", "relative_error"])

        for _ in range(args.npoints):
            kin = random_kinematics()
            exact = closed_form(kin, eps)
            numeric = method1_w_integral(kin, eps)
            err = relative_error(numeric, exact)
            writer.writerow([kin.x, kin.y, mp.re(exact), mp.re(numeric), err])

    print(f"wrote {out}")


if __name__ == "__main__":
    main()
