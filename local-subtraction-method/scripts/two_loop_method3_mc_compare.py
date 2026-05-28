import argparse
import os
import re
import subprocess
import math
from pathlib import Path
from fractions import Fraction

import _path  # noqa: F401
from lsmethod.method3 import method3_prefactor, method3_stripped
from lsmethod.two_loop_points import POINTS, get_point


DEFAULT_EPS_LIST = ["1/5", "3/20", "1/10", "3/40", "1/20"]


def parse_mpf(text):
    return Fraction(str(text))


def parse_float(text):
    return float(Fraction(str(text)))


def delta0_value(point):
    return (1 - point.x) * point.ml2 + point.x * point.Ml2 - point.x * (1 - point.x) * point.p2


def delta1_value(point):
    lam = point.lam
    mu_k_lam = (1 - lam) * point.mk2 + lam * point.Mk2
    return -mu_k_lam + lam * (1 - lam) * point.Ml2 + lam * (1 - point.x) * point.p2


def closed_form_value(point, eps, eta):
    d0 = complex(float(delta0_value(point)), -float(eta))
    d1 = complex(float(delta1_value(point)), -float(eta))
    return math.gamma(float(eps)) ** 2 / float(1 - point.x) * d0 ** (-float(eps)) * d1 ** (-float(eps))


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--point", choices=sorted(POINTS), default="equal_mass_onshell_branch")
    parser.add_argument("--all-points", action="store_true")
    parser.add_argument("--eps", default="1/10")
    parser.add_argument("--eps-list", nargs="*", default=None)
    parser.add_argument("--eta", default="1e-6")
    parser.add_argument("--N", type=int, default=2**18)
    parser.add_argument("--scrambles", type=int, default=4)
    parser.add_argument("--rho", default=None)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--amflow", choices=["fresh", "skip"], default="fresh")
    parser.add_argument("--amflow-eps-order", default="4")
    parser.add_argument("--precision-goal", default="10")
    parser.add_argument("--dps", type=int, default=50)
    return parser.parse_args()


def complex_from_wolfram(text):
    text = text.strip()
    text = text.replace("*^", "e")
    text = text.replace("I", "j")
    text = text.replace(" ", "")
    if text.startswith("Complex[") and text.endswith("]"):
        inner = text[len("Complex[") : -1]
        real, imag = inner.split(",", 1)
        return complex(float(real), float(imag))
    return complex(eval(text, {"__builtins__": {}}, {"j": 1j}))


def run_amflow(point_name, eps, args):
    repo = Path(__file__).resolve().parents[2]
    run_sh = repo / "amflow-project" / "run.sh"
    env = os.environ.copy()
    env["TWO_LOOP_POINT"] = point_name
    env["TWO_LOOP_EPS_VALUE"] = str(float(eps))
    env["AMFLOW_EPS_ORDER"] = args.amflow_eps_order
    env["AMFLOW_PRECISION_GOAL"] = args.precision_goal

    proc = subprocess.run(
        [str(run_sh), "compare-twoloop-fixed-eps"],
        cwd=str(repo / "amflow-project"),
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if proc.returncode != 0:
        return None, proc.stdout

    match = re.search(r"AMFLOW_FIXED_EPS_VALUE=(.+)", proc.stdout)
    if not match:
        return None, proc.stdout
    try:
        return complex_from_wolfram(match.group(1)), proc.stdout
    except Exception:
        return None, proc.stdout


def fmt_complex(z):
    if z is None:
        return "not run"
    z = complex(z)
    return f"{z.real:.12e}{z.imag:+.12e}j"


def fmt_float(x):
    if x is None:
        return "n/a"
    try:
        return f"{float(x):.4e}"
    except Exception:
        return str(x)


def print_rule(widths):
    print("+-" + "-+-".join("-" * width for width in widths) + "-+")


def print_cells(cells, widths):
    print("| " + " | ".join(str(cell).ljust(width) for cell, width in zip(cells, widths)) + " |")


def print_small_table(headers, rows):
    widths = [
        max(len(str(header)), *(len(str(row[i])) for row in rows))
        for i, header in enumerate(headers)
    ]
    print_rule(widths)
    print_cells(headers, widths)
    print_rule(widths)
    for row in rows:
        print_cells(row, widths)
    print_rule(widths)


def print_delta_diagnostics(point, eta):
    d0 = delta0_value(point)
    d1 = delta1_value(point)
    print()
    print(f"== {point.name} ==")
    print(f"support 0 < x, 0 < y, x + y < pPlus: {point.in_support()}")
    print(f"pPlus={point.p_plus}, pMinus={point.p_minus}, pPerp2={point.p_perp2}, p2={point.p2}")
    print(f"x={point.x}, y={point.y}, lambda={point.lam}")
    print(f"Delta0 = {d0}")
    print(f"Delta1 = {d1}")
    print(f"Delta0 - i eta = {complex(float(d0), -float(eta))}")
    print(f"Delta1 - i eta = {complex(float(d1), -float(eta))}")
    print(f"(Delta0 - i eta)^(-eps) and (Delta1 - i eta)^(-eps) are printed in the table branch path.")


def compare_point(point_name, eps, eta, args):
    point = get_point(point_name)
    print_delta_diagnostics(point, eta)

    rho = None if args.rho is None else parse_mpf(args.rho)
    stripped, stripped_err = method3_stripped(
        point=point,
        eps=eps,
        eta=eta,
        n_samples=args.N,
        rho=rho,
        scrambles=args.scrambles,
        seed=args.seed,
    )
    pref = method3_prefactor(point, eps)
    normalized = complex(pref) * stripped
    normalized_err = abs(complex(pref)) * stripped_err

    method12 = closed_form_value(point, eps, Fraction(1, 10**30))
    d0_pow = complex(float(delta0_value(point)), -1e-30) ** (-float(eps))
    d1_pow = complex(float(delta1_value(point)), -1e-30) ** (-float(eps))

    amflow = None
    if args.amflow == "fresh":
        print("Running AMFlow fixed-eps original integral for this point...")
        amflow, output = run_amflow(point_name, eps, args)
        if amflow is None:
            print("AMFlow value unavailable. Last AMFlow output follows:")
            print(output[-4000:])

    return {
        "point_name": point_name,
        "eps": eps,
        "eta": eta,
        "N": args.N,
        "method3_stripped": stripped,
        "method3_stripped_err": stripped_err,
        "method3_normalized": normalized,
        "method3_normalized_err": normalized_err,
        "method12_closed": method12,
        "AMFlow_original": amflow,
        "diff_method3_vs_method12": normalized - method12,
        "diff_method3_vs_AMFlow": None if amflow is None else normalized - amflow,
        "diff_method12_vs_AMFlow": None if amflow is None else method12 - amflow,
        "Delta0_power": d0_pow,
        "Delta1_power": d1_pow,
    }


def print_table(rows):
    print()
    print("== Fixed-eps Comparison Summary ==")
    print_small_table(
        ["point", "eps", "eta", "N"],
        [
            [row["point_name"], str(row["eps"]), str(row["eta"]), str(row["N"])]
            for row in rows
        ],
    )

    for row in rows:
        print()
        print(f"== {row['point_name']} ==")
        print_small_table(
            ["quantity", "value", "error"],
            [
                [
                    "method3_stripped",
                    fmt_complex(row["method3_stripped"]),
                    fmt_float(row["method3_stripped_err"]),
                ],
                [
                    "method3_normalized",
                    fmt_complex(row["method3_normalized"]),
                    fmt_float(row["method3_normalized_err"]),
                ],
                ["method12_closed", fmt_complex(row["method12_closed"]), ""],
                ["AMFlow_original", fmt_complex(row["AMFlow_original"]), ""],
            ],
        )
        print_small_table(
            ["difference", "value"],
            [
                ["method3_normalized - method12_closed", fmt_complex(row["diff_method3_vs_method12"])],
                ["method3_normalized - AMFlow_original", fmt_complex(row["diff_method3_vs_AMFlow"])],
                ["method12_closed - AMFlow_original", fmt_complex(row["diff_method12_vs_AMFlow"])],
            ],
        )

    print()
    print("== Branch Powers ==")
    for row in rows:
        print(
            f"{row['point_name']}: "
            f"(Delta0-i0)^(-eps)={fmt_complex(row['Delta0_power'])}, "
            f"(Delta1-i0)^(-eps)={fmt_complex(row['Delta1_power'])}"
        )


def main():
    args = parse_args()

    if args.eps_list is not None:
        eps_values = [parse_mpf(value) for value in (args.eps_list or DEFAULT_EPS_LIST)]
    else:
        eps_values = [parse_mpf(args.eps)]

    point_names = sorted(POINTS) if args.all_points else [args.point]
    eta = parse_mpf(args.eta)

    rows = []
    for eps in eps_values:
        for point_name in point_names:
            rows.append(compare_point(point_name, eps, eta, args))

    print_table(rows)


if __name__ == "__main__":
    main()
