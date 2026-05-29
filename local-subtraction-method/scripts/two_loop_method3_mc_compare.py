import argparse
import contextlib
import cmath
import os
import re
import subprocess
import math
import sys
from pathlib import Path
from fractions import Fraction

import _path  # noqa: F401
import numpy as np
from lsmethod.method3 import (
    gamma_xy,
    mapped_integrand_probe,
    method3_prefactor,
    method3_prefactor_components,
    method3_stripped,
    method3_stripped_cutoff,
)
from lsmethod.two_loop_points import POINTS, get_point

DEFAULT_EPS_LIST = ["1/5", "3/20", "1/10", "3/40", "1/20"]
EULER_GAMMA = 0.577215664901532860606512090082402431
ZETA3 = 1.20205690315959428539973816151144999


def parse_mpf(text):
    return Fraction(str(text))


def parse_eps_list(values):
    if values is None:
        return [parse_mpf(value) for value in DEFAULT_EPS_LIST]

    items = []
    for value in values:
        items.extend(part.strip() for part in str(value).split(",") if part.strip())
    return [parse_mpf(value) for value in items]


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


def delta_log(delta, branch):
    if delta < 0:
        phase = -math.pi if branch == "-i0" else math.pi
        return math.log(abs(float(delta))) + 1j * phase
    return math.log(float(delta))


def closed_form_value_branch(point, eps, branch="-i0"):
    d0_log = delta_log(delta0_value(point), branch)
    d1_log = delta_log(delta1_value(point), branch)
    return math.gamma(float(eps)) ** 2 / float(1 - point.x) * cmath.exp(-float(eps) * (d0_log + d1_log))


def method12_laurent_coefficients(point, min_power=-2, max_power=2, branch="-i0"):
    if min_power < -2:
        raise ValueError("method12 analytic coefficients are implemented from eps^-2 upward.")
    if max_power > 2:
        raise ValueError("method12 analytic coefficients are implemented through eps^2.")

    s_log = delta_log(delta0_value(point), branch) + delta_log(delta1_value(point), branch)
    q1 = -(2.0 * EULER_GAMMA + s_log)
    q2 = math.pi**2 / 6.0
    q3 = -2.0 * ZETA3 / 3.0
    q4 = math.pi**4 / 180.0

    exp_coeffs = {
        0: 1.0 + 0j,
        1: q1,
        2: q2 + q1**2 / 2.0,
        3: q3 + q1 * q2 + q1**3 / 6.0,
        4: q4 + q1 * q3 + q2**2 / 2.0 + q1**2 * q2 / 2.0 + q1**4 / 24.0,
    }
    xval = float(1 - point.x)
    return {power: exp_coeffs[power + 2] / xval for power in range(min_power, max_power + 1)}


def repo_root():
    return Path(__file__).resolve().parents[2]


def safe_name(value):
    return str(value).replace("/", "_over_").replace(".", "p").replace("-", "m").replace("+", "p")


def parse_args():
    examples = """examples:
  ./run.sh two-loop-method3-mc-compare --all-points --eps 0.1 --amflow skip
  ./run.sh two-loop-method3-mc-compare --all-points --compare-coefficients --amflow skip
  ./run.sh two-loop-method3-mc-compare --point equal_mass_offshell_positive --compare-coefficients --eps-list 0.2,0.15,0.1,0.075,0.05 --amflow skip
"""
    parser = argparse.ArgumentParser(description=__doc__, epilog=examples, formatter_class=argparse.RawDescriptionHelpFormatter)
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
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--debug-include-2p20", action="store_true")
    parser.add_argument("--compare-coefficients", "--eps-fit", action="store_true")
    parser.add_argument("--fit-min-power", type=int, default=-2)
    parser.add_argument("--fit-max-power", type=int, default=2)
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
    repo = repo_root()
    run_sh = repo / "amflow-project" / "run.sh"
    log_dir = repo / "amflow-project" / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"two-loop-method3-mc-compare_amflow_{point_name}_eps_{safe_name(eps)}.log"
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
    log_path.write_text(proc.stdout, encoding="utf-8")
    print(f"AMFlow fixed-eps log written to: {log_path}")
    if proc.returncode != 0:
        return None, proc.stdout

    match = re.search(r"AMFLOW_FIXED_EPS_VALUE=(.+)", proc.stdout)
    if not match:
        return None, proc.stdout
    try:
        return complex_from_wolfram(match.group(1)), proc.stdout
    except Exception:
        return None, proc.stdout


def run_amflow_eps_values(point_name, eps_values, args):
    repo = repo_root()
    run_sh = repo / "amflow-project" / "run.sh"
    log_dir = repo / "amflow-project" / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    eps_label = "_".join(safe_name(eps) for eps in eps_values)
    log_path = log_dir / f"two-loop-method3-mc-compare_amflow_{point_name}_eps_list_{eps_label}.log"
    env = os.environ.copy()
    env["TWO_LOOP_POINT"] = point_name
    env["TWO_LOOP_EPS_VALUE"] = str(float(eps_values[0]))
    env["TWO_LOOP_EPS_LIST"] = ",".join(str(eps) for eps in eps_values)
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
    log_path.write_text(proc.stdout, encoding="utf-8")
    print(f"AMFlow eps-list log written to: {log_path}")
    if proc.returncode != 0:
        return [None for _ in eps_values], proc.stdout

    values = [None for _ in eps_values]
    for match in re.finditer(r"AMFLOW_EPS_LIST_VALUE_INDEX=(\d+) VALUE=(.+)", proc.stdout):
        index = int(match.group(1)) - 1
        if 0 <= index < len(values):
            try:
                values[index] = complex_from_wolfram(match.group(2))
            except Exception:
                values[index] = None

    return values, proc.stdout


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


def ratio_or_none(numerator, denominator):
    if numerator is None or denominator is None or denominator == 0:
        return None
    return numerator / denominator


def method3_supported(point):
    return point.p_perp2 == 0 and point.ml2 == point.Ml2 == point.mk2 == point.Mk2


def method3_unsupported_reason(point):
    if point.p_perp2 != 0:
        return "method3 currently implements p_perp2 = 0"
    if point.ml2 != point.Ml2 or point.ml2 != point.mk2 or point.ml2 != point.Mk2:
        return (
            "method3 currently implements the equal-mass PDF specialization; "
            f"got ml2={point.ml2}, Ml2={point.Ml2}, mk2={point.mk2}, Mk2={point.Mk2}"
        )
    return None


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


class Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for stream in self.streams:
            stream.write(data)

    def flush(self):
        for stream in self.streams:
            stream.flush()


def compare_point(point_name, eps, eta, args):
    point = get_point(point_name)
    print_delta_diagnostics(point, eta)

    rho = None if args.rho is None else parse_mpf(args.rho)
    unsupported_reason = method3_unsupported_reason(point)
    if unsupported_reason is None:
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
    else:
        print(f"Method 3 skipped: {unsupported_reason}")
        stripped = None
        stripped_err = None
        pref = None
        normalized = None
        normalized_err = None

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
        "diff_method3_vs_method12": None if normalized is None else normalized - method12,
        "diff_method3_vs_AMFlow": None if amflow is None or normalized is None else normalized - amflow,
        "diff_method12_vs_AMFlow": None if amflow is None else method12 - amflow,
        "ratio_method12_over_method3_stripped": ratio_or_none(method12, stripped),
        "ratio_method12_over_method3_normalized": ratio_or_none(method12, normalized),
        "Delta0_power": d0_pow,
        "Delta1_power": d1_pow,
    }


def print_method12_debug(point, eps):
    eta_branch = Fraction(1, 10**30)
    x = float(point.x)
    xval = 1.0 - x
    d0 = delta0_value(point)
    d1 = delta1_value(point)
    gamma_eps_sq = math.gamma(float(eps)) ** 2
    d0_minus = complex(float(d0), -float(eta_branch)) ** (-float(eps))
    d1_minus = complex(float(d1), -float(eta_branch)) ** (-float(eps))
    d0_plus = complex(float(d0), float(eta_branch)) ** (-float(eps))
    d1_plus = complex(float(d1), float(eta_branch)) ** (-float(eps))
    method12 = gamma_eps_sq / xval * d0_minus * d1_minus
    method12_plus = gamma_eps_sq / xval * d0_plus * d1_plus

    print()
    print("== Debug: Method 1/2 Closed Form Components ==")
    rows = [
        ["x", point.x],
        ["y", point.y],
        ["lambda = y/(1-x)", point.lam],
        ["pPlus", point.p_plus],
        ["pMinus", point.p_minus],
        ["p2", point.p2],
        ["ml2", point.ml2],
        ["Ml2", point.Ml2],
        ["mk2", point.mk2],
        ["Mk2", point.Mk2],
        ["Delta0", d0],
        ["Delta1", d1],
        ["Xval used", xval],
        ["1-x", 1.0 - x],
        ["Xval - (1-x)", xval - (1.0 - x)],
        ["Gamma[eps]^2", f"{gamma_eps_sq:.16e}"],
        ["Delta0Power = (Delta0 - i0)^(-eps)", fmt_complex(d0_minus)],
        ["Delta1Power = (Delta1 - i0)^(-eps)", fmt_complex(d1_minus)],
        ["Delta0Power opposite = (Delta0 + i0)^(-eps)", fmt_complex(d0_plus)],
        ["Delta1Power opposite = (Delta1 + i0)^(-eps)", fmt_complex(d1_plus)],
        ["method12_closed (-i0)", fmt_complex(method12)],
        ["method12_closed opposite (+i0)", fmt_complex(method12_plus)],
    ]
    print_small_table(["component", "value"], rows)
    if abs(xval - (1.0 - x)) > 0:
        print("LOUD CHECK FAILED: Xval is not 1-x.")


def print_method3_formula_debug(point):
    print()
    print("== Debug: Method 3 Denominator Formula Check ==")
    print("D2(l) = (y - pPlus) * ((l^2 + Ml2)/y - pMinus) - (l^2 + Ml2) + i eta * pPlus/y")
    print("D4 = A(l,k;x,y) - 2 l k cos(theta) + i eta Gamma_xy")
    print("Gamma_xy = (x + y) (pPlus - x - y)/(x y)")
    print("A = ((x - pPlus)/y) l^2 + ((y - pPlus)/x) k^2 + massTerm - (x + y - pPlus) pMinus")
    print("equal-mass massTerm = M^2 * (((x + y)(x + y - pPlus) - x y)/(x y))")
    rows = [
        ["ml2", point.ml2],
        ["Ml2 used in D2", point.Ml2],
        ["mk2", point.mk2],
        ["Mk2", point.Mk2],
        ["M2 used in equal-mass A", point.M2 if method3_supported(point) else "not available"],
        ["all masses equal", point.ml2 == point.Ml2 == point.mk2 == point.Mk2],
        ["x > 0", point.x > 0],
        ["y > 0", point.y > 0],
        ["x + y < pPlus", point.x + point.y < point.p_plus],
        ["Gamma_xy", gamma_xy(point)],
        ["Gamma_xy > 0", gamma_xy(point) > 0],
    ]
    print_small_table(["check", "value"], rows)


def print_method3_prefactor_debug(point, eps, stripped):
    components = method3_prefactor_components(point, eps)
    prefactor = method3_prefactor(point, eps)
    xy = float(point.x * point.y)
    one_minus_x = float(1 - point.x)
    variants = {
        "prefactor": prefactor,
        "prefactor_no_contour_sign = -prefactor": -prefactor,
        "prefactor_without_1_over_xy": prefactor * xy,
        "prefactor_times_xy": prefactor * xy,
        "prefactor_times_1_minus_x": prefactor * one_minus_x,
        "prefactor_divided_by_1_minus_x": prefactor / one_minus_x,
    }

    print()
    print("== Debug: Method 3 Normalization Components ==")
    rows = [
        ["d = 4 - 2 eps", f"{components['d']:.16e}"],
        ["n = d - 2", f"{components['n']:.16e}"],
        ["Omega[n-1]", f"{components['omega_n_minus_1']:.16e}"],
        ["Omega[n-2]", f"{components['omega_n_minus_2']:.16e}"],
        ["(2 Pi I)^2", fmt_complex(components["contour_factor"])],
        ["(I Pi^(d/2))^2", fmt_complex(components["loop_normalization"])],
        ["1/(x y)", f"{components['one_over_xy']:.16e}"],
        ["prefactor", fmt_complex(prefactor)],
    ]
    print_small_table(["component", "value"], rows)

    print()
    print("== Debug: Method 3 Prefactor Variants Applied To Current Stripped Integral ==")
    print_small_table(
        ["variant", "factor", "stripped * factor"],
        [[name, fmt_complex(factor), fmt_complex(stripped * factor)] for name, factor in variants.items()],
    )


def print_integrand_probes(point, eps, eta, rho):
    print()
    print("== Debug: Deterministic Infinite-Map Integrand Probes ==")
    probes = [
        (0.1, 0.1, 0.5),
        (0.5, 0.5, 0.5),
        (0.9, 0.9, 0.5),
        (0.99, 0.99, 0.5),
        (0.999, 0.999, 0.5),
    ]
    for u, v, w in probes:
        probe = mapped_integrand_probe(point, eps, eta, u, v, w, rho=rho)
        print()
        print(f"probe u={u}, v={v}, w={w}")
        print_small_table(
            ["quantity", "value"],
            [
                ["k", f"{probe['k']:.16e}"],
                ["l", f"{probe['ell']:.16e}"],
                ["theta", f"{probe['theta']:.16e}"],
                ["jacobian", f"{probe['jacobian']:.16e}"],
                ["D2", fmt_complex(probe["D2"])],
                ["A", f"{probe['A']:.16e}"],
                ["Gamma_xy", probe["Gamma_xy"]],
                ["D4", fmt_complex(probe["D4"])],
                ["measure", f"{probe['measure']:.16e}"],
                ["integrand_without_jacobian", fmt_complex(probe["integrand_without_jacobian"])],
                ["integrand_with_jacobian", fmt_complex(probe["integrand_with_jacobian"])],
                ["real sign D2", probe["D2_real_sign"]],
                ["real sign D4", probe["D4_real_sign"]],
            ],
        )


def print_cutoff_scan(point, eps, eta, n_samples, scrambles, seed):
    prefactor = method3_prefactor(point, eps)
    print()
    print("== Debug: Finite Radial Cutoff Scan ==")
    rows = []
    for cutoff in [1, 2, 5, 10, 20, 50, 100]:
        stripped, error = method3_stripped_cutoff(
            point=point,
            eps=eps,
            eta=eta,
            n_samples=n_samples,
            cutoff=cutoff,
            scrambles=scrambles,
            seed=seed,
        )
        rows.append([
            cutoff,
            fmt_complex(stripped),
            fmt_float(error),
            fmt_complex(prefactor * stripped),
        ])
    print_small_table(["Kmax=Lmax", "stripped_cutoff_value", "error", "normalized_cutoff_value"], rows)


def print_rho_scan(point, eps, eta, n_samples, scrambles, seed):
    prefactor = method3_prefactor(point, eps)
    print()
    print("== Debug: Infinite-Map Rho Scan ==")
    rows = []
    for rho in [0.25, 0.5, 1, 2, 5, 10, 20]:
        stripped, error = method3_stripped(
            point=point,
            eps=eps,
            eta=eta,
            n_samples=n_samples,
            rho=Fraction(str(rho)),
            scrambles=scrambles,
            seed=seed,
        )
        rows.append([rho, fmt_complex(stripped), fmt_float(error), fmt_complex(prefactor * stripped)])
    print_small_table(["rho", "method3_stripped", "estimated_error", "method3_normalized"], rows)


def print_n_scan(point, eps, eta, rho, scrambles, seed, include_2p20):
    prefactor = method3_prefactor(point, eps)
    ns = [2**14, 2**16, 2**18]
    if include_2p20:
        ns.append(2**20)

    print()
    print("== Debug: N Scan ==")
    rows = []
    for n_samples in ns:
        stripped, error = method3_stripped(
            point=point,
            eps=eps,
            eta=eta,
            n_samples=n_samples,
            rho=rho,
            scrambles=scrambles,
            seed=seed,
        )
        rows.append([n_samples, fmt_complex(stripped), fmt_float(error), fmt_complex(prefactor * stripped)])
    print_small_table(["N", "method3_stripped", "estimated_error", "method3_normalized"], rows)


def print_eta_scan(point, eps, n_samples, rho, scrambles, seed):
    prefactor = method3_prefactor(point, eps)
    print()
    print("== Debug: Eta Scan ==")
    rows = []
    for eta_text in ["1e-2", "1e-3", "1e-4", "1e-5", "1e-6"]:
        eta = Fraction(eta_text)
        stripped, error = method3_stripped(
            point=point,
            eps=eps,
            eta=eta,
            n_samples=n_samples,
            rho=rho,
            scrambles=scrambles,
            seed=seed,
        )
        rows.append([eta_text, fmt_complex(stripped), fmt_float(error), fmt_complex(prefactor * stripped)])
    print_small_table(["eta", "method3_stripped", "estimated_error", "method3_normalized"], rows)


def print_comparison_levels(row):
    point = get_point(row["point_name"])
    prefactor = method3_prefactor(point, row["eps"]) if method3_supported(point) else None
    print()
    print("== Debug: Comparison Levels And Ratios ==")
    print_small_table(
        ["quantity", "value"],
        [
            ["method3_stripped", fmt_complex(row["method3_stripped"])],
            ["method3_prefactor", fmt_complex(prefactor)],
            ["method3_normalized", fmt_complex(row["method3_normalized"])],
            ["method12_closed", fmt_complex(row["method12_closed"])],
            ["method12_closed / method3_stripped", fmt_complex(row["ratio_method12_over_method3_stripped"])],
            ["method12_closed / method3_normalized", fmt_complex(row["ratio_method12_over_method3_normalized"])],
        ],
    )


def print_debug_report(point, eps, eta, args, row):
    rho = None if args.rho is None else parse_mpf(args.rho)
    print()
    print(f"######## Method 3 Debug Report: {point.name}, eps={eps}, eta={eta} ########")
    print_method12_debug(point, eps)
    print_method3_formula_debug(point)
    print_comparison_levels(row)
    unsupported_reason = method3_unsupported_reason(point)
    if unsupported_reason is None:
        print_method3_prefactor_debug(point, eps, row["method3_stripped"])
        print_integrand_probes(point, eps, eta, rho)
        print_cutoff_scan(point, eps, eta, args.N, args.scrambles, args.seed)
        print_rho_scan(point, eps, eta, args.N, args.scrambles, args.seed)
        print_n_scan(point, eps, eta, rho, args.scrambles, args.seed, args.debug_include_2p20)
        print_eta_scan(point, eps, args.N, rho, args.scrambles, args.seed)
    else:
        print(f"Method 3 numerical diagnostics skipped: {unsupported_reason}")
    print("######## End Method 3 Debug Report ########")


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
                ["method12_closed / method3_stripped", fmt_complex(row["ratio_method12_over_method3_stripped"])],
                ["method12_closed / method3_normalized", fmt_complex(row["ratio_method12_over_method3_normalized"])],
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


def fit_laurent_coefficients(eps_values, values, min_power, max_power):
    powers = list(range(min_power, max_power + 1))
    matrix = np.array(
        [[float(eps) ** power for power in powers] for eps in eps_values],
        dtype=np.complex128,
    )
    vector = np.array(values, dtype=np.complex128)

    if len(eps_values) == len(powers):
        coeffs = np.linalg.solve(matrix, vector)
    else:
        coeffs, _, _, _ = np.linalg.lstsq(matrix, vector, rcond=None)

    return {power: coeffs[i] for i, power in enumerate(powers)}


def evaluate_laurent(coefficients, eps):
    eps_f = float(eps)
    return sum(value * eps_f**power for power, value in coefficients.items())


def coefficient_rows(coefficients, min_power, max_power):
    return [[f"c[{power}]", fmt_complex(coefficients[power])] for power in range(min_power, max_power + 1)]


def print_method12_coefficient_debug(point, min_power, max_power):
    d0 = delta0_value(point)
    d1 = delta1_value(point)
    print()
    print("== Method 1/2 Analytic Laurent Coefficients ==")
    branch = "-i0"
    rows = [
        ["point_name", point.name],
        ["Delta0", d0],
        ["Delta1", d1],
        ["branch", branch],
    ]
    if d0 < 0 or d1 < 0:
        rows.append(["negative-Delta phase", "(-a - i0)^(-eps) = a^(-eps) Exp[+ I Pi eps]"])
        rows.append(["fixed-eps expectation", "method12_closed is complex on the default -i0 branch"])
    print_small_table(["quantity", "value"], rows)

    default_coeffs = method12_laurent_coefficients(point, min_power, max_power, branch="-i0")
    print_small_table(["coefficient", "method12_analytic (-i0)"], coefficient_rows(default_coeffs, min_power, max_power))

    if d0 < 0 or d1 < 0:
        plus_coeffs = method12_laurent_coefficients(point, min_power, max_power, branch="+i0")
        print()
        print("== Method 1/2 Alternative Branch Diagnostics ==")
        print_small_table(["coefficient", "method12_analytic (+i0 diagnostic)"], coefficient_rows(plus_coeffs, min_power, max_power))

    return default_coeffs


def print_sample_table(title, eps_values, values, errors=None):
    print()
    print(title)
    rows = []
    for i, eps in enumerate(eps_values):
        error = "" if errors is None else fmt_float(errors[i])
        rows.append([str(eps), fmt_complex(values[i]), error])
    print_small_table(["eps", "value", "error"], rows)


def print_fit_residuals(title, eps_values, values, coefficients):
    print()
    print(title)
    rows = []
    for eps, value in zip(eps_values, values):
        fit_value = evaluate_laurent(coefficients, eps)
        rows.append([str(eps), fmt_complex(value), fmt_complex(fit_value), fmt_complex(value - fit_value)])
    print_small_table(["eps", "value", "fit_value", "residual"], rows)


def compare_coefficients_for_point(point_name, eps_values, eta, args):
    point = get_point(point_name)
    min_power = args.fit_min_power
    max_power = args.fit_max_power
    rho = None if args.rho is None else parse_mpf(args.rho)

    print()
    print(f"######## Coefficient Comparison: {point.name} ########")
    print("epsList = {", ", ".join(str(eps) for eps in eps_values), "}", sep="")
    print(f"fit powers = {min_power}..{max_power}")
    print(f"N = {args.N}, scrambles = {args.scrambles}, eta = {eta}")

    method12_analytic = print_method12_coefficient_debug(point, min_power, max_power)
    method12_values = [closed_form_value_branch(point, eps, branch="-i0") for eps in eps_values]
    method12_fit = fit_laurent_coefficients(eps_values, method12_values, min_power, max_power)

    print_sample_table("== Method 1/2 Samples Used For Fit ==", eps_values, method12_values)
    print()
    print("== Method 1/2 Fit-From-Samples Coefficients ==")
    print_small_table(
        ["coefficient", "analytic", "fit_from_samples", "fit - analytic"],
        [
            [
                f"c[{power}]",
                fmt_complex(method12_analytic[power]),
                fmt_complex(method12_fit[power]),
                fmt_complex(method12_fit[power] - method12_analytic[power]),
            ]
            for power in range(min_power, max_power + 1)
        ],
    )
    print_fit_residuals("== Method 1/2 Fit Residuals ==", eps_values, method12_values, method12_fit)

    method3_values = []
    method3_errors = []
    method3_fit = None
    unsupported_reason = method3_unsupported_reason(point)
    if unsupported_reason is None:
        print()
        print("== Evaluating Method 3 At epsList ==")
        for eps in eps_values:
            stripped, stripped_err = method3_stripped(
                point=point,
                eps=eps,
                eta=eta,
                n_samples=args.N,
                rho=rho,
                scrambles=args.scrambles,
                seed=args.seed,
            )
            prefactor = method3_prefactor(point, eps)
            normalized = prefactor * stripped
            normalized_err = abs(prefactor) * stripped_err
            method3_values.append(normalized)
            method3_errors.append(normalized_err)
            print(f"  eps={eps}: method3_normalized={fmt_complex(normalized)} error={fmt_float(normalized_err)}")

        print_sample_table("== Method 3 Samples Used For Fit ==", eps_values, method3_values, method3_errors)
        method3_fit = fit_laurent_coefficients(eps_values, method3_values, min_power, max_power)
        print()
        print("== Method 3 Fitted Laurent Coefficients ==")
        print_small_table(["coefficient", "method3_fit"], coefficient_rows(method3_fit, min_power, max_power))
        print_fit_residuals("== Method 3 Fit Residuals ==", eps_values, method3_values, method3_fit)
    else:
        print()
        print(f"== Method 3 Fit Skipped: {unsupported_reason} ==")

    amflow_values = None
    amflow_fit = None
    if args.amflow == "fresh":
        print()
        print("== Evaluating AMFlow At epsList ==")
        print("Running AMFlow once for this point and reusing its epsilon expression for all eps values.")
        amflow_values, output = run_amflow_eps_values(point_name, eps_values, args)
        for eps, value in zip(eps_values, amflow_values):
            if value is None:
                print(f"  eps={eps}: AMFlow unavailable")
            else:
                print(f"  eps={eps}: AMFlow_original={fmt_complex(value)}")
        if any(value is None for value in amflow_values):
            print(output[-4000:])

        if all(value is not None for value in amflow_values):
            amflow_fit = fit_laurent_coefficients(eps_values, amflow_values, min_power, max_power)
            print()
            print("== AMFlow Fitted Laurent Coefficients ==")
            print_small_table(["coefficient", "AMFlow_fit"], coefficient_rows(amflow_fit, min_power, max_power))
            print_fit_residuals("== AMFlow Fit Residuals ==", eps_values, amflow_values, amflow_fit)

    print()
    print("== Final Coefficient Comparison ==")
    headers = [
        "coefficient",
        "method12_analytic",
        "method12_fit",
            "method3_fit",
            "method3_fit - method12_analytic",
    ]
    rows = []
    for power in range(min_power, max_power + 1):
        row = [
            f"c[{power}]",
            fmt_complex(method12_analytic[power]),
            fmt_complex(method12_fit[power]),
            fmt_complex(None if method3_fit is None else method3_fit[power]),
            fmt_complex(None if method3_fit is None else method3_fit[power] - method12_analytic[power]),
        ]
        rows.append(row)
    print_small_table(headers, rows)

    if amflow_fit is not None:
        print()
        print("== Final Coefficient Comparison With AMFlow ==")
        print_small_table(
            [
                "coefficient",
                "AMFlow_fit",
                "AMFlow_fit - method12_analytic",
                "method3_fit - AMFlow_fit",
            ],
            [
                [
                    f"c[{power}]",
                    fmt_complex(amflow_fit[power]),
                    fmt_complex(amflow_fit[power] - method12_analytic[power]),
                    fmt_complex(None if method3_fit is None else method3_fit[power] - amflow_fit[power]),
                ]
                for power in range(min_power, max_power + 1)
            ],
        )

    return {
        "point_name": point_name,
        "method12_analytic": method12_analytic,
        "method12_fit": method12_fit,
        "method3_fit": method3_fit,
        "AMFlow_fit": amflow_fit,
    }


def run_coefficient_mode(point_names, eps_values, eta, args):
    output_dir = repo_root() / "amflow-project" / "targets" / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    report_name = "two_loop_method3_coefficients_all_points.txt" if len(point_names) > 1 else f"two_loop_method3_coefficients_{point_names[0]}.txt"
    report_path = output_dir / report_name

    with report_path.open("w", encoding="utf-8") as report:
        tee = Tee(sys.stdout, report)
        with contextlib.redirect_stdout(tee):
            print("Coefficient-mode examples:")
            print("  ./run.sh two-loop-method3-mc-compare --all-points --eps 0.1 --amflow skip")
            print("  ./run.sh two-loop-method3-mc-compare --all-points --compare-coefficients --amflow skip")
            print("  ./run.sh two-loop-method3-mc-compare --point equal_mass_offshell_positive --compare-coefficients --eps-list 0.2,0.15,0.1,0.075,0.05 --amflow skip")
            results = [compare_coefficients_for_point(point_name, eps_values, eta, args) for point_name in point_names]

    print(f"Coefficient report written to: {report_path}")
    return results


def main():
    args = parse_args()

    if args.compare_coefficients:
        eps_values = parse_eps_list(args.eps_list)
    elif args.eps_list is not None:
        eps_values = parse_eps_list(args.eps_list)
    else:
        eps_values = [parse_mpf(args.eps)]

    point_names = sorted(POINTS) if args.all_points else [args.point]
    eta = parse_mpf(args.eta)

    if args.compare_coefficients:
        run_coefficient_mode(point_names, eps_values, eta, args)
        return

    rows = []
    for eps in eps_values:
        for point_name in point_names:
            if args.debug:
                point = get_point(point_name)
                safe_eps = str(eps).replace("/", "_over_").replace(".", "p")
                output_dir = repo_root() / "amflow-project" / "targets" / "output"
                output_dir.mkdir(parents=True, exist_ok=True)
                report_path = output_dir / f"two_loop_method3_debug_{point_name}_eps_{safe_eps}.txt"
                with report_path.open("w", encoding="utf-8") as report:
                    tee = Tee(sys.stdout, report)
                    with contextlib.redirect_stdout(tee):
                        row = compare_point(point_name, eps, eta, args)
                        print_debug_report(point, eps, eta, args, row)
                print(f"Debug report written to: {report_path}")
            else:
                row = compare_point(point_name, eps, eta, args)
            rows.append(row)

    print_table(rows)


if __name__ == "__main__":
    main()
