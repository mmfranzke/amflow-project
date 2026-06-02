import argparse
import contextlib
import cmath
import io
import os
import re
import subprocess
import math
import sys
import shutil
from pathlib import Path
from fractions import Fraction

import _path  # noqa: F401
import numpy as np
import mpmath as mp
from lsmethod.closed_forms import (
    product_closed_form,
    product_laurent_coefficients,
    residue_closed_form,
    residue_ingredients as closed_form_residue_ingredients,
    residue_laurent_coefficients,
)
from lsmethod.kernels import (
    l_residue_A1 as kernel_l_residue_A1,
    l_residue_A1_short_without_extra_term as kernel_l_residue_A1_short,
    l_residue_B1 as kernel_l_residue_B1,
    l_residue_closed_ingredients as kernel_l_residue_closed_ingredients,
    l_residue_delta0 as kernel_l_residue_delta0,
    l_residue_extra_offshell_term as kernel_l_residue_extra_offshell_term,
    l_residue_z as kernel_l_residue_z,
    method1_l_residue_closed as kernel_method1_l_residue_closed,
)
from lsmethod.method3 import (
    gamma4,
    gamma_xy,
    mapped_integrand_probe,
    method3_old_prefactor_diagnostic,
    method3_old_omega_prefactor_diagnostic,
    method3_cutoff_taylor_coefficients,
    method3_cutoff_integrand_expansion_coefficients,
    method3_prefactor,
    method3_prefactor_components,
    method3_stripped,
    method3_stripped_cutoff,
)
from lsmethod.two_loop_points import POINTS, get_point
from lsmethod.two_loop_points import TwoLoopPoint

DEFAULT_EPS_LIST = ["1/5", "3/20", "1/10", "3/40", "1/20"]
DEFAULT_KMAX_LIST = ["1", "2", "5", "10", "20", "50", "100", "200", "500"]
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


def parse_number_list(values, default_values):
    if values is None:
        values = default_values
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
    return delta1_value_current_pdf(point, lam)


def delta1_value_current_pdf(point, lam):
    mu_k_lam = (1 - lam) * point.mk2 + lam * point.Mk2
    return -mu_k_lam + lam * (1 - lam) * point.Ml2 + lam * (1 - lam * point.x) * point.p2


def delta1_value_old_compare(point, lam):
    mu_k_lam = (1 - lam) * point.mk2 + lam * point.Mk2
    return -mu_k_lam + lam * (1 - lam) * point.Ml2 + lam * (1 - point.x) * point.p2


def delta1_value_for_lambda(point, lam):
    return delta1_value_current_pdf(point, lam)


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


def gamma_sq_times_powers_laurent(delta0, delta1, prefactor, min_power=-2, max_power=2, branch="-i0"):
    if min_power < -2:
        raise ValueError("Gamma[eps]^2 variants are implemented from eps^-2 upward.")
    if max_power > 2:
        raise ValueError("Gamma[eps]^2 variants are implemented through eps^2.")

    s_log = delta_log(delta0, branch) + delta_log(delta1, branch)
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
    prefactor_c = complex(float(prefactor), 0.0)
    return {power: prefactor_c * exp_coeffs[power + 2] for power in range(min_power, max_power + 1)}


def gamma_sq_variant_value(delta0, delta1, prefactor, eps, branch="-i0"):
    d0_log = delta_log(delta0, branch)
    d1_log = delta_log(delta1, branch)
    return math.gamma(float(eps)) ** 2 * float(prefactor) * cmath.exp(-float(eps) * (d0_log + d1_log))


def repo_root():
    return Path(__file__).resolve().parents[2]


def safe_name(value):
    return str(value).replace("/", "_over_").replace(".", "p").replace("-", "m").replace("+", "p")


def parse_args():
    examples = """examples:
  ./run.sh two-loop-method3-mc-compare --all-points --eps 0.1 --amflow skip
  ./run.sh two-loop-method3-mc-compare --all-points --compare-coefficients --amflow skip
  ./run.sh two-loop-method3-mc-compare --point equal_mass_offshell_positive --compare-coefficients --eps-list 0.2,0.15,0.1,0.075,0.05 --amflow skip
  ./run.sh two-loop-method3-mc-compare --point equal_mass_offshell_positive --eps 0.1 --amflow skip --debug-normalization
  ./run.sh two-loop-method3-mc-compare --point equal_mass_offshell_positive --radial-cutoff-scan --eps 0.1 --Kmax-list 500,1000,2000 --amflow skip --debug-normalization
  ./run.sh two-loop-method3-mc-compare --point equal_mass_offshell_positive --radial-cutoff-scan --eps-list 0.2,0.15,0.1,0.075,0.05 --Kmax-list 500,1000,2000 --amflow skip --debug-normalization
  ./run.sh two-loop-method3-mc-compare --point equal_mass_offshell_positive --radial-cutoff-scan --eps-list 0.2,0.15,0.1,0.075,0.05 --Kmax-list 100,200,500,1000,2000 --amflow skip --ratio-fit
  ./run.sh two-loop-method3-mc-compare --point equal_mass_offshell_positive --radial-cutoff-scan --tail-extrapolate --eps-list 0.2,0.15,0.1,0.075,0.05 --Kmax-list 100,200,500,1000,2000 --amflow skip
  ./run.sh two-loop-method3-mc-compare --point equal_mass_offshell_positive --method3-coeff-cutoff-log --amflow skip
  ./run.sh two-loop-method3-mc-compare --point equal_mass_offshell_positive --method1-l-residue-direct-check --eps 0.1 --amflow reuse
  ./run.sh two-loop-method3-mc-compare --point equal_mass_offshell_positive --method1-l-residue-closed-check --eps 0.1 --amflow imported
  ./run.sh two-loop-method3-mc-compare --point equal_mass_offshell_positive --compare-closed-form --amflow imported
  ./run.sh two-loop-method3-mc-compare --point equal_mass_offshell_positive --residue-1d-check --eps 0.1 --amflow reuse
  ./run.sh two-loop-method3-mc-compare --show-amflow-import-help
  ./run.sh two-loop-method3-mc-compare --point equal_mass_offshell_positive --compare-coefficients --amflow reuse --diagnose-method12-conventions
"""
    parser = argparse.ArgumentParser(description=__doc__, epilog=examples, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--point", choices=sorted(POINTS), default="equal_mass_onshell_branch")
    parser.add_argument("--points", default=None, help="Comma-separated point names.")
    parser.add_argument("--all-points", action="store_true")
    parser.add_argument("--eps", default="1/10")
    parser.add_argument("--eps-list", nargs="*", default=None)
    parser.add_argument("--eta", default="1e-6")
    parser.add_argument("--N", type=int, default=2**18)
    parser.add_argument("--scrambles", type=int, default=4)
    parser.add_argument("--rho", default=None)
    parser.add_argument("--lc-jacobian", default="1/4")
    parser.add_argument("--method3-label-convention", choices=["pdf", "amflow"], default="pdf")
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--amflow", choices=["fresh", "reuse", "imported", "skip"], default="fresh")
    parser.add_argument("--amflow-import-log", default=None)
    parser.add_argument("--amflow-import-result", default=None)
    parser.add_argument("--amflow-import-trust-metadata", action="store_true")
    parser.add_argument("--show-amflow-import-help", action="store_true")
    parser.add_argument("--amflow-eps-order", default="4")
    parser.add_argument("--precision-goal", default="10")
    parser.add_argument("--dps", type=int, default=50)
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--debug-include-2p20", action="store_true")
    parser.add_argument("--debug-normalization", action="store_true")
    parser.add_argument("--radial-cutoff-scan", action="store_true")
    parser.add_argument("--Kmax-shells", action="store_true")
    parser.add_argument("--Kmax-list", nargs="*", default=None)
    parser.add_argument("--ratio-fit", action="store_true")
    parser.add_argument("--max-relative-error-for-ratio-fit", type=float, default=0.2)
    parser.add_argument("--tail-extrapolate", action="store_true")
    parser.add_argument("--tail-fit-Kmin", default=None)
    parser.add_argument("--tail-fit-Kmax", default=None)
    parser.add_argument("--tail-fit-last-n", type=int, default=None)
    parser.add_argument("--tail-fit-unweighted", action="store_true")
    parser.add_argument("--tail-fit-self-test", action="store_true")
    parser.add_argument("--method3-coeff-cutoff-log", action="store_true")
    parser.add_argument("--method3-integrand-eps-expansion", "--transverse-integrand-expansion-check", action="store_true", dest="method3_integrand_eps_expansion")
    parser.add_argument("--method1-l-residue-direct-check", "--residue-1d-check", action="store_true", dest="method1_l_residue_direct_check")
    parser.add_argument("--method1-l-residue-closed-check", "--compare-closed-form", action="store_true", dest="compare_closed_form")
    parser.add_argument("--transverse-mc-check", action="store_true")
    parser.add_argument("--method3-coeff-order", type=int, default=2)
    parser.add_argument("--method3-coeff-eps-step", default="1e-3")
    parser.add_argument("--log-fit-Kmin", default="20")
    parser.add_argument("--log-fit-degree", type=int, default=None)
    parser.add_argument("--reconstruction-eps", default=None)
    parser.add_argument("--assume-cutoff-dimreg-map", action="store_true")
    parser.add_argument("--compare-coefficients", "--eps-fit", action="store_true")
    parser.add_argument("--diagnose-method12-conventions", "--diagnose-product-form", action="store_true", dest="diagnose_method12_conventions")
    parser.add_argument("--fit-min-power", type=int, default=-2)
    parser.add_argument("--fit-max-power", type=int, default=2)
    return parser.parse_args()


def complex_from_wolfram(text):
    text = text.strip()
    if text.startswith("InputForm[") and text.endswith("]"):
        text = text[len("InputForm[") : -1].strip()
    text = text.replace("\\\n", "")
    text = text.replace("*^", "e")
    text = re.sub(r"(?<=\d)`{1,2}[0-9.]*", "", text)
    text = re.sub(r"(?<=\.)`{1,2}[0-9.]*", "", text)
    text = text.replace("I", "j")
    text = re.sub(r"\s+", "", text)
    if text.startswith("Complex[") and text.endswith("]"):
        inner = text[len("Complex[") : -1]
        real, imag = inner.split(",", 1)
        return complex(float(real), float(imag))
    return complex(eval(text, {"__builtins__": {}}, {"j": 1j}))


def parse_amflow_coefficients(output):
    raw = {}
    normalized = {}
    for match in re.finditer(r"AMFLOW_COEFF_POWER=([-+]?\d+)\s+RAW=(.+?)\s+NORMALIZED=(.+)", output):
        power = int(match.group(1))
        try:
            raw[power] = complex_from_wolfram(match.group(2))
            normalized[power] = complex_from_wolfram(match.group(3))
        except Exception:
            continue
    return raw, normalized


def parse_wolfram_association_coefficients(text, key):
    text = text.replace("\\\n", "")
    match = re.search(rf'"{re.escape(key)}"\s*->\s*<\|(.+?)\|>', text, flags=re.S)
    if not match:
        return {}
    body = match.group(1)
    coeffs = {}
    pattern = re.compile(r"([-+]?\d+)\s*->\s*(.+?)(?=,\s*[-+]?\d+\s*->|\s*$)", flags=re.S)
    for entry in pattern.finditer(body):
        power = int(entry.group(1))
        value_text = entry.group(2).strip().rstrip(",")
        try:
            coeffs[power] = complex_from_wolfram(value_text)
        except Exception:
            continue
    return coeffs


def parse_imported_amflow_coefficients(text):
    raw, normalized = parse_amflow_coefficients(text)
    if normalized:
        return raw, normalized
    raw = parse_wolfram_association_coefficients(text, "rawPieceCombinationCoefficients")
    normalized = parse_wolfram_association_coefficients(text, "doubleDiscontinuityCoefficients")
    if normalized:
        return raw, normalized

    in_raw_amflow_table = False
    for line in text.splitlines():
        if line.strip() == "== Raw AMFlow Laurent Coefficients ==":
            in_raw_amflow_table = True
            continue
        if in_raw_amflow_table and line.startswith("== "):
            break
        if not in_raw_amflow_table:
            continue
        if not line.strip().startswith("| c["):
            continue
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        if len(cells) < 3:
            continue
        power_match = re.fullmatch(r"c\[([-+]?\d+)\]", cells[0])
        if not power_match:
            continue
        power = int(power_match.group(1))
        try:
            if "not run" not in cells[1]:
                raw[power] = complex_from_wolfram(cells[1])
            if "not run" not in cells[2]:
                normalized[power] = complex_from_wolfram(cells[2])
        except Exception:
            continue
    return raw, normalized


def parse_amflow_object_lines(output):
    fields = {}
    for line in output.splitlines():
        if line.startswith("AMFLOW_OBJECT_"):
            key, _, value = line.partition("=")
            fields[key.removeprefix("AMFLOW_OBJECT_")] = value.strip()
    return fields


def parse_amflow_combo_coefficients(output):
    combos = {}
    pattern = re.compile(
        r"AMFLOW_COMBO_COEFF NAME=([A-Za-z0-9_]+)\s+POWER=([-+]?\d+)\s+"
        r"VALUE=(.+?)\s+METHOD12=(.+?)\s+DIFF=(.+?)\s+RATIO=(.+?)\s+SCORE=(.+)"
    )
    for match in pattern.finditer(output):
        name = match.group(1)
        power = int(match.group(2))
        try:
            entry = {
                "value": complex_from_wolfram(match.group(3)),
                "method12": complex_from_wolfram(match.group(4)),
                "diff": complex_from_wolfram(match.group(5)),
                "ratio": complex_from_wolfram(match.group(6)),
                "score": complex_from_wolfram(match.group(7)).real,
            }
        except Exception:
            continue
        combos.setdefault(name, {})[power] = entry
    return combos


def amflow_import_log_dir():
    return repo_root() / "amflow-project" / "logs" / "imported_amflow"


def amflow_import_result_dir():
    return repo_root() / "amflow-project" / "results" / "imported_amflow"


def ensure_amflow_import_dirs():
    amflow_import_log_dir().mkdir(parents=True, exist_ok=True)
    amflow_import_result_dir().mkdir(parents=True, exist_ok=True)


def print_amflow_import_help():
    ensure_amflow_import_dirs()
    log_dir = amflow_import_log_dir()
    result_dir = amflow_import_result_dir()
    print("AMFlow import workflow")
    print()
    print("Put old AMFlow coefficient logs here:")
    print(f"  {log_dir}")
    print()
    print("Put old AMFlow result associations or .res files here:")
    print(f"  {result_dir}")
    print()
    print("Suggested filenames:")
    print("  logs/imported_amflow/<point>_amflow_coefficients.log")
    print("  results/imported_amflow/<point>_amflow_result.wl")
    print("  results/imported_amflow/<point>_amflow_result.res")
    print()
    print("Example import commands:")
    print("  ./run.sh two-loop-method3-mc-compare --amflow-import-log /path/to/amflow-1553698.out --point equal_mass_offshell_positive --amflow skip")
    print("  cp /path/to/amflow-1553698.out logs/imported_amflow/equal_mass_offshell_positive_amflow_coefficients.log")
    print()
    print("Example comparison command:")
    print("  ./run.sh two-loop-method3-mc-compare --point equal_mass_offshell_positive --method1-l-residue-closed-check --eps 0.1 --amflow imported")
    print()
    print("Metadata policy:")
    print("  Imported files should contain AMFLOW_OBJECT_POINT=<point> and AMFLOW_COEFF_POWER lines.")
    print("  If full metadata is missing, use --amflow-import-trust-metadata only after manually checking the file.")


def process_amflow_import_args(args):
    ensure_amflow_import_dirs()
    copied = []
    if args.amflow_import_log:
        source = Path(args.amflow_import_log).expanduser()
        if not source.exists():
            raise SystemExit(f"AMFlow import log does not exist: {source}")
        dest = amflow_import_log_dir() / source.name
        shutil.copy2(source, dest)
        copied.append(dest)
    if args.amflow_import_result:
        source = Path(args.amflow_import_result).expanduser()
        if not source.exists():
            raise SystemExit(f"AMFlow import result does not exist: {source}")
        dest = amflow_import_result_dir() / source.name
        shutil.copy2(source, dest)
        copied.append(dest)
    for dest in copied:
        print(f"Imported AMFlow file copied to: {dest}")


def imported_amflow_candidate_paths(point_name):
    ensure_amflow_import_dirs()
    patterns = [
        f"{point_name}_amflow_coefficients.log",
        f"*{point_name}*coefficients*.log",
        f"*{point_name}*.log",
        f"*{point_name}*.out",
    ]
    paths = []
    for pattern in patterns:
        paths.extend(sorted(amflow_import_log_dir().glob(pattern)))
    paths.extend(sorted(amflow_import_result_dir().glob(f"*{point_name}*")))
    seen = set()
    unique = []
    for path in paths:
        if path not in seen and path.is_file():
            unique.append(path)
            seen.add(path)
    return unique


def imported_metadata_matches_point(text, point_name):
    point = get_point(point_name)
    has_name = (
        f"AMFLOW_OBJECT_POINT={point_name}" in text
        or f'"name" -> "{point_name}"' in text
        or f"name -> {point_name}" in text
        or f"point: {point_name}" in text
        or f"point {point_name}" in text
    )
    expected_fragments = [
        f"pPlus={point.p_plus}",
        f"pMinus={point.p_minus}",
        f"pPerp2={point.p_perp2}",
        f"p2={point.p2}",
        f"x={point.x}",
        f"y={point.y}",
        f"ml2={point.ml2}",
        f"Ml2={point.Ml2}",
        f"mk2={point.mk2}",
        f"Mk2={point.Mk2}",
    ]
    has_full_kinematic_line = all(fragment in text for fragment in expected_fragments)
    has_metadata_association = "pointMetadata" in text and all(str(value) in text for value in [
        point.p_plus,
        point.p_minus,
        point.p2,
        point.p_perp2,
        point.x,
        point.y,
        point.ml2,
        point.Ml2,
        point.mk2,
        point.Mk2,
    ])
    return has_name and (has_full_kinematic_line or has_metadata_association)


def load_imported_amflow_coefficients(point_name, args):
    candidates = imported_amflow_candidate_paths(point_name)
    if not candidates:
        print(f"No imported AMFlow file found for point {point_name}.")
        print_amflow_import_help()
        return {}, {}, {}, {}, ""

    missing_metadata_paths = []
    for path in candidates:
        text = path.read_text(encoding="utf-8", errors="replace")
        raw, normalized = parse_imported_amflow_coefficients(text)
        if not normalized:
            continue
        if not imported_metadata_matches_point(text, point_name):
            missing_metadata_paths.append(path)
            if not args.amflow_import_trust_metadata:
                continue
        print(f"AMFlow imported coefficients loaded from: {path}")
        if missing_metadata_paths and args.amflow_import_trust_metadata:
            print("WARNING: imported AMFlow file lacks full metadata; trusting it because --amflow-import-trust-metadata was passed.")
        return raw, normalized, parse_amflow_object_lines(text), parse_amflow_combo_coefficients(text), text

    if missing_metadata_paths and not args.amflow_import_trust_metadata:
        print("Imported AMFlow file lacks full metadata. Please confirm or rerun AMFlow fresh.")
        print("Candidate files with coefficients but insufficient metadata:")
        for path in missing_metadata_paths:
            print(f"  {path}")
        print("Use --amflow-import-trust-metadata only after manually checking the selected point.")
    else:
        print(f"No parseable imported AMFlow coefficient file found for point {point_name}.")
    return {}, {}, {}, {}, ""


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
    if args.amflow == "reuse":
        env["AMFLOW_REUSE_RESULT"] = "1"

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
    if args.amflow == "reuse":
        env["AMFLOW_REUSE_RESULT"] = "1"

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


def run_amflow_coefficients(point_name, eps, args):
    if args.amflow == "imported":
        return load_imported_amflow_coefficients(point_name, args)

    repo = repo_root()
    run_sh = repo / "amflow-project" / "run.sh"
    log_dir = repo / "amflow-project" / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"two-loop-method3-mc-compare_amflow_{point_name}_coefficients.log"
    env = os.environ.copy()
    env["TWO_LOOP_POINT"] = point_name
    env["TWO_LOOP_EPS_VALUE"] = str(float(eps))
    env["AMFLOW_EPS_ORDER"] = args.amflow_eps_order
    env["AMFLOW_PRECISION_GOAL"] = args.precision_goal
    if args.amflow == "reuse":
        env["AMFLOW_REUSE_RESULT"] = "1"

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
    print(f"AMFlow coefficient log written to: {log_path}")
    if proc.returncode != 0:
        return {}, {}, {}, {}, proc.stdout

    raw, normalized = parse_amflow_coefficients(proc.stdout)
    object_fields = parse_amflow_object_lines(proc.stdout)
    combo_coeffs = parse_amflow_combo_coefficients(proc.stdout)
    return raw, normalized, object_fields, combo_coeffs, proc.stdout


def ensure_amflow_order_for_powers(args, max_power):
    required_order = max_power + 4
    actual_order = int(args.amflow_eps_order)
    if actual_order < required_order:
        print(
            "Requested AMFlow coefficient comparison through "
            f"c[{max_power}], but --amflow-eps-order {args.amflow_eps_order} only provides through "
            f"c[{actual_order - 4}]. AMFlow tables will include only the available powers; "
            f"rerun with --amflow-eps-order {required_order} for higher AMFlow coefficients."
        )


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


def relative_error(error, value):
    if error is None or value is None:
        return None
    try:
        error_f = float(error)
        if math.isnan(error_f) or value == 0:
            return None
        return error_f / abs(value)
    except Exception:
        return None


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


def method3_point_for_convention(point, convention):
    if convention == "pdf":
        return point
    if convention == "amflow":
        return TwoLoopPoint(
            name=point.name + "_method3_amflow_labels",
            p_plus=point.p_plus,
            p_minus=point.p_minus,
            p_perp2=point.p_perp2,
            x=point.y,
            y=point.x,
            ml2=point.ml2,
            Ml2=point.Ml2,
            mk2=point.mk2,
            Mk2=point.Mk2,
        )
    raise ValueError(f"Unknown method3 label convention: {convention}")


def print_method3_label_convention(point, method3_point, convention):
    print(f"Method3 label convention: {convention}")
    if convention == "pdf":
        print("  pdf: method3 uses x=k^+, y=l^+.")
    else:
        print("  amflow: AMFlow uses x=n.l, y=n.k, so method3 receives x_method3=y_point and y_method3=x_point.")
    print(f"  point labels: x={point.x}, y={point.y}; method3 labels: x={method3_point.x}, y={method3_point.y}")
    if convention == "pdf" and point.x != point.y:
        print("WARNING: AMFlow and method3 PDF attach x,y to opposite loops. Consider --method3-label-convention amflow.")


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
    d1_old = delta1_value_old_compare(point, point.lam)
    print()
    print(f"== {point.name} ==")
    print(f"support 0 < x, 0 < y, x + y < pPlus: {point.in_support()}")
    print(f"pPlus={point.p_plus}, pMinus={point.p_minus}, pPerp2={point.p_perp2}, p2={point.p2}")
    print(f"x={point.x}, y={point.y}, lambda={point.lam}")
    print(f"Delta0 = {d0}")
    print(f"Delta1 current PDF = {d1}")
    print(f"Delta1 old compare diagnostic = {d1_old}")
    if point.p2 != 0:
        print("WARNING: off-shell point is sensitive to the Delta1 p2-term convention.")
        print("Current-PDF default uses lambda*(1-lambda*x)*p2. Old compare used lambda*(1-x)*p2.")
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
    method3_point = method3_point_for_convention(point, args.method3_label_convention)
    print_delta_diagnostics(point, eta)
    print_method3_label_convention(point, method3_point, args.method3_label_convention)

    rho = None if args.rho is None else parse_mpf(args.rho)
    unsupported_reason = method3_unsupported_reason(method3_point)
    if unsupported_reason is None:
        stripped, stripped_err = method3_stripped(
            point=method3_point,
            eps=eps,
            eta=eta,
            n_samples=args.N,
            rho=rho,
            scrambles=args.scrambles,
            seed=args.seed,
        )
        lc_jacobian = parse_float(args.lc_jacobian)
        pref = method3_prefactor(method3_point, eps, lc_jacobian=lc_jacobian)
        old_omega_pref = method3_old_omega_prefactor_diagnostic(method3_point, eps, lc_jacobian=lc_jacobian)
        old_pref = method3_old_prefactor_diagnostic(method3_point, eps)
        normalized = complex(pref) * stripped
        normalized_err = abs(complex(pref)) * stripped_err
        normalized_old_omega = complex(old_omega_pref) * stripped
        normalized_old_prefactor = complex(old_pref) * stripped
    else:
        print(f"Method 3 skipped: {unsupported_reason}")
        stripped = None
        stripped_err = None
        pref = None
        old_omega_pref = None
        old_pref = None
        normalized = None
        normalized_err = None
        normalized_old_omega = None
        normalized_old_prefactor = None

    method12 = closed_form_value(point, eps, Fraction(1, 10**30))
    d0_pow = complex(float(delta0_value(point)), -1e-30) ** (-float(eps))
    d1_pow = complex(float(delta1_value(point)), -1e-30) ** (-float(eps))

    amflow = None
    if args.amflow in ("fresh", "reuse", "imported"):
        if args.amflow == "reuse":
            print("Reusing AMFlow fixed-eps original integral expression for this point...")
        else:
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
        "method3_normalized_old_omega_diag": normalized_old_omega,
        "method3_normalized_old_prefactor_diag": normalized_old_prefactor,
        "method12_closed": method12,
        "AMFlow_original": amflow,
        "debug_normalization": args.debug_normalization,
        "lc_jacobian": parse_float(args.lc_jacobian),
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
    d1_old = delta1_value_old_compare(point, point.lam)
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
        ["Delta1 current PDF", d1],
        ["Delta1 old compare diagnostic", d1_old],
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
    print("D4 = A(l,k;x,y) - 2 l k cos(theta) + i eta gamma4")
    print("gamma4 = 1 + (x + y) (pPlus - x - y)/(x y)")
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
        ["old gamma4 without +1", gamma_xy(point)],
        ["gamma4", gamma4(point)],
        ["gamma4 - old", gamma4(point) - float(gamma_xy(point))],
        ["gamma4 > 0", gamma4(point) > 0],
    ]
    print_small_table(["check", "value"], rows)


def method3_prefactor_component_rows(point, eps, lc_jacobian):
    components = method3_prefactor_components(point, eps, lc_jacobian=lc_jacobian)
    p_plus = float(point.p_plus)
    x = float(point.x)
    y = float(point.y)
    one_minus_x = p_plus - x
    lam = y / one_minus_x if one_minus_x != 0 else None
    prefactor = method3_prefactor(point, eps, lc_jacobian=lc_jacobian)
    return [
        ["eps", eps],
        ["d = 4 - 2 eps", f"{components['d']:.16e}"],
        ["n = d - 2", f"{components['n']:.16e}"],
        ["x", point.x],
        ["y", point.y],
        ["pPlus", point.p_plus],
        ["lambda = y/(pPlus-x)", "undefined" if lam is None else f"{lam:.16e}"],
        ["oneMinusX = pPlus - x", f"{one_minus_x:.16e}"],
        ["1 - x diagnostic", f"{1.0 - x:.16e}"],
        ["lcJacobian", f"{components['lc_jacobian']:.16e}"],
        ["contour factor = (2 Pi I)^2", fmt_complex(components["contour_factor"])],
        ["loop normalization denominator = (I Pi^(d/2))^2", fmt_complex(components["loop_normalization"])],
        ["Omega[n]", f"{components['omega_n']:.16e}"],
        ["Omega[n-1]", f"{components['omega_n_minus_1']:.16e}"],
        ["angular prefactor = Omega[n] Omega[n-1]", f"{components['angular_prefactor']:.16e}"],
        ["xy factor = 1/(x y)", f"{components['one_over_xy']:.16e}"],
        ["full method3_prefactor", fmt_complex(prefactor)],
    ]


def print_method3_normalization_definition(point, eps, stripped, lc_jacobian):
    prefactor = method3_prefactor(point, eps, lc_jacobian=lc_jacobian)
    print()
    print("== Method 3 Normalization Definitions ==")
    print("method3_stripped: raw MC integral over dk dl dtheta, including the map Jacobian, excluding the global Method-3 prefactor.")
    print("method3_prefactor: global factor multiplying method3_stripped to obtain method3_normalized.")
    print("method3_normalized: method3_prefactor * method3_stripped.")
    rows = method3_prefactor_component_rows(point, eps, lc_jacobian)
    rows.append(["method3_stripped", fmt_complex(stripped)])
    rows.append(["normalized/stripped ratio", fmt_complex(prefactor)])
    print_small_table(["component", "value"], rows)


def normalization_variant_rows(point, current, method12, lc_jacobian):
    p_plus = float(point.p_plus)
    x = float(point.x)
    one_minus_x = p_plus - x
    variants = [
        ("current = N3 * stripped", current),
        ("flip_sign = -current", None if current is None else -current),
        ("no_lc_jacobian = current / lcJacobian", None if current is None or lc_jacobian == 0 else current / lc_jacobian),
        ("double_lc_jacobian = current * lcJacobian", None if current is None else current * lc_jacobian),
        ("times_2", None if current is None else 2 * current),
        ("times_1_over_2", None if current is None else current / 2),
        ("times_4", None if current is None else 4 * current),
        ("times_1_over_4", None if current is None else current / 4),
        ("times_1_minus_x", None if current is None else current * one_minus_x),
        ("divided_by_1_minus_x", None if current is None or one_minus_x == 0 else current / one_minus_x),
        ("times_y_jacobian", None if current is None else current * one_minus_x),
        ("lambda_density", None if current is None else current * one_minus_x),
        ("y_density_from_lambda", None if current is None or one_minus_x == 0 else current / one_minus_x),
    ]
    return [[name, fmt_complex(value), fmt_complex(ratio_or_none(value, method12))] for name, value in variants]


def print_normalization_variant_audit(point, current, method12, lc_jacobian, title="== Normalization Variant Audit =="):
    print()
    print(title)
    print_small_table(["variant_name", "value", "ratio_to_method12"], normalization_variant_rows(point, current, method12, lc_jacobian))
    print(
        "Note: If method12 is differential in lambda = y/(pPlus-x), but method3 is differential in y, "
        "then compare method3_y_density to method12_lambda_density with dy = (pPlus-x) dlambda. "
        "This convention must be checked against the derivation."
    )


def print_method12_variable_audit(point, eps, method12):
    x = float(point.x)
    y = float(point.y)
    p_plus = float(point.p_plus)
    one_minus_x = p_plus - x
    xval = 1.0 - x
    lam = y / one_minus_x if one_minus_x != 0 else None
    Xval = 1.0 - x
    tol = 1e-14
    print()
    print("== Method 1/2 Variable Convention Audit ==")
    print("implemented formula: Gamma[eps]^2 / Xval * Delta0^(-eps) * Delta1^(-eps)")
    print("current Python implementation sets Xval = 1 - x.")
    rows = [
        ["eps", eps],
        ["Xval", f"{Xval:.16e}"],
        ["pPlus - x", f"{one_minus_x:.16e}"],
        ["1 - x", f"{xval:.16e}"],
        ["lambda = y/(pPlus-x)", "undefined" if lam is None else f"{lam:.16e}"],
        ["y", point.y],
        ["Xval == pPlus - x", abs(Xval - one_minus_x) < tol],
        ["Xval == 1 - x", abs(Xval - xval) < tol],
        ["method12_closed", fmt_complex(method12)],
    ]
    print_small_table(["quantity", "value"], rows)
    print("Diagnostic: the closed form is written in terms of y only through lambda = y/Xval, while the prefactor is 1/Xval.")
    print("Diagnostic: this is consistent with a lambda-to-y Jacobian, but the derivation must decide whether the target density is in y or lambda.")


def print_amflow_delta_convention_audit():
    print()
    print("== AMFlow Delta Convention Audit ==")
    print("Visible in amflow-project/families/TwoLoopKernelUncutFamilies.wl:")
    print("  AMFlowInfo[\"Propagator\"] contains sx (x - n l) and sy (y - n k).")
    print("  The uncut GaugeLink denominators are combined by a double discontinuity.")
    print("Visible in CompareTwoLoopFixedEps.wl:")
    print("  discExpr = (PP + PM + MP + MM)/(2 Pi I)^2.")
    print("Inferred convention from the code: delta constraints correspond to n.l = x and n.k = y.")
    print("No explicit denominator of the form n.k/(pPlus-x) - lambda is visible in the AMFlow family.")


def print_method3_prefactor_debug(point, eps, stripped, lc_jacobian):
    components = method3_prefactor_components(point, eps, lc_jacobian=lc_jacobian)
    prefactor = method3_prefactor(point, eps, lc_jacobian=lc_jacobian)
    old_omega_prefactor = method3_old_omega_prefactor_diagnostic(point, eps, lc_jacobian=lc_jacobian)
    old_prefactor = method3_old_prefactor_diagnostic(point, eps)
    xy = float(point.x * point.y)
    one_minus_x = float(1 - point.x)
    variants = {
        "prefactor": prefactor,
        "old_omega_prefactor_diagnostic": old_omega_prefactor,
        "old_prefactor_diagnostic": old_prefactor,
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
        ["Omega[n]", f"{components['omega_n']:.16e}"],
        ["Omega[n-1]", f"{components['omega_n_minus_1']:.16e}"],
        ["old diagnostic Omega[n-2]", f"{components['omega_n_minus_2']:.16e}"],
        ["new angularPrefactor = Omega[n] Omega[n-1]", f"{components['angular_prefactor']:.16e}"],
        ["old angularPrefactor = Omega[n-1] Omega[n-2]", f"{components['old_angular_prefactor']:.16e}"],
        ["lcJacobian", f"{components['lc_jacobian']:.16e}"],
        ["(2 Pi I)^2", fmt_complex(components["contour_factor"])],
        ["(I Pi^(d/2))^2", fmt_complex(components["loop_normalization"])],
        ["1/(x y)", f"{components['one_over_xy']:.16e}"],
        ["contour sign convention", "current code sign kept"],
        ["prefactor", fmt_complex(prefactor)],
        ["old Omega prefactor diagnostic with lcJacobian", fmt_complex(old_omega_prefactor)],
        ["old prefactor diagnostic", fmt_complex(old_prefactor)],
        ["new / old Omega prefactor", fmt_complex(prefactor / old_omega_prefactor if old_omega_prefactor != 0 else None)],
        ["new / old prefactor", fmt_complex(prefactor / old_prefactor if old_prefactor != 0 else None)],
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
                ["old gamma4 without +1", probe["old_Gamma_xy_without_plus_one"]],
                ["gamma4", probe["gamma4"]],
                ["D4", fmt_complex(probe["D4"])],
                ["measure", f"{probe['measure']:.16e}"],
                ["integrand_without_jacobian", fmt_complex(probe["integrand_without_jacobian"])],
                ["integrand_with_jacobian", fmt_complex(probe["integrand_with_jacobian"])],
                ["real sign D2", probe["D2_real_sign"]],
                ["real sign D4", probe["D4_real_sign"]],
            ],
        )


def print_cutoff_scan(point, eps, eta, n_samples, scrambles, seed, lc_jacobian):
    prefactor = method3_prefactor(point, eps, lc_jacobian=lc_jacobian)
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


def print_rho_scan(point, eps, eta, n_samples, scrambles, seed, lc_jacobian):
    prefactor = method3_prefactor(point, eps, lc_jacobian=lc_jacobian)
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


def print_n_scan(point, eps, eta, rho, scrambles, seed, include_2p20, lc_jacobian):
    prefactor = method3_prefactor(point, eps, lc_jacobian=lc_jacobian)
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


def print_eta_scan(point, eps, n_samples, rho, scrambles, seed, lc_jacobian):
    prefactor = method3_prefactor(point, eps, lc_jacobian=lc_jacobian)
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


def print_comparison_levels(row, lc_jacobian):
    point = get_point(row["point_name"])
    prefactor = method3_prefactor(point, row["eps"], lc_jacobian=lc_jacobian) if method3_supported(point) else None
    print()
    print("== Debug: Comparison Levels And Ratios ==")
    print_small_table(
        ["quantity", "value"],
        [
            ["method3_stripped", fmt_complex(row["method3_stripped"])],
            ["method3_prefactor", fmt_complex(prefactor)],
            ["method3_normalized", fmt_complex(row["method3_normalized"])],
            ["method3_normalized_old_omega_diag", fmt_complex(row["method3_normalized_old_omega_diag"])],
            ["method3_normalized_old_prefactor_diag", fmt_complex(row["method3_normalized_old_prefactor_diag"])],
            ["method12_closed", fmt_complex(row["method12_closed"])],
            ["method12_closed / method3_stripped", fmt_complex(row["ratio_method12_over_method3_stripped"])],
            ["method12_closed / method3_normalized", fmt_complex(row["ratio_method12_over_method3_normalized"])],
        ],
    )


def print_debug_report(point, eps, eta, args, row):
    rho = None if args.rho is None else parse_mpf(args.rho)
    lc_jacobian = parse_float(args.lc_jacobian)
    print()
    print(f"######## Method 3 Debug Report: {point.name}, eps={eps}, eta={eta} ########")
    print_method12_debug(point, eps)
    print_method12_variable_audit(point, eps, row["method12_closed"])
    print_amflow_delta_convention_audit()
    print_method3_formula_debug(point)
    print_comparison_levels(row, lc_jacobian)
    unsupported_reason = method3_unsupported_reason(point)
    if unsupported_reason is None:
        print_method3_normalization_definition(point, eps, row["method3_stripped"], lc_jacobian)
        print_normalization_variant_audit(point, row["method3_normalized"], row["method12_closed"], lc_jacobian)
        print_method3_prefactor_debug(point, eps, row["method3_stripped"], lc_jacobian)
        print_integrand_probes(point, eps, eta, rho)
        print_cutoff_scan(point, eps, eta, args.N, args.scrambles, args.seed, lc_jacobian)
        print_rho_scan(point, eps, eta, args.N, args.scrambles, args.seed, lc_jacobian)
        print_n_scan(point, eps, eta, rho, args.scrambles, args.seed, args.debug_include_2p20, lc_jacobian)
        print_eta_scan(point, eps, args.N, rho, args.scrambles, args.seed, lc_jacobian)
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
                [
                    "method3_old_omega_diag",
                    fmt_complex(row["method3_normalized_old_omega_diag"]),
                    "diagnostic only",
                ],
                [
                    "method3_old_prefactor_diag",
                    fmt_complex(row["method3_normalized_old_prefactor_diag"]),
                    "diagnostic only",
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
        if row.get("debug_normalization"):
            point = get_point(row["point_name"])
            lc_jacobian = row["lc_jacobian"]
            print_method12_variable_audit(point, row["eps"], row["method12_closed"])
            print_amflow_delta_convention_audit()
            if row["method3_stripped"] is not None:
                print_method3_normalization_definition(point, row["eps"], row["method3_stripped"], lc_jacobian)
                print_normalization_variant_audit(point, row["method3_normalized"], row["method12_closed"], lc_jacobian)

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
        ["Delta1 current PDF", d1],
        ["Delta1 old compare diagnostic", delta1_value_old_compare(point, point.lam)],
        ["branch", branch],
    ]
    if point.p2 != 0:
        rows.append(["off-shell warning", "current PDF uses lambda*(1-lambda*x)*p2; old compare used lambda*(1-x)*p2"])
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


def print_amflow_object_report(point, object_fields, raw_coeffs, normalized_coeffs, args):
    print()
    print("== AMFlow Object And Normalization ==")
    one_minus_x = point.p_plus - point.x
    rows = [
        ["point name", point.name],
        ["family / target name", object_fields.get("TARGET", "compare-twoloop-fixed-eps")],
        ["families", object_fields.get("FAMILIES", "not printed")],
        ["selected integrals", object_fields.get("TARGET_INTEGRALS", "j[family,1,1,1,1,1,1,0] for PP,PM,MP,MM")],
        ["selected pieces / prescription", object_fields.get("PRESCRIPTION", "PP+PM+MP+MM double discontinuity")],
        ["propagator exponents", object_fields.get("PROPAGATOR_EXPONENTS", "{1,1,1,1,1,1,0}")],
        ["linear denominators", object_fields.get("LINEAR_DENOMINATORS", "{sx (x - n.l), sy (y - n.k)}")],
        ["pPlus", point.p_plus],
        ["pMinus", point.p_minus],
        ["p2", point.p2],
        ["x", point.x],
        ["y", point.y],
        ["lambda", point.lam],
        ["ml2", point.ml2],
        ["Ml2", point.Ml2],
        ["mk2", point.mk2],
        ["Mk2", point.Mk2],
        ["Delta0", delta0_value(point)],
        ["Delta1 current PDF", delta1_value(point)],
        ["Delta1 old compare diagnostic", delta1_value_old_compare(point, point.lam)],
        ["raw AMFlow object", object_fields.get("RAW_COMBINATION", "sum of PP,PM,MP,MM piece expressions")],
        ["normalization factor applied after AMFlow", object_fields.get("NORMALIZATION_FACTOR", "1/(2 Pi I)^2")],
        ["real part taken", object_fields.get("REAL_PART_TAKEN", "False")],
        ["several pieces summed", object_fields.get("SUMMED_PIECES", "True")],
        ["final normalized object", object_fields.get("DOUBLE_DISCONTINUITY", "raw combination/(2 Pi I)^2")],
        ["AMFlow coefficient powers printed", ", ".join(f"c[{power}]" for power in sorted(normalized_coeffs)) or "none"],
        ["pPlus - x", one_minus_x],
        ["AMFLOW_EPS_ORDER", args.amflow_eps_order],
    ]
    print_small_table(["quantity", "value"], rows)

    print()
    print("== Raw AMFlow Laurent Coefficients ==")
    print_small_table(
        ["coefficient", "AMFlow_raw", "AMFlow_normalized"],
        [
            [f"c[{power}]", fmt_complex(raw_coeffs.get(power)), fmt_complex(normalized_coeffs.get(power))]
            for power in sorted(set(raw_coeffs) | set(normalized_coeffs))
        ],
    )


def print_amflow_method12_coefficient_comparison(point, method12_coeffs, raw_coeffs, normalized_coeffs, min_power, max_power):
    print()
    print("== Direct Method12 vs AMFlow Coefficient Comparison ==")
    rows = []
    for power in range(min_power, max_power + 1):
        method12_value = method12_coeffs.get(power)
        raw_value = raw_coeffs.get(power)
        normalized_value = normalized_coeffs.get(power)
        rows.append(
            [
                f"c[{power}]",
                fmt_complex(method12_value),
                fmt_complex(raw_value),
                fmt_complex(normalized_value),
                fmt_complex(None if normalized_value is None or method12_value is None else normalized_value - method12_value),
                fmt_complex(ratio_or_none(normalized_value, method12_value)),
            ]
        )
    print_small_table(
        [
            "coefficient",
            "method12_analytic",
            "AMFlow_raw",
            "AMFlow_normalized",
            "AMFlow_normalized - method12_analytic",
            "AMFlow_normalized / method12_analytic",
        ],
        rows,
    )

    if delta0_value(point) > 0 and delta1_value(point) > 0:
        complex_coeffs = [
            value for value in normalized_coeffs.values()
            if abs(complex(value).imag) > 1e-8 * max(1.0, abs(complex(value).real))
        ]
        if complex_coeffs:
            print("WARNING: this point has Delta0 > 0 and Delta1 > 0, but AMFlow_normalized has sizeable imaginary coefficients.")
            print("This suggests a different object, prescription, or piece combination is still being compared.")


def print_amflow_fixed_eps_reconstruction(point, eps_values, method12_values, method3_values, normalized_coeffs):
    if not normalized_coeffs:
        return
    print()
    print("== AMFlow Fixed-Epsilon Reconstruction From Laurent Coefficients ==")
    print("No finite-epsilon AMFlow rerun is used here; this evaluates the printed AMFlow Laurent coefficients.")
    rows = []
    for index, eps in enumerate(eps_values):
        amflow_from_coeffs = evaluate_laurent(normalized_coeffs, eps)
        method12 = method12_values[index]
        method3 = None if method3_values is None else method3_values[index]
        rows.append(
            [
                str(eps),
                fmt_complex(amflow_from_coeffs),
                fmt_complex(method12),
                fmt_complex(method3),
                fmt_complex(amflow_from_coeffs - method12),
                fmt_complex(None if method3 is None else amflow_from_coeffs - method3),
            ]
        )
    print_small_table(
        [
            "eps",
            "AMFlow_from_coeffs",
            "method12_closed",
            "method3_fixed_eps",
            "AMFlow_from_coeffs - method12",
            "AMFlow_from_coeffs - method3",
        ],
        rows,
    )


def print_common_factor_diagnostics(point, method12_coeffs, normalized_coeffs, eps_values, method12_values):
    print()
    print("== Common Normalization-Ratio Diagnostics ==")
    one_minus_x = point.p_plus - point.x
    factor_rows = [
        ["1 - x / pPlus convention here pPlus-x", one_minus_x],
        ["1/(pPlus - x)", 1 / one_minus_x],
        ["x", point.x],
        ["y", point.y],
        ["x*y", point.x * point.y],
        ["1/(x*y)", 1 / (point.x * point.y)],
        ["pPlus", point.p_plus],
        ["pPlus^2", point.p_plus**2],
        ["light-cone Jacobian", Fraction(1, 4)],
        ["inverse light-cone Jacobian", 4],
    ]
    print_small_table(["candidate_factor", "value"], factor_rows)

    ratio_rows = []
    for power in sorted(set(method12_coeffs) & set(normalized_coeffs)):
        ratio_rows.append(
            [
                f"c[{power}]",
                fmt_complex(ratio_or_none(method12_coeffs[power], normalized_coeffs[power])),
            ]
        )
    print()
    print("method12 / AMFlow_normalized coefficient-wise:")
    print_small_table(["coefficient", "method12 / AMFlow_normalized"], ratio_rows)

    fixed_rows = []
    for eps, method12 in zip(eps_values, method12_values):
        amflow_from_coeffs = evaluate_laurent(normalized_coeffs, eps)
        fixed_rows.append([str(eps), fmt_complex(ratio_or_none(method12, amflow_from_coeffs))])
    print()
    print("method12 / AMFlow_from_coeffs at fixed eps:")
    print_small_table(["eps", "method12 / AMFlow_from_coeffs"], fixed_rows)


def print_amflow_sign_pattern_diagnostics(combo_coeffs, min_power, max_power):
    if not combo_coeffs:
        print()
        print("== AMFlow Double-Discontinuity Sign Pattern Diagnostic ==")
        print("No signed-combination coefficients were parsed from AMFlow output.")
        return

    print()
    print("== AMFlow Double-Discontinuity Sign Pattern Diagnostic ==")
    print("Default AMFlow_normalized is still combo_all_plus. The signed combinations below are diagnostic only.")
    print("combo_delta_standard = PP - PM - MP + MM corresponds to [P-M]_x [P-M]_y if P/M are opposite i0 sides.")
    summary_rows = []
    for name, coeffs in combo_coeffs.items():
        score_values = [entry["score"] for entry in coeffs.values()]
        score = score_values[0] if score_values else None
        row = [name, fmt_float(score)]
        for power in range(min_power, max_power + 1):
            entry = coeffs.get(power)
            row.extend(
                [
                    fmt_complex(None if entry is None else entry["value"]),
                    fmt_complex(None if entry is None else entry["diff"]),
                    fmt_complex(None if entry is None else entry["ratio"]),
                ]
            )
        summary_rows.append(row)
    summary_rows.sort(key=lambda row: float("inf") if row[1] == "n/a" else float(row[1]))

    headers = ["combination", "score"]
    for power in range(min_power, max_power + 1):
        headers.extend([f"c[{power}]", f"diff c[{power}]", f"ratio c[{power}]"])
    print_small_table(headers, summary_rows)


def method12_convention_variants(point, min_power, max_power):
    def specs_for(label_prefix, label_point, swapped):
        one_minus_x = Fraction(1) - label_point.x
        pplus_minus_x = label_point.p_plus - label_point.x
        lam_default = label_point.y / pplus_minus_x
        lam_direct = label_point.y
        delta0 = delta0_value(label_point)
        base_specs = [
            ("A_current", "Gamma[eps]^2/(1-x) with lambda=y/(pPlus-x)", lam_default, 1 / one_minus_x),
            ("B_no_jacobian", "Gamma[eps]^2 with lambda=y/(pPlus-x)", lam_default, Fraction(1)),
            ("C_half_current", "1/2 * Gamma[eps]^2/(1-x) with lambda=y/(pPlus-x)", lam_default, Fraction(1, 2) / one_minus_x),
            ("D_lambda_direct", "Gamma[eps]^2/(1-x) with lambda=y", lam_direct, 1 / one_minus_x),
            ("E_lambda_direct_no_jacobian", "Gamma[eps]^2 with lambda=y", lam_direct, Fraction(1)),
            ("F_lambda_direct_half", "1/2 * Gamma[eps]^2/(1-x) with lambda=y", lam_direct, Fraction(1, 2) / one_minus_x),
            ("G_pplus_minus_x", "Gamma[eps]^2/(pPlus-x) with lambda=y/(pPlus-x)", lam_default, 1 / pplus_minus_x),
            ("H_y_to_lambda_jacobian", "y-to-lambda diagnostic: J_y_factor=1/(1-x)", lam_default, 1 / one_minus_x),
            ("H_y_density", "y-to-lambda diagnostic: J_y_factor=1", lam_default, Fraction(1)),
            ("H_half_y_to_lambda_jacobian", "y-to-lambda diagnostic: J_y_factor=1/2/(1-x)", lam_default, Fraction(1, 2) / one_minus_x),
        ]
        out = []
        for delta_convention, delta_fn in (
            ("current_pdf_delta1", delta1_value_current_pdf),
            ("old_compare_delta1", delta1_value_old_compare),
        ):
            for name, label, lam, prefactor in base_specs:
                delta1 = delta_fn(label_point, lam)
                out.append(
                    {
                        "name": f"{label_prefix}{name}_{delta_convention}",
                        "label": f"{label}; {delta_convention}; x/y {'swapped' if swapped else 'as printed'}",
                        "lambda": lam,
                        "mu_k_lambda": (1 - lam) * label_point.mk2 + lam * label_point.Mk2,
                        "Delta0": delta0,
                        "Delta1": delta1,
                        "prefactor": prefactor,
                        "branch": "-i0",
                        "delta1_convention": delta_convention,
                        "swapped_xy": swapped,
                        "coefficients": gamma_sq_times_powers_laurent(delta0, delta1, prefactor, min_power, max_power, branch="-i0"),
                    }
                )
        return out

    swapped_point = TwoLoopPoint(
        name=point.name + "_method12_xy_swapped",
        p_plus=point.p_plus,
        p_minus=point.p_minus,
        p_perp2=point.p_perp2,
        x=point.y,
        y=point.x,
        ml2=point.ml2,
        Ml2=point.Ml2,
        mk2=point.mk2,
        Mk2=point.Mk2,
    )
    variants = []
    variants.extend(specs_for("", point, False))
    variants.extend(specs_for("XYswap_", swapped_point, True))
    return variants


def coefficient_score(variant_coeffs, amflow_coeffs, powers):
    diffs = [
        variant_coeffs[power] - amflow_coeffs[power]
        for power in powers
        if power in variant_coeffs and power in amflow_coeffs
    ]
    if not diffs:
        return None
    return math.sqrt(sum(abs(diff) ** 2 for diff in diffs))


def candidate_factor_global_match_rows(method12_coeffs, amflow_coeffs):
    candidates = [
        ("x", None),
    ]
    return candidates


def print_candidate_factor_global_checks(point, current_coeffs, amflow_coeffs):
    one_minus_x = Fraction(1) - point.x
    pplus_minus_x = point.p_plus - point.x
    factors = [
        ("x", point.x),
        ("y", point.y),
        ("x*y", point.x * point.y),
        ("1/(x*y)", 1 / (point.x * point.y)),
        ("1-x", one_minus_x),
        ("1/(1-x)", 1 / one_minus_x),
        ("1/2", Fraction(1, 2)),
        ("2", Fraction(2)),
        ("1/4", Fraction(1, 4)),
        ("4", Fraction(4)),
        ("pPlus", point.p_plus),
        ("pPlus-x", pplus_minus_x),
        ("1/(pPlus-x)", 1 / pplus_minus_x),
    ]
    powers = sorted(set(current_coeffs) & set(amflow_coeffs))
    rows = []
    for name, factor in factors:
        ratios = [
            ratio_or_none(float(factor) * current_coeffs[power], amflow_coeffs[power])
            for power in powers
        ]
        finite_ratios = [ratio for ratio in ratios if ratio is not None]
        if finite_ratios:
            mean = sum(finite_ratios) / len(finite_ratios)
            spread = math.sqrt(sum(abs(ratio - mean) ** 2 for ratio in finite_ratios) / len(finite_ratios))
            score = coefficient_score(
                {power: float(factor) * current_coeffs[power] for power in current_coeffs},
                amflow_coeffs,
                powers,
            )
        else:
            mean = None
            spread = None
            score = None
        rows.append(
            [
                name,
                str(factor),
                fmt_complex(mean),
                fmt_float(spread),
                fmt_float(score),
                "yes" if spread is not None and spread < 1e-8 else "no",
            ]
        )

    print()
    print("== Candidate Global Factor Checks For Current Method12 ==")
    print("A global factor would make (factor * current Method12 coefficient)/AMFlow roughly constant across coefficients.")
    print_small_table(
        ["factor", "value", "mean ratio to AMFlow", "ratio spread", "score", "same ratio all coeffs"],
        rows,
    )


def emit_method12_convention_diagnostic(point, amflow_coeffs, min_power, max_power, fixed_eps):
    variants = method12_convention_variants(point, min_power, max_power)
    powers = [power for power in range(min_power, max_power + 1) if power in amflow_coeffs]
    print()
    print(f"######## Method12 Convention Diagnostic: {point.name} ########")
    print("Default Method12 is unchanged. This section only compares diagnostic variants.")
    print("AMFlow coefficients are reused from the already extracted Laurent object; no finite-epsilon AMFlow scan is run.")

    lam_default = point.y / (point.p_plus - point.x)
    lam_direct = point.y
    print_small_table(
        ["quantity", "value"],
        [
            ["x", point.x],
            ["y", point.y],
            ["pPlus", point.p_plus],
            ["pMinus", point.p_minus],
            ["p2", point.p2],
            ["lambda_default = y/(pPlus-x)", lam_default],
            ["lambda_direct = y", lam_direct],
            ["Delta0", delta0_value(point)],
            ["Delta1_current_pdf(lambda_default)", delta1_value_current_pdf(point, lam_default)],
            ["Delta1_old_compare(lambda_default)", delta1_value_old_compare(point, lam_default)],
            ["Delta1_current_pdf(lambda_direct)", delta1_value_current_pdf(point, lam_direct)],
            ["Delta1_old_compare(lambda_direct)", delta1_value_old_compare(point, lam_direct)],
            ["branch convention", "-i0"],
        ],
    )
    if point.p2 != 0:
        print("WARNING: off-shell point is sensitive to the Delta1 p2-term convention.")
        print("Current-PDF default uses lambda*(1-lambda*x)*p2. Old compare used lambda*(1-x)*p2.")

    detail_rows = []
    for variant in variants:
        d0 = variant["Delta0"]
        d1 = variant["Delta1"]
        detail_rows.append(
            [
                variant["name"],
                variant["label"],
                str(variant["lambda"]),
                str(variant["mu_k_lambda"]),
                str(d0),
                str(d1),
                variant["delta1_convention"],
                str(variant["prefactor"]),
                variant["branch"],
                variant["swapped_xy"],
                "positive" if d0 > 0 else "negative",
                "positive" if d1 > 0 else "negative",
            ]
        )
    print()
    print("== Variant Definitions ==")
    print_small_table(
        [
            "variant",
            "definition",
            "lambda",
            "mu_k(lambda)",
            "Delta0",
            "Delta1",
            "p2-term convention",
            "prefactor before Gamma^2",
            "branch",
            "xy swapped",
            "Delta0 sign",
            "Delta1 sign",
        ],
        detail_rows,
    )

    coefficient_headers = ["coefficient", "AMFlow_normalized"] + [variant["name"] for variant in variants]
    coefficient_rows_out = []
    for power in powers:
        coefficient_rows_out.append(
            [f"c[{power}]", fmt_complex(amflow_coeffs.get(power))]
            + [fmt_complex(variant["coefficients"].get(power)) for variant in variants]
        )
    print()
    print("== Coefficient Values ==")
    print_small_table(coefficient_headers, coefficient_rows_out)

    ratio_rows = []
    diff_rows = []
    for power in powers:
        ratio_rows.append(
            [f"c[{power}]"]
            + [fmt_complex(ratio_or_none(variant["coefficients"].get(power), amflow_coeffs.get(power))) for variant in variants]
        )
        diff_rows.append(
            [f"c[{power}]"]
            + [fmt_complex(variant["coefficients"].get(power) - amflow_coeffs[power]) for variant in variants]
        )
    print()
    print("== Variant / AMFlow_normalized Ratios ==")
    print_small_table(["coefficient"] + [variant["name"] for variant in variants], ratio_rows)
    print()
    print("== Variant - AMFlow_normalized Differences ==")
    print_small_table(["coefficient"] + [variant["name"] for variant in variants], diff_rows)

    ranking = []
    for variant in variants:
        score = coefficient_score(variant["coefficients"], amflow_coeffs, powers)
        ranking.append((float("inf") if score is None else score, variant))
    ranking.sort(key=lambda item: item[0])
    print()
    print("== Variant Ranking By Coefficient Distance ==")
    print("score = sqrt(sum_n |variant_c[n] - AMFlow_c[n]|^2) over available coefficients.")
    print_small_table(
        ["rank", "variant", "score", "definition"],
        [
            [index + 1, variant["name"], fmt_float(score), variant["label"]]
            for index, (score, variant) in enumerate(ranking)
        ],
    )

    amflow_from_coeffs = evaluate_laurent(amflow_coeffs, fixed_eps)
    fixed_rows = []
    for variant in variants:
        value = gamma_sq_variant_value(
            variant["Delta0"],
            variant["Delta1"],
            variant["prefactor"],
            fixed_eps,
            branch=variant["branch"],
        )
        fixed_rows.append(
            [
                variant["name"],
                str(fixed_eps),
                fmt_complex(amflow_from_coeffs),
                fmt_complex(value),
                fmt_complex(value - amflow_from_coeffs),
                fmt_complex(ratio_or_none(value, amflow_from_coeffs)),
            ]
        )
    print()
    print("== Fixed-Epsilon Reconstruction Against AMFlow Coefficients ==")
    print("AMFlow_from_coeffs is reconstructed from Laurent coefficients, not rerun at finite eps.")
    print_small_table(
        ["variant", "eps", "AMFlow_from_coeffs", "variant_value", "variant - AMFlow", "variant / AMFlow"],
        fixed_rows,
    )

    current = variants[0]["coefficients"]
    print_candidate_factor_global_checks(point, current, amflow_coeffs)
    print()
    print("Diagnostic conclusion rule: do not declare a variant correct automatically. A pure global factor would give a stable ratio across c[-2], c[-1], c[0].")


def print_method12_convention_diagnostic(point, amflow_coeffs, min_power, max_power, fixed_eps):
    output = io.StringIO()
    with contextlib.redirect_stdout(output):
        emit_method12_convention_diagnostic(point, amflow_coeffs, min_power, max_power, fixed_eps)
    text = output.getvalue()
    print(text, end="")
    output_dir = repo_root() / "amflow-project" / "targets" / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / f"two_loop_method12_convention_diagnostic_{point.name}.txt"
    path.write_text(text, encoding="utf-8")
    print(f"Method12 convention diagnostic report written to: {path}")


def constant_ratio_fit(method3_values, method12_values):
    y_vec = np.array(method3_values, dtype=np.complex128)
    x_vec = np.array(method12_values, dtype=np.complex128)
    denom = np.vdot(x_vec, x_vec)
    if denom == 0:
        return None
    return np.vdot(x_vec, y_vec) / denom


def ratio_stats(ratios):
    if not ratios:
        return None, None, None
    arr = np.array(ratios, dtype=np.complex128)
    mean = complex(np.mean(arr))
    std = float(np.sqrt(np.mean(np.abs(arr - mean) ** 2)))
    cov = None if mean == 0 else std / abs(mean)
    return mean, std, cov


def print_ratio_vs_eps_summary(point, eps_values, kmax_values, cutoff_values, cutoff_errors, method12_values, max_rel_error):
    print()
    print("== Ratio vs eps Summary By Kmax ==")
    print("Question: is method3_normalized_cutoff / method12_closed approximately independent of eps?")
    headers = ["Kmax"] + [f"ratio(eps={eps})" for eps in eps_values] + ["mean_ratio", "std_ratio", "coefficient_of_variation"]
    rows = []
    for kmax in kmax_values:
        ratios = [ratio_or_none(value, method12) for value, method12 in zip(cutoff_values[kmax], method12_values)]
        mean, std, cov = ratio_stats([ratio for ratio in ratios if ratio is not None])
        rows.append(
            [str(kmax)]
            + [fmt_complex(ratio) for ratio in ratios]
            + [fmt_complex(mean), fmt_float(std), fmt_float(cov)]
        )
    print_small_table(headers, rows)

    print()
    print("== Ratio Fit Inclusion Mask ==")
    rows = []
    for kmax in kmax_values:
        for eps, value, error in zip(eps_values, cutoff_values[kmax], cutoff_errors[kmax]):
            rel = relative_error(error, value)
            included = rel is not None and rel <= max_rel_error
            if rel is None and (error is None or (isinstance(error, float) and math.isnan(error))):
                included = True
            rows.append([str(kmax), str(eps), fmt_float(rel), included])
    print_small_table(["Kmax", "eps", "relative_error", f"included <= {max_rel_error}"], rows)


def print_constant_ratio_fits(eps_values, kmax_values, cutoff_values, cutoff_errors, method12_values, max_rel_error):
    print()
    print("== Constant-Ratio Fits method3_cutoff(eps) ~= C(Kmax) method12(eps) ==")
    print("Interpretation: small residuals after fitting C(Kmax) mean the mismatch is mostly an eps-independent factor at that Kmax.")
    summary_rows = []
    residual_blocks = []
    for kmax in kmax_values:
        included_indices = []
        excluded = []
        for i, (value, error) in enumerate(zip(cutoff_values[kmax], cutoff_errors[kmax])):
            rel = relative_error(error, value)
            include = rel is None or rel <= max_rel_error
            if include:
                included_indices.append(i)
            else:
                excluded.append((eps_values[i], rel))

        if not included_indices:
            summary_rows.append([str(kmax), "not fit", "n/a", "n/a", "n/a", "all eps excluded"])
            continue

        y_vals = [cutoff_values[kmax][i] for i in included_indices]
        x_vals = [method12_values[i] for i in included_indices]
        C = constant_ratio_fit(y_vals, x_vals)
        if C is None:
            summary_rows.append([str(kmax), "not fit", "n/a", "n/a", "n/a", "zero method12 norm"])
            continue

        residuals = [y - C * x for y, x in zip(y_vals, x_vals)]
        residual_norm = float(np.sqrt(np.sum(np.abs(np.array(residuals)) ** 2)))
        rel_residuals = [abs(res) / abs(x) if x != 0 else float("nan") for res, x in zip(residuals, x_vals)]
        max_rel_residual = float(np.nanmax(rel_residuals)) if rel_residuals else float("nan")
        chi_like = float(np.sqrt(np.mean(np.array(rel_residuals) ** 2))) if rel_residuals else float("nan")
        excluded_text = "none" if not excluded else ", ".join(f"{eps} relerr={fmt_float(rel)}" for eps, rel in excluded)
        summary_rows.append([str(kmax), fmt_complex(C), fmt_float(residual_norm), fmt_float(max_rel_residual), fmt_float(chi_like), excluded_text])

        block_rows = []
        for eps, y, x in zip(eps_values, cutoff_values[kmax], method12_values):
            ratio = ratio_or_none(y, x)
            residual = y - C * x
            rel_after = abs(residual) / abs(x) if x != 0 else None
            block_rows.append([str(eps), fmt_complex(ratio), fmt_complex(None if ratio is None else ratio - C), fmt_float(rel_after)])
        residual_blocks.append((kmax, block_rows))

    print_small_table(["Kmax", "C(Kmax)", "residual_norm", "max_relative_residual", "chi_like", "excluded_eps"], summary_rows)
    for kmax, rows in residual_blocks:
        print()
        print(f"== Constant-Ratio Residuals By eps: Kmax={kmax} ==")
        print_small_table(["eps", "ratio_i", "ratio_i - C(Kmax)", "relative_residual_after_C"], rows)


def print_kmax_ratio_trends(eps_values, kmax_values, cutoff_values, cutoff_errors, method12_values, max_rel_error):
    print()
    print("== Ratio Trend In Kmax By eps ==")
    for eps_index, eps in enumerate(eps_values):
        rows = []
        previous_value = None
        previous_error = None
        last_reliable = None
        decreases = 0
        low_significance = False
        for kmax in kmax_values:
            value = cutoff_values[kmax][eps_index]
            error = cutoff_errors[kmax][eps_index]
            method12 = method12_values[eps_index]
            ratio = ratio_or_none(value, method12)
            rel = relative_error(error, value)
            if rel is None or rel <= max_rel_error:
                last_reliable = kmax
            shell_increment = None if previous_value is None else value - previous_value
            shell_error = None
            shell_sigma = None
            if previous_error is not None:
                try:
                    err_i = float(error)
                    err_prev = float(previous_error)
                    if not math.isnan(err_i) and not math.isnan(err_prev):
                        shell_error = math.sqrt(err_i**2 + err_prev**2)
                        if shell_error != 0 and shell_increment is not None:
                            shell_sigma = abs(shell_increment) / shell_error
                            low_significance = low_significance or shell_sigma < 2
                except Exception:
                    pass
            if previous_value is not None and value.real < previous_value.real:
                decreases += 1
            rows.append([str(kmax), fmt_complex(ratio), fmt_float(rel), fmt_complex(shell_increment), fmt_float(shell_sigma)])
            previous_value = value
            previous_error = error

        print()
        print(f"eps = {eps}")
        print_small_table(["Kmax", "ratio", "relative_error", "shell_increment", "|shell|/shell_error"], rows)
        monotonicity = "monotone nondecreasing in real part" if decreases == 0 else f"{decreases} real-part decrease(s)"
        print(f"monotonicity check: {monotonicity}")
        print(f"last reliable Kmax with relative_error <= {max_rel_error}: {last_reliable if last_reliable is not None else 'none'}")
        if low_significance:
            print("WARNING: at least one shell increment has significance < 2 sigma; high-Kmax trend may be noise dominated.")


def tail_model_matrix(model_name, eps, k_values):
    eps_f = float(eps)
    k = np.array([float(value) for value in k_values], dtype=np.float64)
    base = k ** (-2.0 * eps_f)
    if model_name == "model1_K^-2eps":
        # I_K = I_inf - a K^(-2 eps)
        return np.column_stack([np.ones_like(k), -base]), ["I_inf", "a"]
    if model_name == "model2_K^-2eps_K^-4eps":
        # I_K = I_inf - a K^(-2 eps) - b K^(-4 eps)
        return np.column_stack([np.ones_like(k), -base, -(base**2)]), ["I_inf", "a", "b"]
    if model_name == "model3_K^-2eps_logK":
        # I_K = I_inf - K^(-2 eps) (a + b log K)
        return np.column_stack([np.ones_like(k), -base, -(base * np.log(k))]), ["I_inf", "a", "b"]
    raise ValueError(f"Unknown tail model {model_name}")


def weighted_complex_lstsq(matrix, values, errors, force_unweighted=False):
    A = np.array(matrix, dtype=np.complex128)
    y = np.array(values, dtype=np.complex128)
    used_weighted = False
    if not force_unweighted and errors is not None:
        try:
            sigma = np.array([float(error) for error in errors], dtype=np.float64)
            if np.all(np.isfinite(sigma)) and np.all(sigma > 0):
                weights = 1.0 / sigma
                A = A * weights[:, None]
                y = y * weights
                used_weighted = True
        except Exception:
            used_weighted = False
    params, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
    condition = float(np.linalg.cond(A)) if A.size else float("nan")
    return params, condition, used_weighted


def select_tail_fit_indices(kmax_values, cutoff_values_for_eps, cutoff_errors_for_eps, args):
    max_rel = args.max_relative_error_for_ratio_fit
    eligible = []
    excluded = []
    for i, kmax in enumerate(kmax_values):
        rel = relative_error(cutoff_errors_for_eps[i], cutoff_values_for_eps[i])
        if rel is None or rel <= max_rel:
            eligible.append(i)
        else:
            excluded.append((kmax, f"relative_error={fmt_float(rel)}"))

    kmin = None if args.tail_fit_Kmin is None else parse_mpf(args.tail_fit_Kmin)
    kmax_limit = None if args.tail_fit_Kmax is None else parse_mpf(args.tail_fit_Kmax)
    if args.tail_fit_last_n is not None:
        if kmin is not None or kmax_limit is not None:
            print("Tail fit window: --tail-fit-last-n provided, so Kmin/Kmax bounds are ignored.")
        eligible = eligible[-args.tail_fit_last_n :]
        mode = f"last_n={args.tail_fit_last_n}"
    else:
        filtered = []
        for i in eligible:
            keep = True
            if kmin is not None and kmax_values[i] < kmin:
                keep = False
            if kmax_limit is not None and kmax_values[i] > kmax_limit:
                keep = False
            if keep:
                filtered.append(i)
            else:
                excluded.append((kmax_values[i], "outside tail fit K window"))
        eligible = filtered
        mode = "relative-error filtered"
        if kmin is not None or kmax_limit is not None:
            mode += f", Kmin={kmin if kmin is not None else '-inf'}, Kmax={kmax_limit if kmax_limit is not None else '+inf'}"
    return eligible, excluded, mode


def fit_tail_models_for_eps(eps, k_values, values, errors, method12, force_unweighted):
    models = ["model1_K^-2eps", "model2_K^-2eps_K^-4eps", "model3_K^-2eps_logK"]
    results = []
    for model in models:
        matrix, names = tail_model_matrix(model, eps, k_values)
        params, condition, weighted = weighted_complex_lstsq(matrix, values, errors, force_unweighted=force_unweighted)
        fitted = matrix @ params
        residuals = np.array(values, dtype=np.complex128) - fitted
        residual_norm = float(np.sqrt(np.sum(np.abs(residuals) ** 2)))
        rel_residuals = [
            abs(residual) / abs(value) if value != 0 else float("nan")
            for residual, value in zip(residuals, values)
        ]
        max_relative_residual = float(np.nanmax(rel_residuals)) if rel_residuals else float("nan")
        note_parts = []
        if len(k_values) < 3:
            note_parts.append("WARNING: too few K values for reliable tail extrapolation.")
        if len(k_values) < len(names):
            note_parts.append("WARNING: underdetermined tail fit; parameters are not uniquely constrained.")
        if len(k_values) == len(names):
            note_parts.append("WARNING: exactly determined tail fit; residuals are not a convergence test.")
        i_inf = complex(params[0])
        largest_value = values[-1]
        if largest_value.real > 0 and i_inf.real < largest_value.real:
            note_parts.append("WARNING: extrapolated I_inf is below largest cutoff value; fit likely unstable.")
        results.append(
            {
                "model": model,
                "param_names": names,
                "params": params,
                "I_inf": i_inf,
                "method12": method12,
                "ratio": ratio_or_none(i_inf, method12),
                "residual_norm": residual_norm,
                "max_relative_residual": max_relative_residual,
                "condition": condition,
                "weighted": weighted,
                "n_points": len(k_values),
                "fit_quality_note": " ".join(note_parts) if note_parts else "ok",
            }
        )
    return results


def print_tail_fit_self_test():
    print()
    print("== Tail Fit Self Test ==")
    eps = Fraction(1, 10)
    k_values = [Fraction(100), Fraction(200), Fraction(500), Fraction(1000), Fraction(2000)]
    true_i_inf = 171.0 + 0.25j
    a = 300.0 - 0.5j
    b = -80.0 + 0.75j
    values = [
        true_i_inf - a * float(k) ** (-2.0 * float(eps)) - b * float(k) ** (-4.0 * float(eps))
        for k in k_values
    ]
    results = fit_tail_models_for_eps(eps, k_values, values, [1.0 for _ in k_values], true_i_inf, force_unweighted=True)
    rows = [
        [result["model"], fmt_complex(result["I_inf"]), fmt_complex(result["I_inf"] - true_i_inf), fmt_float(result["condition"])]
        for result in results
    ]
    print_small_table(["model", "recovered_I_inf", "recovered - true", "condition_number"], rows)


def print_tail_extrapolation(eps_values, kmax_values, cutoff_values, cutoff_errors, method12_values, args):
    print()
    print("== Fixed-epsilon Tail Extrapolation ==")
    print("Tail extrapolation is a diagnostic. It does not replace a proof of the asymptotic form.")
    print("Its purpose is to test whether the finite-K cutoff deficit is consistent with a slowly decaying UV-regulated tail.")
    if args.tail_fit_self_test:
        print_tail_fit_self_test()

    all_results = {model: [] for model in ["model1_K^-2eps", "model2_K^-2eps_K^-4eps", "model3_K^-2eps_logK"]}
    compact_rows = []
    for eps_index, eps in enumerate(eps_values):
        values_for_eps = [cutoff_values[kmax][eps_index] for kmax in kmax_values]
        errors_for_eps = [cutoff_errors[kmax][eps_index] for kmax in kmax_values]
        indices, excluded, mode = select_tail_fit_indices(kmax_values, values_for_eps, errors_for_eps, args)
        used_k = [kmax_values[i] for i in indices]
        used_values = [values_for_eps[i] for i in indices]
        used_errors = [errors_for_eps[i] for i in indices]
        method12 = method12_values[eps_index]

        print()
        print(f"eps = {eps}; fit window mode: {mode}; K values used = {', '.join(str(k) for k in used_k) if used_k else 'none'}")
        if excluded:
            print("excluded K values: " + ", ".join(f"{k} ({reason})" for k, reason in excluded))
        if len(used_k) < 2:
            print("WARNING: too few K values to fit tail models.")
            continue

        results = fit_tail_models_for_eps(
            eps=eps,
            k_values=used_k,
            values=used_values,
            errors=used_errors,
            method12=method12,
            force_unweighted=args.tail_fit_unweighted,
        )
        rows = []
        compact = {"eps": str(eps), "method12": fmt_complex(method12)}
        ratios_for_model_dependence = []
        for result in results:
            all_results[result["model"]].append(result["ratio"])
            ratios_for_model_dependence.append(result["ratio"])
            params_text = ", ".join(
                f"{name}={fmt_complex(value)}"
                for name, value in zip(result["param_names"], result["params"])
                if name != "I_inf"
            )
            rows.append(
                [
                    str(eps),
                    result["model"],
                    ", ".join(str(k) for k in used_k),
                    fmt_complex(result["I_inf"]),
                    fmt_complex(method12),
                    fmt_complex(result["ratio"]),
                    params_text,
                    fmt_float(result["residual_norm"]),
                    fmt_float(result["max_relative_residual"]),
                    fmt_float(result["condition"]),
                    result["n_points"],
                    "weighted" if result["weighted"] else "unweighted",
                    result["fit_quality_note"],
                ]
            )
            compact[f"I_inf_{result['model']}"] = fmt_complex(result["I_inf"])
            compact[f"ratio_{result['model']}"] = fmt_complex(result["ratio"])
        print_small_table(
            [
                "eps",
                "model",
                "K values used",
                "I_inf",
                "method12_closed",
                "I_inf/method12",
                "tail params",
                "residual_norm",
                "max_relative_residual",
                "condition_number",
                "n_points",
                "fit",
                "fit_quality_note",
            ],
            rows,
        )
        mean_ratio, std_ratio, cov = ratio_stats([ratio for ratio in ratios_for_model_dependence if ratio is not None])
        if cov is not None and cov > 0.25:
            print("WARNING: tail model dependence is large.")
        compact_rows.append(compact)

    if compact_rows:
        print()
        print("== Tail Extrapolation Compact Summary ==")
        headers = [
            "eps",
            "method12_closed",
            "I_inf_model1",
            "ratio_model1",
            "I_inf_model2",
            "ratio_model2",
            "I_inf_model3",
            "ratio_model3",
        ]
        rows = []
        for row in compact_rows:
            rows.append(
                [
                    row["eps"],
                    row["method12"],
                    row.get("I_inf_model1_K^-2eps", "not fit"),
                    row.get("ratio_model1_K^-2eps", "not fit"),
                    row.get("I_inf_model2_K^-2eps_K^-4eps", "not fit"),
                    row.get("ratio_model2_K^-2eps_K^-4eps", "not fit"),
                    row.get("I_inf_model3_K^-2eps_logK", "not fit"),
                    row.get("ratio_model3_K^-2eps_logK", "not fit"),
                ]
            )
        print_small_table(headers, rows)

    print()
    print("== Tail Extrapolated Ratio Stability Over eps ==")
    rows = []
    for model, ratios in all_results.items():
        mean, std, cov = ratio_stats([ratio for ratio in ratios if ratio is not None])
        rows.append([model, fmt_complex(mean), fmt_float(std), fmt_float(cov)])
    print_small_table(["model", "mean_ratio_over_eps", "std_ratio_over_eps", "coefficient_of_variation"], rows)
    print("Interpretation: an eps-independent ratio after Kmax extrapolation suggests a remaining overall convention factor.")
    print("If the ratio remains strongly eps-dependent, Method 3 likely still has a formula or implementation issue beyond the finite-K tail.")


def radial_cutoff_scan_for_point(point_name, eps_values, eta, kmax_values, args):
    point = get_point(point_name)
    method3_point = method3_point_for_convention(point, args.method3_label_convention)
    unsupported_reason = method3_unsupported_reason(method3_point)
    lc_jacobian = parse_float(args.lc_jacobian)
    min_power = args.fit_min_power
    max_power = args.fit_max_power
    coefficient_count = max_power - min_power + 1

    print()
    print(f"######## Radial Cutoff Accumulation Scan: {point.name} ########")
    print("finite-box map: k = Kmax u, l = Kmax v, theta = Pi w")
    print("Jacobian: Kmax^2 Pi")
    print("KmaxList = {", ", ".join(str(kmax) for kmax in kmax_values), "}", sep="")
    print("epsList = {", ", ".join(str(eps) for eps in eps_values), "}", sep="")
    print(f"N = {args.N}, scrambles = {args.scrambles}, eta = {eta}, lcJacobian = {args.lc_jacobian}")
    print_method3_label_convention(point, method3_point, args.method3_label_convention)
    print("Default infinite-map method3 is not changed by this diagnostic.")
    if args.Kmax_shells:
        print("Kmax-shells requested: reporting cumulative shell increments I(K_i)-I(K_{i-1}) from the finite-box data.")

    if unsupported_reason is not None:
        print(f"Radial cutoff scan skipped: {unsupported_reason}")
        return {"point_name": point_name, "cutoff_values": {}}

    if args.debug_normalization:
        first_eps = eps_values[0]
        print_method12_variable_audit(point, first_eps, closed_form_value_branch(point, first_eps, branch="-i0"))
        print_amflow_delta_convention_audit()
        print()
        print("== Method 3 Normalization Definitions For Cutoff Scan ==")
        print("method3_stripped_cutoff: raw finite-box MC integral over dk dl dtheta, including Kmax^2 Pi, excluding the global Method-3 prefactor.")
        print("method3_normalized_cutoff: method3_prefactor * method3_stripped_cutoff.")
        print_small_table(["component", "value"], method3_prefactor_component_rows(method3_point, first_eps, lc_jacobian))

    cutoff_values = {kmax: [] for kmax in kmax_values}
    cutoff_errors = {kmax: [] for kmax in kmax_values}
    method12_values = []

    for eps in eps_values:
        method12 = closed_form_value_branch(point, eps, branch="-i0")
        method12_values.append(method12)
        prefactor = method3_prefactor(method3_point, eps, lc_jacobian=lc_jacobian)
        rows = []
        previous_normalized = None
        previous_error = None
        warnings = []
        print()
        print(f"== eps = {eps} ==")
        for kmax in kmax_values:
            stripped, stripped_err = method3_stripped_cutoff(
                point=method3_point,
                eps=eps,
                eta=eta,
                n_samples=args.N,
                cutoff=kmax,
                scrambles=args.scrambles,
                seed=args.seed,
            )
            normalized = prefactor * stripped
            normalized_err = abs(prefactor) * stripped_err
            cutoff_values[kmax].append(normalized)
            cutoff_errors[kmax].append(normalized_err)
            rel_err = relative_error(normalized_err, normalized)
            shell_increment = None if previous_normalized is None else normalized - previous_normalized
            shell_increment_error = None
            shell_sigma = None
            if previous_error is not None and normalized_err is not None:
                try:
                    err_i = float(normalized_err)
                    err_prev = float(previous_error)
                    if not math.isnan(err_i) and not math.isnan(err_prev):
                        shell_increment_error = math.sqrt(err_i**2 + err_prev**2)
                        if shell_increment_error != 0 and shell_increment is not None:
                            shell_sigma = abs(shell_increment) / shell_increment_error
                except Exception:
                    shell_increment_error = None
            if (
                previous_normalized is not None
                and normalized.real < previous_normalized.real
                and shell_increment_error is not None
                and abs((normalized - previous_normalized).real) <= 2.0 * shell_increment_error
            ):
                warnings.append(
                    "WARNING: high-Kmax cumulative scan is noise dominated; do not infer saturation."
                )
            rows.append(
                [
                    str(kmax),
                    fmt_complex(stripped),
                    fmt_complex(normalized),
                    fmt_float(normalized_err),
                    fmt_float(rel_err),
                    fmt_complex(shell_increment),
                    fmt_float(shell_increment_error),
                    fmt_float(shell_sigma),
                    fmt_complex(method12),
                    fmt_complex(ratio_or_none(normalized, method12)),
                ]
            )
            if args.debug_normalization:
                print_normalization_variant_audit(
                    point,
                    normalized,
                    method12,
                    lc_jacobian,
                    title=f"== Normalization Variant Audit: eps={eps}, Kmax={kmax} ==",
                )
            previous_normalized = normalized
            previous_error = normalized_err

        print_small_table(
            [
                "Kmax",
                "method3_stripped_cutoff",
                "method3_normalized_cutoff",
                "MC error",
                "relative_error",
                "shell_increment",
                "shell_error",
                "|shell|/shell_error",
                "method12_closed",
                "cutoff/method12",
            ],
            rows,
        )
        if warnings:
            print()
            for warning in sorted(set(warnings)):
                print(warning)

    if len(eps_values) > 1:
        print_ratio_vs_eps_summary(
            point=point,
            eps_values=eps_values,
            kmax_values=kmax_values,
            cutoff_values=cutoff_values,
            cutoff_errors=cutoff_errors,
            method12_values=method12_values,
            max_rel_error=args.max_relative_error_for_ratio_fit,
        )
        if args.ratio_fit:
            print_constant_ratio_fits(
                eps_values=eps_values,
                kmax_values=kmax_values,
                cutoff_values=cutoff_values,
                cutoff_errors=cutoff_errors,
                method12_values=method12_values,
                max_rel_error=args.max_relative_error_for_ratio_fit,
            )

    print_kmax_ratio_trends(
        eps_values=eps_values,
        kmax_values=kmax_values,
        cutoff_values=cutoff_values,
        cutoff_errors=cutoff_errors,
        method12_values=method12_values,
        max_rel_error=args.max_relative_error_for_ratio_fit,
    )

    if args.tail_extrapolate:
        print_tail_extrapolation(
            eps_values=eps_values,
            kmax_values=kmax_values,
            cutoff_values=cutoff_values,
            cutoff_errors=cutoff_errors,
            method12_values=method12_values,
            args=args,
        )

    print()
    print("WARNING: At finite Kmax the cutoff integral is finite at eps=0, so Laurent coefficients extracted at fixed Kmax are not physical.")
    print("Use ratio-vs-eps and Kmax extrapolation diagnostics instead.")

    if len(eps_values) >= coefficient_count:
        method12_coeffs = method12_laurent_coefficients(point, min_power, max_power, branch="-i0")
        print()
        print("== Non-Physical Finite-Kmax Laurent Diagnostic ==")
        print(f"fit powers = {min_power}..{max_power}")
        print("This table is only a numerical diagnostic of finite-cutoff samples; it is not the Laurent coefficient of the infinite integral.")
        rows = []
        for kmax in kmax_values:
            fit = fit_laurent_coefficients(eps_values, cutoff_values[kmax], min_power, max_power)
            rows.append(
                [
                    str(kmax),
                    fmt_complex(fit[-2]),
                    fmt_complex(method12_coeffs[-2]),
                    fmt_complex(ratio_or_none(fit[-2], method12_coeffs[-2])),
                ]
            )
        print_small_table(
            ["Kmax", "c[-2]_method3_cutoff", "c[-2]_method12", "ratio"],
            rows,
        )
    else:
        print()
        print(
            "Laurent fit skipped: need at least "
            f"{coefficient_count} eps values for powers {min_power}..{max_power}; got {len(eps_values)}."
        )

    return {
        "point_name": point_name,
        "cutoff_values": cutoff_values,
        "cutoff_errors": cutoff_errors,
        "method12_values": method12_values,
    }


def run_radial_cutoff_scan(point_names, eps_values, eta, args):
    output_dir = repo_root() / "amflow-project" / "targets" / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    kmax_values = parse_number_list(args.Kmax_list, DEFAULT_KMAX_LIST)
    point_label = "all_points" if len(point_names) > 1 else point_names[0]
    eps_label = "_".join(safe_name(eps) for eps in eps_values)
    report_path = output_dir / f"two_loop_method3_radial_cutoff_scan_{point_label}_eps_{eps_label}.txt"

    with report_path.open("w", encoding="utf-8") as report:
        tee = Tee(sys.stdout, report)
        with contextlib.redirect_stdout(tee):
            print("Radial-cutoff-scan examples:")
            print("  ./run.sh two-loop-method3-mc-compare --point equal_mass_offshell_positive --radial-cutoff-scan --eps 0.1 --amflow skip")
            print("  ./run.sh two-loop-method3-mc-compare --point equal_mass_offshell_positive --radial-cutoff-scan --eps-list 0.2,0.15,0.1,0.075,0.05 --amflow skip")
            print("  ./run.sh two-loop-method3-mc-compare --point equal_mass_offshell_positive --radial-cutoff-scan --eps 0.1 --Kmax-list 500,1000,2000 --amflow skip --debug-normalization")
            print("  ./run.sh two-loop-method3-mc-compare --point equal_mass_offshell_positive --radial-cutoff-scan --eps-list 0.2,0.15,0.1,0.075,0.05 --Kmax-list 100,200,500,1000,2000 --amflow skip --ratio-fit")
            print("  ./run.sh two-loop-method3-mc-compare --point equal_mass_offshell_positive --radial-cutoff-scan --tail-extrapolate --eps-list 0.2,0.15,0.1,0.075,0.05 --Kmax-list 100,200,500,1000,2000 --amflow skip")
            results = [
                radial_cutoff_scan_for_point(point_name, eps_values, eta, kmax_values, args)
                for point_name in point_names
            ]

    print(f"Radial cutoff scan report written to: {report_path}")
    return results


def log_fit_degree_for_order(order, args):
    if args.log_fit_degree is not None:
        return args.log_fit_degree
    return order + 2


def fit_polynomial_in_log(k_values, values, degree):
    logs = np.array([math.log(float(k)) for k in k_values], dtype=np.float64)
    matrix = np.column_stack([logs**power for power in range(degree, -1, -1)])
    coeffs, _, _, _ = np.linalg.lstsq(matrix.astype(np.complex128), np.array(values, dtype=np.complex128), rcond=None)
    fitted = matrix @ coeffs
    residuals = np.array(values, dtype=np.complex128) - fitted
    condition = float(np.linalg.cond(matrix)) if matrix.size else float("nan")
    residual_norm = float(np.sqrt(np.sum(np.abs(residuals) ** 2)))
    max_relative_residual = float(
        np.nanmax([
            abs(residual) / abs(value) if value != 0 else float("nan")
            for residual, value in zip(residuals, values)
        ])
    )
    return coeffs, fitted, residuals, condition, residual_norm, max_relative_residual


def print_cutoff_log_dimreg_mapping_note(args):
    print()
    print("== Cutoff-log to DimReg Mapping Status ==")
    if args.assume_cutoff_dimreg_map:
        print("Mapping enabled by --assume-cutoff-dimreg-map.")
        print("No rigorous cutoff-log to dimensional-pole mapping is implemented yet, so no final equality claim is made.")
        print("The log-fit coefficients are printed as diagnostic input for the mapping derivation.")
    else:
        print("No direct equality between cutoff-log coefficients and method12 Laurent coefficients is assumed.")
        print("Use --assume-cutoff-dimreg-map only after the mapping from large cutoff logs to dimensional poles is derived and implemented.")


def run_method3_cutoff_log_for_point(point_name, eta, args):
    point = get_point(point_name)
    method3_point = method3_point_for_convention(point, args.method3_label_convention)
    max_order = args.method3_coeff_order
    eps_step = parse_float(args.method3_coeff_eps_step)
    lc_jacobian = parse_float(args.lc_jacobian)
    kmax_values = parse_number_list(args.Kmax_list, ["5", "10", "20", "50", "100", "200", "500", "1000"])
    log_fit_kmin = parse_mpf(args.log_fit_Kmin)
    reconstruction_eps = parse_mpf(args.reconstruction_eps) if args.reconstruction_eps is not None else parse_mpf(args.eps)
    unsupported_reason = method3_unsupported_reason(method3_point)

    print()
    print(f"######## Method-3 Cutoff-Log Coefficient Diagnostic: {point.name} ########")
    print("This is the expand-first finite-cutoff mode.")
    print("Theta remains numerical. No UV subtraction is introduced.")
    print("The epsilon-expanded coefficient integrals are evaluated only at finite Kmax boxes.")
    print("WARNING: J_r(Kmax) are cutoff-dependent. They are not Laurent coefficients by themselves.")
    print("WARNING: epsilon-expanded Method-3 coefficient integrals are not integrated over k,l in [0,infinity).")
    print(f"KmaxList = {{{', '.join(str(k) for k in kmax_values)}}}")
    print(f"coefficient orders r=0..{max_order}; eps Taylor step = {args.method3_coeff_eps_step}")
    print(f"N = {args.N}, scrambles = {args.scrambles}, eta = {eta}, lcJacobian = {args.lc_jacobian}")
    print(f"log fit uses Kmax >= {log_fit_kmin}")
    print_method3_label_convention(point, method3_point, args.method3_label_convention)

    method12_coeffs = print_method12_coefficient_debug(point, -2, max(2, max_order))
    print_cutoff_log_dimreg_mapping_note(args)

    if args.amflow == "fresh":
        print()
        print("== AMFlow Laurent Coefficients ==")
        print("AMFlow direct Laurent coefficient extraction is not wired into this Python mode yet.")
        print("No finite-eps AMFlow sampling is run here, to avoid reintroducing the old eps-fit workflow.")

    if unsupported_reason is not None:
        print()
        print(f"Method 3 cutoff-log mode skipped: {unsupported_reason}")
        return {"point_name": point_name, "skipped": True}

    coeffs_by_k = {}
    errors_by_k = {}
    print()
    print("== Cutoff-dependent coefficient integrals J_r(Kmax) ==")
    rows = []
    for kmax in kmax_values:
        coeffs, errors = method3_cutoff_taylor_coefficients(
            point=method3_point,
            max_order=max_order,
            eta=eta,
            n_samples=args.N,
            cutoff=kmax,
            scrambles=args.scrambles,
            seed=args.seed,
            lc_jacobian=lc_jacobian,
            eps_step=eps_step,
        )
        coeffs_by_k[kmax] = coeffs
        errors_by_k[kmax] = errors
        row = [str(kmax)]
        for order in range(max_order + 1):
            row.append(fmt_complex(coeffs[order]))
            row.append(fmt_float(errors[order]))
        rows.append(row)
    headers = ["Kmax"]
    for order in range(max_order + 1):
        headers.extend([f"J_{order}(Kmax)", "MC error"])
    print_small_table(headers, rows)

    fit_k_values = [k for k in kmax_values if k >= log_fit_kmin]
    print()
    print("== Large-K Polynomial Fits In L = Log[Kmax] ==")
    print("The default degree is r+2 for J_r, matching the expected growth pattern J_0~L^2, J_1~L^3, J_2~L^4.")
    if len(fit_k_values) < 2:
        print("WARNING: fewer than two Kmax values in the log-fit window.")
    log_fit_results = {}
    for order in range(max_order + 1):
        if not fit_k_values:
            print(f"J_{order}(Kmax) fit skipped: no Kmax values in the log-fit window.")
            continue
        degree = log_fit_degree_for_order(order, args)
        values = [coeffs_by_k[k][order] for k in fit_k_values]
        coeffs, fitted, residuals, condition, residual_norm, max_rel = fit_polynomial_in_log(fit_k_values, values, degree)
        log_fit_results[order] = {
            "degree": degree,
            "coeffs": coeffs,
            "residuals": residuals,
            "condition": condition,
            "residual_norm": residual_norm,
            "max_relative_residual": max_rel,
        }
        print()
        print(f"J_{order}(Kmax) fit, degree {degree}")
        if len(fit_k_values) <= degree + 1:
            print("WARNING: log fit is exactly determined or underdetermined; residuals are not a convergence test.")
        print_small_table(
            ["basis", "coefficient"],
            [[f"L^{power}", fmt_complex(value)] for power, value in zip(range(degree, -1, -1), coeffs)],
        )
        print_small_table(
            ["Kmax", "J_r(Kmax)", "fit_value", "residual"],
            [
                [str(k), fmt_complex(value), fmt_complex(fit), fmt_complex(residual)]
                for k, value, fit, residual in zip(fit_k_values, values, fitted, residuals)
            ],
        )
        print_small_table(
            ["quantity", "value"],
            [
                ["condition_number", fmt_float(condition)],
                ["residual_norm", fmt_float(residual_norm)],
                ["max_relative_residual", fmt_float(max_rel)],
            ],
        )

    print()
    print("== Reconstruction Check At Finite Cutoff ==")
    print(f"eps = {reconstruction_eps}")
    print("Compares direct finite-eps cutoff integral to Sum_r eps^r J_r(Kmax).")
    rec_rows = []
    pref_recon = method3_prefactor(method3_point, reconstruction_eps, lc_jacobian=lc_jacobian)
    for kmax in kmax_values:
        stripped, stripped_err = method3_stripped_cutoff(
            point=method3_point,
            eps=reconstruction_eps,
            eta=eta,
            n_samples=args.N,
            cutoff=kmax,
            scrambles=args.scrambles,
            seed=args.seed,
        )
        direct = pref_recon * stripped
        direct_err = abs(pref_recon) * stripped_err
        eps_f = float(reconstruction_eps)
        reconstructed = sum((eps_f**order) * coeffs_by_k[kmax][order] for order in range(max_order + 1))
        rec_rows.append([
            str(kmax),
            fmt_complex(direct),
            fmt_float(direct_err),
            fmt_complex(reconstructed),
            fmt_complex(reconstructed - direct),
            fmt_complex(ratio_or_none(reconstructed, direct)),
        ])
    print_small_table(
        ["Kmax", "direct I3(eps;Kmax)", "direct error", "Sum eps^r J_r", "recon - direct", "recon/direct"],
        rec_rows,
    )

    print()
    print("== Method12 Laurent Coefficients For Reference ==")
    print("These are analytic dimensional-regularization coefficients. No direct equality to J_r(Kmax) is assumed.")
    print_small_table(
        ["coefficient", "method12_analytic"],
        [[f"c[{power}]", fmt_complex(value)] for power, value in method12_coeffs.items()],
    )

    return {
        "point_name": point_name,
        "coeffs_by_k": coeffs_by_k,
        "errors_by_k": errors_by_k,
        "log_fit_results": log_fit_results,
        "method12_coeffs": method12_coeffs,
    }


def run_method3_cutoff_log_mode(point_names, eta, args):
    output_dir = repo_root() / "amflow-project" / "targets" / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    report_name = (
        "two_loop_method3_cutoff_log_coefficients_all_points.txt"
        if len(point_names) > 1
        else f"two_loop_method3_cutoff_log_coefficients_{point_names[0]}.txt"
    )
    report_path = output_dir / report_name

    with report_path.open("w", encoding="utf-8") as report:
        tee = Tee(sys.stdout, report)
        with contextlib.redirect_stdout(tee):
            print("Cutoff-log coefficient mode examples:")
            print("  ./run.sh two-loop-method3-mc-compare --point equal_mass_offshell_positive --method3-coeff-cutoff-log --amflow skip")
            print("  ./run.sh two-loop-method3-mc-compare --all-points --method3-coeff-cutoff-log --amflow skip")
            results = [run_method3_cutoff_log_for_point(point_name, eta, args) for point_name in point_names]

    print(f"Cutoff-log coefficient report written to: {report_path}")
    return results


def run_method3_integrand_expansion_for_point(point_name, eta, args):
    point = get_point(point_name)
    method3_point = method3_point_for_convention(point, args.method3_label_convention)
    max_order = args.method3_coeff_order
    lc_jacobian = parse_float(args.lc_jacobian)
    kmax_values = parse_number_list(args.Kmax_list, ["5", "10", "20", "50", "100", "200", "500", "1000"])
    reconstruction_eps = parse_mpf(args.reconstruction_eps) if args.reconstruction_eps is not None else parse_mpf(args.eps)
    unsupported_reason = method3_unsupported_reason(method3_point)

    print()
    print(f"######## Method-3 Integrand-Level Epsilon Expansion: {point.name} ########")
    print("This mode analytically expands the finite-cutoff integrand in epsilon.")
    print("Local measure expansion: k^(1-2 eps) l^(1-2 eps) sin(theta)^(-2 eps) = k*l*exp[-2 eps L].")
    print("L = log(k) + log(l) + log(sin(theta)); coefficient r is k*l*(-2 L)^r/r!.")
    print("The Method-3 prefactor is Taylor-expanded numerically and convolved with the local integrand coefficients.")
    print("WARNING: these J_r(Kmax) are cutoff-dependent and UV divergent as Kmax -> infinity. They are not Laurent coefficients.")
    print("Theta remains numerical. No UV subtraction is introduced.")
    print(f"KmaxList = {{{', '.join(str(k) for k in kmax_values)}}}")
    print(f"coefficient orders r=0..{max_order}; N={args.N}; scrambles={args.scrambles}; eta={eta}; lcJacobian={args.lc_jacobian}")
    print_method3_label_convention(point, method3_point, args.method3_label_convention)

    method12_coeffs = print_method12_coefficient_debug(point, -2, max(2, max_order))
    if unsupported_reason is not None:
        print()
        print(f"Method 3 integrand expansion skipped: {unsupported_reason}")
        return {"point_name": point_name, "skipped": True}

    coeffs_by_k = {}
    errors_by_k = {}
    rows = []
    print()
    print("== Integrand-Level Coefficients J_r_integrand(Kmax) ==")
    for kmax in kmax_values:
        coeffs, errors = method3_cutoff_integrand_expansion_coefficients(
            point=method3_point,
            max_order=max_order,
            eta=eta,
            n_samples=args.N,
            cutoff=kmax,
            scrambles=args.scrambles,
            seed=args.seed,
            lc_jacobian=lc_jacobian,
            prefactor_eps_step=parse_float(args.method3_coeff_eps_step),
        )
        coeffs_by_k[kmax] = coeffs
        errors_by_k[kmax] = errors
        row = [str(kmax)]
        for order in range(max_order + 1):
            row.append(fmt_complex(coeffs[order]))
            row.append(fmt_float(errors[order]))
        rows.append(row)

    headers = ["Kmax"]
    for order in range(max_order + 1):
        headers.extend([f"J_{order}_integrand(Kmax)", "MC error"])
    print_small_table(headers, rows)

    print()
    print("== Reconstruction Check At Finite Cutoff ==")
    print(f"eps = {reconstruction_eps}")
    pref_recon = method3_prefactor(method3_point, reconstruction_eps, lc_jacobian=lc_jacobian)
    rec_rows = []
    for kmax in kmax_values:
        stripped, stripped_err = method3_stripped_cutoff(
            point=method3_point,
            eps=reconstruction_eps,
            eta=eta,
            n_samples=args.N,
            cutoff=kmax,
            scrambles=args.scrambles,
            seed=args.seed,
        )
        direct = pref_recon * stripped
        direct_err = abs(pref_recon) * stripped_err
        eps_f = float(reconstruction_eps)
        reconstructed = sum((eps_f**order) * coeffs_by_k[kmax][order] for order in range(max_order + 1))
        rec_rows.append(
            [
                str(kmax),
                fmt_complex(direct),
                fmt_float(direct_err),
                fmt_complex(reconstructed),
                fmt_complex(reconstructed - direct),
                fmt_complex(ratio_or_none(reconstructed, direct)),
            ]
        )
    print_small_table(
        ["Kmax", "direct I3(eps;Kmax)", "direct error", "Sum eps^r J_r", "recon - direct", "recon/direct"],
        rec_rows,
    )

    print()
    print("== Method12 Current-PDF Laurent Coefficients For Reference ==")
    print("No direct equality to finite-cutoff J_r(Kmax) is assumed.")
    print_small_table(
        ["coefficient", "method12_current_pdf"],
        [[f"c[{power}]", fmt_complex(value)] for power, value in method12_coeffs.items()],
    )
    return {"point_name": point_name, "coeffs_by_k": coeffs_by_k, "errors_by_k": errors_by_k}


def run_method3_integrand_expansion_mode(point_names, eta, args):
    output_dir = repo_root() / "amflow-project" / "targets" / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    report_name = (
        "two_loop_transverse_integrand_expansion_all_points.txt"
        if len(point_names) > 1
        else f"two_loop_transverse_integrand_expansion_{point_names[0]}.txt"
    )
    report_path = output_dir / report_name
    with report_path.open("w", encoding="utf-8") as report:
        tee = Tee(sys.stdout, report)
        with contextlib.redirect_stdout(tee):
            print("Integrand-level expansion mode examples:")
            print("  ./run.sh two-loop-method3-mc-compare --point equal_mass_offshell_positive --method3-integrand-eps-expansion --amflow skip")
            results = [run_method3_integrand_expansion_for_point(point_name, eta, args) for point_name in point_names]
    print(f"Method3 integrand expansion report written to: {report_path}")
    return results


def omega_float(m):
    half_m = float(m) / 2.0
    return 2.0 * math.pi**half_m / math.gamma(half_m)


def require_mpmath_for_l_residue():
    try:
        import mpmath  # noqa: F401
    except Exception as exc:
        raise SystemExit(
            "ERROR: --method1-l-residue-direct-check requires mpmath.\n"
            "Install it in the Python environment used by run.sh, for example:\n"
            "  python3 -m pip install mpmath"
        ) from exc


def l_residue_ingredients(point, eps, eta, contour_sign, short_C_diagnostic=False):
    eps_f = float(eps)
    x = float(point.x)
    y = float(point.y)
    p2 = float(point.p2)
    ml2 = float(point.ml2)
    Ml2 = float(point.Ml2)
    mk2 = float(point.mk2)
    Mk2 = float(point.Mk2)
    lam = y / (float(point.p_plus) - x)
    kappa = lam * (1.0 - lam)
    mu_k = (1.0 - lam) * mk2 + lam * Mk2
    delta0 = float(delta0_value(point))
    a1_short = mu_k - kappa * (1.0 - x) * p2 + kappa * (1.0 - x) / x * ml2
    extra_offshell_term = -lam * (1.0 - lam * x) * p2
    a1_full = a1_short + extra_offshell_term
    a1 = a1_short if short_C_diagnostic else a1_full
    b1 = kappa / x
    d = 4.0 - 2.0 * eps_f
    n = 2.0 - 2.0 * eps_f
    omega_n = omega_float(n)
    prefactor = (
        float(contour_sign)
        * math.gamma(eps_f)
        / (1.0 - x)
        * (1.0 / (1j * math.pi ** (d / 2.0)))
        * 0.5
        * (2.0 * math.pi * 1j)
        * omega_n
        / 2.0
    )
    return {
        "eps": eps,
        "eta": eta,
        "x": point.x,
        "y": point.y,
        "pPlus": point.p_plus,
        "pMinus": point.p_minus,
        "p2": point.p2,
        "lambda": lam,
        "kappa": kappa,
        "mu_k_lambda": mu_k,
        "Delta0": delta0,
        "A1": a1,
        "A1_full": a1_full,
        "A1_short_without_extra_term": a1_short,
        "extra_offshell_term": extra_offshell_term,
        "C_full_includes_extra_offshell_term": not short_C_diagnostic,
        "B1": b1,
        "d": d,
        "n": n,
        "Omega_n": omega_n,
        "contour_sign": contour_sign,
        "prefactor": prefactor,
    }


def _adaptive_simpson_complex(f, a, b, tol=1e-9, max_depth=24):
    def simpson(fa, fm, fb, a, b):
        return (b - a) * (fa + 4.0 * fm + fb) / 6.0

    fa = f(a)
    fb = f(b)
    m = (a + b) / 2.0
    fm = f(m)
    whole = simpson(fa, fm, fb, a, b)

    def recurse(a, b, fa, fm, fb, whole, depth):
        m = (a + b) / 2.0
        lm = (a + m) / 2.0
        rm = (m + b) / 2.0
        flm = f(lm)
        frm = f(rm)
        left = simpson(fa, flm, fm, a, m)
        right = simpson(fm, frm, fb, m, b)
        delta = left + right - whole
        if depth <= 0 or abs(delta) < 15.0 * tol:
            return left + right + delta / 15.0
        return (
            recurse(a, m, fa, flm, fm, left, depth - 1)
            + recurse(m, b, fm, frm, fb, right, depth - 1)
        )

    return recurse(a, b, fa, fm, fb, whole, max_depth)


def _integrate_complex_piecewise(f, points, tol=1e-9):
    total = 0j
    for a, b in zip(points[:-1], points[1:]):
        if b <= a:
            continue
        total += _adaptive_simpson_complex(f, float(a), float(b), tol=tol / max(1, len(points) - 1))
    return total


def l_residue_integrand_float(t, ingredients):
    eps_f = float(ingredients["eps"])
    eta_f = float(ingredients["eta"])
    a1 = ingredients["A1"]
    b1 = ingredients["B1"]
    delta0 = ingredients["Delta0"]
    return (
        complex(t, 0.0) ** (-eps_f)
        * complex(a1 + b1 * t, -eta_f) ** (-eps_f)
        / complex(t + delta0, -eta_f)
    )


def l_residue_integral_cutoff_float(ingredients, tmax):
    eps_f = float(ingredients["eps"])
    scale = max(float(tmax), 1.0)
    near_zero = min(1e-10, scale * 1e-12)

    def f(t):
        return l_residue_integrand_float(t, ingredients)

    points = [near_zero]
    for frac in (1e-6, 1e-4, 1e-2, 0.1, 0.3, 0.6, 1.0):
        value = float(tmax) * frac
        if near_zero < value <= float(tmax):
            points.append(value)
    points = sorted(set(points))
    integral = _integrate_complex_piecewise(f, points, tol=1e-8)
    # Analytic small-t estimate of the endpoint omitted by near_zero.
    endpoint_scale = complex(ingredients["A1"], -float(ingredients["eta"])) ** (-eps_f) / complex(ingredients["Delta0"], -float(ingredients["eta"]))
    integral += endpoint_scale * near_zero ** (1.0 - eps_f) / (1.0 - eps_f)
    return ingredients["prefactor"] * integral


def l_residue_integral_infinite_float(ingredients):
    z_end = 1.0 - 1e-12

    def fz(z):
        t = z / (1.0 - z)
        jac = 1.0 / (1.0 - z) ** 2
        return jac * l_residue_integrand_float(t, ingredients)

    z0 = 1e-12
    z_points = [z0, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 0.1, 0.5, 0.9, 0.99, 0.999, 0.9999, 0.99999, 0.999999, z_end]
    integral = _integrate_complex_piecewise(fz, z_points, tol=1e-8)
    eps_f = float(ingredients["eps"])
    endpoint_scale = complex(ingredients["A1"], -float(ingredients["eta"])) ** (-eps_f) / complex(ingredients["Delta0"], -float(ingredients["eta"]))
    integral += endpoint_scale * z0 ** (1.0 - eps_f) / (1.0 - eps_f)
    return ingredients["prefactor"] * integral


def l_residue_integrals_mpmath(point, eps, eta, contour_sign, dps=50, short_C_diagnostic=False):
    try:
        import mpmath as mp  # type: ignore
    except Exception as exc:
        raise RuntimeError(
            "The direct l-residue diagnostic requires mpmath. "
            "Install it in the Python environment used by run.sh, for example: "
            "python3 -m pip install mpmath"
        ) from exc

    mp.mp.dps = int(dps)
    eps_mp = mp.mpf(str(float(eps)))
    eta_mp = mp.mpf(str(float(eta)))
    x = mp.mpf(point.x.numerator) / point.x.denominator
    y = mp.mpf(point.y.numerator) / point.y.denominator
    p_plus = mp.mpf(point.p_plus.numerator) / point.p_plus.denominator
    p2 = mp.mpf(point.p2.numerator) / point.p2.denominator
    ml2 = mp.mpf(point.ml2.numerator) / point.ml2.denominator
    mk2 = mp.mpf(point.mk2.numerator) / point.mk2.denominator
    Mk2 = mp.mpf(point.Mk2.numerator) / point.Mk2.denominator
    lam = y / (p_plus - x)
    kappa = lam * (1 - lam)
    mu_k = (1 - lam) * mk2 + lam * Mk2
    delta0 = mp.mpf(delta0_value(point).numerator) / delta0_value(point).denominator
    a1_short = mu_k - kappa * (1 - x) * p2 + kappa * (1 - x) / x * ml2
    extra_offshell_term = -lam * (1 - lam * x) * p2
    a1 = a1_short if short_C_diagnostic else a1_short + extra_offshell_term
    b1 = kappa / x
    d = 4 - 2 * eps_mp
    n = 2 - 2 * eps_mp
    omega_n = 2 * mp.pi ** (n / 2) / mp.gamma(n / 2)
    pref = (
        contour_sign
        * mp.gamma(eps_mp)
        / (1 - x)
        * (1 / (1j * mp.pi ** (d / 2)))
        * mp.mpf("0.5")
        * (2 * mp.pi * 1j)
        * omega_n
        / 2
    )

    def integrand_t(t):
        return t ** (-eps_mp) * (a1 + b1 * t - 1j * eta_mp) ** (-eps_mp) / (t + delta0 - 1j * eta_mp)

    z0 = mp.mpf("1e-30")
    z_end = 1 - mp.mpf("1e-20")
    z_splits = [z0, mp.mpf("1e-10"), mp.mpf("1e-8"), mp.mpf("1e-6"), mp.mpf("1e-4"), mp.mpf("1e-2"), mp.mpf("0.1"), mp.mpf("0.5"), mp.mpf("0.9"), mp.mpf("0.99"), mp.mpf("0.999"), mp.mpf("0.9999"), mp.mpf("0.99999"), mp.mpf("0.999999"), z_end]
    if a1 < 0 and b1 > 0:
        t_root = -a1 / b1
        z_root = t_root / (1 + t_root)
        if z0 < z_root < z_end:
            z_splits.extend([z_root * (1 - mp.mpf("1e-8")), z_root, z_root * (1 + mp.mpf("1e-8"))])
            z_splits = sorted(set(z_splits))

    def integrand_z(z):
        t = z / (1 - z)
        return integrand_t(t) / (1 - z) ** 2

    # Avoid evaluating the mapped integrand exactly at z=0 or z=1.  The
    # omitted pieces are tiny endpoint estimates for t -> 0 and t -> infinity.
    endpoint_scale = (a1 - 1j * eta_mp) ** (-eps_mp) / (delta0 - 1j * eta_mp)
    small_t_endpoint = endpoint_scale * z0 ** (1 - eps_mp) / (1 - eps_mp)
    t_end = z_end / (1 - z_end)
    large_t_tail = b1 ** (-eps_mp) * t_end ** (-2 * eps_mp) / (2 * eps_mp)
    infinite = pref * (mp.quad(integrand_z, z_splits) + small_t_endpoint + large_t_tail)

    def cutoff_value(tmax):
        tmax = mp.mpf(str(tmax))
        t0 = min(mp.mpf("1e-30"), tmax * mp.mpf("1e-20"))
        splits = [t0, tmax * mp.mpf("1e-8"), tmax * mp.mpf("1e-6"), tmax * mp.mpf("1e-4"), tmax * mp.mpf("1e-2"), tmax * mp.mpf("0.1"), tmax * mp.mpf("0.3"), tmax * mp.mpf("0.6"), tmax]
        if a1 < 0 and b1 > 0:
            t_root = -a1 / b1
            if t0 < t_root < tmax:
                splits.extend([t_root * (1 - mp.mpf("1e-8")), t_root, t_root * (1 + mp.mpf("1e-8"))])
        splits = sorted(set(splits))
        endpoint = endpoint_scale * t0 ** (1 - eps_mp) / (1 - eps_mp)
        return pref * (mp.quad(integrand_t, splits) + endpoint)

    return infinite, cutoff_value


def print_l_residue_prefactor(ingredients):
    print_small_table(
        ["prefactor_piece", "value"],
        [
            ["contour_sign", ingredients["contour_sign"]],
            ["Gamma[eps]/(1-x)", f"{math.gamma(float(ingredients['eps'])) / (1.0 - float(ingredients['x'])):.16e}"],
            ["1/(i*pi^(d/2))", fmt_complex(1.0 / (1j * math.pi ** (ingredients["d"] / 2.0)))],
            ["1/2", "1/2"],
            ["2*pi*i", fmt_complex(2.0 * math.pi * 1j)],
            ["Omega[2-2 eps]/2", f"{ingredients['Omega_n'] / 2.0:.16e}"],
            ["N_l(eps)", fmt_complex(ingredients["prefactor"])],
        ],
    )


def run_method1_l_residue_direct_for_point(point_name, eps, eta, args):
    point = get_point(point_name)
    print()
    print(f"######## Method-1 Direct l^- Residue Check: {point.name} ########")
    print("Diagnostic only. Method12, Method3, and AMFlow defaults are unchanged.")
    print("Reduced representation: integrate t = l_perp^2 after taking the D1 l^- residue.")
    print("The code evaluates both contour_sign=+1 and contour_sign=-1.")
    if point.p_plus != 1 or point.p_perp2 != 0:
        print("This diagnostic currently implements pPlus=1 and pPerp2=0 only.")
        return {"point_name": point_name, "skipped": True}

    method12 = closed_form_value_branch(point, eps, branch="-i0")
    print_small_table(
        ["quantity", "value"],
        [
            ["eps", eps],
            ["eta", eta],
            ["Delta0", delta0_value(point)],
            ["Delta1 current PDF", delta1_value(point)],
            ["lambda", point.lam],
            ["kappa", float(point.lam * (1 - point.lam))],
            ["method12_closed current PDF", fmt_complex(method12)],
        ],
    )
    full_ing = l_residue_ingredients(point, eps, eta, 1, short_C_diagnostic=False)
    short_ing = l_residue_ingredients(point, eps, eta, 1, short_C_diagnostic=True)
    lam_exact = point.lam
    kappa_exact = lam_exact * (1 - lam_exact)
    mu_k_exact = (1 - lam_exact) * point.mk2 + lam_exact * point.Mk2
    extra_exact = -lam_exact * (1 - lam_exact * point.x) * point.p2
    a1_short_exact = mu_k_exact - kappa_exact * (1 - point.x) * point.p2 + kappa_exact * (1 - point.x) / point.x * point.ml2
    a1_full_exact = a1_short_exact + extra_exact
    print()
    print("== Current-PDF C(l) Off-Shell Term Diagnostic ==")
    print_small_table(
        ["quantity", "value"],
        [
            ["C_full_includes_extra_offshell_term", full_ing["C_full_includes_extra_offshell_term"]],
            ["extra_offshell_term exact", extra_exact],
            ["extra_offshell_term = -lambda*(1-lambda*x)*p2", f"{full_ing['extra_offshell_term']:.16e}"],
            ["A1_full exact", a1_full_exact],
            ["A1_full", f"{full_ing['A1_full']:.16e}"],
            ["A1_short_without_extra_term exact", a1_short_exact],
            ["A1_short_without_extra_term", f"{full_ing['A1_short_without_extra_term']:.16e}"],
            ["A1_full - A1_short exact", a1_full_exact - a1_short_exact],
            ["A1_full - A1_short_without_extra_term", f"{full_ing['A1_full'] - full_ing['A1_short_without_extra_term']:.16e}"],
        ],
    )

    amflow_from_coeffs = None
    if args.amflow in ("fresh", "reuse"):
        print()
        print("== AMFlow Coefficient Reuse For l-Residue Diagnostic ==")
        ensure_amflow_order_for_powers(args, 0)
        raw, normalized, _fields, _combos, output = run_amflow_coefficients(point_name, eps, args)
        if normalized:
            amflow_from_coeffs = evaluate_laurent(normalized, eps)
            print(f"AMFlow_from_coeffs({eps}) = {fmt_complex(amflow_from_coeffs)}")
        else:
            print("AMFlow coefficients unavailable. Last AMFlow output follows:")
            print(output[-3000:])

    tmax_values = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]
    results = {}
    for sign in (1, -1):
        ingredients = l_residue_ingredients(point, eps, eta, sign, short_C_diagnostic=False)
        results[sign] = {"ingredients": ingredients}
        print()
        print(f"== Ingredients And Prefactor, contour_sign={sign} ==")
        print_small_table(
            ["ingredient", "value"],
            [
                ["lambda", f"{ingredients['lambda']:.16e}"],
                ["kappa", f"{ingredients['kappa']:.16e}"],
                ["mu_k(lambda)", f"{ingredients['mu_k_lambda']:.16e}"],
                ["C_full_includes_extra_offshell_term", ingredients["C_full_includes_extra_offshell_term"]],
                ["extra_offshell_term", f"{ingredients['extra_offshell_term']:.16e}"],
                ["A1_full", f"{ingredients['A1_full']:.16e}"],
                ["A1_short_without_extra_term", f"{ingredients['A1_short_without_extra_term']:.16e}"],
                ["A1 used", f"{ingredients['A1']:.16e}"],
                ["B1", f"{ingredients['B1']:.16e}"],
                ["d", f"{ingredients['d']:.16e}"],
                ["Omega[2-2 eps]", f"{ingredients['Omega_n']:.16e}"],
            ],
        )
        print_l_residue_prefactor(ingredients)

        infinite_mp, cutoff_mp = l_residue_integrals_mpmath(point, eps, eta, sign, dps=args.dps, short_C_diagnostic=False)
        print(f"mpmath integration backend active with dps={args.dps}.")
        infinite = complex(infinite_mp)
        cutoff_function = lambda tmax: complex(cutoff_mp(tmax))
        results[sign]["infinite"] = infinite
        short_diagnostic = None
        full_closed = None
        short_closed = None
        if sign == 1:
            short_mp, _short_cutoff = l_residue_integrals_mpmath(point, eps, eta, sign, dps=args.dps, short_C_diagnostic=True)
            short_diagnostic = complex(short_mp)
            full_closed = method1_l_residue_closed_value(point, eps, eta, short_C_diagnostic=False)
            short_closed = method1_l_residue_closed_value(point, eps, eta, short_C_diagnostic=True)

        comparison_rows = [
            ["l_residue_full_C_numeric infinite-map", fmt_complex(infinite)],
            ["l_residue_full_C_numeric - Method12", fmt_complex(infinite - method12)],
            ["l_residue_full_C_numeric / Method12", fmt_complex(ratio_or_none(infinite, method12))],
        ]
        if sign == 1:
            comparison_rows.extend(
                [
                    ["l_residue_full_C_closed", fmt_complex(full_closed)],
                    ["l_residue_full_C_closed - full_C_numeric", fmt_complex(full_closed - infinite)],
                    ["l_residue_full_C_closed / full_C_numeric", fmt_complex(ratio_or_none(full_closed, infinite))],
                    ["l_residue_short_C_diagnostic numeric", fmt_complex(short_diagnostic)],
                    ["l_residue_short_C_diagnostic closed", fmt_complex(short_closed)],
                ]
            )
        if amflow_from_coeffs is not None:
            comparison_rows.extend(
                [
                    ["l_residue_full_C_numeric - AMFlow_from_coeffs", fmt_complex(infinite - amflow_from_coeffs)],
                    ["l_residue_full_C_numeric / AMFlow_from_coeffs", fmt_complex(ratio_or_none(infinite, amflow_from_coeffs))],
                    ["l_residue_short_C_diagnostic - AMFlow_from_coeffs", fmt_complex(None if short_diagnostic is None else short_diagnostic - amflow_from_coeffs)],
                    ["l_residue_short_C_diagnostic / AMFlow_from_coeffs", fmt_complex(None if short_diagnostic is None else ratio_or_none(short_diagnostic, amflow_from_coeffs))],
                ]
            )
        print()
        print(f"== Infinite-Domain Result, contour_sign={sign} ==")
        print_small_table(["quantity", "value"], comparison_rows)

        cutoff_rows = []
        for tmax in tmax_values:
            value = cutoff_function(tmax)
            cutoff_rows.append(
                [
                    str(tmax),
                    fmt_complex(value),
                    fmt_complex(value - method12),
                    fmt_complex(ratio_or_none(value, method12)),
                ]
            )
        print()
        print(f"== Cutoff Scan, contour_sign={sign} ==")
        print_small_table(
            ["Tmax", "I_l_residue(Tmax)", "I_l_residue(Tmax) - Method12", "I_l_residue(Tmax)/Method12"],
            cutoff_rows,
        )

    return {"point_name": point_name, "method12": method12, "results": results}


def run_method1_l_residue_direct_mode(point_names, eps, eta, args):
    output_dir = repo_root() / "amflow-project" / "targets" / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    report_name = (
        "two_loop_residue_1d_numeric_all_points.txt"
        if len(point_names) > 1
        else f"two_loop_residue_1d_numeric_{point_names[0]}.txt"
    )
    report_path = output_dir / report_name
    with report_path.open("w", encoding="utf-8") as report:
        tee = Tee(sys.stdout, report)
        with contextlib.redirect_stdout(tee):
            print("Direct l-residue diagnostic example:")
            print("  ./run.sh two-loop-method3-mc-compare --point equal_mass_offshell_positive --method1-l-residue-direct-check --eps 0.1 --amflow reuse")
            results = [run_method1_l_residue_direct_for_point(point_name, eps, eta, args) for point_name in point_names]
    print(f"Method1 l-residue direct report written to: {report_path}")
    return results


def mp_to_complex(value):
    return complex(float(mp.re(value)), float(mp.im(value)))


def method1_l_residue_closed_value(point, eps, eta, short_C_diagnostic=False):
    mp.mp.dps = 80
    return mp_to_complex(
        kernel_method1_l_residue_closed(
            point,
            mp.mpf(str(float(eps))),
            eta=mp.mpf(str(float(eta))),
            short_C_diagnostic=short_C_diagnostic,
        )
    )


def method1_l_residue_closed_coefficients(point, min_power=-2, max_power=2, short_C_diagnostic=False):
    mp.mp.dps = 100
    degree = max_power + 10
    if min_power < -2:
        raise ValueError("l-residue closed coefficients are implemented from c[-2] upward.")
    poly_powers = list(range(0, degree + 1))
    eps_nodes = [
        mp.mpf(1) / n
        for n in (30, 35, 42, 50, 60, 72, 86, 103, 124, 150, 180, 220, 270, 330, 410, 520, 660, 850, 1100, 1450, 1900, 2500)
    ][: degree + 6]
    matrix = mp.matrix([[eps ** power for power in poly_powers] for eps in eps_nodes])
    vector = mp.matrix([
        eps**2 * kernel_method1_l_residue_closed(
            point,
            eps,
            eta=mp.mpf("1e-40"),
            short_C_diagnostic=short_C_diagnostic,
        )
        for eps in eps_nodes
    ])
    if len(eps_nodes) == len(poly_powers):
        solved = mp.lu_solve(matrix, vector)
    else:
        solved = mp.lu_solve(matrix.T * matrix, matrix.T * vector)
    return {power: mp_to_complex(solved[power + 2]) for power in range(min_power, max_power + 1)}


def print_l_residue_closed_ingredients(point, eta, short_C_diagnostic=False):
    ingredients = kernel_l_residue_closed_ingredients(point, short_C_diagnostic=short_C_diagnostic)
    rows = [
        ["X", mp.nstr(ingredients["X"], 30)],
        ["lambda", mp.nstr(ingredients["lambda"], 30)],
        ["kappa", mp.nstr(ingredients["kappa"], 30)],
        ["Delta0", mp.nstr(ingredients["Delta0"], 30)],
        ["C_full_includes_extra_offshell_term", ingredients["C_full_includes_extra_offshell_term"]],
        ["extra_offshell_term", mp.nstr(ingredients["extra_offshell_term"], 30)],
        ["A1_full", mp.nstr(ingredients["A1_full"], 30)],
        ["A1_short_without_extra_term", mp.nstr(ingredients["A1_short_without_extra_term"], 30)],
        ["A1_full - A1_short_without_extra_term", mp.nstr(ingredients["A1_full"] - ingredients["A1_short_without_extra_term"], 30)],
        ["A1 used", mp.nstr(ingredients["A1"], 30)],
        ["B1", mp.nstr(ingredients["B1"], 30)],
        ["z = 1 - B1*Delta0/A1", mp.nstr(ingredients["z"], 30)],
        ["eta", eta],
        ["branch convention", "-i0 in (Delta0-i eta)(A1-i eta)"],
    ]
    print_small_table(["ingredient", "value"], rows)


def score_coefficients(lhs, rhs, powers):
    total = 0.0
    for power in powers:
        if power in lhs and power in rhs:
            total += abs(lhs[power] - rhs[power]) ** 2
    return math.sqrt(total)


def print_l_residue_closed_coefficient_table(point, lres_coeffs, method12_coeffs, amflow_coeffs, min_power, max_power, label="l_residue_full_C_closed"):
    print()
    print(f"== {label} Laurent Coefficients ==")
    expected_leading = 1.0 / (2.0 * float(1 - point.x))
    leading = lres_coeffs.get(-2)
    leading_ok = leading is not None and abs(leading - expected_leading) < 1e-7
    print(f"Expected c[-2] = 1/(2*(1-x)) = {expected_leading:.12e}: {'PASS' if leading_ok else 'FAIL'}")
    rows = []
    for power in range(min_power, max_power + 1):
        amf = amflow_coeffs.get(power) if amflow_coeffs else None
        lres = lres_coeffs.get(power)
        m12 = method12_coeffs.get(power)
        rows.append(
            [
                f"c[{power}]",
                fmt_complex(amf),
                fmt_complex(m12),
                fmt_complex(lres),
                fmt_complex(None if amf is None or lres is None else lres - amf),
                fmt_complex(None if amf is None or lres is None else ratio_or_none(lres, amf)),
                fmt_complex(None if amf is None or m12 is None else m12 - amf),
                fmt_complex(None if amf is None or m12 is None else ratio_or_none(m12, amf)),
            ]
        )
    print_small_table(
        [
            "coefficient",
            "AMFlow_normalized",
            "Method12_current_pdf",
            label,
            f"{label} - AMFlow",
            f"{label} / AMFlow",
            "Method12 - AMFlow",
            "Method12 / AMFlow",
        ],
        rows,
    )


def run_method1_l_residue_closed_for_point(point_name, eps, eta, args):
    point = get_point(point_name)
    min_power = args.fit_min_power
    max_power = args.fit_max_power
    print()
    print(f"######## Method-1 l-Residue Closed Candidate Check: {point.name} ########")
    print("This is a diagnostic candidate expression. Method12 defaults are unchanged.")
    if point.p_plus != 1 or point.p_perp2 != 0:
        print("l-residue closed form currently implemented for pPlus=1, pPerp2=0.")
        return {"point_name": point_name, "skipped": True}

    print()
    print("== l-Residue Closed Ingredients: full current-PDF C(l) ==")
    print_l_residue_closed_ingredients(point, eta, short_C_diagnostic=False)
    print()
    print("== l-Residue Closed Ingredients: short-C diagnostic only ==")
    print_l_residue_closed_ingredients(point, eta, short_C_diagnostic=True)

    method12_value = closed_form_value_branch(point, eps, branch="-i0")
    lres_value = method1_l_residue_closed_value(point, eps, eta, short_C_diagnostic=False)
    lres_short_value = method1_l_residue_closed_value(point, eps, eta, short_C_diagnostic=True)
    method12_coeffs = method12_laurent_coefficients(point, min_power, max_power, branch="-i0")
    lres_coeffs = method1_l_residue_closed_coefficients(point, min_power, max_power, short_C_diagnostic=False)
    lres_short_coeffs = method1_l_residue_closed_coefficients(point, min_power, max_power, short_C_diagnostic=True)

    amflow_coeffs = {}
    amflow_from_coeffs = None
    if args.amflow in ("fresh", "reuse", "imported"):
        print()
        print("== AMFlow Coefficients For l-Residue Closed Candidate ==")
        ensure_amflow_order_for_powers(args, max_power)
        _raw, amflow_coeffs, _fields, _combos, output = run_amflow_coefficients(point_name, eps, args)
        if amflow_coeffs:
            amflow_from_coeffs = evaluate_laurent(amflow_coeffs, eps)
            print(f"AMFlow_from_coeffs({eps}) = {fmt_complex(amflow_from_coeffs)}")
        else:
            print("AMFlow coefficients unavailable. Last AMFlow output follows:")
            print(output[-3000:])

    numeric_value = None
    numeric_short_value = None
    if args.method1_l_residue_direct_check:
        print()
        print("== Numerical Direct l-Residue Cross-Check, contour_sign=+1 ==")
        infinite_mp, _cutoff_mp = l_residue_integrals_mpmath(point, eps, eta, 1, dps=args.dps, short_C_diagnostic=False)
        infinite_short_mp, _cutoff_short_mp = l_residue_integrals_mpmath(point, eps, eta, 1, dps=args.dps, short_C_diagnostic=True)
        numeric_value = complex(infinite_mp)
        numeric_short_value = complex(infinite_short_mp)
        print(f"l_residue_full_C_numeric infinite-map = {fmt_complex(numeric_value)}")
        print(f"l_residue_short_C_diagnostic numeric infinite-map = {fmt_complex(numeric_short_value)}")

    comparison_rows = [
        ["eps", eps],
        ["Method12_closed current PDF", fmt_complex(method12_value)],
        ["l_residue_full_C_closed", fmt_complex(lres_value)],
        ["l_residue_full_C_closed - Method12", fmt_complex(lres_value - method12_value)],
        ["l_residue_full_C_closed / Method12", fmt_complex(ratio_or_none(lres_value, method12_value))],
        ["l_residue_full_C_numeric", fmt_complex(numeric_value)],
        ["l_residue_full_C_closed - numeric", fmt_complex(None if numeric_value is None else lres_value - numeric_value)],
        ["l_residue_full_C_closed / numeric", fmt_complex(None if numeric_value is None else ratio_or_none(lres_value, numeric_value))],
        ["l_residue_short_C_diagnostic closed", fmt_complex(lres_short_value)],
        ["l_residue_short_C_diagnostic numeric", fmt_complex(numeric_short_value)],
        ["AMFlow_from_coeffs", fmt_complex(amflow_from_coeffs)],
        ["l_residue_full_C_closed - AMFlow", fmt_complex(None if amflow_from_coeffs is None else lres_value - amflow_from_coeffs)],
        ["l_residue_full_C_closed / AMFlow", fmt_complex(None if amflow_from_coeffs is None else ratio_or_none(lres_value, amflow_from_coeffs))],
        ["l_residue_short_C_diagnostic - AMFlow", fmt_complex(None if amflow_from_coeffs is None else lres_short_value - amflow_from_coeffs)],
        ["l_residue_short_C_diagnostic / AMFlow", fmt_complex(None if amflow_from_coeffs is None else ratio_or_none(lres_short_value, amflow_from_coeffs))],
    ]
    print()
    print("== Fixed-epsilon Comparison ==")
    print_small_table(["quantity", "value"], comparison_rows)
    print_l_residue_closed_coefficient_table(point, lres_coeffs, method12_coeffs, amflow_coeffs, min_power, max_power, label="l_residue_full_C_closed")
    print()
    print("== Short-C Diagnostic Laurent Coefficients ==")
    print_l_residue_closed_coefficient_table(point, lres_short_coeffs, method12_coeffs, amflow_coeffs, min_power, max_power, label="l_residue_short_C_diagnostic")

    powers = list(range(min_power, max_power + 1))
    lres_score = None if not amflow_coeffs else score_coefficients(lres_coeffs, amflow_coeffs, powers)
    method12_score = None if not amflow_coeffs else score_coefficients(method12_coeffs, amflow_coeffs, powers)
    return {
        "point_name": point_name,
        "Delta0": kernel_l_residue_delta0(point),
        "A1": kernel_l_residue_A1(point),
        "B1": kernel_l_residue_B1(point),
        "z": kernel_l_residue_z(point),
        "method12_coeffs": method12_coeffs,
        "lres_coeffs": lres_coeffs,
        "lres_short_coeffs": lres_short_coeffs,
        "amflow_coeffs": amflow_coeffs,
        "lres_score": lres_score,
        "method12_score": method12_score,
    }


def run_method1_l_residue_closed_mode(point_names, eps, eta, args):
    output_dir = repo_root() / "amflow-project" / "targets" / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    report_name = (
        "two_loop_l_residue_closed_all_points.txt"
        if len(point_names) > 1
        else f"two_loop_l_residue_closed_{point_names[0]}.txt"
    )
    report_path = output_dir / report_name
    with report_path.open("w", encoding="utf-8") as report:
        tee = Tee(sys.stdout, report)
        with contextlib.redirect_stdout(tee):
            print("This is a diagnostic candidate expression. Method12 defaults are unchanged.")
            print("Closed l-residue candidate example:")
            print("  ./run.sh two-loop-method3-mc-compare --point equal_mass_offshell_positive --method1-l-residue-closed-check --eps 0.1 --amflow imported")
            results = [run_method1_l_residue_closed_for_point(point_name, eps, eta, args) for point_name in point_names]
            if len(results) > 1:
                print()
                print("== All-Point Summary ==")
                rows = []
                for result in results:
                    if result.get("skipped"):
                        continue
                    point = get_point(result["point_name"])
                    amf = result["amflow_coeffs"]
                    lres = result["lres_coeffs"]
                    m12 = result["method12_coeffs"]
                    best = "not available"
                    if result["lres_score"] is not None and result["method12_score"] is not None:
                        best = "l-residue" if result["lres_score"] <= result["method12_score"] else "Method12"
                    rows.append(
                        [
                            result["point_name"],
                            point.in_support(),
                            mp.nstr(result["Delta0"], 12),
                            mp.nstr(result["A1"], 12),
                            mp.nstr(result["z"], 12),
                            fmt_complex(amf.get(-2)),
                            fmt_complex(amf.get(-1)),
                            fmt_complex(amf.get(0)),
                            fmt_complex(lres.get(-2)),
                            fmt_complex(lres.get(-1)),
                            fmt_complex(lres.get(0)),
                            fmt_complex(m12.get(-2)),
                            fmt_complex(m12.get(-1)),
                            fmt_complex(m12.get(0)),
                            fmt_float(result["lres_score"]),
                            fmt_float(result["method12_score"]),
                            best,
                        ]
                    )
                print_small_table(
                    [
                        "point",
                        "support ok",
                        "Delta0",
                        "A1",
                        "z",
                        "AMF c[-2]",
                        "AMF c[-1]",
                        "AMF c[0]",
                        "lres c[-2]",
                        "lres c[-1]",
                        "lres c[0]",
                        "M12 c[-2]",
                        "M12 c[-1]",
                        "M12 c[0]",
                        "score lres",
                        "score M12",
                        "best match",
                    ],
                    rows,
                )
    print(f"l-residue closed report written to: {report_path}")
    return results


def coefficient_score_against_amflow(candidate, amflow):
    powers = sorted(set(candidate) & set(amflow))
    if not powers:
        return None
    return math.sqrt(sum(abs(candidate[power] - amflow[power]) ** 2 for power in powers))


def compact_coeff_row(power, residue_coeffs, product_coeffs, amflow_coeffs):
    amf = amflow_coeffs.get(power)
    residue = residue_coeffs.get(power)
    product = product_coeffs.get(power)
    return [
        f"c[{power}]",
        fmt_complex(amf),
        fmt_complex(residue),
        fmt_complex(product),
        fmt_complex(None if amf is None else residue - amf),
        fmt_complex(None if amf is None else product - amf),
        fmt_complex(None if amf is None else ratio_or_none(residue, amf)),
        fmt_complex(None if amf is None else ratio_or_none(product, amf)),
    ]


def compare_closed_form_for_point(point_name, eps, eta, args):
    point = get_point(point_name)
    min_power = args.fit_min_power
    max_power = args.fit_max_power
    print()
    print(f"######## Closed-Form Validation: {point.name} ########")
    print("Objects: residue_closed_form, product_closed_form, AMFlow.")
    print("Historical diagnostics are suppressed here; use --diagnose-product-form for product-form convention variants.")
    if point.p_plus != 1 or point.p_perp2 != 0:
        print("residue_closed_form currently supports P=pPlus=1 and pPerp2=0 only.")
        return {"point_name": point_name, "skipped": True}

    ingredients = closed_form_residue_ingredients(point)
    short_ingredients = closed_form_residue_ingredients(point, short_C_diagnostic=True)
    product_coeffs = product_laurent_coefficients(point, max_order=max_power, min_power=min_power, eta=eta)
    residue_coeffs = residue_laurent_coefficients(point, max_order=max_power, min_power=min_power, eta=Fraction(1, 10**40))
    product_value = product_closed_form(point, eps, eta=eta)
    residue_value = residue_closed_form(point, eps, eta=eta)

    amflow_coeffs = {}
    amflow_from_coeffs = None
    if args.amflow in ("fresh", "reuse", "imported"):
        ensure_amflow_order_for_powers(args, max_power)
        _raw, amflow_coeffs, _fields, _combos, output = run_amflow_coefficients(point_name, eps, args)
        if amflow_coeffs:
            amflow_from_coeffs = evaluate_laurent(amflow_coeffs, eps)
        else:
            print("AMFlow coefficients unavailable. Last AMFlow output follows:")
            print(output[-3000:])

    P = point.p_plus
    X = point.p_plus - point.x
    rows = [
        ["point", point.name],
        ["support", point.in_support()],
        ["P = pPlus", P],
        ["pMinus", point.p_minus],
        ["p2", point.p2],
        ["pPerp2", point.p_perp2],
        ["x", point.x],
        ["y", point.y],
        ["X = P - x", X],
        ["lambda = y/X", point.lam],
        ["kappa", point.lam * (1 - point.lam)],
        ["ml2", point.ml2],
        ["Ml2", point.Ml2],
        ["mk2", point.mk2],
        ["Mk2", point.Mk2],
        ["branch convention", "-i0 in Delta0, A1, and the hypergeometric branch argument"],
    ]
    print()
    print("== Point And Branch ==")
    print_small_table(["quantity", "value"], rows)

    print()
    print("== Residue Ingredients ==")
    print_small_table(
        ["ingredient", "value"],
        [
            ["Delta0", mp.nstr(ingredients["Delta0"], 30)],
            ["A1 full current-PDF", mp.nstr(ingredients["A1"], 30)],
            ["B1", mp.nstr(ingredients["B1"], 30)],
            ["z", mp.nstr(ingredients["z"], 30)],
            ["extra_offshell_term", mp.nstr(ingredients["extra_offshell_term"], 30)],
            ["A1_short_diagnostic", mp.nstr(short_ingredients["A1"], 30)],
            ["z_short_diagnostic", mp.nstr(short_ingredients["z"], 30)],
        ],
    )

    print()
    print("== Fixed-eps Values From Closed Forms ==")
    print_small_table(
        ["quantity", "value"],
        [
            ["eps", eps],
            ["residue_closed_form", fmt_complex(residue_value)],
            ["product_closed_form", fmt_complex(product_value)],
            ["AMFlow_from_coeffs", fmt_complex(amflow_from_coeffs)],
            ["residue_closed_form - AMFlow", fmt_complex(None if amflow_from_coeffs is None else residue_value - amflow_from_coeffs)],
            ["product_closed_form - AMFlow", fmt_complex(None if amflow_from_coeffs is None else product_value - amflow_from_coeffs)],
        ],
    )

    print()
    print("== Laurent Coefficient Comparison ==")
    print_small_table(
        [
            "coefficient",
            "AMFlow",
            "residue_closed_form",
            "product_closed_form",
            "residue-AMFlow",
            "product-AMFlow",
            "residue/AMFlow",
            "product/AMFlow",
        ],
        [compact_coeff_row(power, residue_coeffs, product_coeffs, amflow_coeffs) for power in range(min_power, max_power + 1)],
    )

    residue_score = coefficient_score_against_amflow(residue_coeffs, amflow_coeffs)
    product_score = coefficient_score_against_amflow(product_coeffs, amflow_coeffs)
    best_match = "not available"
    if residue_score is not None and product_score is not None:
        best_match = "residue_closed_form" if residue_score <= product_score else "product_closed_form"
    print()
    print("== Scores ==")
    print_small_table(
        ["quantity", "value"],
        [
            ["score residue_closed_form vs AMFlow", fmt_float(residue_score)],
            ["score product_closed_form vs AMFlow", fmt_float(product_score)],
            ["best match", best_match],
        ],
    )
    return {
        "point_name": point.name,
        "residue_coeffs": residue_coeffs,
        "product_coeffs": product_coeffs,
        "amflow_coeffs": amflow_coeffs,
        "residue_score": residue_score,
        "product_score": product_score,
        "best_match": best_match,
    }


def run_compare_closed_form_mode(point_names, eps, eta, args):
    output_dir = repo_root() / "amflow-project" / "targets" / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    report_name = (
        "two_loop_residue_closed_compare_all_points.txt"
        if len(point_names) > 1
        else f"two_loop_residue_closed_compare_{point_names[0]}.txt"
    )
    report_path = output_dir / report_name
    with report_path.open("w", encoding="utf-8") as report:
        tee = Tee(sys.stdout, report)
        with contextlib.redirect_stdout(tee):
            print("Closed-form validation report.")
            print("Main objects: residue_closed_form, product_closed_form, AMFlow.")
            results = [compare_closed_form_for_point(point_name, eps, eta, args) for point_name in point_names]
            if len(results) > 1:
                print()
                print("== Compact All-Point Summary ==")
                print_small_table(
                    ["point", "score residue", "score product", "best match"],
                    [
                        [
                            result["point_name"],
                            fmt_float(result.get("residue_score")),
                            fmt_float(result.get("product_score")),
                            result.get("best_match", "not available"),
                        ]
                        for result in results
                        if not result.get("skipped")
                    ],
                )
    print(f"Closed-form validation report written to: {report_path}")
    return results


def compare_coefficients_for_point(point_name, eps_values, eta, args):
    point = get_point(point_name)
    method3_point = method3_point_for_convention(point, args.method3_label_convention)
    min_power = args.fit_min_power
    max_power = args.fit_max_power
    rho = None if args.rho is None else parse_mpf(args.rho)

    print()
    print(f"######## Coefficient Comparison: {point.name} ########")
    print("epsList = {", ", ".join(str(eps) for eps in eps_values), "}", sep="")
    print(f"fit powers = {min_power}..{max_power}")
    print(f"N = {args.N}, scrambles = {args.scrambles}, eta = {eta}")
    print_method3_label_convention(point, method3_point, args.method3_label_convention)

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
    unsupported_reason = method3_unsupported_reason(method3_point)
    if unsupported_reason is None:
        print()
        print("== Evaluating Method 3 At epsList ==")
        for eps in eps_values:
            stripped, stripped_err = method3_stripped(
                point=method3_point,
                eps=eps,
                eta=eta,
                n_samples=args.N,
                rho=rho,
                scrambles=args.scrambles,
                seed=args.seed,
            )
            prefactor = method3_prefactor(method3_point, eps, lc_jacobian=parse_float(args.lc_jacobian))
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
    amflow_raw_coeffs = {}
    amflow_normalized_coeffs = {}
    amflow_object_fields = {}
    amflow_combo_coeffs = {}
    if args.amflow in ("fresh", "reuse", "imported"):
        print()
        print("== Extracting AMFlow Laurent Coefficients ==")
        ensure_amflow_order_for_powers(args, max_power)
        if args.amflow == "imported":
            print("Loading imported AMFlow coefficient log/result for this point.")
        elif args.amflow == "reuse":
            print("Reusing cached AMFlow double-discontinuity expression and extracting coefficients.")
        else:
            print("Running one AMFlow solve for this point, then extracting coefficients from the resulting epsilon expression.")
        print("No AMFlow finite-epsilon scan is run in coefficient mode.")
        amflow_raw_coeffs, amflow_normalized_coeffs, amflow_object_fields, amflow_combo_coeffs, output = run_amflow_coefficients(
            point_name,
            eps_values[0],
            args,
        )
        if not amflow_normalized_coeffs:
            print("AMFlow coefficients unavailable. Last AMFlow output follows:")
            print(output[-4000:])
        else:
            print_amflow_object_report(point, amflow_object_fields, amflow_raw_coeffs, amflow_normalized_coeffs, args)
            print_amflow_method12_coefficient_comparison(
                point,
                method12_analytic,
                amflow_raw_coeffs,
                amflow_normalized_coeffs,
                min_power,
                max_power,
            )
            print_amflow_fixed_eps_reconstruction(
                point,
                eps_values,
                method12_values,
                None if not method3_values else method3_values,
                amflow_normalized_coeffs,
            )
            print_common_factor_diagnostics(
                point,
                method12_analytic,
                amflow_normalized_coeffs,
                eps_values,
                method12_values,
            )
            print_amflow_sign_pattern_diagnostics(amflow_combo_coeffs, min_power, max_power)
            if args.diagnose_method12_conventions:
                print_method12_convention_diagnostic(
                    point,
                    amflow_normalized_coeffs,
                    min_power,
                    max_power,
                    parse_mpf(args.eps),
                )
            amflow_fit = amflow_normalized_coeffs
    elif args.diagnose_method12_conventions:
        print()
        print("== Method12 Convention Diagnostic Skipped ==")
        print("Run with --amflow fresh, --amflow reuse, or --amflow imported so AMFlow Laurent coefficients are available.")

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
                "AMFlow_normalized",
                "AMFlow_normalized - method12_analytic",
                "method3_fit - AMFlow_normalized",
            ],
            [
                [
                    f"c[{power}]",
                    fmt_complex(amflow_fit.get(power)),
                    fmt_complex(None if power not in amflow_fit else amflow_fit[power] - method12_analytic[power]),
                    fmt_complex(None if method3_fit is None or power not in amflow_fit else method3_fit[power] - amflow_fit[power]),
                ]
                for power in range(min_power, max_power + 1)
            ],
        )

    return {
        "point_name": point_name,
        "method12_analytic": method12_analytic,
        "method12_fit": method12_fit,
        "method3_fit": method3_fit,
        "AMFlow_raw": amflow_raw_coeffs,
        "AMFlow_normalized": amflow_normalized_coeffs,
        "AMFlow_sign_pattern_diagnostics": amflow_combo_coeffs,
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
    used_args = set(sys.argv[1:])
    if "--method1-l-residue-closed-check" in used_args:
        print("Deprecated alias: --method1-l-residue-closed-check -> --compare-closed-form")
    if "--method1-l-residue-direct-check" in used_args:
        print("Deprecated alias: --method1-l-residue-direct-check -> --residue-1d-check")
    if "--method3-integrand-eps-expansion" in used_args:
        print("Deprecated alias: --method3-integrand-eps-expansion -> --transverse-integrand-expansion-check")
    if "--diagnose-method12-conventions" in used_args:
        print("Deprecated alias: --diagnose-method12-conventions -> --diagnose-product-form")

    if args.show_amflow_import_help:
        print_amflow_import_help()
        return

    process_amflow_import_args(args)

    if args.compare_coefficients:
        eps_values = parse_eps_list(args.eps_list)
    elif args.eps_list is not None:
        eps_values = parse_eps_list(args.eps_list)
    else:
        eps_values = [parse_mpf(args.eps)]

    if args.points is not None:
        point_names = [name.strip() for name in args.points.split(",") if name.strip()]
        unknown = [name for name in point_names if name not in POINTS]
        if unknown:
            raise SystemExit(f"Unknown point(s): {', '.join(unknown)}. Known points: {', '.join(sorted(POINTS))}")
    elif args.all_points:
        point_names = sorted(POINTS)
    else:
        point_names = [args.point]
    eta = parse_mpf(args.eta)

    if args.method3_coeff_cutoff_log:
        run_method3_cutoff_log_mode(point_names, eta, args)
        return

    if args.method3_integrand_eps_expansion:
        run_method3_integrand_expansion_mode(point_names, eta, args)
        return

    if args.compare_closed_form:
        if args.method1_l_residue_direct_check:
            require_mpmath_for_l_residue()
        run_compare_closed_form_mode(point_names, parse_mpf(args.eps), eta, args)
        return

    if args.method1_l_residue_direct_check:
        require_mpmath_for_l_residue()
        run_method1_l_residue_direct_mode(point_names, parse_mpf(args.eps), eta, args)
        return

    if args.compare_coefficients:
        run_coefficient_mode(point_names, eps_values, eta, args)
        return

    if args.radial_cutoff_scan:
        run_radial_cutoff_scan(point_names, eps_values, eta, args)
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
