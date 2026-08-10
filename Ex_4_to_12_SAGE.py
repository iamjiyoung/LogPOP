"""Classify and run the direct-product SAGE comparison for Examples 4--12.

The formulation is applicable when the variables are nonnegative and all
logarithmic weights are positive integers.  ``NL`` means that the finite
expansion exists but is not launched because its size is impractical.
``NA`` means that this direct-product SAGE formulation is not applicable.
"""
# J. Choi, August 5, 2026

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path
import sys
import time
import warnings


SCRIPT_DIR = Path(__file__).resolve().parent
PACKAGE_DIR = SCRIPT_DIR.parent / "python_packages"
sys.path.insert(0, str(PACKAGE_DIR))

MOSEK10_HOME = Path(os.environ.get("MOSEK10_HOME", r"C:\Program Files\Mosek\10.0"))
MOSEK10_PLATFORM = MOSEK10_HOME / "tools" / "platform" / "win64x86"
MOSEK10_PYTHON = MOSEK10_PLATFORM / "purepython" / "3"
MOSEK10_BIN = MOSEK10_PLATFORM / "bin"
if not MOSEK10_PYTHON.is_dir():
    raise RuntimeError("MOSEK 10.0 was not found.")
_mosek_dll_directory = os.add_dll_directory(str(MOSEK10_BIN))
sys.path.insert(0, str(MOSEK10_PYTHON))

import mosek
import numpy as np
import sageopt as so


N_DATA = [
    [77, 23], [63, 37], [49, 40, 11], [83, 2, 15], [63, 17, 20],
    [59, 8, 16, 17], [7, 9, 4, 33, 47], [39, 38, 23], [29, 21, 88, 62],
]
P_DATA = [
    [[0.5, 1], [0.5, 0]],
    [[0.5, 1], [0.5, 0]],
    [[0.5, 0.875], [0.25, 0.125], [0.25, 0]],
    [[0.5, 0.25], [0.5, 0.5], [0, 0.25]],
    [[0.5, 0.25], [0.5, 0.5], [0, 0.25]],
    [[0.25, 0.5], [0.25, 0.5], [0.25, 0], [0.25, 0]],
    [[0.25, 0], [0.25, 0], [0.25, 0], [0.25, 0.5], [0, 0.5]],
    [[0.5, 0, 0.875], [0.25, 0.75, 0.125], [0.25, 0.25, 0]],
    [[0.5, 0, 0], [0.25, 0.75, 0.25], [0.25, 0.25, 0.5], [0, 0, 0.25]],
]
REFERENCE = [
    -53.9276341497, -65.8955680683, -104.3236621937,
    -111.2570075021, -115.1463552694, -132.7325824971,
    -128.8283619893, -107.2934733638, -265.7241438158,
]


CLASSIFICATION = [
    ("4", "NL", "11,163-term integer direct-product expansion"),
    ("5", "NA", "signed variables; the rational weights alone could be rescaled"),
    ("6(i)", "NL", "rational rescaling gives degree above 30,000"),
    ("6(ii)", "NL", "rational rescaling gives a prohibitively large expansion"),
    ("7", "NA", "signed variable; a problem-specific x^2 reformulation would be needed"),
    ("8", "NA", "signed variables"),
    ("9", "NA", "signed variables and a non-log-convex domain"),
    ("10", "NA", "signed variables"),
    ("11(1)--11(7)", "solved", "101 or fewer homogeneous monomials"),
    ("11(8)", "NL", "up to 5,151 homogeneous monomials"),
    ("11(9)", "NL", "up to 20,301 homogeneous monomials"),
    ("12(i)", "NL", "at least 60,192 raw product terms before collection"),
    ("12(ii)", "NL", "more than 700 million raw product terms"),
    ("12(iii)", "NL", "more than 160 million raw product terms"),
    ("12(iv)", "NL", "product expansion is combinatorially prohibitive"),
]


def run_paternity_case(case_index):
    counts = N_DATA[case_index]
    probabilities = P_DATA[case_index]
    dimension = len(probabilities[0])
    z = so.standard_sig_monomials(dimension)
    likelihood_terms = [
        sum(row[j] * z[j] for j in range(dimension))
        for row in probabilities
    ]
    objective = 1
    for polynomial, count in zip(likelihood_terms, counts):
        objective *= polynomial ** int(count)
    number_of_terms = len(objective.c)
    negative_objective = -objective
    simplex = 1-sum(z)
    domain = so.infer_domain(negative_objective, [simplex], [])
    if domain is None:
        raise RuntimeError("Sageopt could not infer the simplex domain.")

    started = time.perf_counter()
    problem = so.sig_relaxation(negative_objective, X=domain, form="dual", ell=0)
    construction_time = time.perf_counter()-started
    started = time.perf_counter()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Argument .* Incorrect array format.*")
        solve_result = problem.solve(solver="MOSEK", verbose=True)
    solver_time = time.perf_counter()-started
    product_upper = -float(problem.value)
    log_upper = math.log(product_upper) if product_upper > 0 else math.nan
    gap = log_upper-REFERENCE[case_index]
    print(f"\n=== Example 11, instance {case_index + 1}: level-zero conditional SAGE ===")
    print(f"Expanded monomial count       : {number_of_terms}")
    print(f"Solver status                 : {problem.status}")
    print(f"Solve-result status           : {solve_result[0]}")
    print(f"SAGE log upper bound          : {log_upper:.10f}")
    print(f"Moment-certified optimum      : {REFERENCE[case_index]:.10f}")
    print(f"SAGE upper-bound gap          : {gap:.10f}")
    print(f"Construction time (seconds)   : {construction_time:.6f}")
    print(f"Solver time (seconds)         : {solver_time:.6f}")
    print(f"Optimizer returned by SAGE    : no")
    print(f"Globality certified by SAGE   : {'yes' if abs(gap) <= 1e-6 else 'no'}")
    return (case_index + 1, number_of_terms, str(problem.status), log_upper, gap,
            construction_time, solver_time)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--case", choices=("classify", "all", *[str(i) for i in range(1, 8)]),
        default="classify",
    )
    args = parser.parse_args()
    version = mosek.Env().getversion()
    print(f"MOSEK Python API version      : {'.'.join(map(str, version))}")
    print(f"sageopt version               : {so.__version__}")
    print("\n=== SAGE applicability classification ===")
    for example, status, reason in CLASSIFICATION:
        print(f"Example {example:<13} {status:<7} {reason}")
    if args.case == "classify":
        return
    selected = range(7) if args.case == "all" else [int(args.case)-1]
    rows = [run_paternity_case(index) for index in selected]
    print("\n=== Example 11 SAGE summary ===")
    print(f"{'case':<6} {'terms':>8} {'status':<12} {'log bound':>14} {'gap':>12} {'build':>10} {'solve':>10}")
    print("-" * 82)
    for case, terms, status, bound, gap, build, solve in rows:
        print(f"{case:<6} {terms:>8} {status:<12} {bound:>14.6f} {gap:>12.6g} {build:>10.3g} {solve:>10.3g}")


if __name__ == "__main__":
    main()
