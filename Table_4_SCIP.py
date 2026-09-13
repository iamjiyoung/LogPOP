"""Compare SCIP with the proposed relaxations on four representative examples.

This is an experimental comparison script.  It reports separately whether
SCIP returned a feasible optimizer candidate and whether its primal--dual gap
certifies global optimality to the requested tolerance.
"""
# J. Choi, September 7, 2026

from __future__ import annotations

import argparse
import math
from pathlib import Path
import sys
import time


SCRIPT_DIR = Path(__file__).resolve().parent
PACKAGE_CANDIDATES = [
    SCRIPT_DIR / "python_packages",
    SCRIPT_DIR.parent / "python_packages",
    Path.home() / "Desktop" / "MATLAB" / "26_LogPOP" / "python_packages",
]
for package_dir in PACKAGE_CANDIDATES:
    if package_dir.is_dir():
        sys.path.insert(0, str(package_dir))
        break

import pyscipopt
from pyscipopt import Model, log, quicksum


EPSILON = 1.0e-9


def add_log_objective(model, polynomials, weights, upper_bounds):
    log_values = []
    for index, (polynomial, upper) in enumerate(
        zip(polynomials, upper_bounds), start=1
    ):
        p_value = model.addVar(lb=EPSILON, ub=upper, name=f"p_{index}")
        model.addCons(p_value == polynomial, name=f"p_link_{index}")
        z_value = model.addVar(
            lb=math.log(EPSILON), ub=math.log(upper), name=f"log_p_{index}"
        )
        model.addCons(z_value <= log(p_value), name=f"log_link_{index}")
        log_values.append(z_value)
    model.setObjective(
        quicksum(weight * z for weight, z in zip(weights, log_values)),
        sense="maximize",
    )


def new_model(name, time_limit, relative_gap):
    model = Model(name)
    model.setRealParam("limits/time", time_limit)
    model.setRealParam("limits/gap", relative_gap)
    model.setRealParam("numerics/feastol", 1.0e-9)
    model.setIntParam("randomization/randomseedshift", 0)
    return model


def simplex_model(name, dimension, time_limit, relative_gap):
    model = new_model(name, time_limit, relative_gap)
    x = [model.addVar(lb=0.0, ub=1.0, name=f"x_{i + 1}") for i in range(dimension)]
    model.addCons(quicksum(x) <= 1.0, name="simplex")
    return model, x


def example_5_3(time_limit, relative_gap):
    model, x = simplex_model("example_5_3", 3, time_limit, relative_gap)
    p = [
        x[0] ** 2 + 2 * x[0] * x[2],
        x[1] ** 2 + 2 * x[1] * x[2],
        2 * x[0] * x[1],
        x[2] ** 2,
    ]
    weights = [182.0, 60.0, 17.0, 176.0]
    add_log_objective(model, p, weights, [1.0] * 4)
    return model, x, lambda v: sum(
        a * math.log(q) for a, q in zip(weights, [
            v[0] ** 2 + 2 * v[0] * v[2],
            v[1] ** 2 + 2 * v[1] * v[2],
            2 * v[0] * v[1],
            v[2] ** 2,
        ])
    ), -492.5353166834


def example_6_1(time_limit, relative_gap):
    model = new_model("example_6_1", time_limit, relative_gap)
    x = [model.addVar(lb=-1.0, ub=1.0, name=f"x_{i + 1}") for i in range(5)]
    model.addCons(quicksum(t * t for t in x) <= 1.0, name="unit_ball")
    p = [
        (x[2] + x[3]) ** 2
        + (1 + 2 * x[0] + 3 * x[1] + 3 * x[3] + 2 * x[4]) ** 2
        + 0.01,
        (x[1] + x[3] - x[4]) ** 2 + (2 * x[1] + 3 * x[4]) ** 2 + 0.02,
        (x[0] + x[2] + x[3]) ** 2 + (x[0] - x[2] + x[3]) ** 2 + 0.03,
    ]
    weights = [30 / 218, 97 / 218, 91 / 218]
    add_log_objective(model, p, weights, [50.0, 30.0, 20.0])

    def objective(v):
        values = [
            (v[2] + v[3]) ** 2
            + (1 + 2 * v[0] + 3 * v[1] + 3 * v[3] + 2 * v[4]) ** 2
            + 0.01,
            (v[1] + v[3] - v[4]) ** 2 + (2 * v[1] + 3 * v[4]) ** 2 + 0.02,
            (v[0] + v[2] + v[3]) ** 2 + (v[0] - v[2] + v[3]) ** 2 + 0.03,
        ]
        return sum(a * math.log(q) for a, q in zip(weights, values))

    return model, x, objective, 1.6292207561



def example_6_4(time_limit, relative_gap):
    model = new_model("example_6_4", time_limit, relative_gap)
    x = [model.addVar(lb=-10.0, ub=10.0, name=f"x_{i + 1}") for i in range(10)]
    p = [
        10 - (2 * (x[0] + 0.5) ** 2 + 3 * (x[1] + 0.5) ** 2
              - 3 * (x[0] - 0.5) * (x[1] + 0.5)),
        12 - (3 * (x[2] - 0.6) ** 2 + 2 * (x[3] + 0.1) ** 2
              - (x[2] - 0.6) * (x[3] - 0.1)),
        11 - (4 * (x[4] + 0.7) ** 2 + (x[5] - 0.2) ** 2),
        15 - (2 * (x[6] - 0.8) ** 2 + 5 * (x[7] + 0.3) ** 2
              - 2 * (x[6] - 0.8) * (x[7] + 0.3)),
        13 - (-3 * (x[8] - 0.9) ** 2 + 2 * (x[9] + 0.4) ** 2),
    ]
    for i, polynomial in enumerate(p, start=1):
        model.addCons(polynomial >= i, name=f"p_lower_{i}")
    model.addCons((x[0] + x[1] + x[2]) ** 2 + (x[3] + x[4]) ** 2 >= 8)
    model.addCons((x[5] - x[6] - x[7]) ** 2 + x[8] ** 2 <= 6)
    model.addCons((x[0] - x[9]) ** 2 + (x[1] - x[8]) ** 2 <= 7)
    model.addCons((x[1] + x[3] - x[5] + x[7] - x[9]) ** 2 >= 9)
    weights = [3.0, 2.0, 1.0, 1.0, 5.0]
    add_log_objective(model, p, weights, [1.0e4] * 5)

    def objective(v):
        values = [
            10 - (2 * (v[0] + 0.5) ** 2 + 3 * (v[1] + 0.5) ** 2
                  - 3 * (v[0] - 0.5) * (v[1] + 0.5)),
            12 - (3 * (v[2] - 0.6) ** 2 + 2 * (v[3] + 0.1) ** 2
                  - (v[2] - 0.6) * (v[3] - 0.1)),
            11 - (4 * (v[4] + 0.7) ** 2 + (v[5] - 0.2) ** 2),
            15 - (2 * (v[6] - 0.8) ** 2 + 5 * (v[7] + 0.3) ** 2
                  - 2 * (v[6] - 0.8) * (v[7] + 0.3)),
            13 - (-3 * (v[8] - 0.9) ** 2 + 2 * (v[9] + 0.4) ** 2),
        ]
        return sum(a * math.log(q) for a, q in zip(weights, values))

    return model, x, objective, 36.3993754433


def example_6_5(time_limit, relative_gap):
    model = new_model("example_6_5", time_limit, relative_gap)
    radius = math.sqrt(20.0)
    x = [model.addVar(lb=-radius, ub=radius, name=f"x_{i + 1}") for i in range(12)]
    h = [
        (x[1] - 0.2) + (x[2] - 1.6),
        (x[3] - 1.0) + (x[4] - 0.4),
        (x[5] + 0.4) + (x[6] - 1.0),
        (x[7] + 0.4) + (x[8] + 0.8),
        (x[9] - 1.1) + (x[10] - 1.5),
        (x[0] - 0.8) - (x[5] - 0.7) + (x[10] - 1.3),
        -(x[1] - 0.3) + x[6] - (x[11] + 0.7),
        (x[2] - 1.7) - (x[7] - 1.0) + (x[0] - 0.7),
        -(x[3] + 1.2) + (x[8] - 1.9) - (x[1] + 0.8),
        (x[4] + 1.7) - (x[9] + 0.9) + (x[2] - 0.5),
    ]
    alpha = [20.0, 25.0, 18.0, 22.0, 30.0, 15.0, 28.0, 19.0, 21.0, 26.0]
    beta = [1.5, 1.2, 2.5, 1.1, 1.6, 1.9, 2.3, 1.7, 2.4, 1.4]
    p = [a - b * h_i**4 for a, b, h_i in zip(alpha, beta, h)]
    model.addCons(quicksum(x[:6]) ** 2 <= 15)
    model.addCons(quicksum(x[6:]) ** 2 <= 15)
    model.addCons((x[0] - x[2] + x[4] - x[6] + x[8] - x[10]) ** 2 <= 8)
    model.addCons((x[1] - x[3] + x[5] - x[7] + x[9] - x[11]) ** 2 <= 8)
    model.addCons((x[0] + x[11]) ** 2 + (x[1] + x[10]) ** 2 <= 9)
    model.addCons((x[2] - x[9]) ** 2 + (x[3] - x[8]) ** 2 + (x[4] - x[7]) ** 2 <= 9)
    model.addCons(quicksum(x) >= 0)
    model.addCons(quicksum(x[0::2]) - quicksum(x[1::2]) <= 5)
    model.addCons(
        2 * x[0] ** 2 + x[1] ** 2 + 2 * x[2] ** 2 + x[3] ** 2
        + 3 * x[4] ** 2 + x[5] ** 2 + 2 * x[6] ** 2 + x[7] ** 2
        + 2 * x[8] ** 2 + x[9] ** 2 + 3 * x[10] ** 2 + x[11] ** 2 <= 20
    )
    model.addCons((x[0] + x[5] + x[11]) ** 2 <= 4)
    model.addCons(x[0] - x[11] <= 3)
    for polynomial in p:
        model.addCons(polynomial >= 0)
    add_log_objective(model, p, [1.0] * 10, alpha)

    def objective(v):
        hv = [
            (v[1] - 0.2) + (v[2] - 1.6),
            (v[3] - 1.0) + (v[4] - 0.4),
            (v[5] + 0.4) + (v[6] - 1.0),
            (v[7] + 0.4) + (v[8] + 0.8),
            (v[9] - 1.1) + (v[10] - 1.5),
            (v[0] - 0.8) - (v[5] - 0.7) + (v[10] - 1.3),
            -(v[1] - 0.3) + v[6] - (v[11] + 0.7),
            (v[2] - 1.7) - (v[7] - 1.0) + (v[0] - 0.7),
            -(v[3] + 1.2) + (v[8] - 1.9) - (v[1] + 0.8),
            (v[4] + 1.7) - (v[9] + 0.9) + (v[2] - 0.5),
        ]
        return sum(math.log(a - b * q**4) for a, b, q in zip(alpha, beta, hv))

    return model, x, objective, 30.2623752132



BUILDERS = {
    "5.3": example_5_3,
    "6.1": example_6_1,
    "6.4": example_6_4,
    "6.5": example_6_5,
}


def run_case(case, time_limit, relative_gap):
    print(f"\n{'=' * 72}\nSCIP comparison for Example {case}")
    model, variables, evaluator, reference = BUILDERS[case](time_limit, relative_gap)
    started = time.perf_counter()
    model.optimize()
    wall_time = time.perf_counter() - started
    status = str(model.getStatus())
    best_solution = model.getBestSol()
    candidate_found = best_solution is not None
    point = None
    evaluated_value = math.nan
    if candidate_found:
        point = [model.getSolVal(best_solution, variable) for variable in variables]
        try:
            evaluated_value = evaluator(point)
        except (ValueError, ZeroDivisionError):
            candidate_found = False
    upper_bound = model.getDualbound()
    reported_gap = model.getGap()
    certified = (
        math.isfinite(upper_bound)
        and math.isfinite(evaluated_value)
        and max(0.0, upper_bound - evaluated_value)
        / max(1.0, abs(evaluated_value))
        <= 1.05 * relative_gap
    )
    matches_reference = None if reference is None or not math.isfinite(evaluated_value) else (
        abs(evaluated_value - reference) <= 1.0e-4 * max(1.0, abs(reference))
    )
    result = {
        "case": case,
        "status": status,
        "candidate": candidate_found,
        "value": evaluated_value,
        "upper": upper_bound,
        "gap": reported_gap,
        "time": model.getSolvingTime(),
        "wall": wall_time,
        "nodes": model.getNTotalNodes(),
        "certified": certified,
        "matches": matches_reference,
        "point": point,
    }
    print(f"Status                         : {status}")
    print(f"Optimizer candidate returned   : {'yes' if candidate_found else 'no'}")
    print(f"Evaluated original objective   : {evaluated_value:.12g}")
    print(f"SCIP certified upper bound     : {upper_bound:.12g}")
    print(f"SCIP reported relative gap     : {reported_gap:.6g}")
    print(f"Globality certified (tolerance): {'yes' if certified else 'no'}")
    if matches_reference is not None:
        print(f"Matches moment-certified value : {'yes' if matches_reference else 'no'}")
    else:
        print("Matches moment-certified value : not available")
    print(f"SCIP solver time               : {result['time']:.6g} seconds")
    print(f"Wall-clock runtime             : {wall_time:.6g} seconds")
    print(f"Branch-and-bound nodes         : {result['nodes']}")
    if point is not None:
        print("Returned point                 : [" + ", ".join(f"{v:.9g}" for v in point) + "]")
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--example",
        choices=(*BUILDERS.keys(), "all"),
        default="all",
    )
    parser.add_argument("--time-limit", type=float, default=600.0)
    parser.add_argument("--relative-gap", type=float, default=1.0e-6)
    args = parser.parse_args()
    if args.example == "all":
        selected = list(BUILDERS)
    else:
        selected = [args.example]
    version_model = Model()
    version = (
        version_model.getMajorVersion(), version_model.getMinorVersion(),
        version_model.getTechVersion(),
    )
    print(f"PySCIPOpt version              : {pyscipopt.__version__}")
    print(f"SCIP version                   : {'.'.join(map(str, version))}")
    print(f"Relative-gap tolerance         : {args.relative_gap:.1e}")
    print(f"Time limit per case            : {args.time_limit:g} seconds")
    results = [run_case(case, args.time_limit, args.relative_gap) for case in selected]
    print("\n=== SCIP comparison summary ===")
    print(f"{'case':<9} {'status':<12} {'candidate':<10} {'value':>14} {'upper':>14} {'cert':<6} {'match':<6} {'elapsed':>9}")
    print("-" * 88)
    for row in results:
        match = "NA" if row["matches"] is None else ("yes" if row["matches"] else "no")
        print(
            f"{row['case']:<9} {row['status']:<12} "
            f"{('yes' if row['candidate'] else 'no'):<10} "
            f"{row['value']:>14.6g} {row['upper']:>14.6g} "
            f"{('yes' if row['certified'] else 'no'):<6} {match:<6} {row['wall']:>9.3g}"
        )


if __name__ == "__main__":
    main()
