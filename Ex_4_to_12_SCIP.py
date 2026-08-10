"""Run SCIP directly on Examples 4--12 from the LogPOP manuscript.

This is an experimental comparison script.  It reports separately whether
SCIP returned a feasible optimizer candidate and whether its primal--dual gap
certifies global optimality to the requested tolerance.
"""
# J. Choi, August 5, 2026

from __future__ import annotations

import argparse
import math
from pathlib import Path
import sys
import time


SCRIPT_DIR = Path(__file__).resolve().parent
PACKAGE_DIR = SCRIPT_DIR.parent / "python_packages"
sys.path.insert(0, str(PACKAGE_DIR))

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


def example_4(time_limit, relative_gap):
    model, x = simplex_model("example_4", 3, time_limit, relative_gap)
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


def example_5(time_limit, relative_gap):
    model = new_model("example_5", time_limit, relative_gap)
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


EXAMPLE_6_WEIGHTS = {
    "6_i": [0.0968, 0.1419, 0.2194, 0.0839, 0.2839, 0.1742],
    "6_ii_1": [0.7691, 0.6389, 0.8931, 0.0607, 0.1758, 0.4163],
    "6_ii_2": [0.1774, 0.3959, 0.4922, 0.4379, 0.6354, 0.1527],
    "6_ii_3": [0.2920, 0.4317, 0.0155, 0.9841, 0.1672, 0.1062],
    "6_ii_4": [0.0835, 0.6260, 0.6609, 0.7298, 0.8908, 0.9823],
    "6_ii_5": [0.3424, 0.7360, 0.7947, 0.5449, 0.6862, 0.8936],
}
EXAMPLE_6_REFERENCE = {
    "6_i": None,
    "6_ii_1": -4.6967648309,
    "6_ii_2": -3.9657931018,
    "6_ii_3": None,
    "6_ii_4": -6.7017803712,
    "6_ii_5": -7.0672207272,
}


def example_6(case, time_limit, relative_gap):
    model, x = simplex_model(f"example_{case}", 3, time_limit, relative_gap)
    model.addCons(quicksum(t * t for t in x) <= 1.0, name="redundant_ball")
    p = [
        x[0] ** 3 + 3 * x[0] ** 2 * x[1] + 3 * x[0] ** 2 * x[2],
        3 * x[0] * x[1] ** 2 + 6 * x[0] * x[1] * x[2],
        3 * x[0] * x[2] ** 2,
        x[1] ** 3 + 3 * x[1] ** 2 * x[2],
        3 * x[1] * x[2] ** 2,
        x[2] ** 3,
    ]
    weights = EXAMPLE_6_WEIGHTS[case]
    add_log_objective(model, p, weights, [1.0] * 6)

    def objective(v):
        values = [
            v[0] ** 3 + 3 * v[0] ** 2 * v[1] + 3 * v[0] ** 2 * v[2],
            3 * v[0] * v[1] ** 2 + 6 * v[0] * v[1] * v[2],
            3 * v[0] * v[2] ** 2,
            v[1] ** 3 + 3 * v[1] ** 2 * v[2],
            3 * v[1] * v[2] ** 2,
            v[2] ** 3,
        ]
        return sum(a * math.log(q) for a, q in zip(weights, values))

    return model, x, objective, EXAMPLE_6_REFERENCE[case]


def example_7(time_limit, relative_gap):
    model = new_model("example_7", time_limit, relative_gap)
    x = [model.addVar(lb=-1.0, ub=1.0, name="x_1")]
    p = [2 + x[0] ** 2, 3 + x[0] ** 2]
    add_log_objective(model, p, [1.0, 1.0], [3.0, 4.0])
    return (
        model,
        x,
        lambda v: math.log(2 + v[0] ** 2) + math.log(3 + v[0] ** 2),
        math.log(12.0),
    )


def example_8(case, time_limit, relative_gap):
    model = new_model(f"example_{case}", time_limit, relative_gap)
    x = [model.addVar(lb=-2.0, ub=2.0, name=f"x_{i + 1}") for i in range(5)]
    model.addCons(quicksum(t * t for t in x) <= 4.0, name="radius_two_ball")
    p1 = (x[2] + x[3] + x[4]) ** 2 - x[0] * x[1]
    delta = x[4] - x[1]
    p2 = 8 - delta**2
    model.addCons(p1 >= 0.0, name="p1_nonnegative")
    model.addCons(p2 >= 0.0, name="p2_nonnegative")
    if case == "8_i":
        p = [p1, p2]
        weights = [20.0, 25.0]
        upper = [14.0, 8.0]

        def objective(v):
            q1 = (v[2] + v[3] + v[4]) ** 2 - v[0] * v[1]
            q2 = 8 - (v[4] - v[1]) ** 2
            return 20 * math.log(q1) + 25 * math.log(q2)

    else:
        p3 = math.sqrt(8.0) - delta
        p4 = math.sqrt(8.0) + delta
        p = [p1, p3, p4]
        weights = [20.0, 25.0, 25.0]
        upper = [14.0, 2 * math.sqrt(8.0), 2 * math.sqrt(8.0)]

        def objective(v):
            q1 = (v[2] + v[3] + v[4]) ** 2 - v[0] * v[1]
            delta_v = v[4] - v[1]
            return (
                20 * math.log(q1)
                + 25 * math.log(math.sqrt(8.0) - delta_v)
                + 25 * math.log(math.sqrt(8.0) + delta_v)
            )

    add_log_objective(model, p, weights, upper)
    return model, x, objective, 99.7251825577


def example_9(time_limit, relative_gap):
    model = new_model("example_9", time_limit, relative_gap)
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


def example_10(time_limit, relative_gap):
    model = new_model("example_10", time_limit, relative_gap)
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


PATERNITY_N = [
    [77, 23], [63, 37], [49, 40, 11], [83, 2, 15], [63, 17, 20],
    [59, 8, 16, 17], [7, 9, 4, 33, 47], [39, 38, 23], [29, 21, 88, 62],
]
PATERNITY_P = [
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
PATERNITY_REF = [
    -53.9276341497, -65.8955680683, -104.3236621937,
    -111.2570075021, -115.1463552694, -132.7325824971,
    -128.8283619893, -107.2934733638, -265.7241438158,
]


def example_11(case_index, time_limit, relative_gap):
    counts = PATERNITY_N[case_index]
    probabilities = PATERNITY_P[case_index]
    dimension = len(probabilities[0])
    model, x = simplex_model(
        f"example_11_{case_index + 1}", dimension, time_limit, relative_gap
    )
    model.addCons(quicksum(x) == 1.0, name="simplex_equality")
    p = [
        quicksum(row[j] * x[j] for j in range(dimension))
        for row in probabilities
    ]
    add_log_objective(model, p, [float(v) for v in counts], [1.0] * len(p))

    def objective(v):
        return sum(
            count * math.log(sum(row[j] * v[j] for j in range(dimension)))
            for count, row in zip(counts, probabilities)
        )

    return model, x, objective, PATERNITY_REF[case_index]


LCM_COUNTS = {
    "12_i": [303, 197],
    "12_ii": [227, 273],
    "12_iii": [131, 174, 49, 146],
    "12_iv": [30, 26, 39, 26, 27, 19, 28, 30, 24, 25, 39, 33, 45, 31, 35, 43],
}
LCM_REF = {
    "12_i": -335.2518754137,
    "12_ii": -344.4545954602,
    "12_iii": -652.6718169656,
    "12_iv": None,
}


def example_12(case, time_limit, relative_gap):
    counts = LCM_COUNTS[case]
    model = new_model(f"example_{case}", time_limit, relative_gap)
    if case == "12_i":
        x = [model.addVar(lb=0, ub=1, name=f"x_{i + 1}") for i in range(6)]
        model.addCons(x[0] + x[1] == 1)
        model.addCons(x[2] + x[3] == 1)
        model.addCons(x[4] + x[5] == 1)
        p = [x[0] * x[2] + x[1] * x[4], x[0] * x[3] + x[1] * x[5]]

        def values(v):
            return [v[0] * v[2] + v[1] * v[4], v[0] * v[3] + v[1] * v[5]]

    elif case == "12_ii":
        x = [model.addVar(lb=0, ub=1, name=f"x_{i + 1}") for i in range(9)]
        model.addCons(x[0] + x[1] + x[2] == 1)
        model.addCons(x[3] + x[4] == 1)
        model.addCons(x[5] + x[6] == 1)
        model.addCons(x[7] + x[8] == 1)
        p = [
            x[0] * x[3] + x[1] * x[5] + x[2] * x[7],
            x[0] * x[4] + x[1] * x[6] + x[2] * x[8],
        ]

        def values(v):
            return [
                v[0] * v[3] + v[1] * v[5] + v[2] * v[7],
                v[0] * v[4] + v[1] * v[6] + v[2] * v[8],
            ]

    elif case == "12_iii":
        x = [model.addVar(lb=0, ub=1, name=f"x_{i + 1}") for i in range(10)]
        for first in [0, 2, 4, 6, 8]:
            model.addCons(x[first] + x[first + 1] == 1)
        p = [
            x[0] * x[2] * x[4] + x[1] * x[6] * x[8],
            x[0] * x[2] * x[5] + x[1] * x[6] * x[9],
            x[0] * x[3] * x[4] + x[1] * x[7] * x[8],
            x[0] * x[3] * x[5] + x[1] * x[7] * x[9],
        ]

        def values(v):
            return [
                v[0] * v[2] * v[4] + v[1] * v[6] * v[8],
                v[0] * v[2] * v[5] + v[1] * v[6] * v[9],
                v[0] * v[3] * v[4] + v[1] * v[7] * v[8],
                v[0] * v[3] * v[5] + v[1] * v[7] * v[9],
            ]

    else:
        x = [model.addVar(lb=0, ub=1, name=f"x_{i + 1}") for i in range(9)]
        patterns = [[(r >> (3 - j)) & 1 for j in range(4)] for r in range(16)]
        p = []
        for pattern in patterns:
            component_1 = 1
            component_2 = 1
            for j, bit in enumerate(pattern):
                component_1 *= x[1 + j] if bit else (1 - x[1 + j])
                component_2 *= x[5 + j] if bit else (1 - x[5 + j])
            p.append(x[0] * component_1 + (1 - x[0]) * component_2)

        def values(v):
            output = []
            for pattern in patterns:
                component_1 = math.prod(v[1 + j] if bit else 1 - v[1 + j] for j, bit in enumerate(pattern))
                component_2 = math.prod(v[5 + j] if bit else 1 - v[5 + j] for j, bit in enumerate(pattern))
                output.append(v[0] * component_1 + (1 - v[0]) * component_2)
            return output

    add_log_objective(model, p, [float(v) for v in counts], [1.0] * len(p))

    def objective(v):
        return sum(count * math.log(q) for count, q in zip(counts, values(v)))

    return model, x, objective, LCM_REF[case]


BUILDERS = {
    "4": example_4,
    "5": example_5,
    **{case: (lambda tl, rg, c=case: example_6(c, tl, rg)) for case in EXAMPLE_6_WEIGHTS},
    "7": example_7,
    "8_i": lambda tl, rg: example_8("8_i", tl, rg),
    "8_ii": lambda tl, rg: example_8("8_ii", tl, rg),
    "9": example_9,
    "10": example_10,
    **{f"11_{i + 1}": (lambda tl, rg, j=i: example_11(j, tl, rg)) for i in range(9)},
    **{case: (lambda tl, rg, c=case: example_12(c, tl, rg)) for case in LCM_COUNTS},
}


def run_case(case, time_limit, relative_gap):
    print(f"\n{'=' * 72}\nSCIP comparison for Example {case.replace('_', '(')}")
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
        choices=(*BUILDERS.keys(), "all", "easy", "hard"),
        default="all",
    )
    parser.add_argument("--time-limit", type=float, default=600.0)
    parser.add_argument("--relative-gap", type=float, default=1.0e-6)
    args = parser.parse_args()
    hard_cases = ["10", "12_iii", "12_iv"]
    if args.example == "all":
        selected = list(BUILDERS)
    elif args.example == "easy":
        selected = [case for case in BUILDERS if case not in hard_cases]
    elif args.example == "hard":
        selected = hard_cases
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
    print(f"{'case':<9} {'status':<12} {'candidate':<10} {'value':>14} {'upper':>14} {'cert':<6} {'match':<6} {'time':>9}")
    print("-" * 88)
    for row in results:
        match = "NA" if row["matches"] is None else ("yes" if row["matches"] else "no")
        print(
            f"{row['case']:<9} {row['status']:<12} "
            f"{('yes' if row['candidate'] else 'no'):<10} "
            f"{row['value']:>14.6g} {row['upper']:>14.6g} "
            f"{('yes' if row['certified'] else 'no'):<6} {match:<6} {row['time']:>9.3g}"
        )


if __name__ == "__main__":
    main()
