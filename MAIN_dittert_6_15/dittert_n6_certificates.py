#!/usr/bin/env python3
"""Exact rational certificates for the dimension-six Dittert proof.

The two-sided part reconstructs 98 four-variable rational polynomials from
formulas in the manuscript and proves their positivity by exact tensor-product
Bernstein subdivision.  The one-sided part checks ten exact univariate
certificates.  SymPy is used only for symbolic expansion and NumPy only as an
object-array container; every Bernstein coefficient and every correctness
decision uses ``fractions.Fraction``.  No verification is implemented with
``assert``, so optimization mode performs the same checks as a normal run.
Importing this module has no verification side effects.
"""
from __future__ import annotations

import argparse
import itertools
import json
from fractions import Fraction as F
from math import comb, factorial
from typing import Any

import numpy as np
import sympy as sp

CERTIFICATE_VERSION = "2.0.0"
MANUSCRIPT_SHA256 = "d3d77fee2f6c164b3a2ada4b37c66384a53688c03fbc03aee7effc15d281ba94"


class VerificationError(RuntimeError):
    """Raised when an exact certificate check fails."""


def require(condition: bool, message: str) -> None:
    """Raise a stable verification error instead of relying on ``assert``."""
    if not condition:
        raise VerificationError(message)


# ---------------------------------------------------------------------------
# Exact constants and elementary gaps
# ---------------------------------------------------------------------------


def v(k: int) -> sp.Rational:
    return sp.Rational(1) if k == 0 else sp.Rational(factorial(k), k**k)


GAMMA = sp.Rational(5, 324)
KAPPA = v(5) * sp.Rational(4**4, 5**4)
E = sp.Rational(11, 50)


def verify_constants() -> dict[str, str]:
    require(KAPPA == sp.Rational(6144, 390625), f"kappa mismatch: {KAPPA}")
    require(KAPPA > GAMMA, "kappa must exceed gamma")

    localization = sp.Rational(15, 319)
    intermediate = sp.Rational(7, 32) ** 2
    require(localization < intermediate, "15/319 < (7/32)^2 failed")
    require(intermediate < E**2, "(7/32)^2 < (11/50)^2 failed")

    energy_gap = (
        (1 - GAMMA) * sp.Rational(25, 32) ** 2 / 2
        - sp.Rational(3, 10)
    )
    require(
        energy_gap == sp.Rational(1547, 3317760) and energy_gap > 0,
        f"row-energy gap mismatch: {energy_gap}",
    )
    large_gap = (1 - GAMMA) ** 2 / 108 - sp.Rational(9, 1024)
    require(
        large_gap == sp.Rational(33853, 181398528) and large_gap > 0,
        f"large-entropy gap mismatch: {large_gap}",
    )
    return {
        "gamma": str(GAMMA),
        "kappa": str(KAPPA),
        "localization_gap": str(intermediate - localization),
        "E_gap": str(E**2 - intermediate),
        "row_energy_gap": str(energy_gap),
        "large_entropy_gap": str(large_gap),
    }


# ---------------------------------------------------------------------------
# Exact tensor-product Bernstein arithmetic on the unit cube
# ---------------------------------------------------------------------------

X, Z, U, V = sp.symbols("x z u v")
VARS = (X, Z, U, V)


def as_fraction(q: sp.Rational) -> F:
    require(bool(q.is_Rational), f"non-rational coefficient encountered: {q}")
    return F(int(q.p), int(q.q))


def power_to_bernstein(poly: sp.Poly) -> tuple[np.ndarray, tuple[int, ...]]:
    """Convert a nonzero rational power-basis polynomial on [0,1]^4."""
    require(not poly.is_zero, "zero polynomial cannot certify strict positivity")
    require(tuple(poly.gens) == VARS, f"unexpected variable order: {poly.gens}")
    require(poly.domain == sp.QQ, f"unexpected coefficient domain: {poly.domain}")
    degrees = tuple(poly.degree(var) for var in VARS)
    require(all(degree >= 0 for degree in degrees), f"invalid multidegree: {degrees}")

    shape = tuple(degree + 1 for degree in degrees)
    arr = np.empty(shape, dtype=object)
    arr.fill(F(0))
    for monomial, coefficient in poly.terms():
        arr[monomial] = as_fraction(coefficient)

    # Tensor-product conversion, one coordinate at a time:
    # x^j = sum_{i=j}^d C(i,j)/C(d,j) B_{i,d}(x).
    for axis, degree in enumerate(degrees):
        if degree == 0:
            continue
        out = np.empty_like(arr)
        other_axes = [j for j in range(arr.ndim) if j != axis]
        other_shape = [arr.shape[j] for j in other_axes]
        denominators = [comb(degree, j) for j in range(degree + 1)]
        for other in itertools.product(*[range(size) for size in other_shape]):
            index: list[Any] = [slice(None)] * arr.ndim
            for j, value in zip(other_axes, other):
                index[j] = value
            line = list(arr[tuple(index)])
            bernstein_line: list[F] = []
            for i in range(degree + 1):
                value = F(0)
                for j in range(i + 1):
                    value += line[j] * F(comb(i, j), denominators[j])
                bernstein_line.append(value)
            out[tuple(index)] = bernstein_line
        arr = out
    return arr, degrees


def split_axis(arr: np.ndarray, axis: int) -> tuple[np.ndarray, np.ndarray]:
    """Exact de Casteljau split at 1/2 along one tensor axis."""
    require(0 <= axis < arr.ndim, f"invalid split axis {axis}")
    degree = arr.shape[axis] - 1
    require(degree >= 1, f"cannot split degree-{degree} axis {axis}")
    left = np.empty_like(arr)
    right = np.empty_like(arr)
    other_axes = [j for j in range(arr.ndim) if j != axis]
    other_shape = [arr.shape[j] for j in other_axes]

    for other in itertools.product(*[range(size) for size in other_shape]):
        index: list[Any] = [slice(None)] * arr.ndim
        for j, value in zip(other_axes, other):
            index[j] = value
        line = list(arr[tuple(index)])
        levels = [line]
        for _ in range(1, degree + 1):
            previous = levels[-1]
            levels.append(
                [(previous[i] + previous[i + 1]) / 2
                 for i in range(len(previous) - 1)]
            )
        left_line = [levels[r][0] for r in range(degree + 1)]
        right_line = [levels[degree - i][i] for i in range(degree + 1)]
        left[tuple(index)] = left_line
        right[tuple(index)] = right_line
    return left, right


def choose_axis(arr: np.ndarray) -> int:
    """Choose the coordinate with the largest exact adjacent variation."""
    best_score = F(-1)
    best_axis: int | None = None
    for axis in range(arr.ndim):
        if arr.shape[axis] <= 1:
            continue
        differences = np.diff(arr, axis=axis)
        score = max(abs(value) for value in differences.flat)
        if score > best_score:
            best_score = score
            best_axis = axis
    require(best_axis is not None, "no positive-degree axis available for subdivision")
    return int(best_axis)


def certify_positive(
    arr: np.ndarray,
    *,
    case: str,
    max_allowed_depth: int = 40,
) -> dict[str, object]:
    """Prove positivity by adaptive exact Bernstein subdivision."""
    stack: list[tuple[np.ndarray, int, str]] = [(arr, 0, "root")]
    nodes = 0
    leaves = 0
    max_depth = 0
    smallest_leaf_coefficient: F | None = None

    while stack:
        box, depth, path = stack.pop()
        nodes += 1
        max_depth = max(max_depth, depth)
        values = list(box.flat)
        require(values, f"{case}: empty Bernstein box at {path}")
        minimum = min(values)
        maximum = max(values)
        if minimum > 0:
            leaves += 1
            if (smallest_leaf_coefficient is None
                    or minimum < smallest_leaf_coefficient):
                smallest_leaf_coefficient = minimum
            continue
        if maximum <= 0:
            raise VerificationError(
                f"{case}: nonpositive Bernstein box at {path}; "
                f"depth={depth}, min={minimum}, max={maximum}"
            )
        if depth >= max_allowed_depth:
            raise VerificationError(
                f"{case}: subdivision depth exceeded at {path}; "
                f"depth={depth}, min={minimum}, max={maximum}"
            )
        axis = choose_axis(box)
        left, right = split_axis(box, axis)
        stack.append((right, depth + 1, f"{path}/x{axis}R"))
        stack.append((left, depth + 1, f"{path}/x{axis}L"))

    require(smallest_leaf_coefficient is not None, f"{case}: no positive leaf")
    return {
        "nodes": nodes,
        "leaves": leaves,
        "depth": max_depth,
        "minimum": smallest_leaf_coefficient,
    }


# ---------------------------------------------------------------------------
# Two-sided Hall blocks
# ---------------------------------------------------------------------------


def sharp_block(a: int, b: int) -> sp.Rational:
    p = 6 - a - b
    require(p >= 1, f"invalid Hall type {(a, b, p)}")
    return v(6 - a) * v(6 - b) / v(p)


def entropy_coefficients(a: int, b: int) -> tuple[sp.Rational, sp.Rational]:
    p = 6 - a - b
    m = 6 - a
    ell = 6 - b
    if p == 1:
        return (
            sp.Rational(1, 2) + sp.Rational(1, 3 * m),
            sp.Rational(1, 2) + sp.Rational(1, 3 * ell),
        )
    coefficient = sp.Rational(p * p, 2 * (a * a + b * b))
    return coefficient, coefficient


def side_expression(
    complement_size: int,
    group_size: int,
    support_size: int,
    excess: sp.Expr,
    spread: sp.Symbol,
    denominator: sp.Expr,
    b_numerator: sp.Expr,
    a0: sp.Expr,
    d: int,
    zeta: int,
    entropy_coefficient: sp.Rational,
) -> sp.Expr:
    """Return (p-t) times one two-level side functional."""
    require(1 <= support_size <= group_size, "invalid two-level support size")
    total = group_size - excess
    prefactor = (1 + excess / complement_size) ** complement_size

    if support_size == group_size:
        high = total / group_size
        product_value = prefactor * high**group_size
        product_over_high = prefactor * high ** (group_size - 1)
        entropy = sp.Rational(0)
    else:
        high = (
            total / group_size
            * (1 + spread * sp.Rational(group_size - support_size,
                                        support_size))
        )
        low = total / group_size * (1 - spread)
        product_value = (
            prefactor
            * high**support_size
            * low ** (group_size - support_size)
        )
        product_over_high = (
            prefactor
            * high ** (support_size - 1)
            * low ** (group_size - support_size)
        )
        entropy = sp.Rational(1, support_size) - sp.Rational(1, group_size)

    return sp.expand(
        product_value * b_numerator
        - denominator * product_over_high
        - entropy_coefficient * a0 * (d + zeta * denominator) * entropy
    )


def two_sided_family(
    a: int, b: int, zeta: int
) -> tuple[sp.Poly, dict[int, sp.Poly], dict[int, sp.Poly]]:
    """Construct the common and one-side pieces for one Hall block type."""
    p = 6 - a - b
    d = min(a, b)
    m = 6 - a
    ell = 6 - b
    require(a >= 1 and b >= a and p >= 1, f"invalid two-sided type {(a, b, p)}")

    total_excess = 2 * E * X
    row_excess = total_excess - E + Z * (2 * E - total_excess)
    column_excess = E - Z * (2 * E - total_excess)
    denominator = p - total_excess

    block_constant = sharp_block(a, b)
    a0 = block_constant * denominator**6 / p**6
    c_x, c_y = entropy_coefficients(a, b)
    b_numerator = p + d + zeta * denominator

    common = (
        (total_excess + d) * (GAMMA - 2)
        - d * a0
        + zeta * denominator * (GAMMA - 2 - a0)
    )
    common_poly = sp.Poly(sp.expand(common), *VARS, domain=sp.QQ)

    row_polys: dict[int, sp.Poly] = {}
    for k in range(1, m + 1):
        expr = side_expression(
            a, m, k, row_excess, U, denominator, b_numerator,
            a0, d, zeta, c_x,
        )
        row_polys[k] = sp.Poly(expr, *VARS, domain=sp.QQ)

    column_polys: dict[int, sp.Poly] = {}
    for ell_support in range(1, ell + 1):
        expr = side_expression(
            b, ell, ell_support, column_excess, V, denominator, b_numerator,
            a0, d, zeta, c_y,
        )
        column_polys[ell_support] = sp.Poly(expr, *VARS, domain=sp.QQ)

    return common_poly, row_polys, column_polys


EXPECTED_TWO_SIDED = {
    (1, 1): (25, 165, 6, F(711760783, 88593750000)),
    (1, 2): (20, 276, 7, F(9468769, 6480000000)),
    (1, 3): (15, 483, 14, F(58951253, 118125000000)),
    (2, 2): (16, 416, 13, F(13905244249, 32256000000000)),
    (1, 4): (10, 352, 14, F(6542213, 472500000000)),
    (2, 3): (
        12, 558, 15,
        F(206812176444269171, 2388787200000000000000),
    ),
}


def verify_two_sided() -> dict[str, Any]:
    total_cases = 0
    output: dict[str, Any] = {}
    for a, b in EXPECTED_TWO_SIDED:
        p = 6 - a - b
        zeta = 14 if p >= 2 else 8
        common, rows, columns = two_sided_family(a, b, zeta)
        m, ell = 6 - a, 6 - b
        require(len(rows) == m and len(columns) == ell, f"type {(a,b)} side count mismatch")

        cases = 0
        nodes = 0
        max_depth = 0
        minimum: F | None = None
        leaves = 0
        degree_set: set[tuple[int, ...]] = set()
        for row_support, row in rows.items():
            for column_support, column in columns.items():
                case = (
                    f"n=6 type=({a},{b}), row_support={row_support}, "
                    f"column_support={column_support}"
                )
                polynomial = -(common + row + column)
                expected_degree = (
                    7,
                    6,
                    0 if row_support == m else m,
                    0 if column_support == ell else ell,
                )
                actual_degree = tuple(polynomial.degree(var) for var in VARS)
                require(
                    actual_degree == expected_degree,
                    f"{case}: multidegree={actual_degree}, expected={expected_degree}",
                )
                degree_set.add(actual_degree)
                bernstein, converted_degree = power_to_bernstein(polynomial)
                require(converted_degree == expected_degree, f"{case}: conversion degree changed")
                result = certify_positive(bernstein, case=case)
                cases += 1
                nodes += int(result["nodes"])
                leaves += int(result["leaves"])
                max_depth = max(max_depth, int(result["depth"]))
                coefficient = result["minimum"]
                require(isinstance(coefficient, F), f"{case}: non-Fraction minimum")
                if minimum is None or coefficient < minimum:
                    minimum = coefficient

        require(minimum is not None, f"type {(a,b)}: no minimum coefficient")
        actual = (cases, nodes, max_depth, minimum)
        expected = EXPECTED_TWO_SIDED[(a, b)]
        require(actual == expected, f"type {(a,b)}: actual={actual}, expected={expected}")
        total_cases += cases
        output[f"{a},{b}"] = {
            "polynomials": cases,
            "nodes": nodes,
            "leaves": leaves,
            "depth": max_depth,
            "minimum": str(minimum),
            "multidegrees": [list(degree) for degree in sorted(degree_set)],
        }

    require(total_cases == 98, f"two-sided total={total_cases}, expected=98")
    return output


# ---------------------------------------------------------------------------
# One-sided Hall blocks
# ---------------------------------------------------------------------------


def univariate_bernstein(
    poly: sp.Poly,
    variable: sp.Symbol,
    lower: sp.Rational,
    upper: sp.Rational,
) -> list[F]:
    require(not poly.is_zero, "zero univariate polynomial")
    require(poly.domain == sp.QQ, f"unexpected univariate domain: {poly.domain}")
    require(lower < upper, f"invalid interval [{lower},{upper}]")
    y = sp.symbols("y")
    transformed = sp.Poly(
        sp.expand(poly.as_expr().subs(variable, lower + (upper - lower) * y)),
        y,
        domain=sp.QQ,
    )
    degree = transformed.degree()
    require(degree >= 0, f"invalid univariate degree: {degree}")
    power = [sp.Rational(0)] * (degree + 1)
    for (j,), coefficient in transformed.terms():
        power[j] = coefficient
    result: list[F] = []
    for i in range(degree + 1):
        value = sp.Rational(0)
        for j in range(i + 1):
            value += power[j] * sp.Rational(comb(i, j), comb(degree, j))
        result.append(as_fraction(value))
    return result


def split_univariate(coefficients: list[F]) -> tuple[list[F], list[F]]:
    degree = len(coefficients) - 1
    require(degree >= 1, "cannot split a constant Bernstein polynomial")
    levels = [coefficients]
    for _ in range(1, degree + 1):
        previous = levels[-1]
        levels.append(
            [(previous[i] + previous[i + 1]) / 2
             for i in range(len(previous) - 1)]
        )
    left = [levels[r][0] for r in range(degree + 1)]
    right = [levels[degree - i][i] for i in range(degree + 1)]
    return left, right


def certify_univariate(
    coefficients: list[F],
    *,
    case: str,
    max_allowed_depth: int = 20,
) -> tuple[int, int, int, F]:
    require(coefficients, f"{case}: empty coefficient list")
    stack: list[tuple[list[F], int, str]] = [(coefficients, 0, "root")]
    nodes = leaves = max_depth = 0
    minimum: F | None = None
    while stack:
        current, depth, path = stack.pop()
        nodes += 1
        max_depth = max(max_depth, depth)
        local_minimum = min(current)
        local_maximum = max(current)
        if local_minimum > 0:
            leaves += 1
            if minimum is None or local_minimum < minimum:
                minimum = local_minimum
            continue
        if local_maximum <= 0:
            raise VerificationError(
                f"{case}: nonpositive interval at {path}; "
                f"depth={depth}, min={local_minimum}, max={local_maximum}"
            )
        if depth >= max_allowed_depth:
            raise VerificationError(
                f"{case}: depth limit at {path}; "
                f"depth={depth}, min={local_minimum}, max={local_maximum}"
            )
        left, right = split_univariate(current)
        stack.extend([
            (right, depth + 1, f"{path}/R"),
            (left, depth + 1, f"{path}/L"),
        ])
    require(minimum is not None, f"{case}: no positive leaf")
    return nodes, leaves, max_depth, minimum


EXPECTED_ONE_SIDED_SMALL = {
    1: (3, 2, 1, F(995279, 28125000000)),
    2: (5, 3, 2, F(2374594200709185001, 36771973056000000000000)),
    3: (7, 4, 3, F(506999743896537994129, 8618431185000000000000000)),
    4: (
        9, 5, 4,
        F(1880953471440335502781485790459,
          40164800436633600000000000000000000),
    ),
    5: (
        11, 6, 5,
        F(60824972467849735194923936599,
          2510300027289600000000000000000000),
    ),
}
EXPECTED_ONE_SIDED_LARGE = (1, 1, 0, F(1, 864))


def verify_one_sided() -> dict[str, Any]:
    tau = sp.symbols("tau")
    energy_constant = (
        sp.Rational(3, 10) * (1 - GAMMA) ** 2 * GAMMA
        / (2 * (KAPPA - GAMMA))
    )
    expected_energy = sp.Rational(39750390625, 5253139008)
    require(energy_constant == expected_energy, f"energy constant={energy_constant}")

    output: dict[str, Any] = {}
    for b in range(1, 6):
        p = 6 - b
        q = 1 - tau
        column_envelope = (1 + sp.Rational(p, b) * tau) ** b * q**p

        small = sp.Poly(
            sp.expand(
                q**2 * (1 - column_envelope + KAPPA * q**6 - GAMMA)
                + energy_constant * tau**2
            ),
            tau,
            domain=sp.QQ,
        )
        require(small.degree() == 8, f"one-sided b={b}: small degree={small.degree()}")
        small_coefficients = univariate_bernstein(
            small, tau, sp.Rational(0), E / p
        )
        small_result = certify_univariate(
            small_coefficients, case=f"n=6 one-sided b={b} small"
        )
        require(
            small_result == EXPECTED_ONE_SIDED_SMALL[b],
            f"one-sided b={b} small: actual={small_result}, "
            f"expected={EXPECTED_ONE_SIDED_SMALL[b]}",
        )

        large_expression = sp.expand(
            1 - column_envelope
            + sp.Rational(3, 32) * tau * q**2
            + GAMMA * (q**6 - 1)
        )
        quotient = sp.cancel(large_expression / tau)
        large = sp.Poly(quotient, tau, domain=sp.QQ)
        require(large.degree() == 5, f"one-sided b={b}: large degree={large.degree()}")
        large_coefficients = univariate_bernstein(
            large, tau, sp.Rational(0), E / p
        )
        large_result = certify_univariate(
            large_coefficients, case=f"n=6 one-sided b={b} large"
        )
        require(
            large_result == EXPECTED_ONE_SIDED_LARGE,
            f"one-sided b={b} large: actual={large_result}, "
            f"expected={EXPECTED_ONE_SIDED_LARGE}",
        )

        output[str(b)] = {
            "small": {
                "degree": 8,
                "nodes": small_result[0],
                "leaves": small_result[1],
                "depth": small_result[2],
                "minimum": str(small_result[3]),
            },
            "large": {
                "degree": 5,
                "nodes": large_result[0],
                "leaves": large_result[1],
                "depth": large_result[2],
                "minimum": str(large_result[3]),
            },
        }
    require(len(output) == 5, f"one-sided family count={len(output)}, expected=5")
    return output


def verify_all() -> dict[str, Any]:
    return {
        "certificate_version": CERTIFICATE_VERSION,
        "manuscript_sha256": MANUSCRIPT_SHA256,
        "constants": verify_constants(),
        "two_sided": verify_two_sided(),
        "one_sided": verify_one_sided(),
        "counts": {
            "four_variable_polynomials": 98,
            "univariate_polynomials": 10,
        },
        "status": "passed",
    }


def print_human(result: dict[str, Any]) -> None:
    for pair, stats in result["two_sided"].items():
        print(
            f"two-sided ({pair}): {stats['polynomials']} polynomials, "
            f"{stats['nodes']} Bernstein nodes, depth {stats['depth']}, "
            f"minimum={stats['minimum']}"
        )
    for b, stats in result["one_sided"].items():
        print(
            f"one-sided b={b}: small depth {stats['small']['depth']}, "
            f"large depth {stats['large']['depth']}"
        )
    print("All exact n=6 Dittert certificates passed (98 + 10 polynomials).")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true", help="print verification result as JSON")
    args = parser.parse_args()
    result = verify_all()
    if args.json:
        print(json.dumps(result, indent=2, sort_keys=True))
    else:
        print_human(result)


if __name__ == "__main__":
    main()
