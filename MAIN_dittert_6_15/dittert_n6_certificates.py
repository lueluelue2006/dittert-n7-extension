#!/usr/bin/env python3
"""Exact rational certificates for the n=6 Dittert proof.

The script verifies the finite inequalities in the accompanying note.  It
uses exact rational arithmetic throughout.  SymPy is used only to expand the
explicit rational polynomials, NumPy only as an object-array container, and
all Bernstein coefficients and de Casteljau subdivisions are fractions.Fraction
objects.  No floating-point comparison occurs.

The two-sided part checks 98 polynomials (six Hall block types and all pairs
of support sizes).  The one-sided part checks ten univariate polynomials.
"""
from __future__ import annotations

from fractions import Fraction as F
from math import comb, factorial
import itertools

import numpy as np
import sympy as sp


# ---------------------------------------------------------------------------
# Exact constants and elementary gaps
# ---------------------------------------------------------------------------

def v(k: int) -> sp.Rational:
    return sp.Rational(1) if k == 0 else sp.Rational(factorial(k), k**k)


GAMMA = sp.Rational(5, 324)
KAPPA = v(5) * sp.Rational(4**4, 5**4)
E = sp.Rational(11, 50)

assert KAPPA == sp.Rational(6144, 390625) > GAMMA
assert sp.Rational(15, 319) < E**2
assert sp.Rational(15, 319) < sp.Rational(7, 32) ** 2
assert (
    (1 - GAMMA) * sp.Rational(25, 32) ** 2 / 2
    - sp.Rational(3, 10)
    == sp.Rational(1547, 3317760)
    > 0
)
assert (
    (1 - GAMMA) ** 2 / 108 - sp.Rational(9, 1024)
    == sp.Rational(33853, 181398528)
    > 0
)


# ---------------------------------------------------------------------------
# Exact tensor-product Bernstein arithmetic on the unit cube
# ---------------------------------------------------------------------------

X, Z, U, V = sp.symbols("x z u v")
VARS = (X, Z, U, V)


def as_fraction(q: sp.Rational) -> F:
    return F(int(q.p), int(q.q))


def power_to_bernstein(poly: sp.Poly) -> tuple[np.ndarray, tuple[int, ...]]:
    """Convert a rational power-basis polynomial on [0,1]^4 to Bernstein form."""
    degrees = tuple(poly.degree(var) for var in VARS)
    shape = tuple(d + 1 for d in degrees)
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
        for other in itertools.product(*[range(s) for s in other_shape]):
            index = [slice(None)] * arr.ndim
            for j, value in zip(other_axes, other):
                index[j] = value
            line = list(arr[tuple(index)])
            bernstein_line = []
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
    degree = arr.shape[axis] - 1
    left = np.empty_like(arr)
    right = np.empty_like(arr)
    other_axes = [j for j in range(arr.ndim) if j != axis]
    other_shape = [arr.shape[j] for j in other_axes]

    for other in itertools.product(*[range(s) for s in other_shape]):
        index = [slice(None)] * arr.ndim
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
    best_axis = 0
    for axis in range(arr.ndim):
        if arr.shape[axis] <= 1:
            continue
        differences = np.diff(arr, axis=axis)
        score = max(abs(value) for value in differences.flat)
        if score > best_score:
            best_score = score
            best_axis = axis
    return best_axis


def certify_positive(arr: np.ndarray) -> dict[str, object]:
    """Prove positivity by adaptive exact Bernstein subdivision."""
    stack: list[tuple[np.ndarray, int]] = [(arr, 0)]
    nodes = 0
    leaves = 0
    max_depth = 0
    smallest_leaf_coefficient: F | None = None

    while stack:
        box, depth = stack.pop()
        nodes += 1
        max_depth = max(max_depth, depth)
        values = list(box.flat)
        minimum = min(values)
        maximum = max(values)
        if minimum > 0:
            leaves += 1
            if (smallest_leaf_coefficient is None
                    or minimum < smallest_leaf_coefficient):
                smallest_leaf_coefficient = minimum
            continue
        if maximum <= 0:
            raise AssertionError("a Bernstein box is nonpositive")
        if depth >= 40:
            raise AssertionError("subdivision depth exceeded")
        axis = choose_axis(box)
        left, right = split_axis(box, axis)
        stack.append((right, depth + 1))
        stack.append((left, depth + 1))

    assert smallest_leaf_coefficient is not None
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
    total = group_size - excess
    prefactor = (1 + excess / complement_size) ** complement_size

    if support_size == group_size:
        high = total / group_size
        product = prefactor * high**group_size
        product_over_high = prefactor * high ** (group_size - 1)
        entropy = sp.Rational(0)
    else:
        high = (
            total / group_size
            * (1 + spread * sp.Rational(group_size - support_size,
                                        support_size))
        )
        low = total / group_size * (1 - spread)
        product = (
            prefactor
            * high**support_size
            * low ** (group_size - support_size)
        )
        product_over_high = (
            prefactor
            * high ** (support_size - 1)
            * low ** (group_size - support_size)
        )
        entropy = (
            sp.Rational(1, support_size)
            - sp.Rational(1, group_size)
        )

    return sp.expand(
        product * b_numerator
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

    row_polys = {}
    for k in range(1, m + 1):
        expr = side_expression(
            a, m, k, row_excess, U, denominator, b_numerator,
            a0, d, zeta, c_x,
        )
        row_polys[k] = sp.Poly(expr, *VARS, domain=sp.QQ)

    column_polys = {}
    for k in range(1, ell + 1):
        expr = side_expression(
            b, ell, k, column_excess, V, denominator, b_numerator,
            a0, d, zeta, c_y,
        )
        column_polys[k] = sp.Poly(expr, *VARS, domain=sp.QQ)

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


def verify_two_sided() -> None:
    total_cases = 0
    for a, b in EXPECTED_TWO_SIDED:
        p = 6 - a - b
        zeta = 14 if p >= 2 else 8
        common, rows, columns = two_sided_family(a, b, zeta)

        cases = 0
        nodes = 0
        max_depth = 0
        minimum: F | None = None
        for k, row in rows.items():
            for ell, column in columns.items():
                polynomial = -(common + row + column)
                bernstein, _ = power_to_bernstein(polynomial)
                result = certify_positive(bernstein)
                cases += 1
                nodes += int(result["nodes"])
                max_depth = max(max_depth, int(result["depth"]))
                coefficient = result["minimum"]
                assert isinstance(coefficient, F)
                if minimum is None or coefficient < minimum:
                    minimum = coefficient

        assert minimum is not None
        actual = (cases, nodes, max_depth, minimum)
        assert actual == EXPECTED_TWO_SIDED[(a, b)], (a, b, actual)
        total_cases += cases
        print(
            f"two-sided {(a, b)}: {cases} polynomials, "
            f"{nodes} Bernstein nodes, depth {max_depth}"
        )

    assert total_cases == 98


# ---------------------------------------------------------------------------
# One-sided Hall blocks
# ---------------------------------------------------------------------------

def univariate_bernstein(
    poly: sp.Poly, variable: sp.Symbol, lower: sp.Rational,
    upper: sp.Rational,
) -> list[F]:
    y = sp.symbols("y")
    transformed = sp.Poly(
        sp.expand(poly.as_expr().subs(variable, lower + (upper - lower) * y)),
        y,
        domain=sp.QQ,
    )
    degree = transformed.degree()
    power = [sp.Rational(0)] * (degree + 1)
    for (j,), coefficient in transformed.terms():
        power[j] = coefficient
    result = []
    for i in range(degree + 1):
        value = sp.Rational(0)
        for j in range(i + 1):
            value += power[j] * sp.Rational(comb(i, j), comb(degree, j))
        result.append(as_fraction(value))
    return result


def split_univariate(coefficients: list[F]) -> tuple[list[F], list[F]]:
    degree = len(coefficients) - 1
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


def certify_univariate(coefficients: list[F]) -> tuple[int, int, int, F]:
    stack = [(coefficients, 0)]
    nodes = leaves = max_depth = 0
    minimum: F | None = None
    while stack:
        current, depth = stack.pop()
        nodes += 1
        max_depth = max(max_depth, depth)
        if min(current) > 0:
            leaves += 1
            local = min(current)
            if minimum is None or local < minimum:
                minimum = local
            continue
        assert max(current) > 0 and depth < 20
        left, right = split_univariate(current)
        stack.extend([(right, depth + 1), (left, depth + 1)])
    assert minimum is not None
    return nodes, leaves, max_depth, minimum


EXPECTED_ONE_SIDED_SMALL = {
    1: (3, 2, 1, F(995279, 28125000000)),
    2: (
        5, 3, 2,
        F(2374594200709185001, 36771973056000000000000),
    ),
    3: (
        7, 4, 3,
        F(506999743896537994129, 8618431185000000000000000),
    ),
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


def verify_one_sided() -> None:
    tau = sp.symbols("tau")
    energy_constant = (
        sp.Rational(3, 10) * (1 - GAMMA) ** 2 * GAMMA
        / (2 * (KAPPA - GAMMA))
    )
    assert energy_constant == sp.Rational(39750390625, 5253139008)

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
        small_coefficients = univariate_bernstein(
            small, tau, sp.Rational(0), E / p
        )
        small_result = certify_univariate(small_coefficients)
        assert small_result == EXPECTED_ONE_SIDED_SMALL[b]

        large_expression = sp.expand(
            1 - column_envelope
            + sp.Rational(3, 32) * tau * q**2
            + GAMMA * (q**6 - 1)
        )
        large = sp.Poly(sp.cancel(large_expression / tau), tau, domain=sp.QQ)
        large_coefficients = univariate_bernstein(
            large, tau, sp.Rational(0), E / p
        )
        large_result = certify_univariate(large_coefficients)
        assert large_result == (1, 1, 0, F(1, 864))

        print(
            f"one-sided b={b}: small depth {small_result[2]}, "
            "large polynomial has globally positive Bernstein coefficients"
        )


if __name__ == "__main__":
    verify_two_sided()
    verify_one_sided()
    print("All exact n=6 Dittert certificates passed (98 + 10 polynomials).")
