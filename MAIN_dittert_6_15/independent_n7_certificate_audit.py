#!/usr/bin/env python3
"""Independent exact reconstruction of the three dimension-seven trees.

This checker rebuilds the three bivariate polynomials directly with SymPy,
compares their exact sparse coefficients with the author generator, and then
uses a separately implemented tensor-Bernstein conversion and subdivision
engine.  No correctness check relies on ``assert``.
"""
from __future__ import annotations

import importlib.util
from fractions import Fraction as F
from math import comb, factorial
from pathlib import Path
from typing import Dict, Tuple

import sympy as sp

HERE = Path(__file__).resolve().parent


class VerificationError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def load_author_module():
    spec = importlib.util.spec_from_file_location(
        "n7sub", HERE / "dittert_n7_n15_certificates.py"
    )
    require(spec is not None and spec.loader is not None, "cannot load author verifier")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def v(k: int) -> sp.Rational:
    return sp.Rational(1) if k == 0 else sp.Rational(factorial(k), k**k)


def direct_bern(
    poly: Dict[Tuple[int, int], F],
    x0: F = F(0),
    x1: F = F(3, 20),
    y0: F = F(0),
    y1: F = F(3, 20),
):
    dx = max((i for i, _ in poly), default=0)
    dy = max((j for _, j in poly), default=0)
    local: Dict[Tuple[int, int], F] = {}
    hx, hy = x1 - x0, y1 - y0
    for (i, j), coefficient in poly.items():
        for r in range(i + 1):
            cx = F(comb(i, r)) * x0 ** (i - r) * hx**r
            for s in range(j + 1):
                cy = F(comb(j, s)) * y0 ** (j - s) * hy**s
                local[(r, s)] = local.get((r, s), F(0)) + coefficient * cx * cy
    grid = {(k, ell): F(0) for k in range(dx + 1) for ell in range(dy + 1)}
    for (r, s), coefficient in local.items():
        for k in range(r, dx + 1):
            fx = F(comb(k, r), comb(dx, r))
            for ell in range(s, dy + 1):
                grid[(k, ell)] += coefficient * fx * F(comb(ell, s), comb(dy, s))
    return grid, (dx, dy)


def split_curve(line):
    degree = len(line) - 1
    levels = [line]
    for _ in range(degree):
        previous = levels[-1]
        levels.append([(previous[i] + previous[i + 1]) / 2 for i in range(len(previous) - 1)])
    return (
        [levels[r][0] for r in range(degree + 1)],
        [levels[degree - i][i] for i in range(degree + 1)],
    )


def split(grid, degree, axis):
    dx, dy = degree
    left, right = {}, {}
    if axis == 0:
        for ell in range(dy + 1):
            first, second = split_curve([grid[(k, ell)] for k in range(dx + 1)])
            for k in range(dx + 1):
                left[(k, ell)] = first[k]
                right[(k, ell)] = second[k]
    else:
        for k in range(dx + 1):
            first, second = split_curve([grid[(k, ell)] for ell in range(dy + 1)])
            for ell in range(dy + 1):
                left[(k, ell)] = first[ell]
                right[(k, ell)] = second[ell]
    return left, right


def certify(pair, polynomials, max_total_depth: int = 12):
    grids, degrees = [], []
    for poly in polynomials:
        grid, degree = direct_bern(poly)
        grids.append(grid)
        degrees.append(degree)
    queue = [(tuple(grids), tuple(degrees), 0, 0, "root")]
    head = 0
    nodes = direct = entropy = max_depth = 0
    while head < len(queue):
        current, current_degrees, xdepth, ydepth, path = queue[head]
        head += 1
        nodes += 1
        max_depth = max(max_depth, xdepth + ydepth)
        if all(value > 0 for value in current[0].values()):
            direct += 1
            continue
        if all(value > 0 for value in current[1].values()) and all(
            value > 0 for value in current[2].values()
        ):
            entropy += 1
            continue
        depth = xdepth + ydepth
        require(
            depth < max_total_depth,
            f"depth limit exceeded for {pair} at {path}: depth={depth}",
        )
        axis = 0 if xdepth <= ydepth else 1
        halves = [split(grid, degree, axis) for grid, degree in zip(current, current_degrees)]
        queue.append(
            (
                tuple(half[0] for half in halves),
                current_degrees,
                xdepth + (axis == 0),
                ydepth + (axis == 1),
                path + ("/xL" if axis == 0 else "/yL"),
            )
        )
        queue.append(
            (
                tuple(half[1] for half in halves),
                current_degrees,
                xdepth + (axis == 0),
                ydepth + (axis == 1),
                path + ("/xR" if axis == 0 else "/yR"),
            )
        )
    return nodes, direct, entropy, max_depth


def main() -> None:
    author = load_author_module()
    e, f = sp.symbols("e f")
    gamma7 = sp.Rational(factorial(7), 7**7)

    def sharp(a, b):
        return v(7 - a) * v(7 - b) / v(1)

    def envelope(k, variable):
        return (
            (1 + variable / sp.Rational(k)) ** k
            * (1 - variable / sp.Rational(7 - k)) ** (7 - k)
        )

    def sparse(poly):
        expanded = sp.Poly(sp.expand(poly), e, f, domain=sp.QQ)
        return {
            (i, j): F(int(coefficient.p), int(coefficient.q))
            for (i, j), coefficient in expanded.terms()
            if coefficient
        }

    polynomials_by_pair = {}
    for a, b in ((1, 5), (2, 4), (3, 3)):
        m, ell, d = 7 - a, 7 - b, min(a, b)
        pa, pb = envelope(a, e), envelope(b, f)
        deficit = 2 - pa - pb
        k_value = gamma7 - deficit
        a0 = sharp(a, b) * (1 - e - f) ** 7
        delta = (m - e) * (ell - f) * (1 - e - f)
        b_value = (
            pa * e / (m - e)
            + pb * f / (ell - f)
            - ((e + f + d) * k_value - d * a0) / (1 - e - f)
        )
        numerator = sp.cancel(delta * b_value)
        f_value = 2 * a0 * numerator**2 - 3 * (k_value - a0) ** 2 * delta**2
        own = (sparse(a0 - k_value), sparse(numerator), sparse(f_value))
        submitted = author.single_bridge_polynomials(a, b)
        require(own == submitted, f"symbolic coefficient mismatch for {(a, b)}")
        polynomials_by_pair[(a, b)] = own
        print(
            "symbolic polynomials matched",
            a,
            b,
            [max(i for i, _ in poly) for poly in own],
            [max(j for _, j in poly) for poly in own],
        )

    for pair, expected in author.EXPECTED_SINGLE_BRIDGE.items():
        actual = certify(pair, polynomials_by_pair[pair])
        require(actual == expected, f"tree mismatch for {pair}: actual={actual}, expected={expected}")
        print("independent tree matched", pair, actual)
    print("ALL N7 INDEPENDENT CHECKS PASSED")


if __name__ == "__main__":
    main()
