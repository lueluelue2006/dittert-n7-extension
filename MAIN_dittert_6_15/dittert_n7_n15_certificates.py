#!/usr/bin/env python3
"""Exact rational certificates for Dittert dimensions 7 through 15.

Every correctness decision uses Python integers or ``fractions.Fraction``.
No validation is implemented with ``assert``; consequently ``python -O``
performs exactly the same certificate checks as an ordinary run.  Importing
this module has no verification side effects.
"""
from __future__ import annotations

import argparse
import json
from collections import deque
from fractions import Fraction as Q
from math import comb, factorial
from typing import Any, Dict, List, Tuple

CERTIFICATE_VERSION = "2.0.0"
MANUSCRIPT_SHA256 = "d3d77fee2f6c164b3a2ada4b37c66384a53688c03fbc03aee7effc15d281ba94"


class VerificationError(RuntimeError):
    """Raised when an exact certificate check fails."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def gamma(n: int) -> Q:
    return Q(factorial(n), n**n)


def G(k: int) -> Q:
    return Q(1) if k <= 1 else Q((k - 1) ** (k - 1), k ** (k - 1))


def v(k: int) -> Q:
    return Q(1) if k == 0 else Q(factorial(k), k**k)


def kappa(n: int) -> Q:
    return v(n - 1) * G(n - 1)


def sharp_block(n: int, a: int, b: int) -> Q:
    p = n - a - b
    require(p >= 1, f"invalid Hall type n={n}, a={a}, b={b}")
    return v(n - a) * v(n - b) / v(p)


def weak_block(n: int, a: int, b: int) -> Q:
    return v(n - a) * G(n - b) ** a


def M1(n: int, lower: Q, c: Q) -> Q:
    g = gamma(n)
    return lower - g - Q(n * n, 4) * c * lower * lower / (1 - g)


def M2(n: int, lower: Q, c: Q, lam: Q) -> Q:
    g = gamma(n)
    denominator = lam * lower + (1 - g) / c
    return lower - g - Q(n * n, 4) * lower * lower / denominator


ETA = {
    7: Q(1, 47),
    8: Q(1, 80),
    9: Q(1, 135),
    10: Q(1, 230),
    11: Q(1, 350),
    12: Q(1, 650),
    13: Q(1, 1100),
    14: Q(1, 1800),
}


def alpha(n: int, k: int) -> Q:
    if k == 0:
        return Q(0)
    p = Q(k, n)
    if 2 * k < n:
        return 2 * n * (p + ETA[n]) * (1 - p - ETA[n])
    return Q(2 * k * (n - k), n)


EXPECTED_LOCALIZATION = {
    7: Q(23263, 1808073127),
    8: Q(4757, 836844800),
    9: Q(242489, 87087962025),
    10: Q(31109, 41313127850),
    11: Q(63007753811, 34945789841847500),
    12: Q(11651173, 90828753405000),
    13: Q(13072911571453, 366471344281458130000),
    14: Q(71438752181, 2511343447066440000),
}
EXPECTED_HALF_GAPS = {
    7: Q(33, 658),
    8: Q(9, 80),
    9: Q(13, 270),
    10: Q(11, 115),
    11: Q(82, 1925),
    12: Q(319, 3900),
    13: Q(537, 14300),
    14: Q(893, 12600),
}


# ---------------------------------------------------------------------------
# A small exact bivariate polynomial/Bernstein package
# ---------------------------------------------------------------------------

Monomial = Tuple[int, int]
Poly = Dict[Monomial, Q]
BernsteinGrid = List[List[Q]]


def clean(poly: Poly) -> Poly:
    return {monomial: coefficient for monomial, coefficient in poly.items() if coefficient}


def pconst(coefficient: Q | int) -> Poly:
    coefficient = Q(coefficient)
    return {} if coefficient == 0 else {(0, 0): coefficient}


def padd(left: Poly, right: Poly) -> Poly:
    result = dict(left)
    for monomial, coefficient in right.items():
        result[monomial] = result.get(monomial, Q(0)) + coefficient
    return clean(result)


def pneg(poly: Poly) -> Poly:
    return {monomial: -coefficient for monomial, coefficient in poly.items()}


def psub(left: Poly, right: Poly) -> Poly:
    return padd(left, pneg(right))


def pscale(poly: Poly, scalar: Q | int) -> Poly:
    scalar = Q(scalar)
    return clean({monomial: scalar * coefficient for monomial, coefficient in poly.items()})


def pmul(left: Poly, right: Poly) -> Poly:
    result: Poly = {}
    for (i, j), a in left.items():
        for (k, ell), b in right.items():
            monomial = (i + k, j + ell)
            result[monomial] = result.get(monomial, Q(0)) + a * b
    return clean(result)


def ppow(poly: Poly, exponent: int) -> Poly:
    if exponent < 0:
        raise VerificationError(f"negative polynomial exponent: {exponent}")
    result = pconst(1)
    base = poly
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = pmul(result, base)
        base = pmul(base, base)
        remaining >>= 1
    return result


E_VAR: Poly = {(1, 0): Q(1)}
F_VAR: Poly = {(0, 1): Q(1)}
ONE = pconst(1)
GAMMA7 = gamma(7)


def affine(c0: Q | int, ce: Q | int = 0, cf: Q | int = 0) -> Poly:
    return padd(pconst(c0), padd(pscale(E_VAR, ce), pscale(F_VAR, cf)))


def product_envelope(k: int, variable: Poly) -> Poly:
    return pmul(
        ppow(padd(ONE, pscale(variable, Q(1, k))), k),
        ppow(psub(ONE, pscale(variable, Q(1, 7 - k))), 7 - k),
    )


def degree_pair(poly: Poly) -> Tuple[int, int]:
    return (
        max((i for i, _ in poly), default=0),
        max((j for _, j in poly), default=0),
    )


def power_to_bernstein_on_box(
    poly: Poly, x0: Q, x1: Q, y0: Q, y1: Q
) -> BernsteinGrid:
    """Return exact tensor-product Bernstein coefficients on a rectangle."""
    require(x0 < x1 and y0 < y1, "Bernstein box must have positive side lengths")
    dx, dy = degree_pair(poly)
    hx, hy = x1 - x0, y1 - y0
    local: Poly = {}
    for (i, j), coefficient in poly.items():
        for r in range(i + 1):
            cx = Q(comb(i, r)) * x0 ** (i - r) * hx**r
            for s in range(j + 1):
                cy = Q(comb(j, s)) * y0 ** (j - s) * hy**s
                local[(r, s)] = local.get((r, s), Q(0)) + coefficient * cx * cy

    grid: BernsteinGrid = [[Q(0) for _ in range(dy + 1)] for _ in range(dx + 1)]
    for k in range(dx + 1):
        for ell in range(dy + 1):
            value = Q(0)
            for r in range(k + 1):
                xr = Q(comb(k, r), comb(dx, r))
                for s in range(ell + 1):
                    coefficient = local.get((r, s), Q(0))
                    if coefficient:
                        value += coefficient * xr * Q(comb(ell, s), comb(dy, s))
            grid[k][ell] = value
    return grid


def split_curve_half(control: List[Q]) -> Tuple[List[Q], List[Q]]:
    degree = len(control) - 1
    temp = list(control)
    left = [Q(0)] * (degree + 1)
    right = [Q(0)] * (degree + 1)
    left[0] = temp[0]
    right[degree] = temp[degree]
    for level in range(1, degree + 1):
        temp = [(temp[i] + temp[i + 1]) / 2 for i in range(degree - level + 1)]
        left[level] = temp[0]
        right[degree - level] = temp[-1]
    return left, right


def split_x(grid: BernsteinGrid) -> Tuple[BernsteinGrid, BernsteinGrid]:
    dx, dy = len(grid) - 1, len(grid[0]) - 1
    left = [[Q(0) for _ in range(dy + 1)] for _ in range(dx + 1)]
    right = [[Q(0) for _ in range(dy + 1)] for _ in range(dx + 1)]
    for ell in range(dy + 1):
        left_column, right_column = split_curve_half([grid[k][ell] for k in range(dx + 1)])
        for k in range(dx + 1):
            left[k][ell] = left_column[k]
            right[k][ell] = right_column[k]
    return left, right


def split_y(grid: BernsteinGrid) -> Tuple[BernsteinGrid, BernsteinGrid]:
    dx, dy = len(grid) - 1, len(grid[0]) - 1
    low = [[Q(0) for _ in range(dy + 1)] for _ in range(dx + 1)]
    high = [[Q(0) for _ in range(dy + 1)] for _ in range(dx + 1)]
    for k in range(dx + 1):
        low[k], high[k] = split_curve_half(grid[k])
    return low, high


def strictly_positive(grid: BernsteinGrid) -> bool:
    return all(coefficient > 0 for row in grid for coefficient in row)


def grid_minimum(grid: BernsteinGrid) -> Q:
    return min(coefficient for row in grid for coefficient in row)


def large_entropy_polynomial() -> Poly:
    """The dimension-seven large-entropy polynomial; no mutable global state."""
    one_minus_t = psub(ONE, E_VAR)
    quadratic = padd(pconst(7), padd(pscale(E_VAR, -21), pscale(ppow(E_VAR, 2), 35)))
    return psub(
        pscale(ppow(one_minus_t, 6), Q(2, 3) * (1 - GAMMA7) ** 2),
        pscale(ppow(quadratic, 2), GAMMA7),
    )


def single_bridge_polynomials(a: int, b: int) -> Tuple[Poly, Poly, Poly]:
    require(a >= 1 and b >= 1 and a + b == 6, f"invalid n=7 bridge type {(a, b)}")
    m, ell, d = 7 - a, 7 - b, min(a, b)
    pa = product_envelope(a, E_VAR)
    pb = product_envelope(b, F_VAR)
    deficit = psub(pconst(2), padd(pa, pb))
    k_value = psub(pconst(GAMMA7), deficit)
    a0 = pscale(ppow(affine(1, -1, -1), 7), sharp_block(7, a, b))

    m_minus_e = affine(m, -1, 0)
    ell_minus_f = affine(ell, 0, -1)
    one_minus_sum = affine(1, -1, -1)
    t_plus_d = affine(d, 1, 1)
    denominator = pmul(pmul(m_minus_e, ell_minus_f), one_minus_sum)
    first = pmul(pmul(pmul(pa, E_VAR), ell_minus_f), one_minus_sum)
    second = pmul(pmul(pmul(pb, F_VAR), m_minus_e), one_minus_sum)
    third = pmul(
        psub(pmul(t_plus_d, k_value), pscale(a0, d)),
        pmul(m_minus_e, ell_minus_f),
    )
    b_numerator = psub(padd(first, second), third)
    g_poly = psub(a0, k_value)
    f_poly = psub(
        pscale(pmul(a0, ppow(b_numerator, 2)), 2),
        pscale(pmul(ppow(psub(k_value, a0), 2), ppow(denominator, 2)), 3),
    )
    require(degree_pair(g_poly) == (7, 7), f"unexpected G bidegree for {(a, b)}")
    require(degree_pair(b_numerator) == (8, 8), f"unexpected N bidegree for {(a, b)}")
    require(degree_pair(f_poly) == (23, 23), f"unexpected F bidegree for {(a, b)}")
    return g_poly, b_numerator, f_poly


def certify_single_bridge(a: int, b: int, max_total_depth: int = 12) -> Dict[str, Any]:
    h = Q(3, 20)
    initial = tuple(
        power_to_bernstein_on_box(poly, Q(0), h, Q(0), h)
        for poly in single_bridge_polynomials(a, b)
    )
    queue = deque([(initial, 0, 0, "")])
    nodes = direct_leaves = entropy_leaves = max_depth = 0
    minimum_direct: Q | None = None
    minimum_entropy: Q | None = None
    while queue:
        grids, xdepth, ydepth, path = queue.popleft()
        nodes += 1
        depth = xdepth + ydepth
        max_depth = max(max_depth, depth)
        if strictly_positive(grids[0]):
            direct_leaves += 1
            local = grid_minimum(grids[0])
            minimum_direct = local if minimum_direct is None else min(minimum_direct, local)
            continue
        if strictly_positive(grids[1]) and strictly_positive(grids[2]):
            entropy_leaves += 1
            local = min(grid_minimum(grids[1]), grid_minimum(grids[2]))
            minimum_entropy = local if minimum_entropy is None else min(minimum_entropy, local)
            continue
        diagnostics = [
            (grid_minimum(grid), max(coefficient for row in grid for coefficient in row))
            for grid in grids
        ]
        require(
            depth < max_total_depth,
            f"bridge {(a,b)} failed at path={path}, depth={depth}, "
            f"coefficient_ranges={diagnostics}",
        )
        if xdepth <= ydepth:
            halves = [split_x(grid) for grid in grids]
            queue.append((tuple(pair[0] for pair in halves), xdepth + 1, ydepth, path + "xL"))
            queue.append((tuple(pair[1] for pair in halves), xdepth + 1, ydepth, path + "xR"))
        else:
            halves = [split_y(grid) for grid in grids]
            queue.append((tuple(pair[0] for pair in halves), xdepth, ydepth + 1, path + "yL"))
            queue.append((tuple(pair[1] for pair in halves), xdepth, ydepth + 1, path + "yR"))
    require(direct_leaves + entropy_leaves > 0, f"bridge {(a,b)} produced no leaves")
    return {
        "nodes": nodes,
        "direct_leaves": direct_leaves,
        "entropy_leaves": entropy_leaves,
        "max_depth": max_depth,
        "minimum_direct": str(minimum_direct) if minimum_direct is not None else None,
        "minimum_entropy": str(minimum_entropy) if minimum_entropy is not None else None,
    }


EXPECTED_SINGLE_BRIDGE = {
    (1, 5): (39, 10, 10, 8),
    (2, 4): (35, 9, 9, 6),
    (3, 3): (47, 11, 13, 8),
}

T_LOW = {8: Q(1, 7), 9: Q(1, 10), 10: Q(1, 10)}
LAM_LOW = {8: Q(20), 9: Q(27), 10: Q(33)}
EXPECTED_TAU_GAPS = {
    8: Q(7277, 6407093),
    9: Q(746489, 477848900),
    10: Q(994933, 156193300),
}
EXPECTED_D0_GAPS = {
    8: Q(4757, 13075700),
    9: Q(2762489, 477848900),
    10: Q(1278433, 156193300),
}
EXPECTED_GAMMA_81_GAPS = {
    8: Q(105557, 10616832),
    9: Q(54569, 4782969),
    10: Q(1516573, 126562500),
}
EXPECTED_LOW_COUNTS = {8: 12, 9: 16, 10: 20}
EXPECTED_LOW_MINIMA = {
    8: (Q(4139711814785404757611563, 300024713494271812137373270016), 1, 1),
    9: (Q(3087705658660821267020074863425, 540713498338770121894081281002569728), 1, 1),
    10: (Q(36184754863040830496579778359572513127, 18255668299923283084495152136002059873437500), 1, 1),
}
EXPECTED_SMALL = {
    8: Q(165702517175104104520946768451687, 215701290564174919711661340674228224),
    9: Q(5586975164420125257678146526492250330079, 16534930139652222569055999993919229392846848),
    10: Q(55944735805258396002952464523577322487640863, 404273588897859606879181224832534790039062500000),
}
EXPECTED_LARGE = {
    8: Q(405685268407642329, 2147483648000000000),
    9: Q(17233957123121, 73811250000000),
    10: Q(7399058328856214668089, 30517578125000000000000),
}
EXPECTED_MID_COUNTS = {11: 35, 12: 42, 13: 49, 14: 56}
EXPECTED_MID_MINIMA = {
    11: (Q(1260916230510486033100687951530284366056569069, 6568905011778983402353515625000000000000000000000000), 0, 9),
    12: (Q(1296568134644485856101383130605756547126510510475, 52823701911706625804697555859966244866493587958274637824), 0, 11),
    13: (Q(2344634926829311765047551592181776557728012393742935099475, 64121715613888291470979410224708337217950361958912179251651280896), 0, 12),
    14: (Q(4746695979316977001752598428103553549748722442212199016699590578325, 277886267525954703535897576147309905022601205806915212282941432903903914488), 0, 13),
}
EXPECTED_N15_INTEGERS = (
    45591579980859375,
    81568443603515625,
    290584966207055614607247802368,
    3939835293455120947239518208,
    18198910159813162803331251792701065161767,
    162574843672041436633515011954046933353332,
)


def verify_localization() -> Dict[str, Any]:
    rows: Dict[str, Any] = {}
    for n in range(7, 15):
        g = gamma(n)
        localization = ETA[n] ** 2 - g / (2 * n * (1 - g))
        half_gap = Q(1, 2) - Q((n - 1) // 2, n) - ETA[n]
        require(localization == EXPECTED_LOCALIZATION[n] and localization > 0, f"n={n} localization mismatch")
        require(half_gap == EXPECTED_HALF_GAPS[n] and half_gap > 0, f"n={n} half-gap mismatch")
        for k in range(n):
            require(alpha(n, k) <= Q(n, 2), f"n={n}, k={k}: alpha exceeds n/2")
        rows[str(n)] = {"localization_gap": str(localization), "half_gap": str(half_gap)}
    return rows


def verify_n7() -> Dict[str, Any]:
    n = 7
    g = GAMMA7
    t7 = Q(5, 24)
    lam7 = Q(329, 24)
    require(t7 * t7 - n * g / (1 - g) == Q(20185, 67351104), "n=7 ordinary tau gap mismatch")
    require(lam7 == comb(n, 2) - comb(n, 3) * t7, "n=7 ordinary lambda mismatch")
    ordinary = []
    for a in range(1, n):
        for b in range(a, n - a):
            p = n - a - b
            if p >= 2:
                c = (alpha(n, a) + alpha(n, b)) / p**2
                ordinary.append((M2(n, sharp_block(n, a, b), c, lam7), a, b))
    expected_ordinary = (Q(217499955444827479495, 21713873952987016007066496), 1, 1)
    require(len(ordinary) == 6, f"n=7 ordinary count={len(ordinary)}")
    require(all(margin > 0 for margin, _, _ in ordinary), "n=7 nonpositive ordinary margin")
    require(min(ordinary) == expected_ordinary, f"n=7 ordinary minimum={min(ordinary)}")

    tmix = Q(3, 20)
    lammix = Q(63, 4)
    require(lammix == comb(n, 2) - comb(n, 3) * tmix, "n=7 mixed lambda mismatch")
    require(max(alpha(n, k) for k in range(1, n)) == alpha(n, 3), "n=7 alpha maximum mismatch")
    require(tmix**2 - alpha(n, 3) * g / (1 - g) == Q(842770143, 723229250800), "n=7 mixed range mismatch")
    mixed = []
    for a, b in ((1, 5), (2, 4), (3, 3)):
        lower = sharp_block(n, a, b)
        mixed.append((M2(n, lower, alpha(n, b), lammix), a, b, "e<=0"))
        mixed.append((M2(n, lower, alpha(n, a), lammix), a, b, "f<=0"))
    expected_mixed = (Q(165035926385, 3897251304526692), 1, 5, "e<=0")
    require(len(mixed) == 6, "n=7 mixed count mismatch")
    require(all(margin > 0 for margin, *_ in mixed), "n=7 nonpositive mixed margin")
    require(min(mixed) == expected_mixed, f"n=7 mixed minimum={min(mixed)}")

    d0sq = n * g / (2 * (1 - g))
    require(Q(9, 400) - d0sq == Q(44361, 46771600), "n=7 d0 gap mismatch")
    energy = (1 - g) * Q(17, 20) ** 2 / 2
    require(energy - Q(1, 3) == Q(7258243, 282357600), "n=7 energy gap mismatch")
    small = Q(1, 3) * (1 - g) ** 2 * g - Q(n * n, 2) * kappa(n) ** 2
    expected_small = Q(22176489486262543850191055, 20672701805255574480144334848)
    require(small == expected_small and small > 0, "n=7 small-entropy mismatch")

    large_grid = power_to_bernstein_on_box(large_entropy_polynomial(), Q(0), Q(3, 20), Q(0), Q(1))
    large_coefficients = [row[0] for row in large_grid]
    expected_large_min = Q(155130361249919329, 1328763571296000000)
    require(all(coefficient > 0 for coefficient in large_coefficients), "n=7 large polynomial has nonpositive coefficient")
    require(min(large_coefficients) == expected_large_min, "n=7 large minimum mismatch")

    trees: Dict[str, Any] = {}
    for pair, expected in EXPECTED_SINGLE_BRIDGE.items():
        result = certify_single_bridge(*pair)
        actual = (result["nodes"], result["direct_leaves"], result["entropy_leaves"], result["max_depth"])
        require(actual == expected, f"n=7 bridge {pair}: actual={actual}, expected={expected}")
        trees[f"{pair[0]},{pair[1]}"] = result
    return {
        "ordinary_cases": len(ordinary),
        "ordinary_minimum": str(expected_ordinary[0]),
        "mixed_cases": len(mixed),
        "mixed_minimum": str(expected_mixed[0]),
        "large_entropy_minimum": str(expected_large_min),
        "single_bridge_trees": trees,
    }


def verify_n8_n10() -> Dict[str, Any]:
    output: Dict[str, Any] = {}
    for n in (8, 9, 10):
        g = gamma(n)
        require(T_LOW[n] ** 2 - n * g / (1 - g) == EXPECTED_TAU_GAPS[n], f"n={n} tau gap mismatch")
        require(Q(1, 100) - Q(n, 2) * g / (1 - g) == EXPECTED_D0_GAPS[n], f"n={n} d0 gap mismatch")
        require(Q(1, 81) - g == EXPECTED_GAMMA_81_GAPS[n], f"n={n} gamma gap mismatch")
        require(LAM_LOW[n] <= comb(n, 2) - comb(n, 3) * T_LOW[n], f"n={n} lambda invalid")
        records = []
        for a in range(1, n):
            for b in range(a, n - a):
                p = n - a - b
                c = (alpha(n, a) + alpha(n, b)) / p**2
                records.append((M2(n, sharp_block(n, a, b), c, LAM_LOW[n]), a, b))
        require(len(records) == EXPECTED_LOW_COUNTS[n], f"n={n} case count mismatch")
        require(all(margin > 0 for margin, _, _ in records), f"n={n} nonpositive second-order margin")
        require(min(records) == EXPECTED_LOW_MINIMA[n], f"n={n} minimum mismatch")
        small = Q(2, 5) * (1 - g) ** 2 * g - Q(n * n, 2) * kappa(n) ** 2
        large = Q(4, 5) * (1 - g) ** 2 * Q(9, 10) ** n - n * n * g
        require(small == EXPECTED_SMALL[n] and small > 0, f"n={n} small entropy mismatch")
        require(large == EXPECTED_LARGE[n] and large > 0, f"n={n} large entropy mismatch")
        output[str(n)] = {
            "cases": len(records),
            "minimum": str(min(records)[0]),
            "small_entropy": str(small),
            "large_entropy": str(large),
        }
    return output


def verify_n11_n14() -> Dict[str, Any]:
    output: Dict[str, Any] = {}
    for n in (11, 12, 13, 14):
        records = []
        for a in range(n):
            for b in range(a, n - a):
                if n == 11 and (a, b) == (0, 10):
                    continue
                p = n - a - b
                c = (alpha(n, a) + alpha(n, b)) / p**2
                lower = kappa(n) if a == 0 else weak_block(n, a, b)
                records.append((M1(n, lower, c), a, b))
        require(len(records) == EXPECTED_MID_COUNTS[n], f"n={n} mid count mismatch")
        require(all(margin > 0 for margin, _, _ in records), f"n={n} nonpositive first-order margin")
        require(min(records) == EXPECTED_MID_MINIMA[n], f"n={n} mid minimum mismatch")
        output[str(n)] = {"cases": len(records), "minimum": str(min(records)[0])}

    require(Q(1, 1000) - gamma(11) == Q(22308624601, 25937424601000), "n=11 gamma bound mismatch")
    g11 = gamma(11)
    cstar = 1 / ((1 - g11) ** 2 / alpha(11, 1) + 1 / alpha(11, 10))
    expected_cstar = Q(1540632960668022501507138780, 1671236290061169784540151329)
    cstar_gap = Q(93, 100) - cstar
    expected_gap = Q(1361678908886539811520195597, 167123629006116978454015132900)
    special = M1(11, kappa(11), Q(93, 100))
    expected_special = Q(82169554295122875121722420238561072941763350843, 656890501177898340235351562500000000000000000000000000)
    require(cstar == expected_cstar, "n=11 cstar mismatch")
    require(cstar_gap == expected_gap and cstar_gap > 0, "n=11 cstar gap mismatch")
    require(special == expected_special and special > 0, "n=11 special margin mismatch")
    output["11_special"] = {"cstar_gap": str(cstar_gap), "margin": str(special)}
    return output


def verify_n15() -> Dict[str, Any]:
    mu1 = v(14) * G(14)
    mu2 = v(13) * G(13) ** 2
    expected_thin = Q(
        int("93713344067889437339457163890050246146212691912980360828364561109008039"),
        int("34960636059289149791385973204315516624092004582151399011826949685248000000000000"),
    )
    expected_thick = Q(
        int("63429645430340668282196872894311697916976706007718554270064778550284048088388617933319028973541376"),
        int("3693029341533257024549996410753592125285633553926582232246120370444607869526685543159813878074388916015625"),
    )
    thin = M1(15, mu1, Q(465, 49))
    thick = M1(15, mu2, Q(15))
    require(thin == expected_thin and thin > 0, "n=15 thin margin mismatch")
    require(thick == expected_thick and thick > 0, "n=15 thick margin mismatch")
    actual_integers = (
        15**15 - 300000 * factorial(15),
        29863 * 15**15 - 10**10 * factorial(15),
        10**10 * factorial(13) * 13**13 - 29937 * 14**26,
        3 * 14**26 - 10**6 * factorial(13) * 13**13,
        10**8 * factorial(13) * 12**24 - 301 * 13**37,
        4 * 13**37 - 10**6 * factorial(13) * 12**24,
    )
    require(actual_integers == EXPECTED_N15_INTEGERS, "n=15 integer tuple mismatch")
    require(all(value > 0 for value in actual_integers), "n=15 integer gap nonpositive")
    coarse1 = 5 * (10**12 * 299999) - 534 * 9 * 300000 * 10**9
    coarse2 = 14 * (10**12 * 299999) - 3375 * 4 * 300000 * 10**9
    require(coarse1 == 58195000000000000, "n=15 first coarse gap mismatch")
    require(coarse2 == 149986000000000000, "n=15 second coarse gap mismatch")
    return {
        "thin_margin": str(thin),
        "thick_margin": str(thick),
        "integer_gaps": list(actual_integers),
        "coarse_gaps": [coarse1, coarse2],
    }


def verify_all() -> Dict[str, Any]:
    return {
        "certificate_version": CERTIFICATE_VERSION,
        "manuscript_sha256": MANUSCRIPT_SHA256,
        "localization": verify_localization(),
        "n7": verify_n7(),
        "n8_n10": verify_n8_n10(),
        "n11_n14": verify_n11_n14(),
        "n15": verify_n15(),
        "status": "passed",
    }


def print_human(result: Dict[str, Any]) -> None:
    n7 = result["n7"]
    print(f"n=7: {n7['ordinary_cases']} ordinary cases, {n7['mixed_cases']} mixed cases")
    for pair, stats in n7["single_bridge_trees"].items():
        print(
            f"n=7 bridge ({pair}): nodes={stats['nodes']}, "
            f"leaves={stats['direct_leaves']}+{stats['entropy_leaves']}, "
            f"depth={stats['max_depth']}"
        )
    for n, stats in result["n8_n10"].items():
        print(f"n={n}: {stats['cases']} second-order Hall cases")
    for n in ("11", "12", "13", "14"):
        print(f"n={n}: {result['n11_n14'][n]['cases']} first-order Hall cases")
    print("n=15: two exact margins and six exact integer inequalities")
    print("All exact rational and Bernstein certificates for dimensions 7 through 15 passed.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true", help="print the verification result as JSON")
    args = parser.parse_args()
    result = verify_all()
    if args.json:
        print(json.dumps(result, indent=2, sort_keys=True))
    else:
        print_human(result)


if __name__ == "__main__":
    main()
