#!/usr/bin/env python3
"""Exact certificates for Dittert's conjecture in dimensions 7 through 15.

The script uses only Python integers and fractions.Fraction.  In particular,
there are no floating-point comparisons.  Besides the scalar certificates for
n=7,...,15, it verifies the three two-variable polynomial alternatives needed
for the single-bridge Hall cuts in dimension seven.  Those alternatives are
certified by dyadic subdivision in tensor-product Bernstein form.
"""
from __future__ import annotations

from collections import deque
from fractions import Fraction as Q
from math import comb, factorial
from typing import Dict, Iterable, List, Tuple


# ---------------------------------------------------------------------------
# Scalar constants
# ---------------------------------------------------------------------------

def gamma(n: int) -> Q:
    return Q(factorial(n), n**n)


def G(k: int) -> Q:
    return Q(1) if k <= 1 else Q((k - 1)**(k - 1), k**(k - 1))


def v(k: int) -> Q:
    return Q(1) if k == 0 else Q(factorial(k), k**k)


def kappa(n: int) -> Q:
    return v(n - 1) * G(n - 1)


def sharp_block(n: int, a: int, b: int) -> Q:
    p = n - a - b
    return v(n - a) * v(n - b) / v(p)


def weak_block(n: int, a: int, b: int) -> Q:
    return v(n - a) * G(n - b)**a


def M1(n: int, L: Q, c: Q) -> Q:
    g = gamma(n)
    return L - g - Q(n*n, 4) * c * L*L / (1 - g)


def M2(n: int, L: Q, c: Q, lam: Q) -> Q:
    g = gamma(n)
    q = lam*L + (1-g)/c
    return L - g - Q(n*n, 4) * L*L/q


ETA = {
    7: Q(1, 47),
    8: Q(1, 80), 9: Q(1, 135), 10: Q(1, 230),
    11: Q(1, 350), 12: Q(1, 650),
    13: Q(1, 1100), 14: Q(1, 1800),
}


def alpha(n: int, k: int) -> Q:
    if k == 0:
        return Q(0)
    p = Q(k, n)
    if 2*k < n:
        return 2*n*(p + ETA[n])*(1 - p - ETA[n])
    return Q(2*k*(n-k), n)


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
    8: Q(9, 80), 9: Q(13, 270), 10: Q(11, 115),
    11: Q(82, 1925), 12: Q(319, 3900),
    13: Q(537, 14300), 14: Q(893, 12600),
}

for n in range(7, 15):
    g = gamma(n)
    localization = ETA[n]**2 - g/(2*n*(1-g))
    half_gap = Q(1, 2) - Q((n-1)//2, n) - ETA[n]
    assert localization == EXPECTED_LOCALIZATION[n] > 0
    assert half_gap == EXPECTED_HALF_GAPS[n] > 0
    for k in range(n):
        assert alpha(n, k) <= Q(n, 2)


# ---------------------------------------------------------------------------
# Dimension seven: scalar and one-variable certificates
# ---------------------------------------------------------------------------

n = 7
g = gamma(n)

# Ordinary two-sided Hall pairs, p >= 2.
T7 = Q(5, 24)
LAM7 = Q(329, 24)
assert T7*T7 - n*g/(1-g) == Q(20185, 67351104) > 0
assert LAM7 == comb(n, 2) - comb(n, 3)*T7

records7 = []
for a in range(1, n):
    for b in range(a, n-a):
        p = n-a-b
        if p >= 2:
            c = (alpha(n, a)+alpha(n, b))/p**2
            records7.append((M2(n, sharp_block(n, a, b), c, LAM7), a, b))
assert len(records7) == 6
assert all(m > 0 for m, _, _ in records7)
assert min(records7) == (
    Q(217499955444827479495, 21713873952987016007066496), 1, 1
)

# Single-bridge pairs for which at least one complementary excess is nonpositive.
TMIX7 = Q(3, 20)
LAMMIX7 = Q(63, 4)
assert LAMMIX7 == comb(n, 2) - comb(n, 3)*TMIX7
assert max(alpha(n, k) for k in range(1, n)) == alpha(n, 3)
assert TMIX7**2 - alpha(n, 3)*g/(1-g) == Q(
    842770143, 723229250800
) > 0

mixed7 = []
for a, b in ((1, 5), (2, 4), (3, 3)):
    L = sharp_block(n, a, b)
    # e <= 0 leaves only the b-column excess; f <= 0 leaves only the a-row excess.
    mixed7.append((M2(n, L, alpha(n, b), LAMMIX7), a, b, "e<=0"))
    mixed7.append((M2(n, L, alpha(n, a), LAMMIX7), a, b, "f<=0"))
assert len(mixed7) == 6
assert all(m > 0 for m, *_ in mixed7)
assert min(mixed7) == (
    Q(165035926385, 3897251304526692), 1, 5, "e<=0"
)

# One-sided Hall cuts.
d0sq7 = n*g/(2*(1-g))
assert Q(9, 400) - d0sq7 == Q(44361, 46771600) > 0
energy7 = (1-g)*Q(17, 20)**2/2
assert energy7 - Q(1, 3) == Q(7258243, 282357600) > 0
small7 = Q(1, 3)*(1-g)**2*g - Q(n*n, 2)*kappa(n)**2
assert small7 == Q(
    22176489486262543850191055,
    20672701805255574480144334848,
) > 0


# ---------------------------------------------------------------------------
# A small exact polynomial package for Bernstein certificates
# ---------------------------------------------------------------------------

Monomial = Tuple[int, int]
Poly = Dict[Monomial, Q]
BernsteinGrid = List[List[Q]]


def clean(p: Poly) -> Poly:
    return {m: c for m, c in p.items() if c}


def pconst(c: Q | int) -> Poly:
    c = Q(c)
    return {} if c == 0 else {(0, 0): c}


def padd(p: Poly, q: Poly) -> Poly:
    out = dict(p)
    for m, c in q.items():
        out[m] = out.get(m, Q(0)) + c
    return clean(out)


def pneg(p: Poly) -> Poly:
    return {m: -c for m, c in p.items()}


def psub(p: Poly, q: Poly) -> Poly:
    return padd(p, pneg(q))


def pscale(p: Poly, c: Q | int) -> Poly:
    c = Q(c)
    return clean({m: c*a for m, a in p.items()})


def pmul(p: Poly, q: Poly) -> Poly:
    out: Poly = {}
    for (i, j), a in p.items():
        for (k, ell), b in q.items():
            m = (i+k, j+ell)
            out[m] = out.get(m, Q(0)) + a*b
    return clean(out)


def ppow(p: Poly, exponent: int) -> Poly:
    assert exponent >= 0
    out = pconst(1)
    base = p
    k = exponent
    while k:
        if k & 1:
            out = pmul(out, base)
        base = pmul(base, base)
        k >>= 1
    return out


E_VAR: Poly = {(1, 0): Q(1)}
F_VAR: Poly = {(0, 1): Q(1)}
ONE = pconst(1)


def affine(c0: Q | int, ce: Q | int = 0, cf: Q | int = 0) -> Poly:
    return padd(pconst(c0), padd(pscale(E_VAR, ce), pscale(F_VAR, cf)))


def product_envelope(k: int, variable: Poly) -> Poly:
    return pmul(
        ppow(padd(ONE, pscale(variable, Q(1, k))), k),
        ppow(psub(ONE, pscale(variable, Q(1, 7-k))), 7-k),
    )


def degree_pair(p: Poly) -> Tuple[int, int]:
    return (
        max((i for i, _ in p), default=0),
        max((j for _, j in p), default=0),
    )


def power_to_bernstein_on_box(
    p: Poly, x0: Q, x1: Q, y0: Q, y1: Q
) -> BernsteinGrid:
    """Exact tensor-product Bernstein coefficients on the given rectangle."""
    dx, dy = degree_pair(p)
    hx, hy = x1-x0, y1-y0
    local: Poly = {}
    for (i, j), c in p.items():
        for r in range(i+1):
            cx = Q(comb(i, r))*x0**(i-r)*hx**r
            for s in range(j+1):
                cy = Q(comb(j, s))*y0**(j-s)*hy**s
                local[(r, s)] = local.get((r, s), Q(0)) + c*cx*cy

    grid: BernsteinGrid = [[Q(0) for _ in range(dy+1)] for _ in range(dx+1)]
    for k in range(dx+1):
        for ell in range(dy+1):
            value = Q(0)
            for r in range(k+1):
                xr = Q(comb(k, r), comb(dx, r))
                for s in range(ell+1):
                    a = local.get((r, s), Q(0))
                    if a:
                        value += a*xr*Q(comb(ell, s), comb(dy, s))
            grid[k][ell] = value
    return grid


def split_curve_half(control: List[Q]) -> Tuple[List[Q], List[Q]]:
    """de Casteljau subdivision of a Bernstein curve at 1/2."""
    d = len(control)-1
    temp = list(control)
    left = [Q(0)]*(d+1)
    right = [Q(0)]*(d+1)
    left[0] = temp[0]
    right[d] = temp[d]
    for r in range(1, d+1):
        temp = [(temp[i]+temp[i+1])/2 for i in range(d-r+1)]
        left[r] = temp[0]
        right[d-r] = temp[-1]
    return left, right


def split_x(grid: BernsteinGrid) -> Tuple[BernsteinGrid, BernsteinGrid]:
    dx, dy = len(grid)-1, len(grid[0])-1
    left = [[Q(0) for _ in range(dy+1)] for _ in range(dx+1)]
    right = [[Q(0) for _ in range(dy+1)] for _ in range(dx+1)]
    for ell in range(dy+1):
        lcol, rcol = split_curve_half([grid[k][ell] for k in range(dx+1)])
        for k in range(dx+1):
            left[k][ell] = lcol[k]
            right[k][ell] = rcol[k]
    return left, right


def split_y(grid: BernsteinGrid) -> Tuple[BernsteinGrid, BernsteinGrid]:
    dx, dy = len(grid)-1, len(grid[0])-1
    low = [[Q(0) for _ in range(dy+1)] for _ in range(dx+1)]
    high = [[Q(0) for _ in range(dy+1)] for _ in range(dx+1)]
    for k in range(dx+1):
        lrow, rrow = split_curve_half(grid[k])
        low[k] = lrow
        high[k] = rrow
    return low, high


def strictly_positive(grid: BernsteinGrid) -> bool:
    return all(c > 0 for row in grid for c in row)


def large_entropy_polynomial() -> Poly:
    # (2/3)(1-gamma_7)^2(1-t)^6 - gamma_7(7-21t+35t^2)^2,
    # represented with t=e and no f dependence.
    one_minus_t = psub(ONE, E_VAR)
    quadratic = padd(pconst(7), padd(pscale(E_VAR, -21), pscale(ppow(E_VAR, 2), 35)))
    return psub(
        pscale(ppow(one_minus_t, 6), Q(2, 3)*(1-g)**2),
        pscale(ppow(quadratic, 2), g),
    )


large_grid = power_to_bernstein_on_box(
    large_entropy_polynomial(), Q(0), Q(3, 20), Q(0), Q(1)
)
large_coefficients = [row[0] for row in large_grid]
assert all(c > 0 for c in large_coefficients)
assert min(large_coefficients) == Q(
    155130361249919329, 1328763571296000000
)


# ---------------------------------------------------------------------------
# Dimension seven: exact two-variable single-bridge certificates
# ---------------------------------------------------------------------------

def single_bridge_polynomials(a: int, b: int) -> Tuple[Poly, Poly, Poly]:
    """Return G=A0-K, the cleared numerator of B, and the final numerator F."""
    assert a >= 1 and b >= 1 and a+b == 6
    m, ell, d = 7-a, 7-b, min(a, b)

    Pa = product_envelope(a, E_VAR)
    Pb = product_envelope(b, F_VAR)
    D = psub(pconst(2), padd(Pa, Pb))
    K = psub(pconst(g), D)
    A0 = pscale(ppow(affine(1, -1, -1), 7), sharp_block(7, a, b))

    m_minus_e = affine(m, -1, 0)
    ell_minus_f = affine(ell, 0, -1)
    one_minus_sum = affine(1, -1, -1)
    t_plus_d = affine(d, 1, 1)

    denominator = pmul(pmul(m_minus_e, ell_minus_f), one_minus_sum)
    first = pmul(pmul(pmul(Pa, E_VAR), ell_minus_f), one_minus_sum)
    second = pmul(pmul(pmul(Pb, F_VAR), m_minus_e), one_minus_sum)
    third = pmul(
        psub(pmul(t_plus_d, K), pscale(A0, d)),
        pmul(m_minus_e, ell_minus_f),
    )
    Bnum = psub(padd(first, second), third)

    Gpoly = psub(A0, K)
    Fpoly = psub(
        pscale(pmul(A0, ppow(Bnum, 2)), 2),
        pscale(pmul(ppow(psub(K, A0), 2), ppow(denominator, 2)), 3),
    )
    return Gpoly, Bnum, Fpoly


def certify_single_bridge(a: int, b: int) -> Tuple[int, int, int, int]:
    """Certify the disjunction G>0 or (Bnum>0 and Fnum>0) on [0,3/20]^2."""
    h = Q(3, 20)
    polys = single_bridge_polynomials(a, b)
    initial = tuple(
        power_to_bernstein_on_box(p, Q(0), h, Q(0), h) for p in polys
    )

    # A node stores the three Bernstein grids and the dyadic depths in x and y.
    queue = deque([(initial, 0, 0)])
    nodes = direct_leaves = entropy_leaves = max_depth = 0
    while queue:
        grids, xdepth, ydepth = queue.popleft()
        nodes += 1
        max_depth = max(max_depth, xdepth+ydepth)

        if strictly_positive(grids[0]):
            direct_leaves += 1
            continue
        if strictly_positive(grids[1]) and strictly_positive(grids[2]):
            entropy_leaves += 1
            continue

        # Split the physically longer side.  Equal sides are split in e first.
        if xdepth <= ydepth:
            split = [split_x(grid) for grid in grids]
            left = tuple(pair[0] for pair in split)
            right = tuple(pair[1] for pair in split)
            queue.append((left, xdepth+1, ydepth))
            queue.append((right, xdepth+1, ydepth))
        else:
            split = [split_y(grid) for grid in grids]
            low = tuple(pair[0] for pair in split)
            high = tuple(pair[1] for pair in split)
            queue.append((low, xdepth, ydepth+1))
            queue.append((high, xdepth, ydepth+1))

        # A failed certificate should stop loudly rather than subdivide forever.
        assert xdepth+ydepth < 12

    return nodes, direct_leaves, entropy_leaves, max_depth


EXPECTED_SINGLE_BRIDGE = {
    (1, 5): (39, 10, 10, 8),
    (2, 4): (35, 9, 9, 6),
    (3, 3): (47, 11, 13, 8),
}
for pair, expected in EXPECTED_SINGLE_BRIDGE.items():
    assert certify_single_bridge(*pair) == expected


# ---------------------------------------------------------------------------
# Dimensions eight through ten
# ---------------------------------------------------------------------------

T = {8: Q(1, 7), 9: Q(1, 10), 10: Q(1, 10)}
LAM = {8: Q(20), 9: Q(27), 10: Q(33)}
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
    8: (Q(4139711814785404757611563,
          300024713494271812137373270016), 1, 1),
    9: (Q(3087705658660821267020074863425,
          540713498338770121894081281002569728), 1, 1),
    10: (Q(36184754863040830496579778359572513127,
           18255668299923283084495152136002059873437500), 1, 1),
}
EXPECTED_SMALL = {
    8: Q(165702517175104104520946768451687,
         215701290564174919711661340674228224),
    9: Q(5586975164420125257678146526492250330079,
         16534930139652222569055999993919229392846848),
    10: Q(55944735805258396002952464523577322487640863,
          404273588897859606879181224832534790039062500000),
}
EXPECTED_LARGE = {
    8: Q(405685268407642329, 2147483648000000000),
    9: Q(17233957123121, 73811250000000),
    10: Q(7399058328856214668089, 30517578125000000000000),
}

for n in (8, 9, 10):
    g = gamma(n)
    assert T[n]**2 - n*g/(1-g) == EXPECTED_TAU_GAPS[n] > 0
    assert Q(1, 100) - Q(n, 2)*g/(1-g) == EXPECTED_D0_GAPS[n] > 0
    assert Q(1, 81) - g == EXPECTED_GAMMA_81_GAPS[n] > 0
    assert LAM[n] <= comb(n, 2) - comb(n, 3)*T[n]

    records = []
    for a in range(1, n):
        for b in range(a, n-a):
            p = n-a-b
            c = (alpha(n, a)+alpha(n, b))/p**2
            records.append((M2(n, sharp_block(n, a, b), c, LAM[n]), a, b))
    assert len(records) == EXPECTED_LOW_COUNTS[n]
    assert all(m > 0 for m, _, _ in records)
    assert min(records) == EXPECTED_LOW_MINIMA[n]

    small = Q(2, 5)*(1-g)**2*g - Q(n*n, 2)*kappa(n)**2
    large = Q(4, 5)*(1-g)**2*Q(9, 10)**n - n*n*g
    assert small == EXPECTED_SMALL[n] > 0
    assert large == EXPECTED_LARGE[n] > 0


# ---------------------------------------------------------------------------
# Dimensions eleven through fourteen
# ---------------------------------------------------------------------------

EXPECTED_MID_COUNTS = {11: 35, 12: 42, 13: 49, 14: 56}
EXPECTED_MID_MINIMA = {
    11: (Q(1260916230510486033100687951530284366056569069,
           6568905011778983402353515625000000000000000000000000), 0, 9),
    12: (Q(1296568134644485856101383130605756547126510510475,
           52823701911706625804697555859966244866493587958274637824), 0, 11),
    13: (Q(2344634926829311765047551592181776557728012393742935099475,
           64121715613888291470979410224708337217950361958912179251651280896),
         0, 12),
    14: (Q(4746695979316977001752598428103553549748722442212199016699590578325,
           277886267525954703535897576147309905022601205806915212282941432903903914488),
         0, 13),
}

for n in (11, 12, 13, 14):
    records = []
    for a in range(n):
        for b in range(a, n-a):
            if n == 11 and (a, b) == (0, 10):
                continue
            p = n-a-b
            c = (alpha(n, a)+alpha(n, b))/p**2
            L = kappa(n) if a == 0 else weak_block(n, a, b)
            records.append((M1(n, L, c), a, b))
    assert len(records) == EXPECTED_MID_COUNTS[n]
    assert all(m > 0 for m, _, _ in records)
    assert min(records) == EXPECTED_MID_MINIMA[n]

assert Q(1, 1000) - gamma(11) == Q(22308624601, 25937424601000) > 0

g = gamma(11)
cstar = 1/((1-g)**2/alpha(11, 1) + 1/alpha(11, 10))
assert cstar == Q(1540632960668022501507138780,
                  1671236290061169784540151329)
assert Q(93, 100) - cstar == Q(
    1361678908886539811520195597,
    167123629006116978454015132900) > 0
assert M1(11, kappa(11), Q(93, 100)) == Q(
    82169554295122875121722420238561072941763350843,
    656890501177898340235351562500000000000000000000000000) > 0


# ---------------------------------------------------------------------------
# Dimension fifteen
# ---------------------------------------------------------------------------

n = 15
g = gamma(n)
mu1 = v(14)*G(14)
mu2 = v(13)*G(13)**2
expected_m15_thin = Q(
    int("93713344067889437339457163890050246146212691912980"
        "360828364561109008039"),
    int("34960636059289149791385973204315516624092004582151"
        "399011826949685248000000000000"),
)
expected_m15_thick = Q(
    int("63429645430340668282196872894311697916976706007718"
        "554270064778550284048088388617933319028973541376"),
    int("36930293415332570245499964107535921252856335539265"
        "82232246120370444607869526685543159813878074388916015625"),
)
assert M1(15, mu1, Q(465, 49)) == expected_m15_thin > 0
assert M1(15, mu2, Q(15)) == expected_m15_thick > 0

EXPECTED_N15_INTEGERS = (
    45591579980859375,
    81568443603515625,
    290584966207055614607247802368,
    3939835293455120947239518208,
    18198910159813162803331251792701065161767,
    162574843672041436633515011954046933353332,
)
actual_n15_integers = (
    15**15 - 300000*factorial(15),
    29863*15**15 - 10**10*factorial(15),
    10**10*factorial(13)*13**13 - 29937*14**26,
    3*14**26 - 10**6*factorial(13)*13**13,
    10**8*factorial(13)*12**24 - 301*13**37,
    4*13**37 - 10**6*factorial(13)*12**24,
)
assert actual_n15_integers == EXPECTED_N15_INTEGERS
assert all(x > 0 for x in actual_n15_integers)
assert 5*(10**12*299999) - 534*9*300000*10**9 == 58195000000000000
assert 14*(10**12*299999) - 3375*4*300000*10**9 == 149986000000000000

print("All exact rational and Bernstein certificates passed.")
