#!/usr/bin/env python3
"""Exact rational certificates for a proposed n=6 Dittert proof.

The script uses only Python's standard library and fractions.Fraction.
It verifies:
  * n=6 relative-entropy localization and the signed excess box;
  * the global and conditional inverse-energy tiers;
  * rational endpoint certificates for the two strengthened core-entropy
    moduli used in codimension one;
  * the one-sided small- and large-entropy inequalities;
  * all two-sided Hall configurations, including all sign patterns, by
    exact bivariate Bernstein subdivision.

The analytic input accompanying this file is the penalized Hall-core
criterion

    4 lambda theta (B + h R)(A B + q R) > R^2,

obtained by retaining both the conditional-product penalty and the core
entropy penalty in the stationarity inequality.  Every polynomial test below
uses exact rational arithmetic; no floating-point comparison occurs.
"""

from collections import Counter
from fractions import Fraction as Q
from math import comb, factorial


def gamma(n):
    return Q(factorial(n), n**n)


def G(k):
    return Q(1) if k <= 1 else Q((k - 1) ** (k - 1), k ** (k - 1))


def v(k):
    return Q(1) if k == 0 else Q(factorial(k), k**k)


def kappa(n):
    return v(n - 1) * G(n - 1)


def sharp_block(n, a, b):
    return v(n - a) * v(n - b) / v(n - a - b)


N = 6
GAMMA = gamma(N)
ETA = Q(1, 27)
T_EXCESS = Q(217, 1000)
THETA_GLOBAL = Q(7, 20)


def alpha(k):
    if k == 0:
        return Q(0)
    p = Q(k, N)
    if 2 * k < N:
        return 2 * N * (p + ETA) * (1 - p - ETA)
    return Q(2 * k * (N - k), N)


# ---------------------------------------------------------------------------
# Localization and the signed excess box.
# ---------------------------------------------------------------------------

assert GAMMA == Q(5, 324)
assert ETA * ETA - GAMMA / (2 * N * (1 - GAMMA)) == Q(61, 930204) > 0
assert Q(1, 2) - Q(2, 6) - ETA == Q(7, 54) > 0
assert tuple(alpha(k) for k in range(N)) == (
    Q(0), Q(473, 243), Q(680, 243), Q(3), Q(8, 3), Q(5, 3)
)
assert all(alpha(k) <= Q(3) for k in range(N))
assert T_EXCESS * T_EXCESS - Q(3) * GAMMA / (1 - GAMMA) == Q(
    21391, 319000000
) > 0
assert kappa(N) - GAMMA == Q(37531, 126562500) > 0


# ---------------------------------------------------------------------------
# Conditional inverse-energy tiers.
#
# If H_0 >= A_* then P,Q >= 1-gamma+A_*.  The first exact inequality in
# each row forces every row and column sum to be at least r_*.  The second
# gives (1-gamma+A_*) r_*^2/2 > theta_*.
# ---------------------------------------------------------------------------

ENERGY_TIERS = (
    (Q(0), Q(169, 200), Q(7, 20), "t350"),
    (Q(3, 800), Q(86533, 100000), Q(37, 100), "t370"),
    (Q(7, 1250), Q(8761, 10000), Q(19, 50), "t380"),
    (Q(73, 10000), Q(8869, 10000), Q(39, 100), "t390"),
    (Q(161, 20000), Q(17843, 20000), Q(79, 200), "t395"),
    (Q(17, 2000), Q(4477, 5000), Q(199, 500), "t398"),
    (Q(4383, 500000), Q(35897, 40000), Q(2, 5), "t400"),
)

for A_star, r_star, theta, _ in ENERGY_TIERS:
    P_star = 1 - GAMMA + A_star
    max_product_with_small_coordinate = r_star * ((N - r_star) / (N - 1)) ** (N - 1)
    assert P_star - max_product_with_small_coordinate > 0
    assert P_star * r_star * r_star / 2 - theta > 0


# ---------------------------------------------------------------------------
# Rational endpoint checks for the pointwise core-entropy moduli.
#
# For h(t)=(1-t)log(1-t), the elementary calculus proof only needs
#   log(m/(m-1)) > 1/m + lambda/m^2.
# The alternating-series lower bounds below are rational.
# ---------------------------------------------------------------------------


def alternating_log_lower(x, even_terms):
    assert even_terms % 2 == 0
    return sum(((-1) ** (k + 1)) * x**k / k for k in range(1, even_terms + 1))


# lambda=23/40 for core marginal dimensions 5 and 2.
assert alternating_log_lower(Q(1, 4), 6) == Q(27419, 122880) > Q(223, 1000)
assert Q(2, 3) > Q(103, 160)  # log 2 > 2/3.

# lambda=301/500 for core marginal dimensions 4 and 3.
assert alternating_log_lower(Q(1, 3), 6) == Q(12581, 43740) > Q(2301, 8000)
assert alternating_log_lower(Q(1, 2), 4) == Q(77, 192) > Q(1801, 4500)


# ---------------------------------------------------------------------------
# One-sided Hall cuts.
# ---------------------------------------------------------------------------

SMALL_ONE_SIDED = (
    THETA_GLOBAL * (1 - GAMMA) ** 2 * GAMMA
    - Q(N * N, 2) * kappa(N) ** 2
)
assert SMALL_ONE_SIDED == Q(
    16250344917716167,
    20759414062500000000,
) > 0


# ---------------------------------------------------------------------------
# Exact bivariate polynomial and Bernstein utilities.
# A polynomial is a dict (i,j) -> coefficient of x^i y^j.
# ---------------------------------------------------------------------------


def clean(p):
    return {ij: c for ij, c in p.items() if c}


def const(c):
    c = Q(c)
    return {} if not c else {(0, 0): c}


X_VAR = {(1, 0): Q(1)}
Y_VAR = {(0, 1): Q(1)}


def add(p, q):
    out = dict(p)
    for ij, c in q.items():
        out[ij] = out.get(ij, Q(0)) + c
    return clean(out)


def neg(p):
    return {ij: -c for ij, c in p.items()}


def sub(p, q):
    return add(p, neg(q))


def scale(p, c):
    c = Q(c)
    return clean({ij: c * a for ij, a in p.items()})


def mul(p, q):
    out = {}
    for (i, j), a in p.items():
        for (k, ell), b in q.items():
            key = (i + k, j + ell)
            out[key] = out.get(key, Q(0)) + a * b
    return clean(out)


def product(*polys):
    out = const(1)
    for p in polys:
        out = mul(out, p)
    return out


def power(p, k):
    out = const(1)
    base = p
    while k:
        if k & 1:
            out = mul(out, base)
        base = mul(base, base)
        k //= 2
    return out


def degree_pair(p):
    return max(i for i, _ in p), max(j for _, j in p)


def power_to_bernstein(p, xmax, ymax):
    dx, dy = degree_pair(p)
    coeff = [[Q(0) for _ in range(dy + 1)] for _ in range(dx + 1)]
    for (i, j), c in p.items():
        coeff[i][j] = c * xmax**i * ymax**j
    ans = [[Q(0) for _ in range(dy + 1)] for _ in range(dx + 1)]
    for k in range(dx + 1):
        for ell in range(dy + 1):
            total = Q(0)
            for i in range(k + 1):
                wx = Q(comb(k, i), comb(dx, i))
                for j in range(ell + 1):
                    wy = Q(comb(ell, j), comb(dy, j))
                    total += coeff[i][j] * wx * wy
            ans[k][ell] = total
    return ans


def split_vector_half(values):
    levels = [values]
    while len(levels[-1]) > 1:
        old = levels[-1]
        levels.append([(old[i] + old[i + 1]) / 2 for i in range(len(old) - 1)])
    d = len(values) - 1
    left = [levels[i][0] for i in range(d + 1)]
    right = [levels[d - i][i] for i in range(d + 1)]
    return left, right


def split_matrix_half(coeffs, axis):
    dx = len(coeffs) - 1
    dy = len(coeffs[0]) - 1
    left = [[Q(0) for _ in range(dy + 1)] for _ in range(dx + 1)]
    right = [[Q(0) for _ in range(dy + 1)] for _ in range(dx + 1)]
    if axis == 0:
        for j in range(dy + 1):
            lo, hi = split_vector_half([coeffs[i][j] for i in range(dx + 1)])
            for i in range(dx + 1):
                left[i][j], right[i][j] = lo[i], hi[i]
    else:
        for i in range(dx + 1):
            left[i], right[i] = split_vector_half(coeffs[i])
    return left, right


def bounds(coeffs):
    flat = [x for row in coeffs for x in row]
    return min(flat), max(flat)


# One-sided large entropy.  It is enough to prove, on 0<=t<=217/1000,
#   (7/10)(1-gamma)^2(1-t)^4
#       > gamma(6-15t+20t^2-15t^3+6t^4-t^5)^2.
T_VAR = X_VAR
SIX_GEOM = const(0)
for k in range(1, 7):
    SIX_GEOM = add(
        SIX_GEOM,
        scale(power(T_VAR, k - 1), Q(((-1) ** (k + 1)) * comb(6, k))),
    )
LARGE_ONE_SIDED_POLY = sub(
    scale(power(sub(const(1), T_VAR), 4), Q(7, 10) * (1 - GAMMA) ** 2),
    scale(power(SIX_GEOM, 2), GAMMA),
)
LARGE_BERNSTEIN = power_to_bernstein(LARGE_ONE_SIDED_POLY, T_EXCESS, Q(1))
LARGE_COEFFS = [row[0] for row in LARGE_BERNSTEIN]
assert min(LARGE_COEFFS) == Q(
    3951217684118237014875556869551,
    64800000000000000000000000000000,
) > 0


# ---------------------------------------------------------------------------
# Penalized two-sided Hall certificates.
# ---------------------------------------------------------------------------


def subset_product(k, signed_variable):
    other = N - k
    return mul(
        power(add(const(1), scale(signed_variable, Q(1, k))), k),
        power(add(const(1), scale(signed_variable, -Q(1, other))), other),
    )


def entropy_lambda(a, b):
    if (a, b) == (1, 4):
        return Q(23, 40)
    if (a, b) == (2, 3):
        return Q(301, 500)
    return Q(1, 2 * (a * a + b * b))


def build_case_polynomials(a, b, quadrant, lam, theta):
    signs = {"++": (1, 1), "-+": (-1, 1), "+-": (1, -1)}
    sx, sy = signs[quadrant]
    e = scale(X_VAR, sx)
    f = scale(Y_VAR, sy)
    e_pos = X_VAR if sx > 0 else const(0)
    f_pos = Y_VAR if sy > 0 else const(0)
    t = add(e_pos, f_pos)

    p = N - a - b
    m = N - a
    ell = N - b
    d = min(a, b)

    pa = subset_product(a, e)
    pb = subset_product(b, f)
    K = add(const(GAMMA - 2), add(pa, pb))
    omega = sub(const(1), scale(t, Q(1, p)))
    A0 = scale(power(omega, N), sharp_block(N, a, b))
    R = sub(K, A0)

    m_minus_e = sub(const(m), e)
    ell_minus_f = sub(const(ell), f)
    denominator = product(m_minus_e, ell_minus_f, omega)

    # Positive-denominator numerator of
    # B = p(P_a e/(m-e)+P_b f/(ell-f))
    #       - ((t+d)K-dA_0)/omega.
    row_term = scale(product(pa, e, ell_minus_f, omega), p)
    col_term = scale(product(pb, f, m_minus_e, omega), p)
    upper_term = mul(
        sub(mul(add(t, const(d)), K), scale(A0, d)),
        mul(m_minus_e, ell_minus_f),
    )
    B_num = sub(add(row_term, col_term), upper_term)

    # q=dA_0/omega and the conservative product-loss penalty
    # h=(t+d)/omega-p(E_+ + F_+).
    q = scale(power(omega, N - 1), d * sharp_block(N, a, b))
    h_num = sub(
        mul(add(t, const(d)), mul(m_minus_e, ell_minus_f)),
        scale(
            add(
                product(e_pos, ell_minus_f, omega),
                product(f_pos, m_minus_e, omega),
            ),
            p,
        ),
    )

    first = add(B_num, mul(h_num, R))
    second = add(mul(A0, B_num), mul(product(q, R), denominator))
    master_num = sub(
        scale(mul(first, second), 4 * lam * theta),
        mul(power(R, 2), power(denominator, 2)),
    )

    return {
        "K": K,
        "A": A0,
        "direct": neg(R),
        "B": B_num,
        "M": master_num,
        "critical": add(e, f),
    }


# On the full signed-excess box the conservative h above is positive.
TWO_SIDED_PAIRS = ((1, 1), (1, 2), (1, 3), (1, 4), (2, 2), (2, 3))
for a, b in TWO_SIDED_PAIRS:
    p = N - a - b
    ell = N - b
    d = min(a, b)
    assert d - Q(2 * p) * T_EXCESS / (ell - T_EXCESS) > 0


EXPECTED_CASES = {
    (1, 1, "++"): (21, 10, Counter({"stationarity": 6, "direct": 5})),
    (1, 1, "-+"): (13, 6, Counter({"direct": 7})),
    (1, 1, "+-"): (11, 5, Counter({"direct": 4, "noncritical": 2})),
    (1, 2, "++"): (21, 9, Counter({"stationarity": 6, "direct": 5})),
    (1, 2, "-+"): (9, 4, Counter({"direct": 5})),
    (1, 2, "+-"): (3, 1, Counter({"direct": 2})),
    (1, 3, "++"): (39, 8, Counter({"direct": 11, "stationarity": 9})),
    (1, 3, "-+"): (11, 5, Counter({"direct": 5, "stationarity": 1})),
    (1, 3, "+-"): (3, 1, Counter({"direct": 2})),
    (1, 4, "++"): (57, 11, Counter({"stationarity": 16, "direct": 13})),
    (1, 4, "-+"): (29, 8, Counter({"direct": 9, "stationarity": 5, "noncritical": 1})),
    (1, 4, "+-"): (65, 14, Counter({"direct": 18, "stationarity": 14, "noncritical": 1})),
    (2, 2, "++"): (23, 8, Counter({"direct": 6, "stationarity": 6})),
    (2, 2, "-+"): (5, 2, Counter({"direct": 3})),
    (2, 2, "+-"): (3, 1, Counter({"direct": 2})),
    (2, 3, "++"): (
        131,
        15,
        Counter({
            "stationarity": 50,
            "direct": 16,
            "t395": 11,
            "t398": 9,
            "t400": 8,
            "t350": 8,
            "t380": 5,
            "t390": 5,
            "t370": 4,
        }),
    ),
    (2, 3, "-+"): (
        23,
        8,
        Counter({"direct": 7, "stationarity": 5, "t400": 4, "t380": 1}),
    ),
    (2, 3, "+-"): (
        47,
        10,
        Counter({
            "direct": 12,
            "stationarity": 11,
            "t400": 9,
            "noncritical": 1,
            "t398": 1,
            "t380": 1,
        }),
    ),
}


def prepare_case(a, b, quadrant):
    lam = entropy_lambda(a, b)
    if (a, b) == (2, 3):
        base = build_case_polynomials(a, b, quadrant, lam, ENERGY_TIERS[0][2])
        polys = {name: base[name] for name in ("K", "A", "direct", "B", "critical")}
        for _, _, theta, name in ENERGY_TIERS:
            polys[name] = build_case_polynomials(a, b, quadrant, lam, theta)["M"]
    else:
        base = build_case_polynomials(a, b, quadrant, lam, THETA_GLOBAL)
        polys = {name: base[name] for name in ("K", "A", "direct", "B", "M", "critical")}
    return {
        name: power_to_bernstein(poly, T_EXCESS, T_EXCESS)
        for name, poly in polys.items()
    }


def certify_case(a, b, quadrant):
    arrays = prepare_case(a, b, quadrant)
    stack = [(arrays, 0)]
    counts = Counter()
    nodes = 0
    maximum_depth = 0

    while stack:
        current, depth = stack.pop()
        nodes += 1
        maximum_depth = max(maximum_depth, depth)

        if quadrant != "++":
            _, critical_hi = bounds(current["critical"])
            if critical_hi <= 0:
                counts["noncritical"] += 1
                continue

        _, K_hi = bounds(current["K"])
        if K_hi < 0:
            counts["K<0"] += 1
            continue

        direct_lo, _ = bounds(current["direct"])
        if direct_lo > 0:
            counts["direct"] += 1
            continue

        B_lo, _ = bounds(current["B"])
        if (a, b) == (2, 3):
            A_lo, _ = bounds(current["A"])
            selected = ENERGY_TIERS[0]
            for tier in ENERGY_TIERS:
                if A_lo >= tier[0]:
                    selected = tier
            master_name = selected[3]
        else:
            master_name = "M"

        master_lo, _ = bounds(current[master_name])
        if B_lo > 0 and master_lo > 0:
            counts["stationarity"] += 1
            if (a, b) == (2, 3):
                counts[master_name] += 1
            continue

        assert depth < 24, (
            a,
            b,
            quadrant,
            depth,
            bounds(current["K"]),
            bounds(current["A"]),
            bounds(current["direct"]),
            bounds(current["B"]),
            bounds(current[master_name]),
        )

        axis = depth % 2
        children = [{}, {}]
        for name, coeffs in current.items():
            left, right = split_matrix_half(coeffs, axis)
            children[0][name] = left
            children[1][name] = right
        stack.append((children[1], depth + 1))
        stack.append((children[0], depth + 1))

    return nodes, maximum_depth, counts


ACTUAL_CASES = {
    key: certify_case(*key)
    for key in EXPECTED_CASES
}
assert ACTUAL_CASES == EXPECTED_CASES

print("All exact n=6 Dittert certificates passed.")
