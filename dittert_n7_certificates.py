#!/usr/bin/env python3
"""Exact rational certificates for a proposed n=7 extension of the Dittert proof.

The script uses only Python's standard library and fractions.Fraction.  It checks:
  * the n=7 localization and Hall-dilation ranges;
  * all p>=2 two-sided Hall pairs by the quadratic scalar criterion;
  * the sign-mixed p=1 cases by a sharper one-sided Hall coefficient;
  * the small- and large-entropy one-sided estimates;
  * the three genuinely hard p=1 pairs via exact bivariate Bernstein certificates.

The last certificate proves, on 0 <= e,f <= 3/20, that either the direct
permanent comparison already contradicts maximality, or the conditional-energy
and core-entropy inequalities contradict one another.
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


def sharp_block(n, a, b):
    return v(n - a) * v(n - b) / v(n - a - b)


N = 7
GAMMA = gamma(N)
ETA = Q(1, 47)
T_GLOBAL = Q(5, 24)
LAMBDA_GLOBAL = Q(329, 24)
T_SIGN = Q(3, 20)
LAMBDA_SIGN = Q(63, 4)
ENERGY = Q(1, 3)


def alpha(k):
    if k == 0:
        return Q(0)
    p = Q(k, N)
    if 2 * k < N:
        return 2 * N * (p + ETA) * (1 - p - ETA)
    return Q(2 * k * (N - k), N)


def M2(L, c, lam):
    q = lam * L + (1 - GAMMA) / c
    return L - GAMMA - Q(N * N, 4) * L * L / q


# Localization, uniform alpha bound, and two useful Hall ranges.
assert ETA * ETA - GAMMA / (2 * N * (1 - GAMMA)) == Q(23263, 1808073127) > 0
assert Q(1, 2) - Q(3, 7) - ETA == Q(33, 658) > 0
assert all(alpha(k) <= Q(7, 2) for k in range(N))
assert T_GLOBAL * T_GLOBAL - N * GAMMA / (1 - GAMMA) == Q(20185, 67351104) > 0
assert T_SIGN * T_SIGN - Q(N, 2) * GAMMA / (1 - GAMMA) == Q(44361, 46771600) > 0

# The row/column energy coefficient is strictly larger than 1/3.
assert Q(1, 160) - GAMMA == Q(2449, 18823840) > 0
assert Q(9, 400) - Q(N, 2) * GAMMA / (1 - GAMMA) == Q(44361, 46771600) > 0
assert Q(159, 160) * Q(17, 20) ** 2 / 2 - ENERGY == Q(9853, 384000) > 0


# All two-sided pairs with p >= 2.
EXPECTED_P_GE_2 = {
    (1, 1): Q(217499955444827479495, 21713873952987016007066496),
    (1, 2): Q(34682878961789264, 639687598598410092675),
    (1, 3): Q(95714954803771495, 1019889385414932689664),
    (1, 4): Q(7349321695167505, 145541465008251155487),
    (2, 2): Q(174030801342732607248, 983554287703870421328125),
    (2, 3): Q(4247735832777429, 22030864916302726750),
}

actual_p_ge_2 = {}
for a in range(1, N):
    for b in range(a, N - a):
        p = N - a - b
        if p < 2:
            continue
        c = (alpha(a) + alpha(b)) / (p * p)
        actual_p_ge_2[(a, b)] = M2(sharp_block(N, a, b), c, LAMBDA_GLOBAL)
assert actual_p_ge_2 == EXPECTED_P_GE_2
assert all(x > 0 for x in actual_p_ge_2.values())


# p=1 with one of the two subset excesses nonpositive.  Only the positive
# excess enters the Hall bound, so c is alpha(a) or alpha(b), not their sum.
EXPECTED_SIGN_MIXED = {
    (1, 5, "e<=0"): Q(165035926385, 3897251304526692),
    (1, 5, "f<=0"): Q(132694230562295, 291795945546131964),
    (2, 4, "e<=0"): Q(7338502353824, 23595287183461875),
    (2, 4, "f<=0"): Q(8455201033462752, 16651481270863560625),
    (3, 3, "e<=0"): Q(5244714741663513, 11810357798135411200),
    (3, 3, "f<=0"): Q(5244714741663513, 11810357798135411200),
}

actual_sign_mixed = {}
for a, b in ((1, 5), (2, 4), (3, 3)):
    L = sharp_block(N, a, b)
    actual_sign_mixed[(a, b, "e<=0")] = M2(L, alpha(b), LAMBDA_SIGN)
    actual_sign_mixed[(a, b, "f<=0")] = M2(L, alpha(a), LAMBDA_SIGN)
assert actual_sign_mixed == EXPECTED_SIGN_MIXED
assert all(x > 0 for x in actual_sign_mixed.values())


# One-sided cuts: the small-entropy certificate with energy coefficient 1/3.
KAPPA = v(6) * G(6)
SMALL_ENTROPY = ENERGY * (1 - GAMMA) ** 2 * GAMMA - Q(49, 2) * KAPPA**2
assert KAPPA == Q(15625, 2519424)
assert SMALL_ENTROPY == Q(
    22176489486262543850191055,
    20672701805255574480144334848,
) > 0


# ---------------------------------------------------------------------------
# Exact polynomial and Bernstein utilities.
# A bivariate polynomial is a dict (i,j) -> coefficient of e^i f^j.
# ---------------------------------------------------------------------------


def clean(p):
    return {ij: c for ij, c in p.items() if c}


def const(c):
    return {} if not c else {(0, 0): Q(c)}


E_VAR = {(1, 0): Q(1)}
F_VAR = {(0, 1): Q(1)}


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


def subset_product(k, variable):
    other = N - k
    return mul(
        power(add(const(1), scale(variable, Q(1, k))), k),
        power(add(const(1), scale(variable, -Q(1, other))), other),
    )


def power_to_bernstein(p, xmax, ymax):
    dx, dy = degree_pair(p)
    a = [[Q(0) for _ in range(dy + 1)] for _ in range(dx + 1)]
    for (i, j), c in p.items():
        a[i][j] = c * xmax**i * ymax**j
    ans = [[Q(0) for _ in range(dy + 1)] for _ in range(dx + 1)]
    for k in range(dx + 1):
        for ell in range(dy + 1):
            total = Q(0)
            for i in range(k + 1):
                wx = Q(comb(k, i), comb(dx, i))
                for j in range(ell + 1):
                    wy = Q(comb(ell, j), comb(dy, j))
                    total += a[i][j] * wx * wy
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


# One-sided large entropy.  For 0 <= t <= 3/20, it is enough to prove
#   (2/3)(1-gamma)^2(1-t)^6
#       > gamma(7-21t+35t^2)^2.
# Every Bernstein coefficient on [0,3/20] is positive.
T_VAR = E_VAR
large_poly = sub(
    scale(power(sub(const(1), T_VAR), 6), Q(2, 3) * (1 - GAMMA) ** 2),
    scale(power(add(add(const(7), scale(T_VAR, -21)), scale(power(T_VAR, 2), 35)), 2), GAMMA),
)
large_bernstein = power_to_bernstein(large_poly, T_SIGN, Q(1))
large_coeffs = [row[0] for row in large_bernstein]
assert min(large_coeffs) == Q(
    155130361249919329,
    1328763571296000000,
) > 0


# The hard p=1 certificate.  Put t=e+f, D=(1-P_a)+(1-P_b),
# K=gamma-D, A=L(1-t)^7, and
#
#   B = P_a e/(m-e) + P_b f/(ell-f)
#       - ((t+d)K-dA)/(1-t).
#
# Conditional product energy plus critical-core stationarity requires
#       s u >= (1/3) B^2,
# while u+A exp(s/2)<=K implies
#       s u <= (K-A)^2/(2A).
# Hence 2 A B^2 > 3(K-A)^2 is impossible.  B_NUM and CERT below are the
# positive-denominator cleared versions of B and that strict inequality.


def hard_pair_polynomials(a, b):
    m = N - a
    ell = N - b
    d = min(a, b)
    pa = subset_product(a, E_VAR)
    pb = subset_product(b, F_VAR)
    deficit = sub(const(2), add(pa, pb))
    K = sub(const(GAMMA), deficit)
    t = add(E_VAR, F_VAR)
    A = scale(power(sub(const(1), t), N), sharp_block(N, a, b))
    direct = sub(A, K)
    denominator = mul(mul(sub(const(m), E_VAR), sub(const(ell), F_VAR)), sub(const(1), t))
    b_num = sub(
        add(
            mul(mul(mul(pa, E_VAR), sub(const(ell), F_VAR)), sub(const(1), t)),
            mul(mul(mul(pb, F_VAR), sub(const(m), E_VAR)), sub(const(1), t)),
        ),
        mul(
            sub(mul(add(t, const(d)), K), scale(A, d)),
            mul(sub(const(m), E_VAR), sub(const(ell), F_VAR)),
        ),
    )
    cert = sub(
        scale(mul(A, power(b_num, 2)), 2),
        scale(mul(power(sub(K, A), 2), power(denominator, 2)), 3),
    )
    return {"K": K, "direct": direct, "B": b_num, "cert": cert}


def certify_hard_pair(a, b):
    polys = hard_pair_polynomials(a, b)
    initial = {
        name: power_to_bernstein(poly, T_SIGN, T_SIGN)
        for name, poly in polys.items()
    }
    stack = [(initial, 0)]
    counts = Counter()
    nodes = 0
    while stack:
        arrays, depth = stack.pop()
        nodes += 1
        _, k_hi = bounds(arrays["K"])
        if k_hi < 0:
            counts["deficit"] += 1
            continue
        direct_lo, _ = bounds(arrays["direct"])
        if direct_lo > 0:
            counts["direct"] += 1
            continue
        b_lo, _ = bounds(arrays["B"])
        cert_lo, _ = bounds(arrays["cert"])
        if b_lo > 0 and cert_lo > 0:
            counts["stationarity"] += 1
            continue
        assert depth < 24, (a, b, depth, bounds(arrays["direct"]), bounds(arrays["B"]), bounds(arrays["cert"]))
        axis = depth % 2
        children = [{}, {}]
        for name, coeffs in arrays.items():
            left, right = split_matrix_half(coeffs, axis)
            children[0][name] = left
            children[1][name] = right
        stack.append((children[1], depth + 1))
        stack.append((children[0], depth + 1))
    return nodes, counts


EXPECTED_BERNSTEIN = {
    (1, 5): (39, Counter({"direct": 10, "stationarity": 10})),
    (2, 4): (35, Counter({"direct": 9, "stationarity": 9})),
    (3, 3): (47, Counter({"stationarity": 13, "direct": 11})),
}
actual_bernstein = {pair: certify_hard_pair(*pair) for pair in EXPECTED_BERNSTEIN}
assert actual_bernstein == EXPECTED_BERNSTEIN

print("All exact n=7 extension certificates passed.")
