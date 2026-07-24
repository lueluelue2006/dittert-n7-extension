#!/usr/bin/env python3
"""Exact certificates for the dimension-six Dittert extension.

The script uses only Python integers and fractions.Fraction.  It verifies:

* the row/column floor and all subset-excess localization constants;
* the conditional inverse-energy coefficient C_* = 3523/10000;
* the strengthened m=3,4 cofactor-entropy constants used for p=1;
* every p>=2 two-sided Hall face by exact bivariate Bernstein subdivision;
* every p=1 sign branch by exact trivariate Bernstein subdivision and the
  quadratic tangent-duality criterion with lambda = 1/6;
* the small- and large-entropy one-sided Hall estimates.

Bernstein coefficients are initially computed as exact Fractions.  Before
subdivision, each coefficient tensor is cleared to integers.  A midpoint
(de Casteljau) split is then performed by integer additions and powers of two;
all sign decisions therefore remain exact and no floating-point computation
occurs.
"""

from collections import Counter
from fractions import Fraction as Q
from itertools import product
from math import comb, factorial, lcm


# ---------------------------------------------------------------------------
# Basic constants.
# ---------------------------------------------------------------------------

N = 6
GAMMA = Q(factorial(N), N**N)          # 5/324
ROW_FLOOR = Q(423, 500)
ENERGY = Q(3523, 10000)
T_ONE_SIDED = Q(31, 200)
TANGENT_LAMBDA = Q(1, 6)
TANGENT_SLOPE = 1 + TANGENT_LAMBDA     # 7/6

H = {
    1: Q(17, 100),
    2: Q(21, 100),
    3: Q(11, 50),
    4: Q(1, 5),
    5: Q(31, 200),
}


def v(k):
    return Q(1) if k == 0 else Q(factorial(k), k**k)


def G(k):
    return Q(1) if k <= 1 else Q((k - 1) ** (k - 1), k ** (k - 1))


def sharp_block(a, b):
    return v(N - a) * v(N - b) / v(N - a - b)


assert GAMMA == Q(5, 324)

# If one row sum were <= 423/500, AM--GM on the other five rows would force
# the row product below 1-gamma.
ROW_FLOOR_GAP = 1 - GAMMA - ROW_FLOOR * ((N - ROW_FLOOR) / 5) ** 5
assert ROW_FLOOR_GAP == Q(
    29233981969641209,
    3955078125000000000000,
) > 0


def subset_envelope(k, t):
    return (1 + t / k) ** k * (1 - t / (N - k)) ** (N - k)


EXPECTED_SUBSET_GAPS = {
    1: Q(99932137192289, 253125000000000000),
    2: Q(38396337601399, 82944000000000000),
    3: Q(7001004061, 11390625000000),
    4: Q(10159, 1296000000),
    5: Q(3512017710706961, 16200000000000000000),
}
for k in range(1, N):
    gap = 1 - GAMMA - subset_envelope(k, H[k])
    assert gap == EXPECTED_SUBSET_GAPS[k] > 0

# The raw conditional-energy coefficient is slightly larger than C_*.
assert (1 - GAMMA) * ROW_FLOOR**2 / 2 - ENERGY == Q(71, 2000000) > 0

# For every one-sided cut, tau < h_b/(6-b) <= 31/200.
assert max(H[b] / (N - b) for b in range(1, N)) == T_ONE_SIDED


# ---------------------------------------------------------------------------
# Sparse exact polynomial arithmetic in (e,f,x).
# ---------------------------------------------------------------------------

DIMS = 3
ZERO_EXP = (0, 0, 0)


def clean(poly):
    return {mon: coeff for mon, coeff in poly.items() if coeff}


def const(value):
    return {} if not value else {ZERO_EXP: Q(value)}


def variable(axis):
    exp = [0] * DIMS
    exp[axis] = 1
    return {tuple(exp): Q(1)}


E_VAR = variable(0)
F_VAR = variable(1)
X_VAR = variable(2)


def add(p, q):
    out = dict(p)
    for mon, coeff in q.items():
        out[mon] = out.get(mon, Q(0)) + coeff
    return clean(out)


def neg(p):
    return {mon: -coeff for mon, coeff in p.items()}


def sub(p, q):
    return add(p, neg(q))


def scale(p, scalar):
    scalar = Q(scalar)
    return clean({mon: scalar * coeff for mon, coeff in p.items()})


def mul(p, q):
    out = {}
    for alpha, a in p.items():
        for beta, b in q.items():
            mon = tuple(i + j for i, j in zip(alpha, beta))
            out[mon] = out.get(mon, Q(0)) + a * b
    return clean(out)


def power(p, exponent):
    out = const(1)
    base = p
    while exponent:
        if exponent & 1:
            out = mul(out, base)
        base = mul(base, base)
        exponent //= 2
    return out


def degree_tuple(p):
    return tuple(max((mon[i] for mon in p), default=0) for i in range(DIMS))


def subset_product(k, variable_poly):
    return mul(
        power(add(const(1), scale(variable_poly, Q(1, k))), k),
        power(
            add(const(1), scale(variable_poly, -Q(1, N - k))),
            N - k,
        ),
    )


# ---------------------------------------------------------------------------
# Exact Bernstein conversion and integer midpoint subdivision.
# ---------------------------------------------------------------------------


def affine_axis(coeffs, axis, lower, upper):
    """Substitute x_axis = lower + (upper-lower) u_axis."""
    width = upper - lower
    exponents = {mon[axis] for mon in coeffs}
    cache = {
        exponent: [
            Q(comb(exponent, j)) * lower ** (exponent - j) * width**j
            for j in range(exponent + 1)
        ]
        for exponent in exponents
    }
    out = {}
    for mon, coeff in coeffs.items():
        exponent = mon[axis]
        for j, factor in enumerate(cache[exponent]):
            beta = list(mon)
            beta[axis] = j
            beta = tuple(beta)
            out[beta] = out.get(beta, Q(0)) + coeff * factor
    return clean(out)


def power_to_bernstein_axis(coeffs, axis, degree):
    cache = {
        exponent: [
            Q(comb(k, exponent), comb(degree, exponent))
            for k in range(exponent, degree + 1)
        ]
        for exponent in range(degree + 1)
    }
    out = {}
    for mon, coeff in coeffs.items():
        exponent = mon[axis]
        for offset, factor in enumerate(cache[exponent]):
            beta = list(mon)
            beta[axis] = exponent + offset
            beta = tuple(beta)
            out[beta] = out.get(beta, Q(0)) + coeff * factor
    return clean(out)


def bernstein_fraction_tensor(poly, box):
    """Return (multidegree, exact Bernstein coefficient dictionary)."""
    degrees = degree_tuple(poly)
    coeffs = poly
    for axis, (lower, upper) in enumerate(box):
        coeffs = affine_axis(coeffs, axis, lower, upper)
    for axis, degree in enumerate(degrees):
        coeffs = power_to_bernstein_axis(coeffs, axis, degree)
    return degrees, coeffs


def integer_tensor(fraction_tensor):
    """Clear all coefficient denominators, preserving every sign."""
    degrees, coeffs = fraction_tensor
    denominator = 1
    for value in coeffs.values():
        denominator = lcm(denominator, value.denominator)
    dense = {
        index: int(coeffs.get(index, Q(0)) * denominator)
        for index in product(*[range(degree + 1) for degree in degrees])
    }
    return degrees, dense


def tensor_bounds(tensor):
    values = tensor[1].values()
    return min(values), max(values)


def split_integer_tensor(tensor, axis):
    """Exact midpoint de Casteljau split, kept integral by common scaling."""
    degrees, coeffs = tensor
    dimension = len(degrees)
    degree = degrees[axis]
    other_axes = [i for i in range(dimension) if i != axis]
    left = {}
    right = {}

    for other_index in product(*[range(degrees[i] + 1) for i in other_axes]):
        index = [0] * dimension
        for i, value in zip(other_axes, other_index):
            index[i] = value

        line = []
        for k in range(degree + 1):
            index[axis] = k
            line.append(coeffs[tuple(index)])

        levels = [line]
        for _ in range(degree):
            previous = levels[-1]
            levels.append(
                [previous[i] + previous[i + 1] for i in range(len(previous) - 1)]
            )

        # Actual midpoint values at level j have denominator 2^j.  Multiplying
        # every child coefficient by 2^degree gives a common positive scale.
        for k in range(degree + 1):
            index[axis] = k
            left[tuple(index)] = levels[k][0] << (degree - k)
            right[tuple(index)] = levels[degree - k][k] << k

    return (degrees, left), (degrees, right)


# ---------------------------------------------------------------------------
# The m=3,4 cofactor-entropy coefficient 5/8.
# ---------------------------------------------------------------------------

# Boundary gap in the m=4 reduction:
# log(256/243)-5/96 > 0.  The four-term alternating lower sum is exact.
LOG_X = Q(13, 243)
LOG_LOWER_GAP = (
    LOG_X
    - LOG_X**2 / 2
    + LOG_X**3 / 3
    - LOG_X**4 / 4
    - Q(5, 96)
)
assert LOG_LOWER_GAP == Q(3635617, 111577100832) > 0

# At an interior minimizer for m=4, the only nonuniform candidate is
# (t,(1-t)/3,(1-t)/3,(1-t)/3).  Retaining power-sum terms k=3,...,7 gives
# the following exact positive factorization.
T_VAR = E_VAR
BETA = scale(sub(const(1), T_VAR), Q(1, 3))
finite_tail = const(0)
for k in range(3, 8):
    pk = add(power(T_VAR, k), scale(power(BETA, k), 3))
    uniform_pk = Q(4, 4**k)
    finite_tail = add(
        finite_tail,
        scale(sub(pk, const(uniform_pk)), Q(1, k * (k - 1))),
    )
p2 = add(power(T_VAR, 2), scale(power(BETA, 2), 3))
finite_tail = sub(finite_tail, scale(sub(p2, const(Q(1, 4))), Q(1, 8)))
positive_factor = add(
    add(
        add(
            add(
                add(scale(power(T_VAR, 5), 133120), scale(power(T_VAR, 4), 255232)),
                scale(power(T_VAR, 3), 387328),
            ),
            scale(power(T_VAR, 2), 696752),
        ),
        scale(T_VAR, 1028168),
    ),
    const(2057),
)
factorized_tail = scale(
    mul(power(sub(scale(T_VAR, 4), const(1)), 2), positive_factor),
    Q(1, 89579520),
)
assert finite_tail == factorized_tail


# ---------------------------------------------------------------------------
# p >= 2 two-sided Hall faces.
# ---------------------------------------------------------------------------


def p_ge_2_polynomials(a, b):
    p = N - a - b
    m = N - a
    ell = N - b
    d = min(a, b)

    pa = subset_product(a, E_VAR)
    pb = subset_product(b, F_VAR)
    deficit_envelope = sub(const(2), add(pa, pb))
    K = sub(const(GAMMA), deficit_envelope)
    t = add(E_VAR, F_VAR)
    z0 = add(const(1), scale(t, -Q(1, p)))
    A0 = scale(power(z0, N), sharp_block(a, b))
    direct = sub(A0, K)

    m_minus_e = sub(const(m), E_VAR)
    ell_minus_f = sub(const(ell), F_VAR)
    denominator = mul(mul(m_minus_e, ell_minus_f), z0)

    b_numerator = sub(
        add(
            mul(mul(mul(scale(pa, p), E_VAR), ell_minus_f), z0),
            mul(mul(mul(scale(pb, p), F_VAR), m_minus_e), z0),
        ),
        mul(
            sub(mul(add(t, const(d)), K), scale(A0, d)),
            mul(m_minus_e, ell_minus_f),
        ),
    )

    certificate = sub(
        scale(mul(A0, power(b_numerator, 2)), 2 * ENERGY),
        scale(
            mul(power(sub(K, A0), 2), power(denominator, 2)),
            a * a + b * b,
        ),
    )
    return {
        "direct": direct,
        "B": b_numerator,
        "certificate": certificate,
    }


def certify_p_ge_2(a, b, max_depth=24):
    box = (
        (-Q(77 * a, 500), H[a]),
        (-Q(77 * b, 500), H[b]),
        (Q(0), Q(0)),
    )
    arrays = {
        name: integer_tensor(bernstein_fraction_tensor(poly, box))
        for name, poly in p_ge_2_polynomials(a, b).items()
    }

    stack = [(arrays, 0)]
    counts = Counter()
    nodes = 0
    attained_depth = 0
    while stack:
        current, depth = stack.pop()
        nodes += 1
        attained_depth = max(attained_depth, depth)

        if tensor_bounds(current["direct"])[0] > 0:
            counts["direct"] += 1
            continue
        if (
            tensor_bounds(current["B"])[0] > 0
            and tensor_bounds(current["certificate"])[0] > 0
        ):
            counts["stationarity"] += 1
            continue

        assert depth < max_depth, (
            a,
            b,
            depth,
            {name: tensor_bounds(value) for name, value in current.items()},
        )
        axis = depth % 2
        children = [{}, {}]
        for name, tensor in current.items():
            left, right = split_integer_tensor(tensor, axis)
            children[0][name] = left
            children[1][name] = right
        stack.append((children[1], depth + 1))
        stack.append((children[0], depth + 1))

    return nodes, attained_depth, counts


EXPECTED_P_GE_2 = {
    (1, 1): (65, 12, Counter({"direct": 26, "stationarity": 7})),
    (1, 2): (51, 12, Counter({"direct": 19, "stationarity": 7})),
    (1, 3): (89, 13, Counter({"direct": 30, "stationarity": 15})),
    (2, 2): (83, 12, Counter({"direct": 27, "stationarity": 15})),
}
ACTUAL_P_GE_2 = {
    pair: certify_p_ge_2(*pair) for pair in EXPECTED_P_GE_2
}
assert ACTUAL_P_GE_2 == EXPECTED_P_GE_2


# ---------------------------------------------------------------------------
# p = 1: second-order tangent duality in (e,f,x).
# ---------------------------------------------------------------------------


def p_equals_1_polynomials(a, b, signs):
    m = N - a
    ell = N - b
    d = min(a, b)
    entropy_kappa = Q(1, 2) if (a, b) == (1, 4) else Q(5, 8)

    pa = subset_product(a, E_VAR)
    pb = subset_product(b, F_VAR)
    deficit_envelope = sub(const(2), add(pa, pb))
    K = sub(const(GAMMA), deficit_envelope)

    t = add(E_VAR, F_VAR)
    tau = sub(t, X_VAR)
    z = sub(const(1), tau)
    A0 = scale(power(z, N), v(m) * v(ell))
    direct = sub(A0, K)

    # R_lambda = K-A0+A0*lambda^2/2.
    r_lambda = add(
        sub(K, A0),
        scale(A0, TANGENT_LAMBDA**2 / 2),
    )

    m_minus_e = sub(const(m), E_VAR)
    ell_minus_f = sub(const(ell), F_VAR)
    ef_denominator = mul(m_minus_e, ell_minus_f)
    common_denominator = mul(ef_denominator, z)

    b_numerator = sub(
        add(
            mul(mul(mul(pa, E_VAR), ell_minus_f), z),
            mul(mul(mul(pb, F_VAR), m_minus_e), z),
        ),
        mul(
            sub(mul(add(tau, const(d)), K), scale(A0, d)),
            ef_denominator,
        ),
    )

    # B_lambda = B - (d/z) A0 lambda^2/2.
    b_lambda_numerator = sub(
        b_numerator,
        scale(
            mul(A0, ef_denominator),
            d * TANGENT_LAMBDA**2 / 2,
        ),
    )

    # F1 = B_lambda + (d/z) R_lambda.
    f1_numerator = add(
        b_lambda_numerator,
        scale(mul(r_lambda, ef_denominator), d),
    )

    # W = E_+ + F_+ - (tau+d)/z, cleared by the same denominator.
    w_numerator = neg(mul(add(tau, const(d)), ef_denominator))
    if signs[0] == "+":
        w_numerator = add(
            w_numerator,
            mul(mul(E_VAR, ell_minus_f), z),
        )
    if signs[1] == "+":
        w_numerator = add(
            w_numerator,
            mul(mul(F_VAR, m_minus_e), z),
        )

    # F2 = B_lambda - W R_lambda.
    f2_numerator = sub(
        b_lambda_numerator,
        mul(w_numerator, r_lambda),
    )

    tangent_certificate = sub(
        scale(
            mul(mul(A0, f1_numerator), f2_numerator),
            4 * TANGENT_SLOPE * entropy_kappa * ENERGY,
        ),
        mul(power(r_lambda, 2), power(common_denominator, 2)),
    )

    return {
        "direct": direct,
        "B_lambda": b_lambda_numerator,
        "F1": f1_numerator,
        "F2": f2_numerator,
        "certificate": tangent_certificate,
    }


def p_equals_1_box(a, b, signs):
    if (a, b) == (1, 4):
        if signs == "++":
            e_range = (Q(0), H[1])
            f_range = (Q(0), H[4])
        elif signs == "-+":
            e_range = (-Q(77, 500), Q(0))
            f_range = (Q(0), H[4])
        else:
            e_range = (Q(0), H[1])
            f_range = (-H[1], Q(0))
    else:
        if signs == "++":
            e_range = (Q(0), H[2])
            f_range = (Q(0), H[3])
        elif signs == "-+":
            e_range = (-H[3], Q(0))
            f_range = (Q(0), H[3])
        else:
            e_range = (Q(0), H[2])
            f_range = (-H[2], Q(0))

    # Criticality gives 0 <= x < e+f.  The rectangular upper bound below is
    # deliberately larger, so the certificate proves a stronger statement.
    x_range = (Q(0), e_range[1] + f_range[1])
    return e_range, f_range, x_range


def certify_p_equals_1(a, b, signs, max_depth=24):
    box = p_equals_1_box(a, b, signs)
    arrays = {
        name: integer_tensor(bernstein_fraction_tensor(poly, box))
        for name, poly in p_equals_1_polynomials(a, b, signs).items()
    }

    stack = [(arrays, 0)]
    counts = Counter()
    nodes = 0
    attained_depth = 0
    while stack:
        current, depth = stack.pop()
        nodes += 1
        attained_depth = max(attained_depth, depth)

        if tensor_bounds(current["direct"])[0] > 0:
            counts["direct"] += 1
            continue

        tangent_names = ("B_lambda", "F1", "F2", "certificate")
        if all(tensor_bounds(current[name])[0] > 0 for name in tangent_names):
            counts["tangent"] += 1
            continue

        assert depth < max_depth, (
            a,
            b,
            signs,
            depth,
            {name: tensor_bounds(value) for name, value in current.items()},
        )
        axis = depth % 3
        children = [{}, {}]
        for name, tensor in current.items():
            left, right = split_integer_tensor(tensor, axis)
            children[0][name] = left
            children[1][name] = right
        stack.append((children[1], depth + 1))
        stack.append((children[0], depth + 1))

    return nodes, attained_depth, counts


EXPECTED_P_EQUALS_1 = {
    (1, 4, "++"): (485, 18, Counter({"direct": 130, "tangent": 113})),
    (1, 4, "-+"): (59, 13, Counter({"direct": 23, "tangent": 7})),
    (1, 4, "+-"): (89, 15, Counter({"direct": 38, "tangent": 7})),
    (2, 3, "++"): (427, 17, Counter({"direct": 116, "tangent": 98})),
    (2, 3, "-+"): (51, 14, Counter({"direct": 21, "tangent": 5})),
    (2, 3, "+-"): (75, 14, Counter({"direct": 31, "tangent": 7})),
}
ACTUAL_P_EQUALS_1 = {
    key: certify_p_equals_1(*key) for key in EXPECTED_P_EQUALS_1
}
assert ACTUAL_P_EQUALS_1 == EXPECTED_P_EQUALS_1


# ---------------------------------------------------------------------------
# One-sided Hall cuts.
# ---------------------------------------------------------------------------

KAPPA_6 = v(5) * G(5)
assert KAPPA_6 == Q(6144, 390625) > GAMMA

SMALL_ENTROPY = (
    ENERGY * (1 - GAMMA) ** 2 * GAMMA
    - Q(N * N, 2) * KAPPA_6**2
)
assert SMALL_ENTROPY == Q(
    67858438997036543,
    83037656250000000000,
) > 0

# Exact large-entropy polynomial on 0 <= t <= 31/200:
# 2 C_* (1-gamma)^2 (1-t)^4
#   - gamma (6-15t+20t^2-15t^3+6t^4-t^5)^2.
one_minus_t = sub(const(1), T_VAR)
ratio_polynomial = add(
    add(
        add(
            add(
                add(const(6), scale(T_VAR, -15)),
                scale(power(T_VAR, 2), 20),
            ),
            scale(power(T_VAR, 3), -15),
        ),
        scale(power(T_VAR, 4), 6),
    ),
    scale(power(T_VAR, 5), -1),
)
large_entropy_poly = sub(
    scale(power(one_minus_t, 4), 2 * ENERGY * (1 - GAMMA) ** 2),
    scale(power(ratio_polynomial, 2), GAMMA),
)
large_tensor = bernstein_fraction_tensor(
    large_entropy_poly,
    ((Q(0), T_ONE_SIDED), (Q(0), Q(0)), (Q(0), Q(0))),
)
large_coefficients = list(large_tensor[1].values())
assert min(large_coefficients) == Q(
    47529302914242746885557519,
    537477120000000000000000000,
) > 0


print("All exact n=6 extension certificates passed.")
print("p>=2 subdivision audit:", ACTUAL_P_GE_2)
print("p=1 subdivision audit:", ACTUAL_P_EQUALS_1)
