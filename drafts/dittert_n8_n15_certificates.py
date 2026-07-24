#!/usr/bin/env python3
"""Exact rational certificates for Dittert dimensions 8 through 15.

The script uses only Python integers and fractions.Fraction.  Every displayed
finite certificate in the accompanying manuscript is compared with a literal
expected rational number; no floating-point computation occurs.
"""
from fractions import Fraction as Q
from math import comb, factorial


def gamma(n):
    return Q(factorial(n), n**n)


def G(k):
    return Q(1) if k <= 1 else Q((k - 1)**(k - 1), k**(k - 1))


def v(k):
    return Q(1) if k == 0 else Q(factorial(k), k**k)


def kappa(n):
    return v(n - 1) * G(n - 1)


def M1(n, L, c):
    g = gamma(n)
    return L - g - Q(n*n, 4) * c * L*L / (1 - g)


ETA = {
    8: Q(1, 80), 9: Q(1, 135), 10: Q(1, 230),
    11: Q(1, 350), 12: Q(1, 650),
    13: Q(1, 1100), 14: Q(1, 1800),
}


def alpha(n, k):
    if k == 0:
        return Q(0)
    p = Q(k, n)
    if 2*k < n:
        return 2*n*(p + ETA[n])*(1 - p - ETA[n])
    return Q(2*k*(n-k), n)


EXPECTED_LOCALIZATION = {
    8: Q(4757, 836844800),
    9: Q(242489, 87087962025),
    10: Q(31109, 41313127850),
    11: Q(63007753811, 34945789841847500),
    12: Q(11651173, 90828753405000),
    13: Q(13072911571453, 366471344281458130000),
    14: Q(71438752181, 2511343447066440000),
}
EXPECTED_HALF_GAPS = {
    8: Q(9, 80), 9: Q(13, 270), 10: Q(11, 115),
    11: Q(82, 1925), 12: Q(319, 3900),
    13: Q(537, 14300), 14: Q(893, 12600),
}

for n in range(8, 15):
    g = gamma(n)
    localization = ETA[n]**2 - g / (2*n*(1-g))
    half_gap = Q(1, 2) - Q((n-1)//2, n) - ETA[n]
    assert localization == EXPECTED_LOCALIZATION[n] > 0
    assert half_gap == EXPECTED_HALF_GAPS[n] > 0
    for k in range(n):
        assert alpha(n, k) <= Q(n, 2)


# Dimensions 8, 9, 10.
T = {8: Q(1, 7), 9: Q(1, 10), 10: Q(1, 10)}
LAM = {8: Q(20), 9: Q(27), 10: Q(33)}


def sharp_block(n, a, b):
    p = n-a-b
    return v(n-a)*v(n-b)/v(p)


def M2(n, L, c):
    g = gamma(n)
    q = LAM[n]*L + (1-g)/c
    return L-g-Q(n*n, 4)*L*L/q


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
    tau_gap = T[n]**2 - n*g/(1-g)
    d0_gap = Q(1, 100) - Q(n, 2)*g/(1-g)
    assert tau_gap == EXPECTED_TAU_GAPS[n] > 0
    assert d0_gap == EXPECTED_D0_GAPS[n] > 0
    assert Q(1, 81) - g == EXPECTED_GAMMA_81_GAPS[n] > 0
    assert LAM[n] <= comb(n, 2) - comb(n, 3)*T[n]

    records = []
    for a in range(1, n):
        for b in range(a, n-a):
            p = n-a-b
            c = (alpha(n, a)+alpha(n, b))/(p*p)
            records.append((M2(n, sharp_block(n, a, b), c), a, b))
    assert len(records) == EXPECTED_LOW_COUNTS[n]
    assert all(m > 0 for m, _, _ in records)
    assert min(records) == EXPECTED_LOW_MINIMA[n]

    small = Q(2, 5)*(1-g)**2*g - Q(n*n, 2)*kappa(n)**2
    large = Q(4, 5)*(1-g)**2*Q(9, 10)**n - n*n*g
    assert small == EXPECTED_SMALL[n] > 0
    assert large == EXPECTED_LARGE[n] > 0


# Dimensions 11, 12, 13, 14.
def weak_block(n, a, b):
    return v(n-a)*G(n-b)**a


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
            c = (alpha(n, a)+alpha(n, b))/(p*p)
            L = kappa(n) if a == 0 else weak_block(n, a, b)
            records.append((M1(n, L, c), a, b))
    assert len(records) == EXPECTED_MID_COUNTS[n]
    assert all(m > 0 for m, _, _ in records)
    assert min(records) == EXPECTED_MID_MINIMA[n]

assert Q(1, 1000) - gamma(11) == Q(22308624601, 25937424601000) > 0

g = gamma(11)
cstar = 1 / ((1-g)**2/alpha(11, 1) + 1/alpha(11, 10))
assert cstar == Q(1540632960668022501507138780,
                  1671236290061169784540151329)
assert Q(93, 100) - cstar == Q(
    1361678908886539811520195597,
    167123629006116978454015132900) > 0
assert M1(11, kappa(11), Q(93, 100)) == Q(
    82169554295122875121722420238561072941763350843,
    656890501177898340235351562500000000000000000000000000) > 0


# Dimension 15.
g = gamma(15)
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

# Cleared-denominator gaps in the two coarse n=15 estimates.
assert 5*(10**12*299999) - 534*9*300000*10**9 == 58195000000000000
assert 14*(10**12*299999) - 3375*4*300000*10**9 == 149986000000000000

print("All exact rational certificates passed.")
