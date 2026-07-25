#!/usr/bin/env python3
"""Independent exact audit of the ten n=6 one-sided certificates.

This checker uses only Python's standard library and Fraction.  It rebuilds
all five small-entropy and five large-entropy polynomials directly from the
formulas, converts them to Bernstein form on the required intervals, and
performs exact dyadic de Casteljau subdivision.
"""
from fractions import Fraction as F
from math import comb

# Coefficients are stored in increasing powers.
def trim(a):
    a=list(a)
    while len(a)>1 and a[-1]==0: a.pop()
    return a

def add(a,b):
    n=max(len(a),len(b)); out=[F(0)]*n
    for i,x in enumerate(a): out[i]+=x
    for i,x in enumerate(b): out[i]+=x
    return trim(out)

def neg(a): return [-x for x in a]
def sub(a,b): return add(a,neg(b))
def scale(a,c):
    c=F(c); return trim([c*x for x in a])

def mul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): out[i+j]+=x*y
    return trim(out)

def power(a,n):
    out=[F(1)]; base=a
    while n:
        if n&1: out=mul(out,base)
        base=mul(base,base); n//=2
    return out

def divide_by_x(a):
    if not a or a[0] != 0: raise RuntimeError(f'not divisible by x: constant={a[0] if a else None}')
    return trim(a[1:])

def affine_substitute(a, lo, hi):
    # p(lo+(hi-lo)y)
    lo, h = F(lo), F(hi)-F(lo)
    out=[F(0)]
    affine=[lo,h]
    for i,c in enumerate(a):
        out=add(out,scale(power(affine,i),c))
    return trim(out)

def power_to_bernstein_on_interval(a,lo,hi):
    p=affine_substitute(a,lo,hi)
    d=len(p)-1
    out=[]
    for i in range(d+1):
        out.append(sum((p[j]*F(comb(i,j),comb(d,j)) for j in range(i+1)),F(0)))
    return out

def split_half(c):
    d=len(c)-1; levels=[list(c)]
    for _ in range(d):
        prev=levels[-1]
        levels.append([(prev[i]+prev[i+1])/2 for i in range(len(prev)-1)])
    left=[levels[r][0] for r in range(d+1)]
    right=[levels[d-i][i] for i in range(d+1)]
    return left,right

def certify(c):
    stack=[(c,0)]; nodes=leaves=depth=0; minimum=None
    while stack:
        cur,d=stack.pop(); nodes+=1; depth=max(depth,d)
        mn,mx=min(cur),max(cur)
        if mn>0:
            leaves+=1; minimum=mn if minimum is None or mn<minimum else minimum
            continue
        if mx<=0 or d>=20:
            raise RuntimeError(f'failed at depth {d}: min={mn}, max={mx}')
        l,r=split_half(cur); stack.extend(((r,d+1),(l,d+1)))
    if minimum is None: raise RuntimeError('no positive leaf')
    return nodes,leaves,depth,minimum

GAMMA=F(5,324)
KAPPA=F(6144,390625)
E=F(11,50)
ENERGY=F(39750390625,5253139008)
T=[F(0),F(1)]
ONE_MINUS_T=[F(1),F(-1)]

EXPECTED_SMALL={
 1:(3,2,1,F(995279,28125000000)),
 2:(5,3,2,F(2374594200709185001,36771973056000000000000)),
 3:(7,4,3,F(506999743896537994129,8618431185000000000000000)),
 4:(9,5,4,F(1880953471440335502781485790459,40164800436633600000000000000000000)),
 5:(11,6,5,F(60824972467849735194923936599,2510300027289600000000000000000000)),
}
EXPECTED_LARGE=(1,1,0,F(1,864))

def main():
    for b in range(1,6):
        p=6-b
        q=ONE_MINUS_T
        column=mul(power(add([F(1)],scale(T,F(p,b))),b),power(q,p))
        small=add(mul(power(q,2),add(sub([F(1)],column),add(scale(power(q,6),KAPPA),[-GAMMA]))),scale(power(T,2),ENERGY))
        sr=certify(power_to_bernstein_on_interval(small,F(0),E/p))
        if sr!=EXPECTED_SMALL[b]: raise RuntimeError(('small',b,sr,EXPECTED_SMALL[b]))

        numerator=add(sub([F(1)],column),add(scale(mul(T,power(q,2)),F(3,32)),scale(sub(power(q,6),[F(1)]),GAMMA)))
        large=divide_by_x(numerator)
        lr=certify(power_to_bernstein_on_interval(large,F(0),E/p))
        if lr!=EXPECTED_LARGE: raise RuntimeError(('large',b,lr,EXPECTED_LARGE))
        print(f'independent n=6 one-sided b={b}: small={sr[:3]}, large={lr[:3]}')
    print('ALL 10 N=6 UNIVARIATE CERTIFICATES INDEPENDENTLY MATCHED')

if __name__=='__main__': main()
