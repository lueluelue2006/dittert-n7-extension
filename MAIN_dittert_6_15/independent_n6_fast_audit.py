#!/usr/bin/env python3
"""Independent exact audit of all n=6 four-variable certificates.

This script reconstructs the polynomials from the formulas using a small
sparse Fraction polynomial algebra (no SymPy and no submitted generator),
then performs its own dense power-to-Bernstein conversion and dyadic
subdivision.  It checks the published family counts, node counts, depths,
and minimum terminal coefficients.
"""
from fractions import Fraction as F
from itertools import product
from math import comb, factorial
import numpy as np

NVAR=4
ZERO=(0,0,0,0)

def clean(p): return {m:c for m,c in p.items() if c}
def const(c):
    c=F(c); return {} if c==0 else {ZERO:c}
def var(i):
    m=[0]*NVAR; m[i]=1; return {tuple(m):F(1)}
def add(p,q):
    r=dict(p)
    for m,c in q.items(): r[m]=r.get(m,F(0))+c
    return clean(r)
def neg(p): return {m:-c for m,c in p.items()}
def sub(p,q): return add(p,neg(q))
def scale(p,c):
    c=F(c); return clean({m:c*a for m,a in p.items()})
def mul(p,q):
    r={}
    for m,a in p.items():
        for n,b in q.items():
            k=tuple(x+y for x,y in zip(m,n))
            r[k]=r.get(k,F(0))+a*b
    return clean(r)
def power(p,n):
    r=const(1); b=p
    while n:
        if n&1: r=mul(r,b)
        b=mul(b,b); n//=2
    return r

def v(k): return F(1) if k==0 else F(factorial(k),k**k)
def sharp(a,b):
    p=6-a-b; return v(6-a)*v(6-b)/v(p)

X,Z,U,V=(var(i) for i in range(4))
E=F(11,50); GAMMA=F(5,324)

def family(a,b,zeta):
    p=6-a-b; d=min(a,b); m=6-a; ell=6-b
    t=scale(X,2*E)
    e=add(sub(t,const(E)),mul(Z,sub(const(2*E),t)))
    f=sub(const(E),mul(Z,sub(const(2*E),t)))
    den=sub(const(p),t)
    a0=scale(power(den,6),sharp(a,b)/F(p**6))
    if p==1:
        cx=F(1,2)+F(1,3*m); cy=F(1,2)+F(1,3*ell)
    else:
        cx=cy=F(p*p,2*(a*a+b*b))
    bnum=add(const(p+d),scale(den,zeta))
    common=add(scale(add(t,const(d)),GAMMA-2),scale(a0,-d))
    common=add(common,scale(mul(den,sub(const(GAMMA-2),a0)),zeta))
    def side(comp,gsize,k,excess,spread,c):
        total=sub(const(gsize),excess)
        pref=power(add(const(1),scale(excess,F(1,comp))),comp)
        if k==gsize:
            high=scale(total,F(1,gsize))
            prod=mul(pref,power(high,gsize))
            prod_over=mul(pref,power(high,gsize-1))
            ent=F(0)
        else:
            high=mul(scale(total,F(1,gsize)),add(const(1),scale(spread,F(gsize-k,k))))
            low=mul(scale(total,F(1,gsize)),sub(const(1),spread))
            prod=mul(pref,mul(power(high,k),power(low,gsize-k)))
            prod_over=mul(pref,mul(power(high,k-1),power(low,gsize-k)))
            ent=F(1,k)-F(1,gsize)
        out=add(mul(prod,bnum),scale(mul(den,prod_over),-1))
        out=add(out,scale(mul(a0,add(const(d),scale(den,zeta))),-c*ent))
        return out
    rows={k:side(a,m,k,e,U,cx) for k in range(1,m+1)}
    cols={k:side(b,ell,k,f,V,cy) for k in range(1,ell+1)}
    return common,rows,cols

EXPECTED={
 (1,1):(25,165,6,F(711760783,88593750000)),
 (1,2):(20,276,7,F(9468769,6480000000)),
 (1,3):(15,483,14,F(58951253,118125000000)),
 (2,2):(16,416,13,F(13905244249,32256000000000)),
 (1,4):(10,352,14,F(6542213,472500000000)),
 (2,3):(12,558,15,F(206812176444269171,2388787200000000000000)),
}

def bernstein_grid(poly):
    deg=tuple(max((m[a] for m in poly),default=0) for a in range(NVAR))
    shape=tuple(d+1 for d in deg)
    arr=np.empty(shape,dtype=object); arr.fill(F(0))
    for m,c in poly.items(): arr[m]=c
    # Convert each coordinate line independently. x^j = sum_{i>=j} C(i,j)/C(d,j) B_i^d.
    for axis,d in enumerate(deg):
        if d==0: continue
        moved=np.moveaxis(arr,axis,-1)
        flat=moved.reshape((-1,d+1))
        out=np.empty_like(flat)
        coeff=[[F(comb(i,j),comb(d,j)) for j in range(i+1)] for i in range(d+1)]
        for r in range(flat.shape[0]):
            line=flat[r]
            for i in range(d+1):
                s=F(0)
                for j,c in enumerate(coeff[i]): s += line[j]*c
                out[r,i]=s
        arr=np.moveaxis(out.reshape(moved.shape),-1,axis)
    return arr

def split(arr,axis):
    moved=np.moveaxis(arr,axis,-1)
    d=moved.shape[-1]-1
    flat=moved.reshape((-1,d+1))
    lf=np.empty_like(flat); rf=np.empty_like(flat)
    for r in range(flat.shape[0]):
        temp=list(flat[r]); left=[F(0)]*(d+1); right=[F(0)]*(d+1)
        left[0]=temp[0]; right[d]=temp[d]
        for level in range(1,d+1):
            temp=[(temp[i]+temp[i+1])/2 for i in range(len(temp)-1)]
            left[level]=temp[0]; right[d-level]=temp[-1]
        lf[r]=left; rf[r]=right
    l=np.moveaxis(lf.reshape(moved.shape),-1,axis)
    r=np.moveaxis(rf.reshape(moved.shape),-1,axis)
    return l,r

def choose_axis(arr):
    best=F(-1); axis_best=0
    for a in range(arr.ndim):
        if arr.shape[a]<=1: continue
        score=max(abs(x) for x in np.diff(arr,axis=a).flat)
        if score>best: best=score; axis_best=a
    return axis_best

def certify(arr):
    stack=[(arr,0)]; nodes=0; depth=0; minimum=None
    while stack:
        box,d=stack.pop(); nodes+=1; depth=max(depth,d)
        vals=list(box.flat); mn=min(vals); mx=max(vals)
        if mn>0:
            minimum=mn if minimum is None or mn<minimum else minimum
            continue
        if mx<=0 or d>=40:
            raise RuntimeError(f'failed box at depth {d}: min={mn} max={mx}')
        a=choose_axis(box); left,right=split(box,a)
        stack.append((right,d+1)); stack.append((left,d+1))
    if minimum is None: raise RuntimeError('no positive leaf')
    return nodes,depth,minimum

def main():
    total=0
    for a,b in EXPECTED:
        p=6-a-b; zeta=14 if p>=2 else 8
        common,rows,cols=family(a,b,zeta)
        cases=nodes=maxdepth=0; minimum=None
        for row in rows.values():
            for col in cols.values():
                poly=neg(add(common,add(row,col)))
                result=certify(bernstein_grid(poly))
                cases+=1; nodes+=result[0]; maxdepth=max(maxdepth,result[1])
                minimum=result[2] if minimum is None or result[2]<minimum else minimum
        actual=(cases,nodes,maxdepth,minimum)
        if actual!=EXPECTED[(a,b)]:
            raise RuntimeError(((a,b),actual,EXPECTED[(a,b)]))
        print('independent exact tree matched', (a,b), actual[:3])
        total+=cases
    if total!=98: raise RuntimeError(total)
    print('ALL 98 N=6 FOUR-VARIABLE CERTIFICATES INDEPENDENTLY MATCHED')

if __name__=='__main__': main()
