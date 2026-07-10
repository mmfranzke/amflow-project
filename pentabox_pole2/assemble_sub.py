#!/usr/bin/env python3
"""Assemble the SUB-SECTOR contribution to poles(I[1]) = Sum_{j != top} c_j(d) M_j, in pure Python
(sympy for exact rational coeffs + Gamma-function bubble-chains; pySecDec floats for the rest).
No Mathematica needed. Compares [1/eps^2] to c_-2 = -J0/2.

Master Laurents:
  - bubble-chains (23 sectors): exact, from the massless bubble G(a,b) (computed here);
  - everything else: pySecDec, from master_laurents.json;
  - the 4 top-sector (255) masters: SKIPPED (this is the decisive isolation test).
"""
import re, json
import sympy as sp
from sympy import gamma, Rational, series, sympify

d, eps = sp.symbols('d eps')

# ---- massless bubble G(a,b) and the 23 exact bubble-chain sectors --------------------------------
def G(a,b): return gamma(a+b-d/2)*gamma(d/2-a)*gamma(d/2-b)/(gamma(a)*gamma(b)*gamma(d-a-b))
# scales are -invariant (>0 at the Euclidean point): from momentum routing (see pole2_assembly.wl)
mm2,s12,m23,s123,mm5,s15,mm3,mm4,s34 = 5,7,14,13,24,26,45,3,33
BC = {  # sector -> ('B1',sa,sb) or ('B2',sa)
 53:('B1',s12,mm5),54:('B1',mm2,mm5),57:('B1',s123,mm5),58:('B1',m23,mm5),60:('B1',mm3,mm5),
 85:('B1',s12,s123),86:('B1',mm2,s123),89:('B1',s123,s123),90:('B1',m23,s123),92:('B1',mm3,s123),
 101:('B1',s12,mm4),102:('B1',mm2,mm4),105:('B1',s123,mm4),106:('B1',m23,mm4),108:('B1',mm3,mm4),
 148:('B2',s12),152:('B2',s123),161:('B2',mm5),162:('B2',s15),164:('B2',s34),168:('B2',mm4),
 194:('B2',m23),196:('B2',mm3)}
def bubblechain_laurent(sec):
    t=BC[sec]
    if t[0]=='B1': expr=G(1,1)**2 * t[1]**(-eps) * t[2]**(-eps)
    else:          expr=G(1,1)*G(eps,1) * t[1]**(1-2*eps)
    s=series(expr.subs(d,4-2*eps),eps,0,2).removeO()
    return {p: complex(s.coeff(eps,p)) for p in (-2,-1,0)}

def sector_of(idx): return sum(1<<i for i,a in enumerate(idx[:8]) if a>0)

# ---- parse pySecDec results ---------------------------------------------------------------------
term=re.compile(r'\(\(([-\d.eE+]+),([-\d.eE+]+)\)\s*\+/-\s*\([^)]*\)\)(\*eps\^(-?\d+))?')
def psd_laurent(s):
    d_={}
    for m in term.finditer(s):
        p=int(m.group(4)) if m.group(4) else 0
        d_[p]=float(m.group(1))
    return d_
psd=json.load(open('master_laurents.json'))
psd={k:v for k,v in psd.items() if v!='FAILED'}

# ---- parse kira_target.m: I[1] -> sum of pentabox[idx]*(coeff(d)) --------------------------------
def load_reduction():
    txt=open('kira_target.m').read()
    terms=[]
    for m in re.finditer(r'pentabox\[([-\d,]+)\]\*\((.*?)\)\s*(?=\n \+|\n\}|\Z)', txt, re.S):
        idx=tuple(int(x) for x in m.group(1).split(','))
        coeff=sympify(m.group(2).replace('^','**'))
        terms.append((idx,coeff))
    return terms

def as_series(coeff_or_laur, kind):
    """Return dict power->complex for a coeff (sympy in d) or a master Laurent dict."""
    if kind=='coeff':
        s=series(coeff_or_laur.subs(d,4-2*eps),eps,0,2).removeO()
        return {p: complex(s.coeff(eps,p)) for p in range(-2,2)}
    return coeff_or_laur

def main():
    terms=load_reduction()
    # accumulate I_sub as dict power(eps)->complex, keeping -3..1
    acc={p:0j for p in range(-3,2)}
    top=0; missing=[]
    for idx,coeff in terms:
        s=sector_of(idx)
        if s==255: top+=1; continue
        if s in BC: M=bubblechain_laurent(s)
        else:
            k=','.join(map(str,idx))
            if k not in psd: missing.append(idx); continue
            M=psd_laurent(psd[k])
        C={p:complex(v) for p,v in as_series(coeff,'coeff').items()}
        for a,cv in C.items():
            for b,mv in M.items():
                if a+b in acc: acc[a+b]+=cv*mv
    print(f"top masters skipped: {top}; missing (unevaluated): {len(missing)}")
    if missing: print("  first missing:", missing[:6])
    P3,P2,P1 = acc[-3].real, acc[-2].real, acc[-1].real
    cm2,cm1 = -1.914807547604452136e-5, 2.197205645392600987e-4
    print(f"[1/eps^3] sub sum (expect ~0): {P3:.6e}")
    print(f"[1/eps^2] = P2_sub          : {P2:.10e}")
    print(f"[1/eps]   = P1_sub          : {P1:.10e}")
    print(f"c_-2 = {cm2:.10e}   P2_sub - c_-2 = {P2-cm2:.4e}   |ratio|={abs((P2-cm2)/cm2):.4e}")
    print("=> P2_sub ~ c_-2 : poles isolated to sub-sectors (top net 0).  O(1): top masters contribute.")

if __name__=='__main__': main()
