#!/usr/bin/env python3
"""Parse master_laurents.json (pySecDec result strings) -> masterLaurents.wl (Mathematica assoc
   {a1,..,a11} -> {m[-2],m[-1],m[0],m[1]} using the central REAL parts). High precision preserved."""
import json, re
res=json.load(open('master_laurents.json'))
# a pySecDec term: ((RE,IM) +/- (dRE,dIM))*eps^N   or   ((RE,IM) +/- (...))  [eps^0]
term=re.compile(r'\(\(([-\d.eE+]+),([-\d.eE+]+)\)\s*\+/-\s*\([^)]*\)\)(\*eps\^(-?\d+))?')
def laurent(s):
    d={}
    for m in term.finditer(s):
        re_val=m.group(1); power=int(m.group(4)) if m.group(4) else 0
        d[power]=re_val
    return d
out=[]
for k,v in res.items():
    if v=='FAILED': continue
    idx=k.split(',')
    L=laurent(v)
    m2=L.get(-2,'0'); m1=L.get(-1,'0'); m0=L.get(0,'0'); mp1=L.get(1,'0')
    # Mathematica: use `30 precision markers not needed; give as strings -> ToExpression
    out.append(f"{{{','.join(idx)}}} -> {{{m2},{m1},{m0},{mp1}}}")
open('masterLaurents.wl','w').write("<|\n"+",\n".join(out)+"\n|>\n")
print(f"wrote masterLaurents.wl: {len(out)} masters parsed")
