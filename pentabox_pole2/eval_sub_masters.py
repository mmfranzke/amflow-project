import re,os,json,subprocess,time
from pySecDec import LoopIntegralFromPropagators, loop_package
from pySecDec.integral_interface import IntegralLibrary
props=['l2**2','(l2+p1)**2','(l2+p1+p2)**2','(l2+p1+p2+p3)**2','l1**2',
 '(l1-p1-p2-p3-p4)**2','(l1-p1-p2-p3)**2','(l1+l2)**2','(l1+p2)**2','(l1+p3)**2','(l2+p4)**2']
rr=[('p1*p1',0),('p2*p2',-5),('p3*p3',-45),('p4*p4',-3),('p1*p2',-1),('p1*p3','3/2'),
    ('p1*p4','1/2'),('p2*p3',18),('p2*p4',-12),('p3*p4','15/2')]
BC={53,54,57,58,60,85,86,89,90,92,101,102,105,106,108,148,152,161,162,164,168,194,196}
def sec(idx): return sum(1<<i for i,a in enumerate(idx[:8]) if a>0)
def nden(idx): return sum(1 for a in idx[:8] if a>0)
def masters():
    seen=set(); out=[]
    for ln in open('kira_target.m'):
        m=re.search(r'pentabox\[([-\d,]+)\]',ln)
        if m and '->' not in ln:
            idx=tuple(int(x) for x in m.group(1).split(','))
            if idx not in seen: seen.add(idx); out.append(idx)
    return out
allm=[x for x in masters() if sec(x)!=255 and sec(x) not in BC]  # skip 4 top + 23 exact bubble-chains
allm.sort(key=nden)   # fast (few-line) first
RES='master_laurents.json'; res=json.load(open(RES)) if os.path.exists(RES) else {}
def ev(idx):
    name='M'+''.join(('m' if x<0 else '')+str(abs(x)) for x in idx)
    subprocess.run(['rm','-rf',name])
    li=LoopIntegralFromPropagators(loop_momenta=['l1','l2'],external_momenta=['p1','p2','p3','p4'],
        propagators=props,powerlist=list(idx),replacement_rules=rr)
    loop_package(name=name,loop_integral=li,requested_orders=[3],real_parameters=[],processes=1)
    subprocess.run(['make','-C',name,'-j4'],check=True,capture_output=True)
    lib=IntegralLibrary(f'{name}/{name}_pylink.so')
    lib.use_Qmc(verbosity=0,transform='korobov3',minn=10**6,maxeval=10**7)
    _,_,s=lib(); subprocess.run(['rm','-rf',name]); return str(s)
print(f"{len(allm)} sub-sector masters to evaluate (top+bubblechains skipped)",flush=True)
for i,idx in enumerate(allm):
    k=','.join(map(str,idx))
    if k in res and res[k]!='FAILED': continue
    t0=time.time()
    try: res[k]=ev(idx); json.dump(res,open(RES,'w')); print(f"[{i+1}/{len(allm)}] {k} nden={nden(idx)} ({time.time()-t0:.0f}s) OK",flush=True)
    except Exception as e: res[k]='FAILED'; json.dump(res,open(RES,'w')); print(f"[{i+1}/{len(allm)}] {k} FAILED {str(e)[:80]}",flush=True)
print("DONE SUB MASTERS",flush=True)
