import os
os.environ["OMP_NUM_THREADS"]="4"; os.environ["OPENBLAS_NUM_THREADS"]="4"; os.environ["MKL_NUM_THREADS"]="4"
os.environ["NUMBA_NUM_THREADS"]="4"
import sys, time, json, warnings
import numpy as np, pandas as pd
warnings.simplefilter("ignore")
WS=os.environ["WS"]; HO=os.path.join(WS,"handoff"); sys.path.insert(0, os.path.join(WS,"dca-actin/src"))
import sasa_dca as sd
from pydca.meanfield_dca import msa_numerics as mn
from Bio.PDB import MMCIFParser
from Bio.PDB.SASA import ShrakeRupley
def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)

fi=np.load(f"{HO}/cache_fi.npy"); reg_fij=np.load(f"{HO}/cache_regfij.npy")
meta=json.load(open(f"{HO}/cache_q.json")); q=meta["q"]; neff=meta["neff"]; L=410; qm=q-1
iu,ju=np.triu_indices(L,1); P=len(iu)
log(f"cache loaded, P={P}")

def fn_apc(theta_vec, tag):
    t=time.time()
    reg_fi=sd._reg_single_per_site(fi, np.asarray(theta_vec,float), q)
    corr=mn.construct_corr_mat(reg_fi, reg_fij, L, q)
    J=mn.compute_couplings(corr).reshape(L,qm,L,qm)
    B=J[iu,:,ju,:]
    B=B - B.mean(1,keepdims=True) - B.mean(2,keepdims=True) + B.mean((1,2),keepdims=True)
    FN=np.sqrt((B**2).sum((1,2)))
    S=np.zeros((L,L)); S[iu,ju]=FN; S[ju,iu]=FN
    col=S.sum(0)/(L-1); tot=S.sum()/(L*(L-1))
    apc=FN - col[iu]*col[ju]/tot
    log(f"[{tag}] FN-APC {time.time()-t:.0f}s")
    return np.column_stack([iu,ju,apc])

# monomer-exposure theta
col2res=pd.read_csv(f"{WS}/dca-actin/results/col2res.csv")
c2r={int(r.dca_col_1based)-1:int(r.actin_resnum) for _,r in col2res.iterrows() if pd.notna(r.actin_resnum) and str(r.actin_resnum)!=""}
s=MMCIFParser(QUIET=True).get_structure("s",f"{WS}/dca-actin/data/pdb/8d13.cif")[0]
for cid in [c.id for c in s]:
    if cid!="A": s.detach_child(cid)
ShrakeRupley().compute(s, level="R")
mono={r.id[1]:r.sasa for r in s["A"] if r.id[0]==" "}
theta_sasa=sd.theta_from_exposure(sd.exposure_percentile({c:mono[rn] for c,rn in c2r.items() if rn in mono},L),0.3,0.7)

np.save(f"{HO}/fnapc_plain.npy", fn_apc(np.full(L,0.5),"PLAIN"))
np.save(f"{HO}/fnapc_sasa.npy",  fn_apc(theta_sasa,"SASA"))
json.dump({"neff":neff}, open(f"{HO}/fnapc_meta.json","w"))
log("SAVED both FN-APC arrays")
