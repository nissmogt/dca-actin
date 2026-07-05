import os
os.environ["OMP_NUM_THREADS"]="1"; os.environ["OPENBLAS_NUM_THREADS"]="1"
os.environ["MKL_NUM_THREADS"]="1"; os.environ["NUMBA_NUM_THREADS"]="4"
import sys,time,json,warnings; import numpy as np
warnings.simplefilter("ignore")
WS=os.environ["WS"]
from pydca.meanfield_dca import meanfield_dca as mfd
def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}",flush=True)
fasta=f"{WS}/dca-actin/data/aln/actin_profilin_paired.fasta"
t=time.time()
mf=mfd.MeanFieldDCA(fasta,"protein",pseudocount=0.5,seqid=0.8)
log(f"init; neff={float(mf.effective_num_sequences):.0f}")
res=mf.compute_sorted_DI()   # list of ((i,j), di) sorted; standard uniform-pseudocount DI
log(f"DI done {time.time()-t:.0f}s, {len(res)} pairs")
# APC correct
L=537
di={}
for (i,j),v in res: di[(i,j)]=v
# build symmetric, APC
S=np.zeros((L,L))
for (i,j),v in di.items(): S[i,j]=v; S[j,i]=v
col=S.sum(0)/(L-1); tot=S.sum()/(L*(L-1))
out=[]
for (i,j),v in di.items():
    out.append((i,j,v, v-col[i]*col[j]/tot))
arr=np.array(out)  # i,j,di,diapc
np.save(f"{WS}/handoff/paired_di.npy", arr)
json.dump({"neff":float(mf.effective_num_sequences),"L":L,"npairs":len(out)},
          open(f"{WS}/handoff/paired_di_meta.json","w"))
log(f"SAVED paired_di.npy ({len(out)} pairs)")
