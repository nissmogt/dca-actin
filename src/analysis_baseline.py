#!/usr/bin/env python3
"""
analysis_baseline.py — reproduce the Chapter 4 actin baseline end to end.

Runnable pipeline (each step writes into results/):
  1. Parse the PF00022 Stockholm alignment -> match-state FASTA (L columns).
  2. Register DCA columns to 8D13 chain-A actin numbering (global alignment).
  3. Run mean-field DCA (src/dca_engine.py) -> DI table ranked by diapc.
  4. Monomer PPV vs 8D13 chain-A contacts (heavy-atom min-dist <=8 A, |i-j|>=5).
  5. Inter-protomer distances from the 3-protomer 8D13 assembly.
  6. Z-score (diapc background) at the thesis ladder [12,10,9,8,5.6,4.5,4,3.5,2.5,1].

Conventions: rank by APC-corrected DI (diapc); seqid 0.8; pseudocount 0.5;
contacts heavy-atom min-dist <=8 A, |i-j|>=5. Thread-pin BEFORE numpy import.

Usage:
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=4 \
        python src/analysis_baseline.py \
        --sto data/aln/PF00022_full.sto --cif data/pdb/8d13.cif --outdir results

Reference values (docs/BASELINE_RESULTS.md): NEFF=4564, Neff/L=11.13;
monomer PPV@25=0.92,@50=0.90; 0 inter-protomer contacts in top-410;
GLU195-GLY366 inter-dist 15.81 A (thesis 16.5). Inference is ~20 min
single-process on an M1 (L=410, ~35k seqs, reweighting-dominated).
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import sys, csv, json, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

CONTACT_THRESH = 8.0
MIN_SEP = 5
Z_LADDER = [12, 10, 9, 8, 5.6, 4.5, 4, 3.5, 2.5, 1]


# ---------- 1. Stockholm -> match-state FASTA ----------
def parse_stockholm_to_fasta(sto_path, fasta_out, max_gap=0.5):
    seqs, order = {}, []
    with open(sto_path) as f:
        for line in f:
            if line.startswith("#") or line.startswith("//") or not line.strip():
                continue
            name, seq = line.rstrip("\n").split(None, 1)
            seqs.setdefault(name, []).append(seq)
            if name not in order:
                order.append(name)
    for k in seqs:
        seqs[k] = "".join(seqs[k])
    L_aln = len(seqs[order[0]])
    arr = np.frombuffer("".join(seqs[n] for n in order).encode("ascii"),
                        dtype=np.uint8).reshape(len(order), L_aln)
    # Pfam convention: insert columns contain '.' (46) or lowercase (97-122)
    col_has_insert = ((arr == 46) | ((arr >= 97) & (arr <= 122))).any(axis=0)
    match_cols = np.where(~col_has_insert)[0]
    sub = arr[:, match_cols].copy()
    low = (sub >= 97) & (sub <= 122); sub[low] -= 32          # uppercase
    sub[sub == 46] = 45                                        # '.' -> '-'
    keep = (sub == 45).mean(axis=1) <= max_gap                 # drop >50%-gap fragments
    with open(fasta_out, "w") as fo:
        for i, n in enumerate(order):
            if keep[i]:
                fo.write(">" + n + "\n" + sub[i].tobytes().decode("ascii") + "\n")
    L = len(match_cols)
    kept_idx = np.where(keep)[0]
    subk = sub[kept_idx]
    # consensus over kept sequences (most common non-gap residue per column)
    cons = []
    for j in range(L):
        col = subk[:, j]; col = col[col != 45]
        if len(col) == 0:
            cons.append("X"); continue
        vals, counts = np.unique(col, return_counts=True)
        cons.append(chr(vals[counts.argmax()]))
    return L, int(keep.sum()), "".join(cons)


# ---------- 2. structure parsing + column->residue registration ----------
def _heavy_by_res(chain):
    out = {}
    for r in chain:
        if r.id[0] != " ":
            continue
        coords = [a.coord for a in r if a.element != "H"]
        if coords:
            out[r.id[1]] = np.array(coords)
    return out


def _min_dist(a, b):
    return np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(-1)).min()


def load_structure(cif_path):
    from Bio.PDB import MMCIFParser
    from Bio.PDB.Polypeptide import is_aa
    import warnings; warnings.simplefilter("ignore")
    model = MMCIFParser(QUIET=True).get_structure("s", cif_path)[0]
    return model


def chain_sequence(model, cid):
    from Bio.Data.IUPACData import protein_letters_3to1
    res = []
    for r in model[cid]:
        if r.id[0] != " ":
            continue
        aa = protein_letters_3to1.get(r.resname.capitalize())
        if aa:
            res.append((r.id[1], aa))
    return [n for n, _ in res], "".join(a for _, a in res)


def register(consensus, seqA, numsA, out_csv, struct_map):
    from Bio import pairwise2
    aln = pairwise2.align.globalms(consensus, seqA, 2, -1, -8, -0.5,
                                   one_alignment_only=True)[0]
    ci = si = -1; col2res = {}; ident = naligned = 0
    for ca, cb in zip(aln.seqA, aln.seqB):
        if ca != "-": ci += 1
        if cb != "-": si += 1
        if ca != "-" and cb != "-":
            col2res[ci] = numsA[si]; naligned += 1
            if ca == cb: ident += 1
    with open(out_csv, "w", newline="") as fo:
        w = csv.writer(fo); w.writerow(["dca_col_1based", "actin_resnum", "consensus_aa", "struct_aa"])
        for c in range(len(consensus)):
            rn = col2res.get(c)
            w.writerow([c + 1, rn if rn is not None else "", consensus[c],
                        struct_map.get(rn, "") if rn else ""])
    return col2res, ident / max(naligned, 1)


# ---------- 4/5. distance matrices ----------
def intra_distmap(heavy):
    rn = sorted(heavy); d = {}
    for x in range(len(rn)):
        for y in range(x + 1, len(rn)):
            d[(rn[x], rn[y])] = _min_dist(heavy[rn[x]], heavy[rn[y]])
    return d


def inter_distmap(chains):
    d = {}
    for X, Y in [(a, b) for a in chains for b in chains if a is not b]:
        for rx, Ax in X.items():
            for ry, Ay in Y.items():
                k = (min(rx, ry), max(rx, ry)); v = _min_dist(Ax, Ay)
                if k not in d or v < d[k]:
                    d[k] = v
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sto", default="data/aln/PF00022_full.sto")
    ap.add_argument("--cif", default="data/pdb/8d13.cif")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--fasta", default="data/aln/actin_PF00022.fasta")
    ap.add_argument("--skip-inference", action="store_true",
                    help="reuse an existing DI table instead of recomputing")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    import pandas as pd

    # 1. alignment
    L, nseq, consensus = parse_stockholm_to_fasta(a.sto, a.fasta)
    print(f"[1] L={L} match-state cols, {nseq} seqs -> {a.fasta}")

    # 2. registration
    model = load_structure(a.cif)
    numsA, seqA = chain_sequence(model, "A")
    struct_map = dict(zip(numsA, seqA))
    col2res, ident = register(consensus, seqA, numsA,
                              os.path.join(a.outdir, "col2res.csv"), struct_map)
    print(f"[2] registration identity {ident:.3f}; GLU195->col {[c+1 for c,r in col2res.items() if r==195]}, "
          f"GLY366->col {[c+1 for c,r in col2res.items() if r==366]}")

    # 3. inference
    import dca_engine
    di_path = os.path.join(a.outdir, f"DI_actin_n{L}.txt")
    if not (a.skip_inference and os.path.exists(di_path)):
        neff = dca_engine.inference("actin", "", a.fasta, L)
        os.replace(os.path.join(os.path.dirname(a.fasta), f"DI_actin_n{L}.txt"), di_path)
        print(f"[3] NEFF={neff}, Neff/L={neff/L:.2f} -> {di_path}")
    di = pd.read_csv(di_path)

    # 4. monomer PPV
    chA = _heavy_by_res(model["A"]); dmat = intra_distmap(chA)
    c2r = {c + 1: r for c, r in col2res.items()}
    di["sep_col"] = (di.j - di.i).abs()
    dc = di[di.sep_col >= MIN_SEP].copy()

    def rp(i, j):
        ri, rj = c2r.get(int(i)), c2r.get(int(j))
        return None if ri is None or rj is None else (min(ri, rj), max(ri, rj))

    dc["respair"] = [rp(i, j) for i, j in zip(dc.i, dc.j)]
    dc = dc[dc.respair.notna()].copy()
    dc = dc[dc.respair.map(lambda p: p[1] - p[0]) >= MIN_SEP].copy()
    dc["dmin_intra"] = dc.respair.map(lambda p: dmat.get(p, np.nan))
    dc = dc.sort_values("diapc", ascending=False).reset_index(drop=True)
    dc["intra_contact"] = dc.dmin_intra <= CONTACT_THRESH
    ppv = {n: round(dc.head(n).intra_contact.mean(), 3) for n in (25, 50, 100, 200)}
    pd.DataFrame([{"top_n": n, "ppv": round(dc.head(n).intra_contact.mean(), 4),
                   "n_contacts": int(dc.head(n).intra_contact.sum())}
                  for n in range(10, L + 1, 10)]).to_csv(
        os.path.join(a.outdir, "monomer_ppv.csv"), index=False)
    print(f"[4] monomer PPV {ppv}")

    # 5. inter-protomer
    chains = [_heavy_by_res(model[c]) for c in ("A", "B", "C")]
    inter = inter_distmap(chains)
    dc["dmin_inter"] = dc.respair.map(lambda p: inter.get(p, np.nan))
    dc["inter_contact"] = dc.dmin_inter <= CONTACT_THRESH
    dc["filament_specific"] = dc.inter_contact & (~dc.intra_contact)
    top410 = dc.head(L)
    print(f"[5] top-{L}: {int(top410.intra_contact.sum())} intra-contacts, "
          f"{int(top410.filament_specific.sum())} filament-specific inter-contacts")

    # 6. Z-score
    mu, sd = float(di.diapc.mean()), float(di.diapc.std())
    dc["Z"] = (dc.diapc - mu) / sd
    zrows = []
    for z in Z_LADDER:
        sel = dc[dc.Z >= z]
        zrows.append({"Z": z, "n_pairs": len(sel),
                      "intra_contacts": int(sel.intra_contact.sum()),
                      "inter_contacts": int(sel.inter_contact.sum()),
                      "filament_specific": int(sel.filament_specific.sum())})
    pd.DataFrame(zrows).to_csv(os.path.join(a.outdir, "interface_zscores.csv"), index=False)
    dc.drop(columns=["respair"]).assign(
        res_i=dc.respair.map(lambda p: p[0]), res_j=dc.respair.map(lambda p: p[1])
    ).to_csv(os.path.join(a.outdir, "interface_pairs.csv"), index=False)
    g = dc[(dc.i == [c+1 for c,r in col2res.items() if r==195][0]) &
           (dc.j == [c+1 for c,r in col2res.items() if r==366][0])]
    if len(g):
        print(f"[6] GLU195-GLY366: Z={g.Z.iloc[0]:.2f}, inter-dist {g.dmin_inter.iloc[0]:.2f} A")
    print("done ->", a.outdir)


if __name__ == "__main__":
    main()
