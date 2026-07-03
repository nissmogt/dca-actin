#!/usr/bin/env python3
"""
pdb_register.py — fetch PDB structures, build heavy-atom distance matrices, and
register MSA alignment columns to PDB residue numbering by global alignment.

Replaces the repo's hardcoded per-system `shift` (e.g. shift=4 for 1cc8A) with a
robust sequence-registration step, so contact validation works for homolog
structures across all 138 systems.

Output: a PDB distance DataFrame whose i,j columns are 1-indexed *DCA alignment
column* indices (not raw PDB residue ordinals), so it merges directly with the
DCA dataframe on ['i','j'] as analysis/analysis_pipeline.py expects.
"""
import os
import sys
import urllib.request
import numpy as np
import pandas as pd
import Bio.PDB
from Bio import Align

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import data.tools.pdb as dpdb


def fetch_pdb(pdbid, outdir="pdb"):
    """Download a PDB file from RCSB (cached). Returns local path or None on failure."""
    pdbid = pdbid.lower()
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, f"{pdbid}.pdb")
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return path
    for base in ("https://files.rcsb.org/download/",):
        try:
            urllib.request.urlretrieve(f"{base}{pdbid}.pdb", path)
            if os.path.getsize(path) > 0:
                return path
        except Exception as e:
            print(f"  fetch {pdbid} failed: {e}")
    if os.path.exists(path):
        os.remove(path)
    return None


def sysid_to_pdb_chain(sysid):
    """1ctfA -> ('1ctf','A'); 3dqgA -> ('3dqg','A'). Chain = trailing letter."""
    return sysid[:4].lower(), sysid[4:]


def get_chain_residues(pdb_path, chain_id):
    """Return (residue_list, one_letter_seq) for standard AA residues with a CA."""
    parser = Bio.PDB.PDBParser(QUIET=True)
    struct = parser.get_structure("s", pdb_path)
    model = struct[0]
    if chain_id not in [c.id for c in model]:
        # fall back to first chain
        chain_id = list(model.get_chains())[0].id
    chain = model[chain_id]
    res = [r for r in chain.get_residues()
           if Bio.PDB.is_aa(r) and r.has_id("CA")]
    seq = dpdb.three2one([r.get_resname() for r in res])
    return res, seq, chain_id


def _aligner():
    a = Align.PairwiseAligner()
    a.mode = "global"
    a.open_gap_score = -10
    a.extend_gap_score = -0.5
    a.match_score = 2
    a.mismatch_score = -1
    return a


def register_columns_to_pdb(ref_aln_seq, pdb_seq):
    """Map DCA alignment column index (1..L) -> PDB residue ordinal (1..Npdb).

    ref_aln_seq : the reference (first) sequence of the *filtered* MSA, WITH gaps
                  (length L = number of DCA columns).
    pdb_seq     : one-letter PDB chain sequence (length Npdb).

    Returns (col2pdb, n_mapped, identity) where col2pdb[col] = pdb_ordinal (both
    1-indexed) for columns whose reference residue aligns to a PDB residue.
    """
    ref_ungapped = ref_aln_seq.replace("-", "").upper()
    # column (1..L) -> reference residue ordinal (1..len(ref_ungapped)) or None
    col2ref = {}
    rpos = 0
    for k, ch in enumerate(ref_aln_seq, start=1):
        if ch != "-":
            rpos += 1
            col2ref[k] = rpos

    aln = _aligner().align(ref_ungapped, pdb_seq)[0]
    # ref residue ordinal -> pdb residue ordinal via aligned blocks
    ref2pdb = {}
    ai = aln.aligned  # ((ref_blocks),(pdb_blocks)) 0-indexed half-open
    n_ident = 0
    for (rs, re), (ps, pe) in zip(ai[0], ai[1]):
        for off in range(re - rs):
            r_ord = rs + off + 1
            p_ord = ps + off + 1
            ref2pdb[r_ord] = p_ord
            if ref_ungapped[rs + off] == pdb_seq[ps + off]:
                n_ident += 1

    col2pdb = {c: ref2pdb[r] for c, r in col2ref.items() if r in ref2pdb}
    identity = n_ident / max(1, len(ref_ungapped))
    return col2pdb, len(col2pdb), identity


def build_registered_distance_df(sysid, ref_aln_seq, pdb_dir="pdb",
                                 heavy_atom=True):
    """Fetch PDB, compute heavy-atom min-distance matrix, and reindex it by DCA
    alignment column. Returns (df_pdb_cols, info) where df_pdb_cols has columns
    i,j (DCA cols, i<j), d, si, sj (PDB resnums), or (None, info) on failure.
    """
    pdbid, chain = sysid_to_pdb_chain(sysid)
    info = {"sysid": sysid, "pdb": pdbid, "chain": chain,
            "status": "ok", "n_mapped": 0, "identity": 0.0, "L": 0, "Npdb": 0}
    path = fetch_pdb(pdbid, pdb_dir)
    if path is None:
        info["status"] = "fetch_failed"
        return None, info

    res_list, pdb_seq, used_chain = get_chain_residues(path, chain)
    info["chain"] = used_chain
    info["Npdb"] = len(res_list)
    if len(res_list) == 0:
        info["status"] = "no_residues"
        return None, info

    col2pdb, n_mapped, identity = register_columns_to_pdb(ref_aln_seq, pdb_seq)
    info["n_mapped"] = n_mapped
    info["identity"] = round(identity, 3)
    info["L"] = len(ref_aln_seq)
    if n_mapped < 5:
        info["status"] = "registration_failed"
        return None, info

    # pdb ordinal (1..Npdb) -> DCA column
    pdb2col = {p: c for c, p in col2pdb.items()}

    # heavy-atom min distances between mapped residue pairs
    rows = []
    mapped_ord = sorted(pdb2col.keys())
    for a_i in range(len(mapped_ord)):
        for b_i in range(a_i + 1, len(mapped_ord)):
            oa, ob = mapped_ord[a_i], mapped_ord[b_i]
            ra, rb = res_list[oa - 1], res_list[ob - 1]
            if heavy_atom:
                d, _ = dpdb.calc_min_dist(ra, rb)
            else:
                d = dpdb.calc_ca_distance(ra, rb)
            ca, cb = pdb2col[oa], pdb2col[ob]
            i, j = (ca, cb) if ca < cb else (cb, ca)
            rows.append((i, j, d, ra.id[1], rb.id[1]))
    df = pd.DataFrame(rows, columns=["i", "j", "d", "si", "sj"])
    return df, info


if __name__ == "__main__":
    sysid = sys.argv[1] if len(sys.argv) > 1 else "1ctfA"
    with open(f"aln/{sysid}.fa") as f:
        f.readline(); ref = f.readline().strip()
    df, info = build_registered_distance_df(sysid, ref)
    print(info)
    if df is not None:
        print(df.head())
        print("contacts d<8, |i-j|>=5:",
              int(((df.d < 8) & ((df.j - df.i) >= 5)).sum()))
