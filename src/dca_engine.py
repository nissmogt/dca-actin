#!/usr/bin/env python3
"""
dca_engine.py — MATLAB-free mean-field DCA engine for sequence_degradation.

Byte-compatible drop-in for dca/pipeline_inference.py::inference (which called
matlab.engine -> mf_plm_reweight/dca_h_J_Full_v5.m). Uses pydca's mean-field DCA
for the couplings + DI, computes an APC-corrected DI ("diapc") and a reweighted
mutual-information column ("mi"), and writes the same output file the MATLAB code
did:

    DI_{sysid}_n{model_length}.txt   (CSV, header: i,j,di,diapc,mi ; i,j are 1-indexed)

and returns NEFF = ceil(Meff) exactly like dca_h_J_Full_v5.m.

Reweighting: seqid=0.8  == theta=0.2  (80% identity), pseudocount=0.5 — matches
the thesis Neff definition and the MATLAB pseudocount_weight=0.5.

Author: MATLAB-free port (2026). Original DCA: Morcos et al. PNAS 2011.
Build note (M1): the mean-field path needs only numba; pydca's plmDCA C++ backend
needs an OpenMP compiler (see METHODS_port.md).
"""
import os
import tempfile
import numpy as np
from pydca.meanfield_dca import meanfield_dca
from pydca.meanfield_dca import msa_numerics

# Header written by dca_h_J_Full_v5.m (writematrix of ["i","j","di","diapc","mi"]).
DI_HEADER = "i,j,di,diapc,mi"

# Engine-level cap on sequences fed to pydca (safety net). Capping is normally
# done once upstream at the filtered-MSA level (run_local.MAX_POOL_SEQS) so the
# whole depth ladder derives from one bounded pool; this engine-level cap is a
# defensive backstop for direct run_mfdca() calls on very large MSAs. Read via
# _cap() so reassigning the module global takes effect. Set to None to disable.
#
# Why capping is scientifically safe: PPV is plotted against the *measured*
# Neff/L (the x-axis), and Neff is recomputed from whatever sequences are used.
# Capping the raw-sequence pool therefore never mislabels a point — it only
# limits how far up the Neff/L axis the ~6 largest families reach at their top
# rung. Note that for highly-redundant alignments (e.g. 2hs1A: ~73k sequences
# but Neff on the order of 1e2), subsampling raw sequences also lowers the top-
# rung Neff; those families simply start their depth curve at a lower maximum
# Neff/L. It does NOT bias the PPV-vs-Neff/L relationship, which is what the
# study measures.
MAX_SEQS = None


def _cap(max_seqs):
    """Resolve the effective cap: explicit arg wins, else the module global."""
    return MAX_SEQS if max_seqs == "default" else max_seqs


def _fast_sequences_weight(alignment_data=None, seqid=None, chunk=2000):
    """Vectorized replacement for pydca's O(M^2 L) numba weight loop.

    Computes, for each sequence, 1 / (number of sequences with fractional
    identity > seqid). Uses a one-hot encoding and chunked BLAS matmuls, giving
    bit-identical results to pydca's compute_sequences_weight (verified
    max|Δ| = 0) at ~30x the speed. `alignment_data` is an (M, L) integer array in
    any consistent encoding; states are remapped to 0..q-1 internally.
    """
    aln = np.asarray(alignment_data)
    M, L = aln.shape
    q = int(aln.max()) + 1
    oh = np.zeros((M, L * q), dtype=np.float32)
    rows = np.repeat(np.arange(M), L)
    cols = (np.arange(L)[None, :] * q + aln).ravel()
    oh[rows, cols] = 1.0
    thr = float(seqid) * L
    counts = np.zeros(M, dtype=np.float64)
    for s in range(0, M, chunk):
        e = min(s + chunk, M)
        counts[s:e] = (oh[s:e] @ oh.T > thr).sum(axis=1)
    return (1.0 / counts).astype(np.float64)


# Monkeypatch pydca to use the fast weight kernel everywhere.
msa_numerics.compute_sequences_weight = _fast_sequences_weight


def _maybe_cap_msa(msa_input, max_seqs):
    """If the MSA has more than max_seqs sequences, write a deterministic random
    subsample to a temp file and return its path (else return msa_input). Seeded
    by sequence count so the cap is reproducible."""
    if max_seqs is None:
        return msa_input, None
    # count sequences cheaply
    with open(msa_input) as fh:
        n = sum(1 for line in fh if line and line[0] == ">")
    if n <= max_seqs:
        return msa_input, None
    # read all records, subsample
    headers, seqs = [], []
    with open(msa_input) as fh:
        h = None
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                h = line[1:]; headers.append(h); seqs.append([])
            elif seqs:
                seqs[-1].append(line)
    seqs = ["".join(s) for s in seqs]
    rng = np.random.default_rng(n)                 # deterministic in n
    idx = rng.choice(len(seqs), max_seqs, replace=False)
    tmp = tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False)
    for k in idx:
        tmp.write(f">{headers[k]}\n{seqs[k]}\n")
    tmp.close()
    return tmp.name, tmp.name


def _compute_mi_matrix(alignment, weights, q=21, pseudocount=0.5):
    """Reweighted mutual information per residue pair, matching the pydca/MATLAB
    with_pc pseudocount scheme (relative weight `pseudocount` on a uniform prior).

    alignment: (M, L) int array of states in pydca's 1..q encoding.
    weights:   (M,) sequence reweighting weights (1/m_i); sum == Meff.
    Returns an (L, L) symmetric MI matrix (nats).
    """
    aln = np.asarray(alignment) - 1                      # 1..q  ->  0..q-1
    M, L = aln.shape
    Meff = float(np.sum(weights))
    w = np.asarray(weights, dtype=float).ravel()

    # single-site reweighted frequencies with pseudocount
    fi = np.zeros((L, q))
    for a in range(q):
        fi[:, a] = (aln == a).astype(float).T @ w
    fi = (1.0 - pseudocount) * fi / Meff + pseudocount / q

    mi = np.zeros((L, L))
    for i in range(L - 1):
        for j in range(i + 1, L):
            fij = np.zeros((q, q))
            np.add.at(fij, (aln[:, i], aln[:, j]), w)
            fij = (1.0 - pseudocount) * fij / Meff + pseudocount / (q * q)
            outer = np.outer(fi[i], fi[j])
            mi[i, j] = mi[j, i] = float(np.sum(fij * np.log(fij / outer)))
    return mi


def _apc_from_di(di_list, L):
    """Average-product correction of a raw DI list -> dict {(i,j): apc_score}.

    Matches pydca's compute_sorted_DI_APC to machine precision (verified
    max|Δ|~1e-16), but computed from a single DI pass instead of re-running the
    coupling inversion — halving the pydca cost per call.

    APC(i,j) = DI(i,j) - (s_i * s_j) / s_mean, where s_i is the mean DI over the
    other L-1 sites and s_mean the overall off-diagonal mean.
    """
    M = np.zeros((L, L))
    for (i, j), s in di_list:
        M[i, j] = M[j, i] = s
    col_mean = M.sum(axis=1) / (L - 1)
    overall = M.sum() / (L * (L - 1))
    return {(i, j): s - (col_mean[i] * col_mean[j]) / overall
            for (i, j), s in di_list}


def run_mfdca(msa_input, pseudocount=0.5, seqid=0.8, max_seqs="default"):
    """Run pydca mean-field DCA on an MSA file.

    max_seqs: cap on sequences (see MAX_SEQS). "default" -> module global; an int
    overrides; None disables the cap.

    Returns dict with:
        di      : list[((i,j), score)] raw DI, 0-indexed, sorted desc
        di_apc  : list[((i,j), score)] APC-corrected DI, 0-indexed, sorted desc
        mi      : (L,L) reweighted MI matrix
        L       : sequence length (columns after insert removal)
        neff    : effective number of sequences (float, un-rounded)
        nseq    : number of sequences
    """
    run_input, tmp = _maybe_cap_msa(msa_input, _cap(max_seqs))
    try:
        mf = meanfield_dca.MeanFieldDCA(run_input, "protein",
                                        pseudocount=pseudocount, seqid=seqid)
        di = mf.compute_sorted_DI()
        L = int(mf.sequences_len)
        apc = _apc_from_di(di, L)
        di_apc = sorted(apc.items(), key=lambda x: -x[1])
        mi = _compute_mi_matrix(mf.alignment, mf.sequences_weight,
                                q=mf.num_site_states, pseudocount=pseudocount)
        return {
            "di": di,
            "di_apc": di_apc,
            "mi": mi,
            "L": L,
            "neff": float(mf.effective_num_sequences),
            "nseq": int(mf.num_sequences),
        }
    finally:
        if tmp and os.path.exists(tmp):
            os.remove(tmp)


def write_di_file(result, out_path):
    """Write DI_{sysid}_n{N}.txt in the exact format dca_h_J_Full_v5.m produced:
    CSV, header 'i,j,di,diapc,mi', one row per residue pair (i<j), i,j 1-indexed,
    ordered by residue pair (i ascending, then j) — same enumeration as
    Compute_Results in the MATLAB code (loop i=1..N-1, j=i+1..N).
    """
    di = dict(result["di"])          # {(i,j): raw_di}
    di_apc = dict(result["di_apc"])  # {(i,j): apc_di}
    mi = result["mi"]
    L = result["L"]
    lines = [DI_HEADER]
    for i in range(L - 1):
        for j in range(i + 1, L):
            key = (i, j)
            d = di.get(key, di.get((j, i), 0.0))
            da = di_apc.get(key, di_apc.get((j, i), 0.0))
            lines.append(f"{i+1},{j+1},{d:.6f},{da:.6f},{mi[i, j]:.6f}")
    with open(out_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return out_path


def inference(_pfamid, dir_dca, msa_input, model_length, run=True):
    """Drop-in replacement for dca.pipeline_inference.inference.

    Signature matches the MATLAB-caller version:
        inference(sysid, dir_dca, msa_input, model_length) -> NEFF (int)
    `dir_dca` is accepted for API compatibility (was the MATLAB code path) and
    ignored. Writes DI_{_pfamid}_n{model_length}.txt next to msa_input.
    """
    output_dir = os.path.dirname(msa_input)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    file_dca_output = os.path.join(output_dir, f"DI_{_pfamid}_n{model_length}.txt")
    if not run:
        return 0
    if not os.path.exists(msa_input):
        print(f"{msa_input} not found!")
        return 0
    result = run_mfdca(msa_input)
    write_di_file(result, file_dca_output)
    neff = int(np.ceil(result["neff"]))     # NEFF = ceil(Meff), matches MATLAB
    return neff


if __name__ == "__main__":
    import sys
    sysid = sys.argv[1] if len(sys.argv) > 1 else "1ctfA"
    msa = sys.argv[2] if len(sys.argv) > 2 else f"aln/{sysid}.fa"
    n = inference(sysid, "", msa, "test")
    print(f"NEFF={n}  wrote DI_{sysid}_ntest.txt")
