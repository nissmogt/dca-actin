#!/usr/bin/env python3
"""
sasa_dca.py — SASA-conditioned mean-field DCA for the actin filament interface.

Two ways to bring solvent-accessible surface area (SASA) into the coevolutionary
model, run as a contrast (see docs/SASA_DCA.md):

  (A) FIELD-PRIOR RE-FIT  (the non-circular test)
      pydca's mean-field regularization shrinks each column's single-site
      statistics toward uniform with ONE scalar pseudocount theta:
          reg_fi[i,a] = theta/q + (1-theta) * fi[i,a]
      We replace theta with a PER-COLUMN theta_i tied to MONOMER surface
      exposure e_i (SASA of the isolated protomer, knowable without the
      filament). Buried core columns -> larger theta_i (down-weight the
      dominant monomer-fold coevolution); exposed surface columns -> smaller
      theta_i (trust the data where the interface can live). The modified
      reg_fi reshapes the correlation matrix C and therefore every coupling
      J = -C^-1 through the inversion — a genuine global re-fit, then DI is
      recomputed. Conditioning on monomer exposure (not on assembly burial)
      keeps the test honest: we never feed the interface into the prior.

  (B) DI RE-RANKING BASELINE  (partly circular)
      Keep the couplings from the plain fit; re-score each pair as
          DI'_ij = DI_ij * w(dSASA_i, dSASA_j)
      with w increasing in assembly burial. Buried inter-protomer pairs are
      boosted by construction, so this cannot test whether coevolution ENCODES
      the interface — only whether structure-informed re-weighting improves
      recovery. It is the comparison baseline for (A).

The engine's own DI path is reproduced exactly by the manual pipeline here: on
a 600-sequence subsample of the actin MSA (same L=410 columns), the manual
reg_fi -> construct_corr_mat -> compute_couplings -> compute_two_site_model_fields
-> compute_direct_info chain matched MeanFieldDCA.compute_sorted_DI to max|Δ|=0
over the first 2000 ranked pairs. The check is algebraic (independent of
sequence count), so the per-column theta is a legitimate re-fit of the same
model, not an approximation.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
from pydca.meanfield_dca import meanfield_dca as mfd, msa_numerics as mn


def exposure_percentile(mono_sasa_by_col, L):
    """Map a {col: monomer_SASA} dict to a per-column exposure percentile in
    [0,1] (0 = most buried, 1 = most exposed). Unmapped columns get 0.5."""
    vals = np.array([mono_sasa_by_col.get(i, np.nan) for i in range(L)], float)
    fin = np.isfinite(vals)
    pr = np.full(L, 0.5)
    order = np.argsort(np.argsort(vals[fin]))
    pr[np.where(fin)[0]] = order / max(order.max(), 1)
    return pr


def theta_from_exposure(expo_pct, theta_min=0.3, theta_max=0.7):
    """Per-column pseudocount: buried (expo=0) -> theta_max (shrink),
    exposed (expo=1) -> theta_min (trust). Reduces to a constant when
    theta_min == theta_max."""
    return theta_max - (theta_max - theta_min) * expo_pct


def _reg_single_per_site(fi, theta_vec, q):
    r = fi.copy()
    for i in range(r.shape[0]):
        t = theta_vec[i]
        r[i, :] = t / q + (1.0 - t) * fi[i, :]
    return r


def refit_di(msa_path, theta_vec, pseudocount=0.5, seqid=0.8):
    """Run the full mean-field DI pipeline with a PER-COLUMN single-site
    pseudocount theta_vec. Pair-frequency regularization keeps the scalar
    pseudocount (kept fixed so the change is isolated to the single-site prior).
    Returns dict {(i,j): DI} 0-indexed, plus L."""
    mf = mfd.MeanFieldDCA(msa_path, "protein", pseudocount=pseudocount, seqid=seqid)
    L = int(mf.sequences_len); q = int(mf.num_site_states)
    fi = mf.get_single_site_freqs()
    reg_fi = _reg_single_per_site(fi, np.asarray(theta_vec, float), q)
    reg_fij = mf.get_reg_pair_site_freqs()
    corr = mn.construct_corr_mat(reg_fi, reg_fij, L, q)
    J = mn.compute_couplings(corr)
    fields_ij = mn.compute_two_site_model_fields(J, reg_fi, L, q)
    di_flat = mn.compute_direct_info(couplings=J, fields_ij=fields_ij,
                                     reg_fi=reg_fi, seqs_len=L, num_site_states=q)
    di = {}; k = 0
    for i in range(L - 1):
        for j in range(i + 1, L):
            di[(i, j)] = di_flat[k]; k += 1
    return {"di": di, "L": L, "neff": float(mf.effective_num_sequences)}


def apc_correct(di, L):
    """APC-correct a {(i,j): score} DI dict; returns {(i,j): diapc}."""
    S = np.zeros((L, L))
    for (i, j), v in di.items():
        S[i, j] = S[j, i] = v
    col = S.mean(0); tot = S.mean()
    apc = {}
    for (i, j) in di:
        apc[(i, j)] = S[i, j] - (col[i] * col[j] / tot if tot else 0.0)
    return apc


def rerank_by_sasa(diapc, dsasa_by_col, alpha=1.0, scale=30.0):
    """Baseline (B): multiply diapc by a burial weight
    w = (1 + alpha * min(dSASA_i,dSASA_j)/scale). dsasa_by_col maps col->ΔSASA."""
    out = {}
    for (i, j), v in diapc.items():
        b = min(dsasa_by_col.get(i, 0.0), dsasa_by_col.get(j, 0.0))
        out[(i, j)] = v * (1.0 + alpha * b / scale)
    return out
