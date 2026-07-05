#!/usr/bin/env python3
"""
zscore.py — Z-score significance filter for DCA couplings.

Corrected reimplementation of nissmogt/dca-interface/zscore.py::zscore_calc
(thesis Ch. 2, Mehrabiani 2022). The thesis defines the Z of a DCA pair as its
coupling score standardized against a reference distribution:

    z = (score - mean(reference)) / std(reference)

and subsets the scored pairs by pair class (`which`): monomer (intra-chain,
chain_1 == chain_2), interface (inter-chain, chain_1 != chain_2), or all.

Original (dca-interface/zscore.py, lines ~94-108):

    def zscore_calc(data, n_pairs, score='fn_apc', which='monomer', stats=False):
        if which == 'monomer':   _df = data[data.chain_1 == data.chain_2]
        if which == 'interface': _df = data[data.chain_1 != data.chain_2]
        else:                    _df = data
        mean = np.mean(data[score]); std = np.std(data[score])
        _df['zscore'] = (_df[score] - mean) / std
        ...

NOTE — behavior change vs the original: the two `if`s above are NOT chained with
`elif`, so for which='monomer' the second `if` is False and its `else` branch
fires, overwriting `_df` back to the FULL data — i.e. which='monomer' is a
latent no-op in the original (Z ends up computed on all pairs, not the monomer
subset). This port fixes that: `_subset()` isolates the requested pair class via
proper elif/masking, so which='monomer'/'interface' return the intended subsets.
What carries over faithfully is the Z definition itself: mean/std are taken over
the FULL reference passed in, and the reference (score column + row set) is the
single knob that sets the Z scale. (The actin analysis only exercises
which='all', where original and port agree.)

SCORE-TYPE NOTE. The thesis computes Z on plmDCA Frobenius-norm (FN / CN-APC)
scores pooled across many systems; the local actin engine produces mean-field
APC-corrected DI (`diapc`) for one family. The METHOD is identical; absolute Z
values are NOT comparable across score types (see docs/ZSCORE_PORT.md).
"""
import numpy as np
import pandas as pd


def zscore_calc(data, score="diapc", which="all", reference=None):
    """Standardize DCA scores into Z against a reference distribution.

    Parameters
    ----------
    data : DataFrame with at least the `score` column. For monomer/interface
        subsetting it must also carry a pair-class label (see `which`).
    score : score column to standardize (default 'diapc' for the mean-field
        engine; the thesis used 'fn_apc'/'fn').
    which : 'all' | 'monomer' | 'interface' — which subset of `data` to return Z
        for. 'monomer' keeps intra-chain pairs, 'interface' inter-chain pairs;
        requires a boolean/2-value 'is_interface' column or 'chain_1'/'chain_2'.
    reference : optional DataFrame or 1-D array giving the reference score
        distribution. If None, the reference is `data[score]` itself (the thesis
        default — mean/std over the full passed-in table).

    Returns
    -------
    (df_out, mean, std) — df_out is the `which` subset with a 'zscore' column.
    """
    if reference is None:
        ref = np.asarray(data[score], dtype=float)
    elif isinstance(reference, pd.DataFrame):
        ref = np.asarray(reference[score], dtype=float)
    else:
        ref = np.asarray(reference, dtype=float)
    mean = float(np.mean(ref))
    std = float(np.std(ref))

    if which == "monomer":
        _df = _subset(data, interface=False)
    elif which == "interface":
        _df = _subset(data, interface=True)
    else:
        _df = data.copy()
    _df = _df.copy()
    _df["zscore"] = (_df[score].astype(float) - mean) / std
    return _df, mean, std


def _subset(data, interface):
    """Return intra- (interface=False) or inter-chain (interface=True) pairs."""
    if "is_interface" in data.columns:
        mask = data["is_interface"].astype(bool)
    elif {"chain_1", "chain_2"}.issubset(data.columns):
        mask = data["chain_1"] != data["chain_2"]
    else:
        raise ValueError("need 'is_interface' or 'chain_1'/'chain_2' to subset by pair class")
    return data[mask if interface else ~mask]


# --- reference builders (name the convention explicitly) ---------------------

def reference_within_family(di_df, score="diapc"):
    """Ch.3/current convention: reference = the full within-family score table."""
    return np.asarray(di_df[score], dtype=float)


def reference_monomeric_contacts(pairs_df, score="diapc", contact_col="intra_contact"):
    """Ch.2 prose convention: reference = scores of monomeric (intra-protomer)
    CONTACT pairs. Requires a boolean contact column."""
    return np.asarray(pairs_df.loc[pairs_df[contact_col].astype(bool), score], dtype=float)


def reference_monomeric_background(pairs_df, score="diapc", intra_dist_col="dmin_intra",
                                   contact_thresh=8.0):
    """Reference = scores of monomeric NON-contact pairs (intra-protomer,
    heavy-atom min-dist > contact_thresh) — a null/background distribution."""
    m = pairs_df[intra_dist_col].astype(float) > contact_thresh
    return np.asarray(pairs_df.loc[m, score], dtype=float)


if __name__ == "__main__":
    import json, sys
    df = pd.read_csv(sys.argv[1] if len(sys.argv) > 1 else "results/interface_pairs.csv")
    out, mu, sd = zscore_calc(df, score="diapc", which="all")
    print(json.dumps({"n": len(out), "mean": mu, "std": sd,
                      "max_Z": float(out.zscore.max())}, indent=2))
