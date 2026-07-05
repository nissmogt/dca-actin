# Z-reference port — Chapter 2 method into dca-actin

`src/zscore.py` ports `nissmogt/dca-interface/zscore.py::zscore_calc` — the
thesis Chapter 2 Z-score significance filter — into this repo.

## Method (Z definition unchanged from the thesis)
For a reference score distribution, the Z of a DCA pair is
`z = (score − mean(reference)) / std(reference)`, evaluated on a pair-class
subset (`which` = monomer / interface / all). The Z definition is transcribed
faithfully. One corrected behavior: in the original `zscore_calc` the `monomer`
and `interface` branches are separate `if`s (not `elif`), so `which='monomer'`
falls through to the `else` and silently reverts to the full table — a latent
no-op. The port fixes this so `which='monomer'/'interface'` return the intended
subsets; the actin analysis uses only `which='all'`, where the two agree.
Reference builders name each convention:
- `reference_within_family` — full within-family score table (Ch.3 / prior actin run).
- `reference_monomeric_contacts` — scores of monomeric (intra-protomer) contact pairs (Ch.2 prose).
- `reference_monomeric_background` — scores of monomeric non-contact pairs (a null background).

## Score-type caveat (why absolute Z ≠ the thesis 2.09)
The thesis computes Z on **plmDCA Frobenius-norm (FN / CN-APC)** scores pooled
across many systems. The local actin engine produces **mean-field APC-corrected
DI (`diapc`)** for a single family. The standardization method is identical, but
the score type and reference population differ, so absolute Z values are not
directly comparable. The port reproduces the *procedure*, not the *number*.

## GLU195–GLY366 across reference conventions
| reference | n(ref) | mean | std | Z(GLU195–GLY366) | inter-protomer contacts with Z ≥ 2.5 |
|---|---|---|---|---|---|
| within-family, full DI | 83,845 | 0.000000 | 0.01026 | 0.15 | 0 |
| within-family, contact-eligible | 64,265 | −0.000404 | 0.003699 | 0.53 | 1 |
| monomeric contacts (Ch.2 prose) | 3,114 | 0.004908 | 0.011564 | −0.29 | 0 |
| monomeric non-contact background | 61,151 | −0.000675 | 0.002462 | **0.91** | **3** |
| — thesis (plmDCA FN, multi-system) | — | — | — | 2.09 | — |

## Takeaway
Under every reference convention the local mean-field Z for GLU195–GLY366 is
below the thesis 2.09 (score-type difference), but the qualitative conclusion is
reference-independent: **the filament interface does not register above the
monomer signal.** Even under the most favorable (background) reference, only 3
inter-protomer contacts reach Z ≥ 2.5, versus 100+ intra-protomer contacts —
the monomer-dominance finding stands.

Outputs: `results/interface_zscores_ported.csv` (Z per pair under each
reference), `results/zscore_reference_comparison.csv` (this table).
