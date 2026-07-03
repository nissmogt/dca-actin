# dca-actin

Mean-field **Direct Coupling Analysis (DCA)** applied to the actin filament —
thesis Chapter 4 (Mehrabiani 2022, Rice University / Onuchic lab, CTBP).

The open question: coevolution recovers the actin **monomer fold** almost
perfectly, but the **filament-interface** contacts stay buried under the
dominant monomer signal. This repo migrates the Chapter 4 mean-field DCA code
off legacy MATLAB/HPC infrastructure and onto the local `pydca` engine shared
with `sequence_degradation`, then attacks monomer-signal dominance.

## Baseline to reproduce
Mean-field DCA on the actin Pfam family (**31,592 sequences**) gave:
- near-perfect monomer-fold recovery;
- essentially **no** filament-interface contacts;
- a single standout inter-protomer pair **GLU195–GLY366** (16.5 Å in the
  filament model, **Z = 2.09**);
- two residue-range clusters (100–150, 250–300) hypothesized to relate to
  nucleation-phase contacts.

See [`docs/BASELINE.md`](docs/BASELINE.md) for the reproduction target and
method conventions.

## Layout
```
src/            shared DCA engine (pydca mean-field, drop-in) + PDB registration
data/aln/       actin Pfam MSA(s)          (not tracked; see data/README.md)
data/pdb/       F-actin reference model     (not tracked)
results/        DI tables, contact maps, figures
notebooks/      exploratory analysis
docs/BASELINE.md   Chapter-4 baseline + method conventions
```

## Method conventions (shared across the DCA program)
- Rank pairs by **APC-corrected DI** (`diapc`).
- Sequence reweighting identity threshold **0.8** (theta 0.2); pseudocount **0.5**.
- Contacts: heavy-atom minimum distance **≤ 8 Å**, sequence separation **|i−j| ≥ 5**.
  (Filament interface contacts are *inter*-protomer, so this |i−j| rule is
  applied within a protomer; cross-protomer pairs are treated separately.)
- Z-score significance ladder: `[12, 10, 9, 8, 5.6, 4.5, 4, 3.5, 2.5, 1]`.
- Register MSA columns to PDB residues by global alignment; low identity
  (< 0.20) flags unreliable ground truth, not a DCA failure.

## Status
Scaffold. Engine and PDB-registration modules migrated from
`sequence_degradation`; actin data assembly and baseline reproduction are the
first work items (see issues / `docs/BASELINE.md`).
