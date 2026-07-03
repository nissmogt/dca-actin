# Chapter 4 baseline — actin filament DCA

Target to reproduce before attempting any new methodology.

## Input
- **MSA:** actin Pfam family, **PF00022 (Actin)**. Thesis used **31,592
  sequences**. Gap-filter, then run mean-field DCA.
- **Reference structure:** an F-actin filament model providing *inter-protomer*
  heavy-atom distances (the thesis reports a filament-model distance for the
  standout pair). Record the exact accession used.

## Engine
`src/dca_engine.py` — pydca mean-field DCA, byte-compatible drop-in for the
thesis MATLAB `dca_h_J_Full_v5.m`. Writes `DI_{sysid}_n{L}.txt`
(`i,j,di,diapc,mi`, 1-indexed) and returns `NEFF = ceil(Meff)`.
Reweighting seqid = 0.8 (theta 0.2), pseudocount 0.5.

## Expected result (the baseline)
1. **Monomer fold recovered** — top `diapc` pairs reproduce intra-protomer
   contacts (validate against a single actin monomer, e.g. an F-actin protomer
   or G-actin such as PDB 1ATN chain A).
2. **Filament interface essentially absent** — inter-protomer contacts do *not*
   surface among high-ranked pairs.
3. **Single standout inter-protomer pair:** **GLU195–GLY366**, 16.5 Å in the
   filament model, **Z = 2.09**.
4. **Two residue-range clusters:** **100–150** and **250–300**, hypothesized to
   relate to nucleation-phase contacts.

## Z-score
Significance ladder (shared program convention):
`[12, 10, 9, 8, 5.6, 4.5, 4, 3.5, 2.5, 1]`.
Z is computed against the monomeric-contact reference distribution
(as in Chapter 2 / `dca-interface`).

## Open problem (why this repo exists)
Monomer coevolution signal dominates and masks the filament-interface signal.
Reaching the interface likely needs explicit monomer-signal subtraction,
phylogenetic deconvolution, or a coarse-grained model — see README and the
thesis Ch. 4 perspectives.
