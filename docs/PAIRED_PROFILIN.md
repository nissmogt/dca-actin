# Actin–profilin paired-MSA DCA

First heterodimer partner for the ABP interface program. Tests whether
inter-protein coevolution (mean-field DCA on a species-paired actin+profilin
alignment) predicts the known actin–profilin interface.

## Inputs
- **Actin**: Pfam PF00022, 410 match-state columns (same as the monomer baseline).
- **Profilin**: Pfam PF00235 (9,504 proteins; 6,699 in the full alignment),
  127 match-state columns.
- **Validation structure**: PDB **2PAV** (1.8 Å X-ray, ternary profilin–actin +
  VASP poly-Pro). Actin = chain A, profilin-1 = chain P.

## Pairing
1,424 UniProt species mnemonics are shared between the two families. Actin
carries 3–7 paralogs (isoforms) per species and profilin 1–6, so all-vs-all
pairing would give ~1.9M mostly-mispaired combinations. We take **one
least-gapped representative per species per family** (clade-level "9XXXX" codes
dropped) and concatenate → **1,312 paired rows, length 537 (410 + 127)**.
Registration of each block to 2PAV: actin identity 0.475 (358 cols), profilin
0.452 (124 cols) — both well above the 0.20 unreliable-ground-truth flag. Ground
truth: **286 actin–profilin heavy-atom contacts ≤8 Å**, 261 mappable to columns.

## Result — inter-block signal does not recover the interface
Mean-field DCA (uniform pseudocount 0.5, seqid 0.8), APC-corrected DI. The
inter-block (actin-column × profilin-column) submatrix is 410×127 = 52,070 pairs.

| top-N inter-block | PPV | true contacts |
|---|---|---|
| 5 | 0.00 | 0 |
| 50 | 0.00 | 0 |
| 100 | 0.00 | 0 |
| 500 | 0.002 | 1 |

- **AUC = 0.441** (below chance); best true contact ranks **#462 / 52,070**;
  median true-contact rank 30,561.
- Top inter-block DCA scores (max diapc 0.037) are an order of magnitude below
  the intra-actin fold signal (max 0.227) and do not land on the profilin-binding
  surface (Fig. paired_interface.png).

## Cause — depth, not alignment
**Neff = 397, Neff/L = 0.74** — about 15× shallower per position than the monomer
run (Neff/L = 11.13). Mean-field DCA needs roughly Neff/L ≳ 1–2 for reliable
*inter*-protein contacts; 0.74 is well below that floor. The interface columns
themselves are well-occupied (median 100% non-gap), so this is a sequence-
**diversity** limit created by the one-representative-per-species pairing, not a
gap or registration artifact.

## Interpretation and next step
Consistent with the thesis Ch.4 finding that coevolutionary interface signal for
actin assemblies is weak, now extended to a heterodimeric partner: at this depth,
mean-field DCA does not resolve the actin–profilin interface. The clean 1:1
pairing traded depth for pairing accuracy. To raise depth without inflating
mispairs, the next attempt should **augment pairs by allowing multiple
actin/profilin paralogs per species matched by best reciprocal sequence identity**
(potentially 2–3× the rows), and/or expand the taxonomic sampling. A plmDCA
engine (higher information efficiency at low Neff) would also help; the mean-field
engine is the current constraint.

## Files
- `data/aln/actin_profilin_paired.fasta` — 1,312 paired sequences
- `results/paired_msa_meta.json`, `results/paired_register.json`
- `results/paired_interface_pairs.csv` — inter-block pairs ranked, residue-annotated
- `results/paired_ppv.csv`, `results/paired_interface.png`
- `src/run_paired_dca.py`

## Control: paralog augmentation does not rescue depth
To test whether the clean 1:1 pairing was throwing away usable diversity, we
re-ran with **all profilin paralogs per species** (kept clade-level codes too),
paired to each species' representative actin: 6,447 rows (5× more) at the same
length 537. Reweighting (seqid 0.8) collapsed the added paralogs almost entirely
— **Neff rose only 397 → 481 (Neff/L 0.74 → 0.90)** — and the interface still did
not emerge (top-100 PPV = 0, AUC = 0.460, best true contact rank #467). Within-
species profilin paralogs are too similar to their orthologs to add *effective*
diversity. This confirms the limit is genuine sequence diversity (deep taxonomic
sampling), not the pairing rule. The real levers are a **plmDCA engine** (more
information-efficient at low Neff) or a **partner with intrinsically more
diversity and a larger interface** (Arp2/3, formin FH2). Files:
`data/aln/actin_profilin_paired_aug.fasta`, `src/run_paired_dca_aug.py`.
