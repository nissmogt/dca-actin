# SASA, hydrophobicity, and the actin filament interface

Testing the hypothesis (KMM): *the filament interface forms through hydrophobic /
surface interactions, and coevolution encodes these as higher-order, long-range
signal.* Method: (i) measure the interface by solvent-accessible surface area
(SASA) burial on the assembled filament (8D13, middle protomer B); (ii) bring
SASA into the DCA Hamiltonian and re-derive the couplings; (iii) compare against
the mechanosensitive-crack model of Zsolnay, Kovar & Voth 2024
(Biophys. J., DOI 10.1016/j.bpj.2024.06.014).

## 1. The interface is patchy and polar, not hydrophobic
Middle protomer B, Shrake–Rupley SASA (monomer vs full ABC filament):
- 36 / 370 residues bury >5 Å² on assembly; total buried ≈1510 Å² per protomer.
- Interface burial clusters in windows 100–150, 265–295, 320–355 — the same
  ranges as the DCA coevolution clusters and the thesis Ch.4 clusters.
- Interface residues are **not** more hydrophobic than the rest of the surface
  (Kyte–Doolittle median −0.70 vs −0.70; Mann–Whitney interface>rest p=0.31;
  only 42% have KD>0). Actin's filament interface is built from polar/charged and
  shape-complementary contacts. **A hydrophobicity-weighted prior would be
  mis-specified; burial (ΔSASA) is the right structural variable.**

## 2. Plain DCA does not rank the interface
Rigorous test (FN-APC ranking filament-specific inter-protomer contacts,
338 positives / 63,907 contact-eligible pairs):
- **AUC = 0.416; zero filament-specific contacts in the top 500.**
- Filament-specific contacts are much longer-range in sequence than intra
  contacts (median |i−j| 104 vs 36, MWU p≈2e−46) — consistent with the
  "higher-order / long-range" idea about *where* the interface lives, but the
  plain model does not surface it.
(An earlier Spearman ρ=0.156, p=0.004 was a within-contact-subset correlation
only; it does not translate into above-chance ranking. Corrected here.)

## 3. The D-loop↔W-loop crack IS coevolutionarily encoded, weakly
Zsolnay et al. 2024 (all-atom MD, starting structure PDB 6DJO, 500 pN tension)
identify the strain-sensitive contact as the **D-loop (residues 40–50, hydrophobic
V43/M44/V45/M47) to W-loop (Y169)** connection between longitudinal subunits
i and i−2, which "buries a large fraction of the surface area between subunits."
Our independent measurements agree and add coevolution:
- **W-loop Y169 buries 150 Å²** — the single largest burial of any residue in
  the protomer, matching the paper's identification of Y169 as the flip residue
  that initiates the crack.
- The top crack-touching inter-protomer DCA pairs are exactly the D-loop×W-loop
  contact: 42–165, 44–351, 45–351, 41–172, 43–171, 40–166, and M44–Y169 directly
  (5.4 Å, positive DI). So coevolution **does** encode the specific crack contact.
- But it is weak: best crack inter-contact ranks #6278 / 64,265 (90th percentile),
  DI ≈ 0.002 vs top monomer DI 0.21. Encoded but far below the fold signal —
  consistent with a strain-activated cryptic site rather than a strong constitutive
  contact. (Fig. crack_dloop_wloop_dca.png)

## 4. SASA-conditioned re-fit: a negative methods result
Two formulations (both run; the contrast is the result):
- **(A) field-prior re-fit** — replace pydca's single scalar pseudocount θ with a
  per-column θ_i tied to monomer surface exposure (buried→θ 0.7, exposed→θ 0.3),
  re-invert the covariance for new couplings, re-score (FN-APC). Conditions on
  *monomer* exposure (knowable without the filament), not assembly burial, so it
  is non-circular.
- **(B) re-ranking baseline** — keep plain couplings, multiply FN by a burial
  weight w(ΔSASA). Circular by construction.

Result (sasa_dca_comparison.csv / .png):

| method | AUC | top-500 filament | max FN-APC |
|---|---|---|---|
| plain FN-APC | 0.416 | 0 | 3.6 |
| re-rank (B) | 0.416 | 0 | 3.6 |
| re-fit θ[.45,.55] | 0.463 | 0 | 501 |
| re-fit θ[.4,.6] | 0.530 | 1 | 482 |
| re-fit θ[.3,.7] | 0.454 | 9 | 1706 |

**Conditioning mean-field DCA on monomer burial does not recover the interface.**
The per-column θ is numerically unstable: even a mild [0.45,0.55] range inflates
max FN-APC from 3.6 to ~500, and the top ranks become dominated by a few blown-up
columns (e.g. col 57 in 330 of the top-500 pairs). Any apparent AUC>0.5 (med θ)
and the "9 top-500 hits" at aggressive θ are byproducts of this inflation, not
interface signal. The re-ranking baseline (B) does nothing because the interface
residues lack plain-FN to boost.

**Interpretation.** The mean-field engine's uniform-pseudocount assumption is
load-bearing — the covariance inversion J=−C⁻¹ is ill-conditioned under per-column
regularization variation. A genuine SASA-conditioned model would need a
regularization that preserves column-scale conditioning (e.g. plmDCA with per-site
Gaussian priors on the fields, or a shrinkage that renormalizes column norms after
re-weighting), not a per-column pseudocount on the mean-field marginals.

## Files
- `results/actin_sasa_burial.csv`, `results/sasa_hydrophobicity.png` — burial + KD
- `results/crack_interface_pairs.csv`, `results/crack_dloop_wloop_dca.png` — crack
- `results/sasa_dca_comparison.csv`, `results/sasa_dca_comparison.png` — A vs B
- `src/sasa_dca.py`, `src/run_fnapc.py` — engine + runner
