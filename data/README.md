# Data

Large inputs are not tracked in git. To reproduce the Chapter-4 baseline:

- **`aln/`** — the actin Pfam family MSA (thesis: 31,592 sequences). Actin is
  Pfam **Actin (PF00022)**. Fetch the full alignment from InterPro/Pfam and
  place it here (FASTA or Stockholm). Gap-filter before inference.
- **`pdb/`** — an F-actin filament reference model for contact validation
  (e.g. a cryo-EM filament such as PDB 6BNO / 3J8I, or the model used in the
  thesis). Registration maps MSA columns to protomer residue numbering.

Record the exact accession and any filtering used when you add a file.
