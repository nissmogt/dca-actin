#!/usr/bin/env python3
"""
run_actin_dca.py — run mean-field DCA on the actin Pfam MSA using the shared
engine, and rank pairs by APC-corrected DI (diapc).

Thread-pinning MUST happen before numpy/pydca import (each worker one thread, or
the numba weight loop loses ~10x). This script is single-process; for sweeps use
a batch driver that sets these per worker.

The shared engine writes DI_{sysid}_n{L}.txt into the directory that holds the
MSA (it ignores dir_dca and uses model_length only as a filename tag). Point the
MSA at data/aln/ and the DI table lands beside it; copy into results/ after.

Usage:
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
        python src/run_actin_dca.py data/aln/actin_PF00022.fasta actin [L]

Prints NEFF and the DI-file path.
"""
import os
# thread-pin before any numpy/numba import
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMBA_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import sys
sys.path.insert(0, os.path.dirname(__file__))
from dca_engine import inference  # noqa: E402


def main():
    if len(sys.argv) < 3:
        sys.exit("usage: run_actin_dca.py <msa.fasta> <sysid> [model_length_tag]")
    msa, sysid = sys.argv[1], sys.argv[2]
    L = sys.argv[3] if len(sys.argv) > 3 else "0"
    neff = inference(sysid, "", msa, L)
    outdir = os.path.dirname(msa) or "."
    print(f"{sysid}: NEFF={neff}  -> {outdir}/DI_{sysid}_n{L}.txt")


if __name__ == "__main__":
    main()
