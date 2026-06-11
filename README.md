# ER-57bpDNA-cryoEM-3DVA

Analysis code and 3D variability analysis (3DVA) latent coordinates for the
single-particle cryo-EM study of an estrogen receptor alpha (ERα) DNA-binding
domain dimer bound to a 57-bp duplex containing the consensus vitellogenin A2
estrogen response element (ERE) with native flanking sequence.

This repository contains the code to reproduce the conformational-landscape
analysis showing that the flanking DNA samples a **single continuous bending
landscape** around an invariant protein–ERE interface, with the C-, S-, and
I-shape states being representative geometries rather than discrete conformers.

## Contents

- `make_fig4B.py` — reproduces the latent-landscape panel (main-text Fig. 4B):
  black scatter of all consensus particles, per-class translucent density
  fields for the C/S/I states, and class-centroid labels.
- `requirements.txt` — Python dependencies.

## Data

The per-particle 3DVA latent coordinates are deposited on Zenodo
(DOI: `10.5281/zenodo.20649656`), not in this repository, due to size:

| File | Contents |
|------|----------|
| `3DVA_J2538_8A_latent_with_class.csv` | Primary run (8 Å filter). Columns: `uid, PC1, PC2, PC3, class` (class = C/S/I or empty). Reproduces Fig. 4B and Supplementary Fig. S4A. |
| `3DVA_J2541_10A_latent.csv` | Robustness run, 10 Å filter. Columns: `uid, PC1, PC2, PC3`. |
| `3DVA_J2581_loosermask_latent.csv` | Robustness run, more permissive mask. |
| `3DVA_J2587_6A_diffmask_latent.csv` | Robustness run, 6 Å filter, different mask. |

The three robustness runs share the same particle `uid` set as the primary
run; class assignments can be joined from the primary file by `uid`.

## Reproducing Fig. 4B

```bash
pip install -r requirements.txt
# download 3DVA_J2538_8A_latent_with_class.csv from Zenodo into this folder
python make_fig4B.py
# -> Fig4B_latent_landscape.pdf / .png / .svg
```

To run on a different file: `python make_fig4B.py path/to/latent.csv`.

## 3DVA parameters

3DVA (cryoSPARC v4.7.1) was run on the full consensus particle set
(n = 1,176,651), aligned by homogeneous refinement followed by local
refinement on the invariant protein–ERE core, prior to heterogeneous
classification. Settings: 3 variability components, 8 Å filter resolution,
C1 symmetry, soft mask enclosing the whole complex including the flanking DNA.
Equivalent results were obtained at 6 and 10 Å and with alternative masks
(the robustness files above).

## Citation

If you use this code or data, please cite:

> [Authors]. [Title]. [Journal] (Year). DOI: [paper DOI]

and the Zenodo archive: DOI `10.5281/zenodo.20649656`.

## License

MIT (see `LICENSE`).
