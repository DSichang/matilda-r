# matilda tutorial

**`matilda-tutorial.Rmd`** — the complete Matilda workflow in **R** on TEA-seq. (A parallel,
identically-structured **Python** tutorial — pure Python, `matilda-sc` function API — lives
alongside as a Jupyter notebook; R users use this one, Python users use that one.)

1. **Read your data** into a `SingleCellExperiment`, shown for four formats — native `.h5`
   (`rhdf5`), `.h5ad` (`zellkonverter`), 10x (`Seurat::Read10X`), and a `Seurat` `.rds` — each
   verified to give the same object.
2. **Train** (`matilda_train`).
3. **Classification** of held-out query cells (+ per-cell-type accuracy plot).
4. **Dimension reduction** (+ latent-space UMAP).
5. **Feature selection** (+ marker-importance heatmap).
6. **Simulation** (+ real-vs-synthetic UMAP).
7. **The data types Matilda supports** (`rna_only`/`CITEseq`/`SHAREseq`/`TEAseq`, + accuracy-by-modality bar).
8. **Reproducibility**.

## Requirements
- the `matilda` package (Python is provisioned automatically by basilisk);
- the TEA-seq demo data — point the tutorial at it with
  `options(matilda.demo = "...", matilda.demo_formats = "...")` (the Rmd has
  sensible defaults), or edit the two paths near the top;
- `zellkonverter` (`.h5ad`) and `Seurat` (10x + `.rds`) for the alternative loaders;
  `scater` for the UMAPs; `ggplot2` for the plots.

## Run
```r
rmarkdown::render("matilda-tutorial.Rmd")
```
Takes a few minutes on a CPU laptop, seconds on a GPU.

## Reproducibility note
The R wrapper reproduces the original Matilda Python **exactly on the same
hardware** (verified bit-identical on GPU; see `inst/scripts/parity_check.R`).
Exact figures differ slightly between GPU and CPU (a floating-point property of
PyTorch, not the wrapper), so the accuracy you see on a CPU laptop will be very
close to — but not identical to — the GPU number. The tutorial computes and
prints the accuracy on your machine.
```
