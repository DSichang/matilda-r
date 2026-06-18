# matilda tutorial

**`matilda-tutorial.Rmd`** — a **Python ⇄ R mirror** on TEA-seq: for every step it shows
the original Matilda **Python CLI** command and the **R** call that runs the same code, then
the output and a plot.

1. **Read the data** into a `SingleCellExperiment`, shown explicitly (native `.h5` via
   `rhdf5`, plus a `Seurat`/`.h5ad` cross-check verified element-identical).
2. **Train** (`matilda_train` ⇄ `main_matilda_train.py`).
3. **Multi-task on the training data** — simulation, dimension reduction, feature selection,
   each with a plot mirroring the official `qc/visualize_simulated_data.Rmd` (real-vs-synthetic
   UMAP) and `qc/visulize_latent_space.Rmd` (latent UMAP), plus a marker-importance heatmap.
4. **Multi-task on the query data** — classification (+ per-cell-type accuracy plot),
   dimension reduction, feature selection (the same verbs with `reference =`, i.e. `--query True`).

## Requirements
- the `matilda` package (Python is provisioned automatically by basilisk);
- the TEA-seq demo data — point the tutorial at it with
  `options(matilda.demo = "...", matilda.demo_formats = "...")` (the Rmd has
  sensible defaults), or edit the two paths near the top;
- `Seurat` + `zellkonverter` for the alternative loaders; `scater` for the UMAPs;
  `ggplot2` for the plots.

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
