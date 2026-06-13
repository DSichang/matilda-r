# matilda tutorials

Runnable, per-task tutorials that reproduce every Matilda function from R and
validate that the result matches the original Matilda Python implementation.

| # | File | Task |
|---|------|------|
| 0 | `00-data-preprocessing.Rmd`  | Get your data into an accepted format (h5ad / h5mu / 10x h5 / Seurat / matrices) |
| 1 | `01-training.Rmd`            | Train the multimodal VAE + classifier |
| 2 | `02-classification.Rmd`     | Predict cell types for query cells |
| 3 | `03-dimension-reduction.Rmd`| Integrated multimodal latent space |
| 4 | `04-feature-selection.Rmd`  | Per-cell-type marker importance |
| 5 | `05-simulation.Rmd`         | Generate synthetic cells |

`_setup.R` holds shared helpers (`load_teaseq_sce()`); each tutorial sources it.

## Requirements
- The `matilda` package installed.
- The TEA-seq demo data — `matilda_example_teaseq()` resolves it; otherwise set
  `options(matilda.demo = "/path/to/TEAseq")` or the `MATILDA_DEMO` env var.
- The first call builds the bundled Python environment (one-off); a GPU is used
  if present, otherwise CPU.

## Run
From this directory, knit any tutorial, e.g.
```r
rmarkdown::render("01-training.Rmd")
```
Run `01-training.Rmd` first — it caches the trained model to
`matilda_teaseq_ref.rds`, which tutorials 2–5 reuse (they retrain automatically
if it is missing).

## Reproduce the parity check
`inst/scripts/parity_check.R` runs matilda (R) and the original Matilda Python on
the **same checkpoint** and reports the per-task difference. Confirmed on TEA-seq:

| Task | matildaR vs Matilda Python |
|------|----------------------------|
| Classification | identical predictions (accuracy 0.8092 = 0.8092) |
| Dimension reduction | bit-identical (max\|Δ\| = 0) |
| Feature selection | bit-identical (max\|Δ\| = 0, all cell types) |
| Simulation | bit-identical (max\|Δ\| = 0) |
