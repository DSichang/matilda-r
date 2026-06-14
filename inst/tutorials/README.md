# matilda tutorial

**`matilda-tutorial.Rmd`** — a single, end-to-end walkthrough:

1. **Load your data** into a `SingleCellExperiment`, shown explicitly for every
   common format — native Matilda `.h5` (via `rhdf5`), a `Seurat` object, `.h5ad`
   (AnnData/scanpy, via `zellkonverter`), or plain matrices. The loaders are
   verified to give the *same* data, so downstream results don't depend on format.
2. **Run the whole workflow** — `matilda_train` → `matilda_reduce` →
   `matilda_classify` → `matilda_markers` → `matilda_simulate`, plus query
   classification.

## Requirements
- the `matilda` package (Python is provisioned automatically by basilisk);
- the TEA-seq demo data — point the tutorial at it with
  `options(matilda.demo = "...", matilda.demo_formats = "...")` (the Rmd has
  sensible defaults), or edit the two paths near the top;
- `Seurat` and `zellkonverter` for the Seurat / `.h5ad` loading options.

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
