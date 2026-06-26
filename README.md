# matilda

**Multi-task learning from single-cell multimodal omics (RNA + ADT + ATAC), in R.**

`matilda` is an R interface to [Matilda](https://github.com/PYangLab/Matilda). One multimodal
variational autoencoder + classifier is trained once and reused for five tasks —
**classification, dimension reduction, feature selection, data simulation**, and training-time
**augmentation**. The R package runs Matilda's **unchanged** PyTorch code through `basilisk`
+ `reticulate`, so you never install or manage Python, and results are **bit-identical** to the
original on the same hardware (`inst/scripts/parity_check.R`).

The API is Seurat-style: object in → object out. You pass a `SingleCellExperiment` /
`MultiAssayExperiment` / `Seurat`; the trained model is stored in `metadata()`, and results are
written back (`colData$matilda_pred`, `reducedDim "MATILDA"`).

```r
sce <- matilda_train(sce, label = "cell_type")   # train (model stored in the object)
sce <- matilda_reduce(sce)                        # reducedDim "MATILDA"
sce <- matilda_classify(query, reference = sce)   # colData$matilda_pred / $matilda_prob
mk  <- matilda_markers(sce)                        # per-cell-type feature importance
sim <- matilda_simulate(sce, celltype = "B.Naive", n = 200)
```

---

## How it works (data flow)

```
your object ─▶ convert.R ─▶ io_write.R ─▶ bridge.R ─▶ [ vendored Python in basilisk env ] ─▶ io_read.R ─▶ your object
 (SCE/MAE/         (→ rna/adt/atac      (→ Matilda      (run the unchanged          (parse output      (pred / latent /
  Seurat)           matrices + labels)   .h5 + .csv)     main_matilda_*.py)          CSV/txt files)      markers / sim SCE)
```

The R layer is a thin, faithful driver; all the modeling math lives in the vendored Python.

---

## Repository structure

### R package — the binding (`R/`)

| File | Responsibility |
|------|----------------|
| `matilda-package.R` | package-level roxygen / imports |
| `train.R`   | `matilda_train()`, `matilda_train_files()` — train; build the `matilda_model` |
| `tasks.R`   | `matilda_classify/reduce/markers/simulate()`, `matilda_task_files()`, shared `.run_task()`, result write-backs |
| `model.R`   | the `matilda_model` S3 class: constructor, `print`, accessor, store/resolve |
| `convert.R` | `.as_modalities()` — SCE/MAE/Seurat/list → `rna/adt/atac` matrices + labels + mode |
| `io_write.R`| `.write_h5_matilda()`, `.write_cty_csv()` — R objects → Matilda `.h5`/`.csv` |
| `io_read.R` | parse Python outputs (classification / latent / markers / simulation) |
| `bridge.R`  | `.matilda_run()` — run a vendored script in the env (sets `sys.argv` + cwd + `py_run_file`) |
| `basilisk.R`| the bundled Python environment (torch 2.1.2, captum, scanpy, pandas, …) |
| `rundir.R`  | stage the temporary `../` tree the scripts expect; seed the trained checkpoint |
| `utils.R`   | `.mode_of()`, `.ncells()`, `.h5_features()`, `.pkg_py()` |
| `data.R`    | `matilda_example_sce()`, `matilda_example_teaseq()` |

### Python — the engine (`inst/python/matilda/`, vendored byte-identical to upstream)

| File | Role | Triggered by |
|------|------|--------------|
| `main_matilda_train.py`     | multimodal training | `matilda_train(_files)` (CITEseq/SHAREseq/TEAseq) |
| `main_matilda_task.py`      | multimodal tasks (classify / dim_reduce / fs / simulation) | `matilda_classify/reduce/markers/simulate`, `matilda_task_files` |
| `main_matilda_rna_train.py` | RNA-only training | `matilda_train(_files)` when only RNA |
| `main_matilda_rna_task.py`  | RNA-only tasks | the task verbs when `mode == "rna_only"` |
| `util.py`                   | h5/csv IO, log2 + z-score, Dataset, train/infer helpers | imported by all scripts |
| `learn/model.py`            | VAE + classifier: CITEseq / SHAREseq / TEAseq | imported by the multimodal scripts |
| `learn/model_rna.py`        | VAE + classifier: rna_only | imported by the rna scripts |
| `learn/train.py`            | `train_model()` (training loop) | imported by the train scripts |
| `learn/predict.py`          | `test_model()` (inference / accuracy) | imported by the task scripts |

### Supporting files

| Path | Purpose |
|------|---------|
| `inst/tutorials/matilda-tutorial.Rmd` | **the tutorial** — a Python ⇄ R mirror on TEA-seq (see below) |
| `inst/scripts/parity_check.R` | R-vs-Python bit-parity across all tasks |
| `tests/testthat/` | unit tests (+ `helper-skip.R`, `helper-toy.R`) |
| `man/` | generated documentation |

---

## API reference

Every R verb maps one-to-one to a Python CLI invocation (the R wrapper builds the argv and runs
the script). `device = c("auto","cpu","cuda")` on every call. Supplying `reference=` (object API)
or `query = TRUE` (file API) ⇒ the `--query True` branch.

### Train

| R | Python | Effect |
|---|--------|--------|
| `matilda_train(x, label, assay="counts", adt_exp="ADT", atac_exp="ATAC", batch_size=64, epochs=30, lr=0.02, z_dim=100, hidden_rna=185, hidden_adt=30, hidden_atac=185, augmentation=TRUE, seed=1, device)` | `main_matilda_train.py --rna --adt --atac --cty [--batch_size --epochs --lr --z_dim --hidden_* --seed --augmentation]` | train on an SCE/MAE/Seurat (or matrices); store model in `metadata(x)$matilda`; return `x` |
| `matilda_train_files(rna, adt=NULL, atac=NULL, cty, …same hyperparams…, device)` | same script | train from file paths; return a `matilda_model` |

### Tasks — object API (model carried in the object, results written back)

| R | Python | Returns / writes back |
|---|--------|-----------------------|
| `matilda_classify(x, reference=NULL, label=NULL, assay, adt_exp, atac_exp, device)` | `main_matilda_task.py --classification True [--query True]` | `colData$matilda_pred`, `$matilda_prob` |
| `matilda_reduce(x, reference=NULL, …)` | `--dim_reduce True [--query True]` | `reducedDim(x, "MATILDA")` |
| `matilda_markers(x, reference=NULL, method=c("IntegratedGradient","Saliency"), …)` | `--fs True --fs_method <method> [--query True]` | `data.frame(celltype, feature, importance)` |
| `matilda_simulate(x, reference=NULL, celltype=NULL, n=100, label=NULL, …)` | `--simulation True --simulation_ct <celltype \| -1> --simulation_num <n> [--query True]` | a simulated `SingleCellExperiment` (`celltype=NULL` ⇒ `-1` = all cells) |

### Tasks — file-path API (power-user / validation)

| R | Python | Effect |
|---|--------|--------|
| `matilda_task_files(model, rna, adt=NULL, atac=NULL, cty, classification=FALSE, fs=FALSE, dim_reduce=FALSE, simulation=FALSE, query=FALSE, fs_method="IntegratedGradient", simulation_ct=-1, simulation_num=100, outdir=".", device)` | combinable `--classification/--fs/--dim_reduce/--simulation [--query True]` | runs the task(s) and copies the `output/` tree to `outdir` |

### Model handle & example data

| R | Effect |
|---|--------|
| `matilda_model(object)` | extract the stored `matilda_model` from an SCE/MAE (or pass one through) |
| `print(<matilda_model>)` | summary: mode, cell types, latent dim, per-modality feature counts, hyperparameters |
| `matilda_example_sce(n_cells=60)` | a tiny synthetic multimodal SCE for examples / quick trials |
| `matilda_example_teaseq(dir=NULL)` | resolve the 8 TEA-seq demo file paths (`options(matilda.demo=)` / `MATILDA_DEMO`) |
| `matilda_download_example()` | download + cache the TEA-seq demo dataset (~75 MB) used by the tutorial; returns the local dir |

**Modes** are auto-detected from which modalities are present: `TEAseq` (RNA+ADT+ATAC),
`CITEseq` (RNA+ADT), `SHAREseq` (RNA+ATAC), `rna_only` (RNA). The mode selects both the Python
script and the model architecture.

---

## Tutorial

[`inst/tutorials/matilda-tutorial.Rmd`](inst/tutorials/matilda-tutorial.Rmd) — the complete
workflow in **R** on real TEA-seq. (A separate, identically-structured **Python** tutorial — pure
Python, `matilda-sc` function API — is provided as a Jupyter notebook for Python users.)

1. **Read your data** — load TEA-seq from `.h5`, `.h5ad`, 10x, or a `Seurat` `.rds` (each verified
   to give the same `SingleCellExperiment`).
2. **Train** (`matilda_train`).
3. **Classification** of held-out query cells (+ per-cell-type accuracy plot).
4. **Dimension reduction** (+ latent-space UMAP).
5. **Feature selection** (+ marker heatmap).
6. **Simulation** (+ real-vs-synthetic UMAP).
7. **The data types Matilda supports** — `rna_only` / `CITEseq` / `SHAREseq` / `TEAseq` on the
   *same* TEA-seq cells (a modality ablation; adding the ADT panel helps most).

```r
rmarkdown::render("inst/tutorials/matilda-tutorial.Rmd")
```

---

## Reproducibility

On the **same hardware** the R wrapper is bit-identical to the original Python (verified across
all tasks by `inst/scripts/parity_check.R`, `max|Δ| = 0`): TEA-seq query accuracy **0.8092** on
GPU, **0.8148** on CPU. The only difference between machines is PyTorch GPU↔CPU floating point —
not the wrapper. UMAP layouts vary run-to-run by construction; cluster structure is stable.
