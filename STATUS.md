# matilda R package — build status

Autonomous build log (updated as work proceeds). Canonical copy on cmri:
`/media/disk2/Sichang/matilda-r/matilda/`. Dev location: local
`/Users/sdin3712/Documents/matilda-r/matilda/` → rsynced to cmri.

## Design / plan
- Spec:  `project_matilda/docs/superpowers/specs/2026-06-10-matilda-r-package-v1-design.md`
- Plan:  `project_matilda/docs/superpowers/plans/2026-06-11-matilda-r-package-v1.md`
- API: Seurat-style — object in / object out, model stored in `metadata(obj)$matilda`,
  results written back (`colData$matilda_pred`, `reducedDim "MATILDA"`), pipeable.

## Progress
- [x] Vendored upstream python verbatim → `inst/python/matilda/` (+ empty `__init__.py`)
- [x] DESCRIPTION / NAMESPACE / package doc
- [x] basilisk env (modern torch 2.1.2 for RTX 4090; CPU fallback)
- [x] bridge `.matilda_run` (argv + cwd + py_run_file), rundir staging
- [x] IO: `.write_h5_matilda` (HDF5Array, matches data_to_h5.R), `.write_cty_csv`
- [x] parsers: classification / latent / markers / sim
- [x] convert `.as_modalities` (SCE/MAE/Seurat/matrix)
- [x] model S3 + `matilda_model()` accessor + store/resolve
- [x] verbs: `matilda_train` / `matilda_classify` / `matilda_reduce` /
      `matilda_markers` / `matilda_simulate` + `*_files` path APIs
- [x] tests: pure-R (io/convert/model/io_read) + integration (bridge/train/tasks)
- [x] real-data validation harness `inst/scripts/validate_vs_python.R`
- [x] deps install on cmri (DONE — Bioc stack + devtools/roxygen2/BiocCheck)
- [x] HDF5Array 1.38.0 (DONE)
- [x] basilisk.utils 1.22.0 — was MISSING, installed (it's what lets basilisk build
      its own conda; without it reticulate 1.46 wrongly bootstrapped pyenv)
- [x] R CMD INSTALL + pure-R tests pass (io_write/convert/model/io_read green)
- [~] basilisk env build (IN PROGRESS via nohup → `logs/run_all.log`, phase [1]):
      downloads torch 2.1.2 + scanpy; then verifies torch/cuda + h5 transpose
- [ ] integration tests (queued, phase [2])
- [ ] TEA-seq end-to-end train→classify→accuracy (queued, phase [3])
- [ ] roxygen docs / vignette / BiocCheck (Tasks 16–18)

## Known notes
- torch pinned to 2.1.2 (not upstream 1.9.1) so it drives the RTX 4090 (sm_89);
  "same model quality" = statistical parity, validated in-env.
- benign "stack imbalance" warnings from Bioc S4 dispatch under test_dir (not failures).

## How to check / resume
```bash
ssh kserver 'tail -20 /media/disk2/Sichang/matilda-r/logs/deps_install.log'
ssh kserver 'tail -20 /media/disk2/Sichang/matilda-r/logs/validate.log'   # once started
```
