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
- [ ] deps install on cmri (running via nohup → `logs/deps_install.log`)
- [ ] HDF5Array (chained → `logs/hdf5array.log`)
- [ ] R CMD INSTALL + pure-R tests pass
- [ ] basilisk env build (first use) + integration tests
- [ ] TEA-seq end-to-end validation (nohup → `logs/validate.log`)
- [ ] roxygen docs / vignette / BiocCheck (Tasks 16–18)

## How to check / resume
```bash
ssh kserver 'tail -20 /media/disk2/Sichang/matilda-r/logs/deps_install.log'
ssh kserver 'tail -20 /media/disk2/Sichang/matilda-r/logs/validate.log'   # once started
```
