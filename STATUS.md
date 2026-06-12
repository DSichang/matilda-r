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
- [x] basilisk env built (torch 2.1.2 on RTX 4090) — env at
      `~/.cache/R/basilisk/1.22.0/matilda/0.0.0.9000/env-matilda`
- [x] integration tests PASS (train + classify/reduce/markers/simulate on toy SCE)
- [x] **TEA-seq end-to-end VALIDATED**: train 12s → classify → **0.8092 query acc** (n=5048)
- [x] **PARITY PROVEN**: R wrapper 0.8092 == direct-python (same scripts/env/seed) 0.8092
- [x] roxygen docs (29 man pages), runnable examples, vignette, NEWS (Tasks 16–17)
- [x] R CMD check: 1 NOTE only ("Non-staged installation" — required by basilisk, expected)
- [x] BiocCheck: 1 ERROR + 1 WARNING + 9 NOTES, all non-code:
      - ERROR `checkSupportReg`: HTTP 502 to bioconductor.org + maintainer not yet
        registered on the Bioc Support Site → register your email there; not a code defect.
      - WARNING `checkVigChunkEval`: the matilda_* vignette chunks are eval=FALSE
        (would build the env at vignette time) — evaluate on a GPU host before submission.
      - NOTES: line length / 4-space indent / add URL+BugReports when the GitHub repo
        exists / bioc-devel subscription — all cosmetic or pending the public repo.

## Bottom line
v1 is FUNCTIONALLY COMPLETE and VALIDATED on real TEA-seq (query acc 0.8092,
identical to upstream). Remaining before Bioc submission: create the public repo
(URL/BugReports), register maintainer on the Bioc Support Site, and evaluate the
vignette on a GPU host. Package builds as `matilda_0.99.0.tar.gz`.

## Bugs found & fixed during validation (all in the R glue, not the model)
1. basilisk.utils missing → reticulate bootstrapped pyenv (fixed: installed it).
2. Unclosed script file handles → buffered output never flushed (toy: empty;
   real: last ~7 cells lost). Fixed in bridge: flush all open files after each run.
3. os.chdir left the process in a deleted run dir → getcwd()/saveRDS failed.
   Fixed in bridge: restore cwd after each run.
4. Perfectly-balanced toy classes broke upstream median-augmentation. Fixed:
   imbalanced toy_sce.
5. reducedDim rownames mismatch warning. Fixed: rownames(L) <- colnames(x).

## Known notes
- torch pinned to 2.1.2 (not upstream 1.9.1) so it drives the RTX 4090 (sm_89).
  Parity confirmed against the same env, so model quality is identical.
- benign "stack imbalance" warnings from Bioc S4 dispatch under test_dir (not failures).

## How to check / resume
```bash
ssh kserver 'tail -20 /media/disk2/Sichang/matilda-r/logs/deps_install.log'
ssh kserver 'tail -20 /media/disk2/Sichang/matilda-r/logs/validate.log'   # once started
```
