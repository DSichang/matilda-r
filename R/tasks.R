# Downstream tasks: data-first, object-carries-model, results written back.

#' Shared task runner: stage a run dir, seed the model, write inputs, run the script.
#' @keywords internal
.run_task <- function(model, x, label, flags, query, assay, adt_exp, atac_exp, device) {
  n <- .ncells(x)
  lab <- if (is.null(label)) rep(model$label_levels[1], n) else label
  mods <- .as_modalities(x, label = lab, rna = assay, adt = adt_exp, atac = atac_exp)
  rundir <- .stage_rundir(); .seed_rundir_model(rundir, model)
  td <- file.path(rundir, "data")
  pth <- function(nm) file.path(td, paste0(nm, ".h5"))
  .write_h5_matilda(mods$rna, pth("rna"))
  if (!is.null(mods$adt))  .write_h5_matilda(mods$adt,  pth("adt"))
  if (!is.null(mods$atac)) .write_h5_matilda(mods$atac, pth("atac"))
  ctyf <- file.path(td, "cty.csv"); .write_cty_csv(mods$cty, ctyf)
  args <- c("--rna", pth("rna"), "--cty", ctyf,
            if (!is.null(mods$adt))  c("--adt",  pth("adt")),
            if (!is.null(mods$atac)) c("--atac", pth("atac")),
            if (query) c("--query", "True"),
            "--z_dim", model$dims$z_dim, "--hidden_rna", model$dims$hidden_rna,
            "--hidden_adt", model$dims$hidden_adt, "--hidden_atac", model$dims$hidden_atac,
            "--seed", model$hyper$seed, flags)
  script <- if (model$mode == "RNAseq") "main_matilda_rna_task.py" else "main_matilda_task.py"
  .matilda_run(script, as.character(args), rundir, device = device)
  list(rundir = rundir, out = file.path(rundir, "output"),
       sub = if (query) "query" else "reference")
}

#' @keywords internal
.write_back_pred <- function(x, pred, prob) {
  if (methods::is(x, "SummarizedExperiment")) {
    SummarizedExperiment::colData(x)$matilda_pred <- pred
    SummarizedExperiment::colData(x)$matilda_prob <- prob
    return(x)
  }
  if (methods::is(x, "MultiAssayExperiment")) {
    MultiAssayExperiment::colData(x)$matilda_pred <- pred
    MultiAssayExperiment::colData(x)$matilda_prob <- prob
    return(x)
  }
  list(pred = pred, prob = prob)
}

#' @keywords internal
.write_back_latent <- function(x, L) {
  if (methods::is(x, "SingleCellExperiment")) {
    SingleCellExperiment::reducedDim(x, "MATILDA") <- L
    return(x)
  }
  if (methods::is(x, "SummarizedExperiment") || methods::is(x, "MultiAssayExperiment")) {
    S4Vectors::metadata(x)$MATILDA <- L
    return(x)
  }
  list(latent = L)
}

#' @keywords internal
.build_sim <- function(s, mode) {
  if (is.null(s$rna)) return(s)
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = s$rna))
  if (!is.null(s$adt)) {
    SingleCellExperiment::altExp(sce, "ADT") <-
      SummarizedExperiment::SummarizedExperiment(list(counts = s$adt))
  }
  if (!is.null(s$atac)) {
    SingleCellExperiment::altExp(sce, "ATAC") <-
      SummarizedExperiment::SummarizedExperiment(list(counts = s$atac))
  }
  if (!is.null(s$label)) SummarizedExperiment::colData(sce)$label <- s$label
  sce
}

#' Classify cells with a trained Matilda model.
#'
#' @param x SCE/MAE (with a model, or a query) / matrix list.
#' @param reference a trained object/model to use; \code{NULL} = use \code{x}'s own.
#' @param label optional cell-type labels for the ground-truth column of the report.
#' @param assay,adt_exp,atac_exp assay/altExp selectors.
#' @param device "auto"/"cpu"/"cuda".
#' @return \code{x} with \code{colData$matilda_pred}/\code{$matilda_prob}, or a list for matrices.
#' @export
matilda_classify <- function(x, reference = NULL, label = NULL,
                             assay = "counts", adt_exp = "ADT", atac_exp = "ATAC",
                             device = c("auto", "cpu", "cuda")) {
  device <- match.arg(device)
  model <- .resolve_model(x, reference)
  r <- .run_task(model, x, label, c("--classification", "True"),
                 query = !is.null(reference), assay, adt_exp, atac_exp, device)
  on.exit(unlink(r$rundir, recursive = TRUE), add = TRUE)
  f <- file.path(r$out, "classification", model$mode, r$sub, "accuracy_each_cell.txt")
  df <- .parse_classification(f)
  .write_back_pred(x, df$predicted, df$prob)
}

#' Project cells into the Matilda integrated latent space.
#'
#' @inheritParams matilda_classify
#' @return \code{x} with \code{reducedDim "MATILDA"}, or a list for matrices.
#' @export
matilda_reduce <- function(x, reference = NULL, label = NULL,
                           assay = "counts", adt_exp = "ADT", atac_exp = "ATAC",
                           device = c("auto", "cpu", "cuda")) {
  device <- match.arg(device)
  model <- .resolve_model(x, reference)
  r <- .run_task(model, x, label, c("--dim_reduce", "True"),
                 query = !is.null(reference), assay, adt_exp, atac_exp, device)
  on.exit(unlink(r$rundir, recursive = TRUE), add = TRUE)
  L <- .read_latent(file.path(r$out, "dim_reduce", model$mode, r$sub, "latent_space.csv"))
  .write_back_latent(x, L)
}

#' Per-cell-type feature importance (markers) via integrated gradients / saliency.
#'
#' @inheritParams matilda_classify
#' @param method "IntegratedGradient" (default) or "Saliency".
#' @return data.frame(celltype, feature, importance).
#' @export
matilda_markers <- function(x, reference = NULL, label = NULL,
                            method = c("IntegratedGradient", "Saliency"),
                            assay = "counts", adt_exp = "ADT", atac_exp = "ATAC",
                            device = c("auto", "cpu", "cuda")) {
  method <- match.arg(method); device <- match.arg(device)
  model <- .resolve_model(x, reference)
  if (is.null(label)) label <- model$label_col
  if (is.null(label)) {
    stop("matilda_markers needs cell-type labels; pass label= (a colData column or a vector).")
  }
  r <- .run_task(model, x, label, c("--fs", "True", "--fs_method", method),
                 query = !is.null(reference), assay, adt_exp, atac_exp, device)
  on.exit(unlink(r$rundir, recursive = TRUE), add = TRUE)
  .read_markers(file.path(r$out, "marker", model$mode, r$sub))
}

#' Simulate cells for a cell type (or all types) from a trained model.
#'
#' @inheritParams matilda_classify
#' @param celltype cell-type name to simulate; \code{NULL} = all types.
#' @param n number of cells to simulate.
#' @return a SingleCellExperiment of simulated cells.
#' @export
matilda_simulate <- function(x, reference = NULL, celltype = NULL, n = 100L, label = NULL,
                             assay = "counts", adt_exp = "ADT", atac_exp = "ATAC",
                             device = c("auto", "cpu", "cuda")) {
  device <- match.arg(device)
  model <- .resolve_model(x, reference)
  if (is.null(label)) label <- model$label_col
  if (is.null(label)) stop("matilda_simulate needs cell-type labels; pass label=.")
  ct <- if (is.null(celltype)) "-1" else as.character(celltype)
  r <- .run_task(model, x, label,
                 c("--simulation", "True", "--simulation_ct", ct,
                   "--simulation_num", as.character(as.integer(n))),
                 query = !is.null(reference), assay, adt_exp, atac_exp, device)
  on.exit(unlink(r$rundir, recursive = TRUE), add = TRUE)
  sim <- .read_sim(file.path(r$out, "simulation_result", model$mode, r$sub))
  .build_sim(sim$sim, model$mode)
}

#' Run Matilda tasks from file paths (mirrors main_matilda_task.py); writes to outdir.
#'
#' Power-user / validation entry point: drives the unchanged task script on raw
#' .h5/.csv inputs and copies the produced \code{output/} tree to \code{outdir}.
#'
#' @param model a \code{matilda_model}.
#' @param rna,adt,atac,cty input file paths.
#' @param classification,fs,dim_reduce,simulation,query task flags.
#' @param fs_method,simulation_ct,simulation_num task options.
#' @param outdir directory to copy results into.
#' @param device "auto"/"cpu"/"cuda".
#' @return the output directory, invisibly.
#' @export
matilda_task_files <- function(model, rna, adt = NULL, atac = NULL, cty,
                               classification = FALSE, fs = FALSE, dim_reduce = FALSE,
                               simulation = FALSE, query = FALSE,
                               fs_method = "IntegratedGradient", simulation_ct = -1,
                               simulation_num = 100L, outdir = ".",
                               device = c("auto", "cpu", "cuda")) {
  device <- match.arg(device)
  rundir <- .stage_rundir(); on.exit(unlink(rundir, recursive = TRUE), add = TRUE)
  .seed_rundir_model(rundir, model)
  args <- c("--rna", rna, "--cty", cty,
            if (!is.null(adt))  c("--adt",  adt),
            if (!is.null(atac)) c("--atac", atac),
            if (classification) c("--classification", "True"),
            if (fs) c("--fs", "True", "--fs_method", fs_method),
            if (dim_reduce) c("--dim_reduce", "True"),
            if (simulation) c("--simulation", "True", "--simulation_ct", as.character(simulation_ct),
                              "--simulation_num", as.character(as.integer(simulation_num))),
            if (query) c("--query", "True"),
            "--z_dim", model$dims$z_dim, "--hidden_rna", model$dims$hidden_rna,
            "--hidden_adt", model$dims$hidden_adt, "--hidden_atac", model$dims$hidden_atac,
            "--seed", model$hyper$seed)
  script <- if (model$mode == "RNAseq") "main_matilda_rna_task.py" else "main_matilda_task.py"
  .matilda_run(script, as.character(args), rundir, device = device)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  files <- list.files(file.path(rundir, "output"), full.names = TRUE)
  if (length(files)) file.copy(files, outdir, recursive = TRUE)
  invisible(normalizePath(outdir))
}
