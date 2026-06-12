#' Locate the TEA-seq demo dataset (train/test .h5 + label .csv).
#'
#' For v1 this resolves a local copy via \code{options(matilda.demo=)},
#' the \code{MATILDA_DEMO} environment variable, or the development server path.
#' (A \pkg{BiocFileCache} download for public distribution is a planned upgrade.)
#'
#' @param dir directory holding the demo files; \code{NULL} auto-resolves.
#' @return a named character vector of the eight file paths.
#' @export
matilda_example_teaseq <- function(dir = NULL) {
  if (is.null(dir)) {
    dir <- getOption("matilda.demo",
                     Sys.getenv("MATILDA_DEMO",
                                "/media/disk2/Sichang/matilda_test/Matilda/data/TEAseq"))
  }
  files <- c(train_rna = "train_rna.h5", train_adt = "train_adt.h5",
             train_atac = "train_atac.h5", train_cty = "train_cty.csv",
             test_rna = "test_rna.h5", test_adt = "test_adt.h5",
             test_atac = "test_atac.h5", test_cty = "test_cty.csv")
  paths <- stats::setNames(file.path(dir, files), names(files))
  miss <- !file.exists(paths)
  if (any(miss)) {
    stop("TEA-seq demo files not found in '", dir, "': ",
         paste(files[miss], collapse = ", "),
         ". Set options(matilda.demo=) or MATILDA_DEMO.")
  }
  paths
}
