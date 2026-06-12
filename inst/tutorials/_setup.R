# Shared helpers for the matilda tutorials.
# source("_setup.R") at the top of each tutorial.
suppressMessages({
  library(matilda)
  library(SingleCellExperiment)
})

#' Read a Matilda .h5 (matrix/data, matrix/features, matrix/barcodes) into a
#' genes x cells matrix, robust to the on-disk orientation.
read_matilda_h5 <- function(path) {
  feats <- as.character(rhdf5::h5read(path, "matrix/features"))
  cells <- as.character(rhdf5::h5read(path, "matrix/barcodes"))
  d <- rhdf5::h5read(path, "matrix/data")
  if (nrow(d) != length(feats)) d <- t(d)      # ensure features x cells
  dimnames(d) <- list(feats, cells)
  d
}

#' Build a SingleCellExperiment for a TEA-seq split ("train" or "test"):
#' RNA in the main assay, ADT and ATAC as altExps, labels in colData$cell_type.
load_teaseq_sce <- function(split = c("train", "test"),
                            files = matilda_example_teaseq()) {
  split <- match.arg(split)
  g <- function(mod) files[[paste0(split, "_", mod)]]
  sce <- SingleCellExperiment(assays = list(counts = read_matilda_h5(g("rna"))))
  altExp(sce, "ADT")  <- SummarizedExperiment(list(counts = read_matilda_h5(g("adt"))))
  altExp(sce, "ATAC") <- SummarizedExperiment(list(counts = read_matilda_h5(g("atac"))))
  sce$cell_type <- as.character(utils::read.csv(g("cty"), header = FALSE)[[2]][-1])
  sce
}
