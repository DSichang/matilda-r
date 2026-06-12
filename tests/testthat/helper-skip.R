skip_if_no_matilda_env <- function() {
  testthat::skip_if_not(
    isTRUE(tryCatch(matilda:::.have_env(), error = function(e) FALSE)),
    "basilisk matilda env not available"
  )
}

demo_dir <- function() {
  Sys.getenv("MATILDA_DEMO", "/media/disk2/Sichang/matilda_test/Matilda/data/TEAseq")
}

skip_if_no_demo <- function() {
  testthat::skip_if_not(
    file.exists(file.path(demo_dir(), "train_rna.h5")),
    "TEAseq demo not found"
  )
}
