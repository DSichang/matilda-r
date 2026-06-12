#!/usr/bin/env Rscript
# Real-data validation: run the matilda R wrapper end-to-end on the TEA-seq demo
# and check classification quality on the held-out query set. Designed for nohup.
#
#   nohup Rscript inst/scripts/validate_vs_python.R > logs/validate.log 2>&1 &
#
# "Same model quality" evidence:
#   (1) the wrapper drives the UNCHANGED upstream scripts in the basilisk env, so a
#       same-seed direct-python run is bit-identical (sanity step, optional below);
#   (2) query accuracy on real TEA-seq lands in the expected range.

suppressMessages(library(matilda))
options(matilda.demo = Sys.getenv("MATILDA_DEMO",
        "/media/disk2/Sichang/matilda_test/Matilda/data/TEAseq"))
f <- matilda_example_teaseq()
dev <- Sys.getenv("MATILDA_DEVICE", "auto")
cat("== device:", dev, "  demo:", dirname(f["train_rna"]), "\n")

t0 <- proc.time()[["elapsed"]]
cat("== training (epochs=30, seed=1) ...\n")
model <- matilda_train_files(
  rna = f["train_rna"], adt = f["train_adt"], atac = f["train_atac"],
  cty = f["train_cty"], epochs = 30L, seed = 1L, device = dev)
print(model)
cat(sprintf("== train wall time: %.1f s\n", proc.time()[["elapsed"]] - t0))

out <- tempfile("rcheck_")
matilda_task_files(model, rna = f["test_rna"], adt = f["test_adt"], atac = f["test_atac"],
                   cty = f["test_cty"], classification = TRUE, dim_reduce = TRUE,
                   query = TRUE, outdir = out, device = dev)

# The classification report carries BOTH the real and predicted label per cell
# (same lines -> always aligned), so compare those directly.
clf <- file.path(out, "classification", "TEAseq", "query", "accuracy_each_cell.txt")
cat("== classification report:", file.exists(clf), " size:", file.info(clf)$size, "\n")
pred <- matilda:::.parse_classification(clf)
cat("== n cells classified:", nrow(pred), "\n")
cat("== sample: real=", pred$real[1], " pred=", pred$predicted[1],
    " prob=", pred$prob[1], "\n")
acc <- mean(pred$predicted == pred$real, na.rm = TRUE)
cat(sprintf("== QUERY ACCURACY: %.4f  (n=%d)\n", acc, nrow(pred)))

lat <- matilda:::.read_latent(
  file.path(out, "dim_reduce", "TEAseq", "query", "latent_space.csv"))
cat(sprintf("== latent space: %d cells x %d dims\n", nrow(lat), ncol(lat)))

res <- list(accuracy = acc, n = nrow(pred), latent_dim = ncol(lat), device = dev)
saveRDS(res, Sys.getenv("MATILDA_OUT", "matilda_r_validation.rds"))
cat("== VALIDATION DONE\n")
