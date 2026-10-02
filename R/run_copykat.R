## Run CopyKAT on one exported panel dataset, in both reference modes.
##   reffree - no normal cells declared (CopyKAT picks its own baseline)
##   refimm  - immune cells declared as the normal reference
## the workspace R library is not durable across session restarts; reinstall if needed
if (!requireNamespace("copykat", quietly = TRUE)) {
  tb <- file.path(if (exists("ROOT")) ROOT else ".", "copykat.tar.gz")
  if (file.exists(tb)) install.packages(tb, repos = NULL, type = "source")
}
suppressPackageStartupMessages({ library(Matrix); library(copykat); library(jsonlite) })

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) && exists("TARGET")) args <- TARGET
name <- args[1]
max_cells <- if (length(args) > 1) as.integer(args[2]) else 1200L
seed <- 0L

ROOT <- if (exists("ROOT")) ROOT else normalizePath(".")
d <- file.path(ROOT, "rmat", name)
m <- readMM(file.path(d, "counts.mtx"))
genes <- readLines(file.path(d, "genes.txt"))
cells <- readLines(file.path(d, "cells.txt"))
meta <- read.csv(file.path(d, "meta.csv"), stringsAsFactors = FALSE)
rownames(m) <- genes; colnames(m) <- cells

## CopyKAT is single-threaded here (the sandbox blocks the sockets its cluster needs),
## so datasets are capped at max_cells, sampled proportionally across cell types.
set.seed(seed)
if (ncol(m) > max_cells) {
  f <- max_cells / ncol(m)
  idx <- unlist(lapply(split(seq_len(nrow(meta)), meta$cell_type), function(ii)
    if (length(ii) <= 3) ii else sample(ii, max(3, round(length(ii) * f)))))
  idx <- sort(idx)
  m <- m[, idx, drop = FALSE]; meta <- meta[idx, , drop = FALSE]
}
mat <- as.matrix(m)
storage.mode(mat) <- "numeric"

dir.create(file.path(ROOT, "calls"), showWarnings = FALSE)
res <- list(dataset = name, n_cells = ncol(mat), n_genes = nrow(mat), modes = list())
preds <- list()

for (mode in c("reffree", "refimm")) {
  norm_names <- if (mode == "refimm") meta$cell[meta$is_immune %in% c(TRUE, "True", "true")] else ""
  if (mode == "refimm" && length(norm_names) < 50) {
    res$modes[[mode]] <- list(skipped = "fewer than 50 immune cells"); next
  }
  t0 <- Sys.time()
  ck <- tryCatch(
    copykat(rawmat = mat, id.type = "S", ngene.chr = 5, win.size = 25, KS.cut = 0.1,
            sam.name = paste0(name, "_", mode), distance = "euclidean",
            norm.cell.names = norm_names, output.seg = "FALSE",
            plot.genes = "FALSE", genome = "hg20", n.cores = 1),
    error = function(e) list(error = conditionMessage(e)))
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (!is.null(ck$error)) {
    res$modes[[mode]] <- list(error = substr(ck$error, 1, 300), runtime_s = round(el, 1))
    next
  }
  p <- as.data.frame(ck$prediction, stringsAsFactors = FALSE)
  p$mode <- mode
  preds[[mode]] <- p
  tab <- table(p$copykat.pred)
  res$modes[[mode]] <- list(runtime_s = round(el, 1),
                            counts = as.list(setNames(as.integer(tab), names(tab))),
                            fpr_aneuploid = unname(sum(p$copykat.pred == "aneuploid") / nrow(p)))
}

if (length(preds)) {
  out <- do.call(rbind, preds)
  out <- merge(out, meta, by.x = "cell.names", by.y = "cell", all.x = TRUE)
  write.csv(out, file.path(ROOT, "calls", paste0(name, "_copykat.csv")), row.names = FALSE)
}
write(toJSON(res, auto_unbox = TRUE, null = "null"), file.path(ROOT, "calls", paste0(name, "_copykat.json")))
cat(name, "DONE\n")
