## Install the malignant-cell callers audited in this study.
LOG <- "install_callers_log.txt"
say <- function(...) { cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)),
                           file = LOG, append = TRUE) }
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 1800,
        Ncpus = max(1, parallel::detectCores() - 2))

## Compilers live in a SEPARATE conda env: putting m2w64-toolchain in the R env
## itself makes the R process fail at startup with a system-DLL relocation error,
## so we add the toolchain to PATH at build time only.
## Optional: path to a conda env holding m2w64-toolchain (Windows only).
BT <- Sys.getenv("R_BUILDTOOLS", "")
tc <- c(file.path(BT, "Library/mingw-w64/bin"), file.path(BT, "Library/bin"))
if (all(dir.exists(tc))) {
  Sys.setenv(PATH = paste(c(tc, Sys.getenv("PATH")), collapse = .Platform$path.sep))
  Sys.setenv(BINPREF = paste0(BT, "/Library/mingw-w64/bin/"))
  say("toolchain on PATH: ", paste(tc, collapse = ";"))
} else {
  say("WARNING: buildtools toolchain not found; packages needing compilation will fail")
}
say("R ", R.version.string, " | libpath ", .libPaths()[1])
say("gcc: ", Sys.which("gcc"), " | make: ", Sys.which("make"))

try_install <- function(label, expr) {
  ok <- tryCatch({ eval(expr); TRUE },
                 error = function(e) { say(label, " ERROR: ", conditionMessage(e)); FALSE },
                 warning = function(w) { say(label, " WARN: ", conditionMessage(w)); TRUE })
  say(label, if (ok) " done" else " failed")
  ok
}

## --- CRAN dependencies (Windows binaries; no compilation needed) -------------
cran <- c("parallelDist", "dlm", "gplots", "RColorBrewer", "mixtools", "cluster",
          "MCMCpack", "transport", "pheatmap", "forcats", "dplyr", "ggplot2",
          "Rtsne", "ape", "ggrepel", "doParallel", "beanplot", "squash",
          "tidytree", "fastcluster", "Matrix", "remotes", "coda", "futile.logger",
          "argparse", "reshape")
for (p in cran) {
  if (!requireNamespace(p, quietly = TRUE)) {
    try_install(paste0("CRAN:", p), bquote(install.packages(.(p))))
  }
}
say("CRAN present: ", paste(cran[vapply(cran, function(p) requireNamespace(p, quietly = TRUE), TRUE)], collapse = ","))
say("CRAN MISSING: ", paste(cran[!vapply(cran, function(p) requireNamespace(p, quietly = TRUE), TRUE)], collapse = ","))

## --- Bioconductor ------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
say("Bioc version target: ", as.character(BiocManager::version()))
bioc <- c("infercnv", "scran", "fgsea", "SingleCellExperiment", "SummarizedExperiment",
          "edgeR", "limma", "ggtree", "biomaRt")
for (p in bioc) {
  if (!requireNamespace(p, quietly = TRUE)) {
    try_install(paste0("BIOC:", p),
                bquote(BiocManager::install(.(p), ask = FALSE, update = FALSE, force = FALSE)))
  }
}

## --- GitHub-only callers -----------------------------------------------------
gh <- list(
  copykat  = list(repo = "navinlabcode/copykat",   subdir = NULL),
  yaGST    = list(repo = "miccec/yaGST",           subdir = NULL),
  SCEVAN   = list(repo = "AntonioDeFalco/SCEVAN",  subdir = NULL),
  CONICSmat= list(repo = "diazlab/CONICS",         subdir = "CONICSmat")
)
for (nm in names(gh)) {
  if (!requireNamespace(nm, quietly = TRUE)) {
    g <- gh[[nm]]
    try_install(paste0("GH:", nm), bquote(remotes::install_github(
      .(g$repo), subdir = .(g$subdir), upgrade = "never", dependencies = TRUE,
      build_vignettes = FALSE, quiet = FALSE)))
  }
}

## --- Final status ------------------------------------------------------------
callers <- c("infercnv", "copykat", "SCEVAN", "CONICSmat")
st <- vapply(callers, function(p) {
  v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
  if (is.na(v)) "MISSING" else v
}, "")
say("CALLER STATUS: ", paste(names(st), st, sep = "=", collapse = " | "))
say("rjags available: ", requireNamespace("rjags", quietly = TRUE))
say("ALL DONE")
writeLines(jsonlite::toJSON(as.list(st), auto_unbox = TRUE), file.path("logs", "caller_status.json"))
