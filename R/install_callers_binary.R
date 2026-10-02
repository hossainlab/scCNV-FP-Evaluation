## Install the malignant-cell callers. On this Windows host there is no in-env
## compiler, so everything that can come as a pre-built Windows binary must:
## forcing type="binary" stops R from preferring newer *source* releases.
LOG <- "install2_log.txt"
say <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)),
                         file = LOG, append = TRUE)

## Optional: path to a conda env holding m2w64-toolchain (Windows only).
BT <- Sys.getenv("R_BUILDTOOLS", "")
tc <- c(file.path(BT, "Library/mingw-w64/bin"), file.path(BT, "Library/bin"))
if (all(dir.exists(tc))) {
  Sys.setenv(PATH = paste(c(tc, Sys.getenv("PATH")), collapse = .Platform$path.sep))
  Sys.setenv(BINPREF = paste0(BT, "/Library/mingw-w64/bin/"))
}
options(repos = c(CRAN = "https://cloud.r-project.org"),
        pkgType = "binary",
        install.packages.check.source = "no",
        install.packages.compile.from.source = "never",
        timeout = 3600, Ncpus = 4)
say("R ", getRversion(), " | pkgType=", getOption("pkgType"))

have <- function(p) requireNamespace(p, quietly = TRUE)
bin_install <- function(pkgs, repos = getOption("repos"), label = "CRAN") {
  todo <- pkgs[!vapply(pkgs, have, TRUE)]
  if (!length(todo)) return(invisible(NULL))
  for (p in todo) {
    ok <- tryCatch({
      suppressWarnings(install.packages(p, repos = repos, type = "binary",
                                        dependencies = c("Depends", "Imports", "LinkingTo"),
                                        quiet = TRUE))
      have(p)
    }, error = function(e) { say(label, ":", p, " ERROR ", conditionMessage(e)); FALSE })
    say(label, ":", p, if (ok) " ok" else " FAILED")
  }
}

## --- recommended/base packages conda's r-base omits -------------------------
bin_install(c("MASS", "cluster", "survival", "nlme", "KernSmooth", "Matrix", "lattice",
              "mgcv", "codetools", "foreign", "nnet", "rpart", "spatial", "class"))

## --- compiled CRAN dependencies of the callers ------------------------------
bin_install(c("Rcpp", "RcppArmadillo", "RcppParallel", "RcppEigen", "bitops", "caTools",
              "gtools", "gplots", "parallelDist", "dlm", "mixtools", "MCMCpack",
              "pheatmap", "Rtsne", "ape", "fastcluster", "igraph", "data.table"))

## --- pure-R CRAN dependencies ------------------------------------------------
bin_install(c("dplyr", "ggplot2", "forcats", "ggrepel", "doParallel", "foreach",
              "beanplot", "squash", "coda", "futile.logger", "argparse", "reshape",
              "reshape2", "tidyr", "tibble", "purrr", "stringr", "yulab.utils",
              "tidytree", "ggpubr", "cowplot", "gridExtra", "plyr", "e1071",
              "digest", "jsonlite", "httr", "curl", "XML", "xml2", "rjson"))

## --- Bioconductor (Windows binaries; avoid the archive mirror BiocManager probes)
if (!have("BiocManager")) install.packages("BiocManager", type = "binary", quiet = TRUE)
bioc_repos <- c(
  BioCsoft = sprintf("https://bioconductor.org/packages/%s/bioc", BiocManager::version()),
  BioCann  = sprintf("https://bioconductor.org/packages/%s/data/annotation", BiocManager::version()),
  BioCexp  = sprintf("https://bioconductor.org/packages/%s/data/experiment", BiocManager::version()),
  CRAN     = "https://cloud.r-project.org")
say("bioc repos: ", bioc_repos[["BioCsoft"]])
bin_install(c("BiocGenerics", "S4Vectors", "IRanges", "GenomeInfoDb", "GenomicRanges",
              "Biobase", "SummarizedExperiment", "SingleCellExperiment", "DelayedArray",
              "BiocParallel", "beachmat", "scuttle", "scran", "edgeR", "limma",
              "fgsea", "biomaRt", "ggtree", "phylogram", "rjags", "infercnv"),
            repos = bioc_repos, label = "BIOC")

## --- GitHub-only callers -----------------------------------------------------
gh <- list(copykat   = c("navinlabcode/copykat", NA),
           yaGST     = c("miccec/yaGST", NA),
           SCEVAN    = c("AntonioDeFalco/SCEVAN", NA),
           CONICSmat = c("diazlab/CONICS", "CONICSmat"))
for (nm in names(gh)) {
  if (have(nm)) { say("GH:", nm, " already present"); next }
  g <- gh[[nm]]
  ok <- tryCatch({
    remotes::install_github(g[1], subdir = if (is.na(g[2])) NULL else g[2],
                            upgrade = "never", dependencies = TRUE,
                            build = FALSE, quiet = TRUE, type = "binary")
    have(nm)
  }, error = function(e) { say("GH:", nm, " ERROR ", substr(conditionMessage(e), 1, 300)); FALSE })
  say("GH:", nm, if (ok) " ok" else " FAILED")
}

callers <- c("infercnv", "copykat", "SCEVAN", "CONICSmat")
st <- vapply(callers, function(p)
  tryCatch(as.character(utils::packageVersion(p)), error = function(e) "MISSING"), "")
say("CALLER STATUS: ", paste(names(st), st, sep = "=", collapse = " | "))
say("rjags: ", have("rjags"), " (inferCNV HMM mode requires it)")
writeLines(jsonlite::toJSON(as.list(st), auto_unbox = TRUE), file.path("logs", "caller_status.json"))
say("ALL DONE")
