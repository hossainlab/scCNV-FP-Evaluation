## Batch CopyKAT over the panel. CopyKAT is single-threaded here and emits a large
## amount of progress text plus per-run plots, so each dataset runs inside its own
## scratch directory with output captured to a log.
ROOT <- normalizePath(".")
sets <- sort(list.dirs(file.path(ROOT, "rmat"), recursive = FALSE, full.names = FALSE))
CAP <- "600"
logf <- file.path(ROOT, "copykat_log.txt")
if (!file.exists(logf)) cat("", file = logf)

for (s in sets) {
  if (file.exists(file.path(ROOT, "calls", paste0(s, "_copykat.json")))) {
    cat(sprintf("[%s] %s SKIP\n", format(Sys.time(), "%H:%M:%S"), s), file = logf, append = TRUE)
    next
  }
  scratch <- file.path(ROOT, "ck_scratch", s)
  dir.create(scratch, recursive = TRUE, showWarnings = FALSE)
  t0 <- Sys.time()
  ok <- FALSE
  con <- file(file.path(scratch, "stdout.txt"), open = "wt")
  setwd(scratch)
  try({
    sink(con, type = "output"); sink(con, type = "message")
    TARGET <- c(s, CAP)
    source(file.path(ROOT, "R", "run_copykat.R"), local = TRUE)
    ok <- TRUE
  }, silent = TRUE)
  try(sink(type = "message"), silent = TRUE); try(sink(type = "output"), silent = TRUE)
  close(con); setwd(ROOT)
  el <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  cat(sprintf("[%s] %s %s %s min\n", format(Sys.time(), "%H:%M:%S"), s,
              if (ok) "OK" else "FAIL", el), file = logf, append = TRUE)
}
cat("ALL DONE\n", file = logf, append = TRUE)
