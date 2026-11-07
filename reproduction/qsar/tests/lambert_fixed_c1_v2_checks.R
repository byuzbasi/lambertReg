#!/usr/bin/env Rscript
driver <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L]))
root <- dirname(dirname(driver))
for (path in c("R/penalties.R", "R/lambert_family.R", "R/performance_models.R",
    "R/performance_run.R", "R/qsar_aquatic_realdata_v1.R", "R/diabetes_realdata_v1.R",
    "R/hitters_realdata_v1.R", "R/lambert_fixed_c1_realdata_v1.R",
    "R/lambert_fixed_c1_realdata_v2.R")) source(file.path(root, path))
e <- lf2_environment(root); count <- 0L
gate <- function(value) { stopifnot(isTRUE(value)); count <<- count + 1L }
default <- formals(ncvreg::ncvreg)$gamma
gate(eval(default, list(penalty = "MCP")) == 3)
gate(eval(default, list(penalty = "SCAD")) == 3.7)
for (dataset in c("diabetes", "qsar", "hitters")) {
  old <- lf1_environment(root)$rcv1_config(root, dataset, "production")
  cfg <- e$rcv1_config(root, dataset, "production")
  gate(identical(cfg$base$lambert_shapes, 1) && identical(cfg$base$mcp_gamma, 3) && identical(cfg$base$scad_a, 3.7))
  gate(cfg$total_units == 500L && identical(old$seeds, cfg$seeds))
  kept <- setdiff(names(cfg$base), c("schema", "mcp_gamma", "scad_a", "comparator_shapes"))
  gate(identical(old$base[kept], cfg$base[kept]))
  settings <- perf_settings(cfg$base)
  gate(all(table(settings$family)[c("lambert", "mcp", "scad")] == 1L))
}
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) {
  stopifnot(length(args) == 2L, args[1] == "--smoke")
  out <- file.path(root, args[2]); state <- e$rcv1_validate_complete(out, root, "qsar")
  gate(length(state$completed) == 4L)
  for (i in 1:4) {
    shard <- e$rcv1_read_shard(out, i, state); cv <- shard$fit$fit$cv
    gate(all(cv$shape[cv$family == "lambert"] == 1) &&
      all(cv$shape[cv$family == "mcp"] == 3) && all(cv$shape[cv$family == "scad"] == 3.7))
  }
  gate(any(grepl("interrupted_smoke", readLines(file.path(out, "progress.tsv")))))
  wrong <- shard; wrong$signature <- "wrong"
  gate(inherits(try(e$rcv1_validate_shard(wrong, state$tasks[4, ], state), silent = TRUE), "try-error"))
}
cat("PASS:", count, "default-shape checks; production fits: 0.\n")
