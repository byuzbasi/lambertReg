#!/usr/bin/env Rscript
# Portable runner for the manuscript's frozen simulation.
main <- function(args) {
  if (!length(args) || !args[1] %in% c("design", "smoke", "start", "resume", "verify"))
    stop("Usage: Rscript --vanilla run.R design | smoke/start/resume/verify OUTPUT [WORKERS]")
  action <- args[1]
  if (action != "design" && !length(args) %in% 2:3) stop("Specify the output directory")
  workers <- if (length(args) == 3L) as.integer(args[3]) else 1L
  if (is.na(workers) || workers < 1L || workers > 56L) stop("WORKERS must be between 1 and 56")
  source("R/lambert_simulation_main_v2.R")
  stage <- if (action == "verify") readRDS(file.path(args[2], "config.rds"))$study_stage
    else if (action == "smoke") "smoke" else "main"
  root <- normalizePath(".")
  e <- lms_environment(root, stage, "fast")
  for (f in c("lambert_bic_comparison_v1.R", "lambert_ic_pilot_v1.R", "ic_scoring.R"))
    sys.source(file.path("R", f), e)
  cfg <- lms_config(root, stage, "fast")
  spec <- e$lbc_spec(root)
  tasks <- lms_tasks(cfg)
  bic_set <- with(tasks, profile == spec$profile & p %in% spec$p &
                  rho %in% spec$rho & snr %in% spec$snr)
  if (stage == "main")
    stopifnot(nrow(tasks) == 15600L, length(unique(tasks$scenario)) == 52L, sum(bic_set) == 4800L)
  if (action == "design") {
    cat("52 CV scenarios x 300 replications, 16 paired BIC scenarios x 300,",
        "known-support Oracle in all 52 scenarios; model fits: 0.\n")
    return(invisible(NULL))
  }

  # Only the dispatcher, paths, and run receipts are portable adaptations.
  # The generator, tuning, estimates and numerical gates remain in the engines.
  e$perf_source_paths <- e$perf_truba_source_paths <- function()
    c("run.R", list.files(c("R", "src", "scripts", "protocol"), recursive = TRUE, full.names = TRUE))
  e$perf_truba_release_id <- function() "lambert_simulation_public_v1"
  e$perf_truba_validate_release <- function(root, quiet = FALSE) {
    stopifnot(all(file.exists(file.path(root, e$perf_source_paths()))))
    invisible(TRUE)
  }
  config <- e$perf_config
  e$perf_config <- function(root, stage) {
    result <- config(root, stage)
    result$bic <- spec
    result$ic <- list(reproduction_tolerance = spec$reproduction_tolerance)
    result
  }
  base_task <- e$perf_run_task
  e$perf_run_task <- function(task, cfg, seeds, pulse = function(...) NULL) {
    result <- base_task(task, cfg, seeds, pulse)
    train <- e$perf_training(task, cfg, seeds)
    test <- e$perf_test_data(task, cfg, seeds, train$truth)
    prep <- e$perf_prepare(train$X, train$y, cfg)
    oracle <- e$lbc_oracle(prep, cfg)
    e$lbc_oracle_check(oracle, prep, cfg)
    stopifnot(oracle$ok)
    ref <- e$lic_metrics(oracle$beta, prep, test, train$truth, cfg, e)
    result$oracle <- cbind(task, ref[c("mse", "normalized_mse", "sigma2")], kkt = oracle$kkt)
    if (stage == "smoke" || bic_set[task$unit]) {
      settings <- e$perf_settings(cfg)
      paths <- lapply(seq_len(nrow(settings)), function(j)
        e$perf_path(prep, settings[j, ], cfg, pulse))
      cache <- new.env(parent = emptyenv())
      candidates <- lapply(seq_along(paths), function(j)
        e$lbc_candidates(paths[[j]], prep, settings[j, ], cfg, e, cache))
      verified <- e$lbc_reproduction(paths, result$rows, prep, test, train$truth, cfg, e)
      stopifnot(all(verified$passed))
      task$source_unit <- task$unit
      result$bic <- e$lbc_rows(task, cfg, paths, candidates, result$rows, oracle,
                              prep, test, train$truth, e)
      result$information_criteria <- e$score_ic(task, cfg, paths, candidates,
        result$bic, prep, test, train$truth, e)
    }
    result
  }
  out <- args[2]
  e$perf_preflight(TRUE, root)
  if (action != "verify") {
    completed <- e$perf_truba_run(root, out, cfg$stage, workers,
      max_seconds = if (stage == "smoke") 300 else 252000, interval = 30,
      resume = action == "resume")
    if (!completed$complete) return(invisible(NULL))
  }
  state <- e$perf_truba_scan(out, root, complete = TRUE)
  shards <- file.path(out, "shards", sprintf("unit_%05d", tasks$unit), "result.rds")
  extras <- lapply(shards, function(path) readRDS(path)[c("bic", "oracle", "information_criteria")])
  # Write derived exports once and verify any existing copies.
  for (kind in c("oracle", "bic", "information_criteria")) {
    value <- do.call(rbind, lapply(extras, function(x) x[[kind]]))
    path <- paste0(out, "-", kind, ".csv")
    if (file.exists(path)) {
      expected <- capture.output(write.csv(value, row.names = FALSE))
      stopifnot(identical(readLines(path), expected))
    } else e$perf_atomic(path, function(tmp) write.csv(value, tmp, row.names = FALSE))
  }
}
threads <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")
do.call(Sys.setenv, setNames(as.list(rep("1", length(threads))), threads))
main(commandArgs(TRUE))
