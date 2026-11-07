# Execution-only overlay. The fixed scientific sources remain unchanged.
lft_load <- function(root) {
  for (f in c("R/penalties.R", "R/lambert_family.R", "R/performance_run.R",
      "R/lambert_fixed_metrics_v1.R", "R/lambert_fixed_simulation_v1.R"))
    sys.source(file.path(root, f), envir = .GlobalEnv)
  lft_environment(root)
}

lft_environment <- function(root) {
  e <- lfb_sim_environment(root)
  sys.source(file.path(root, "R/performance_truba_v1.R"), envir = e)
  for (name in c("lft_sources", "lft_guard", "lft_preflight", "lft_release",
      "lft_scan", "lft_validation_pulse", "lft_options", "lft_main")) {
    e[[name]] <- get(name, envir = environment(lft_environment))
    environment(e[[name]]) <- e
  }
  e$perf_truba_release_id <- function() "lambert_fixed_truba_v1_a01"
  e$perf_truba_source_paths <- e$lft_sources
  e$perf_truba_operational_guard <- e$lft_guard
  e$perf_truba_runtime_preflight <- e$lft_preflight
  e$perf_truba_validate_release <- e$lft_release
  e$perf_truba_scan <- e$lft_scan
  e
}

lft_sources <- function() unique(c(perf_source_paths(),
  "R/performance_truba_v1.R", "R/lambert_fixed_truba_v1.R",
  "scripts/lambert_fixed_truba_v1.R", "scripts/validate_lambert_fixed_truba_v1.R",
  "scripts/analyze_lambert_fixed_truba_v1.R", "scripts/analyze_lambert_fixed_benchmark_v1.R",
  "R/lambert_fixed_c1_realdata_v1.R", "R/lambert_fixed_c1_realdata_v2.R",
  "tests/lambert_fixed_truba_v1_checks.R", "truba/run_lambert_fixed_v1.slurm",
  "truba/build_lambert_fixed_v1_bundle.R", "truba/README_TRUBA_LAMBERT_FIXED_V1.md"))

lft_guard <- function(stage, workers, max_seconds, task_limit, interval,
                      require_slurm = stage == "confirmation") {
  stopifnot(stage %in% c("smoke", "confirmation"), .Platform$OS.type == "unix")
  if (length(workers) != 1L || !is.finite(workers) || workers != floor(workers) ||
      !workers %in% if (stage == "smoke") 1:2 else 56L)
    stop("Use 1-2 smoke workers or exactly 56 confirmation workers")
  ceiling <- if (stage == "smoke") 300 else 70 * 3600
  if (length(max_seconds) != 1L || !is.finite(max_seconds) ||
      max_seconds <= 0 || max_seconds > ceiling) stop("Invalid soft dispatch budget")
  if (length(task_limit) != 1L || is.na(task_limit) || task_limit < 1 ||
      (!is.infinite(task_limit) && task_limit != floor(task_limit)) ||
      (stage == "confirmation" && !identical(task_limit, Inf)))
    stop("Task limits are positive smoke-only controls; confirmation requires Inf")
  if (length(interval) != 1L || !is.finite(interval) || interval <= 0 || interval > 60)
    stop("Heartbeat interval must be positive and at most 60 seconds")
  if (require_slurm && (!nzchar(Sys.getenv("SLURM_JOB_ID")) ||
      Sys.getenv("SLURM_CPUS_PER_TASK") != "56"))
    stop("Confirmation requires sbatch and exactly 56 allocated CPUs")
  invisible(TRUE)
}

lft_preflight <- function(root, stage, require_release = FALSE,
                          require_truba_runtime = FALSE) {
  perf_preflight(compile = FALSE, root = root)
  expected <- c(glmnet = "5.0", ncvreg = "3.16.0", mltools = "0.3.5")
  for (p in names(expected)) {
    if (!requireNamespace(p, quietly = TRUE)) stop("Required installed package missing: ", p)
    if (as.character(utils::packageVersion(p)) != expected[[p]])
      stop("Required package version: ", p, " ", expected[[p]], "; no automatic installation")
  }
  if (require_truba_runtime) {
    lib <- Sys.getenv("R_LIBS_USER")
    if (as.character(getRversion()) != "4.3.0" || !grepl("linux", R.version$platform) ||
        Sys.getenv("R_LIBS_USER") != lib || !lib %in% .libPaths())
      stop("Use the configured TRUBA Linux R 4.3.0 module and user library")
  }
  if (require_release) perf_truba_validate_release(root)
  s <- perf_software()
  cat("Lightweight fixed-Lambert preflight accepted; model fits: 0.\n", s$R, "\n", s$platform, "\n")
  print(unlist(s$packages)); cat("Library paths:\n", paste(.libPaths(), collapse = "\n"), "\n")
  invisible(s)
}

lft_release <- function(root, quiet = FALSE) {
  p <- file.path(root, "release", perf_truba_release_id(), "LOCAL_VALIDATED.rds")
  if (!file.exists(p)) stop("Missing local validation receipt")
  x <- readRDS(p)
  if (!identical(x$schema, "lambert_fixed_truba_release_v1") || !isTRUE(x$accepted) ||
      x$production_fits != 0L || !nrow(x$checks) || !all(x$checks$passed))
    stop("Invalid local validation receipt")
  if (!identical(x$sources, perf_truba_sources(root))) stop("Frozen source inventory changed")
  if (!identical(x$identity, perf_truba_portable_identity(root, "confirmation")))
    stop("Frozen scientific design or seeds changed")
  if (!identical(x$evidence, perf_truba_inventory(root, x$evidence$file)))
    stop("Local validation evidence changed")
  if (!quiet) cat("FROZEN_LOCAL_VALIDATION=PASS; production fits during validation: 0\n")
  invisible(x)
}

lft_validation_pulse <- function(out, done, total, now = as.numeric(Sys.time())) {
  # Only the active aggregation writer may update an unfinished run. Portable
  # or completed-output verification must never modify archived inputs.
  owner_path <- file.path(paste0(out, ".run_lock"), "owner.json")
  path <- file.path(out, "progress.json")
  if (!file.exists(owner_path) || !file.exists(path) ||
      file.exists(file.path(out, "manifest.csv")) || file.exists(file.path(out, "COMPLETED")))
    return(invisible(FALSE))
  owner <- jsonlite::fromJSON(owner_path)
  if (!identical(as.integer(owner$pid), Sys.getpid()) ||
      !identical(owner$host, unname(Sys.info()[["nodename"]]))) return(invisible(FALSE))
  p <- jsonlite::fromJSON(path)
  if (!identical(p$phase, "aggregation_and_validation")) return(invisible(FALSE))
  delta <- max(0, now - as.numeric(as.POSIXct(p$last_heartbeat_utc,
    format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")))
  p$elapsed_seconds <- p$elapsed_seconds + delta
  p$invocation_elapsed_seconds <- p$invocation_elapsed_seconds + delta
  p$throughput_units_per_second <- p$completed_this_invocation / p$invocation_elapsed_seconds
  deadline <- suppressWarnings(as.numeric(Sys.getenv("SLURM_JOB_END_TIME", NA_character_)))
  p$slurm_remaining_seconds <- if (is.finite(deadline)) max(0, deadline - now) else NA_real_
  p$last_heartbeat_utc <- perf_truba_utc(now)
  p$detail <- sprintf("aggregation validation: %d/%d stored units checked", done, total)
  perf_atomic_json(p, path)
  p <- lapply(p, function(x) if (is.null(x)) NA else x)
  write.table(as.data.frame(p), file.path(out, "progress.tsv"), sep = "\t",
    row.names = FALSE, col.names = FALSE, append = TRUE, quote = TRUE, na = "NA")
  invisible(TRUE)
}

# Same shard, numerical and manifest gates as the shared scanner, with bounded
# model memory: retain only small method rows, not all fitted-model/CV objects.
lft_scan <- function(out, root, complete = FALSE, portable = FALSE, quiet = FALSE) {
  state <- if (portable) perf_truba_verify_metadata_portable(out, root) else perf_verify_metadata(out, root)
  if (portable) {
    pkgs <- state$stored_software$packages
    if (!identical(pkgs$mltools, "0.3.5") || !identical(pkgs$glmnet, "5.0") ||
        !identical(pkgs$ncvreg, "3.16.0")) stop("Archived scientific package contract differs")
  }
  state$execution_signature <- perf_truba_verify_execution_metadata(out, root, state$cfg$stage)
  state$attempts <- perf_truba_read_attempts(out, state)
  shard_root <- file.path(out, "shards")
  if (!dir.exists(shard_root)) stop("Missing shard directory")
  actual <- list.dirs(shard_root, recursive = FALSE, full.names = FALSE)
  if (any(!grepl("^unit_[0-9]{5}$", actual))) stop("Unexpected or partial shard directory")
  present <- vapply(state$tasks$unit, function(i) dir.exists(perf_unit_path(out, i)), logical(1))
  if (!setequal(actual, sprintf("unit_%05d", which(present)))) stop("Unexpected shard identities")
  rows <- vector("list", sum(present)); k <- 0L; last <- proc.time()[[3]]
  for (i in which(present)) {
    k <- k + 1L
    rows[[k]] <- perf_truba_read_unit(out, i, state)$rows
    if (proc.time()[[3]] - last >= 30) {
      cat(sprintf("[%s] validation | %d/%d stored units checked | ETA unknown; no model fits\n",
        perf_truba_utc(as.numeric(Sys.time())), k, sum(present)))
      flush.console(); lft_validation_pulse(out, k, sum(present)); last <- proc.time()[[3]]
    }
  }
  state$completed <- which(present)
  state$rows <- if (length(rows)) do.call(rbind, rows) else NULL
  if (complete) {
    if (!all(present) || !file.exists(file.path(out, "COMPLETED"))) stop("Run incomplete")
    if (length(list.files(file.path(out, "failures"), all.files = TRUE, no.. = TRUE)))
      stop("Unresolved failure evidence")
    perf_verify_manifest(out)
    r <- readRDS(file.path(out, "results.rds"))
    if (!identical(r, state$rows)) stop("Aggregate differs from validated shards")
    if (!isTRUE(all.equal(read.csv(file.path(out, "results.csv")), r,
        tolerance = 1e-12, check.attributes = FALSE))) stop("CSV differs from aggregate")
    p <- jsonlite::fromJSON(file.path(out, "progress.json"))
    if (p$phase != "completed" || p$validated_completed != nrow(state$tasks) ||
        p$running != 0 || p$failed != 0 || p$pending != 0 ||
        !identical(p$scientific_signature, state$signature) ||
        !identical(p$execution_signature, state$execution_signature)) stop("Final progress inconsistent")
    if (!identical(readLines(file.path(out, "COMPLETED")), state$signature))
      stop("Completion signature mismatch")
    if (!quiet) cat(sprintf(paste0("COMPLETE: %d validated units, %d method rows; ",
      "CV/seed/source/execution/shard/manifest gates verified.\n",
      "Failed method fits retained: %d; partial searches: %d.\n"),
      nrow(state$tasks), nrow(r), sum(!r$ok), sum(r$partial_search)))
  }
  invisible(state)
}

lft_options <- function(args) {
  if (any(!grepl("^--[a-z-]+=.+$", args))) stop("Use --key=value arguments")
  keys <- sub("^--([^=]+)=.*$", "\\1", args)
  if (anyDuplicated(keys) || any(!keys %in% c("stage", "action", "out", "workers",
      "max-seconds", "task-limit", "interval"))) stop("Unknown or duplicate option")
  x <- setNames(as.list(sub("^--[^=]+=", "", args)), keys)
  defaults <- list(stage = "confirmation", action = "design", workers = "56",
    `max-seconds` = "252000", `task-limit` = "Inf", interval = "30",
    out = "outputs/performance/lambert_fixed_simulation_truba_v1_a01")
  for (k in setdiff(names(defaults), keys)) x[[k]] <- defaults[[k]]
  x$stage <- match.arg(x$stage, c("smoke", "confirmation"))
  x$action <- match.arg(x$action, c("design", "preflight", "start", "resume", "verify", "portable-verify"))
  if (x$stage == "smoke") {
    if (!"workers" %in% keys) x$workers <- "2"
    if (!"max-seconds" %in% keys) x[["max-seconds"]] <- "300"
    if (!"out" %in% keys) stop("Smoke requires an explicit versioned --out")
  }
  x
}

lft_main <- function(root, args) {
  x <- lft_options(args)
  out <- if (startsWith(x$out, "/")) x$out else file.path(root, x$out)
  if (x$action == "design") {
    cfg <- perf_config(root, x$stage); tasks <- perf_tasks(cfg)
    cat("DESIGN ONLY:", nrow(tasks), "units;", nrow(tasks) * length(cfg$methods), "method rows; model fits: 0.\n")
    return(invisible(tasks))
  }
  if (x$action == "preflight") return(perf_truba_runtime_preflight(root, x$stage,
    require_release = x$stage == "confirmation",
    require_truba_runtime = Sys.getenv("LFT_REQUIRE_TRUBA_RUNTIME") == "1"))
  if (x$action %in% c("verify", "portable-verify")) return(perf_truba_scan(out, root,
    complete = TRUE, portable = x$action == "portable-verify"))
  workers <- as.numeric(x$workers); seconds <- as.numeric(x[["max-seconds"]])
  limit <- as.numeric(x[["task-limit"]]); interval <- as.numeric(x$interval)
  perf_truba_operational_guard(x$stage, workers, seconds, limit, interval)
  if (x$action == "start" && file.exists(out)) stop("Fresh start refuses existing output")
  if (x$action == "resume" && !dir.exists(out)) stop("Resume requires an existing TRUBA output")
  perf_truba_runtime_preflight(root, x$stage, x$stage == "confirmation", x$stage == "confirmation")
  # Compile the existing Rcpp solver once; forked workers inherit it.
  perf_preflight(compile = TRUE, root = root)
  z <- perf_truba_run(root, out, x$stage, workers, seconds, limit, interval,
    resume = x$action == "resume")
  if (!isTRUE(z$complete)) {
    cat("CHECKPOINTED: valid shards preserved; resume the identical frozen bundle.\n")
    return(invisible(3L))
  }
  invisible(0L)
}
