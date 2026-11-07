# Repeated low-dimensional nested CV. The frozen estimator and single-split
# dataset adapters are used without modification.

rcv1_engine <- function(root, dataset) {
  if (!dataset %in% c("diabetes", "qsar", "hitters"))
    stop("Unknown repeated-CV dataset")
  if (dataset == "diabetes") return(diabetes_v1_engine(root))
  if (dataset == "hitters") return(hitters_v1_engine(root))
  engine <- new.env(parent = globalenv())
  sys.source(file.path(root, "R/qsar_aquatic_realdata_v1.R"), envir = engine)
  engine
}

rcv1_seed_table <- function(base, repeats) {
  stopifnot(repeats %in% c(2L, 10L))
  offsets <- as.integer(1000L * (seq_len(repeats) - 1L))
  seeds <- data.frame(repeat_id = seq_len(repeats),
    outer_seed = as.integer(base$outer_seed + offsets))
  for (fold in seq_len(5L))
    seeds[[paste0("inner_seed_", fold)]] <-
      as.integer(base$inner_seeds[[fold]] + offsets)
  if (anyDuplicated(unlist(seeds[-1L])) || any(unlist(seeds[-1L]) <= 0L))
    stop("Repeated-CV seed table is invalid")
  seeds
}

rcv1_config <- function(root, dataset, mode, engine = rcv1_engine(root, dataset)) {
  if (!mode %in% c("smoke", "production")) stop("Invalid repeated-CV mode")
  base <- engine$qsar_v1_config(root, mode)
  if (base$outer_folds != 5L || base$folds != 5L ||
      !identical(engine$qsar_v1_methods(), c("ols", base$methods)))
    stop("Frozen five-fold method contract changed")
  repeats <- if (mode == "smoke") 2L else 10L
  active_folds <- if (mode == "smoke") 2L else 5L
  if (base$lambda_points != if (mode == "smoke") 4L else 60L)
    stop("Frozen penalty grid size changed")
  list(schema = "repeated_lowdim_realdata_v1", dataset = dataset,
    mode = mode, repeats = repeats, outer_folds = 5L,
    active_folds = active_folds, inner_folds = 5L,
    total_units = repeats * active_folds,
    base = base, seeds = rcv1_seed_table(base, repeats),
    methods = engine$qsar_v1_methods(),
    primary_summary = "mean_of_complete_repeat_out_of_fold_rmse",
    supportive_summaries = c("median_repeat_rmse", "rmse_iqr",
      "repeat_mae", "repeat_medae", "paired_repeat_differences"))
}

rcv1_rep_config <- function(cfg, repeat_id) {
  if (repeat_id < 1L || repeat_id > cfg$repeats)
    stop("Invalid repetition identity")
  base <- cfg$base
  seeds <- cfg$seeds[cfg$seeds$repeat_id == repeat_id, , drop = FALSE]
  if (nrow(seeds) != 1L) stop("Missing repetition seeds")
  base$outer_seed <- as.integer(seeds$outer_seed)
  base$inner_seeds <- as.integer(unlist(seeds[paste0("inner_seed_", 1:5)],
    use.names = FALSE))
  base
}

rcv1_tasks <- function(cfg) {
  data.frame(unit = seq_len(cfg$total_units),
    repeat_id = rep(seq_len(cfg$repeats), each = cfg$active_folds),
    outer_fold = rep(seq_len(cfg$active_folds), times = cfg$repeats))
}

rcv1_folds <- function(cfg, engine) {
  lapply(seq_len(cfg$repeats), function(repeat_id)
    engine$qsar_v1_folds(rcv1_rep_config(cfg, repeat_id)))
}

rcv1_validate_design <- function(data, folds, cfg) {
  if (!identical(dim(data$X), c(cfg$base$n, cfg$base$p)) ||
      length(data$y) != cfg$base$n || length(folds) != cfg$repeats)
    stop("Repeated-CV data dimensions changed")
  for (repeat_id in seq_len(cfg$repeats)) {
    partition <- folds[[repeat_id]]
    if (length(partition$outer) != cfg$base$n ||
        !identical(sort(unique(partition$outer)), seq_len(5L)) ||
        length(partition$inner) != 5L)
      stop("Repeated-CV outer partition changed")
    for (outer_fold in seq_len(5L)) {
      train <- which(partition$outer != outer_fold)
      inner <- partition$inner[[outer_fold]]
      if (length(inner) != length(train) ||
          !identical(sort(unique(inner)), seq_len(5L)))
        stop("Repeated-CV inner partition changed")
      for (inner_fold in seq_len(5L)) {
        design <- data$X[train[inner != inner_fold], , drop = FALSE]
        if (any(apply(design, 2L, stats::sd) <= 0) ||
            qr(cbind(1, scale(design)), tol = 1e-10)$rank != cfg$base$p + 1L)
          stop(sprintf("Rank/variance gate failed at repeat %d outer %d inner %d",
            repeat_id, outer_fold, inner_fold))
      }
    }
  }
  invisible(TRUE)
}

rcv1_source_paths <- function(engine) unique(c(engine$qsar_v1_source_paths(),
  "R/repeated_lowdim_realdata_v1.R",
  "scripts/repeated_lowdim_realdata_v1.R",
  "protocol/repeated_lowdim_realdata_v1.md"))

rcv1_sources <- function(root, engine) {
  files <- rcv1_source_paths(engine)
  paths <- file.path(root, files)
  if (!all(file.exists(paths))) stop("Repeated-CV source file is missing")
  data.frame(file = files, bytes = as.numeric(file.info(paths)$size),
    sha256 = unname(vapply(paths, perf_hash, character(1))))
}

rcv1_raw_paths <- function(root, cfg) {
  base <- cfg$base
  if (cfg$dataset == "diabetes")
    return(file.path(root, base$raw_dir, diabetes_v1_raw_names()))
  if (cfg$dataset == "hitters")
    return(file.path(root, base$raw_dir, hitters_v1_raw_names()))
  file.path(root, base$raw_dir, c("qsar_aquatic_toxicity.zip",
    "qsar_aquatic_toxicity.csv", "ACQUISITION_RECEIPT.json", "manifest.csv"))
}

rcv1_raw_inputs <- function(root, cfg) {
  paths <- rcv1_raw_paths(root, cfg)
  if (!all(file.exists(paths)) || any(dir.exists(paths)))
    stop("Repeated-CV raw input is missing")
  data.frame(file = basename(paths), bytes = as.numeric(file.info(paths)$size),
    sha256 = unname(vapply(paths, perf_hash, character(1))))
}

rcv1_signature <- function(out) {
  files <- c("config.json", "sources.csv", "raw_inputs.csv", "tasks.csv",
    "folds.rds", "software.json")
  digest::digest(paste(files,
    vapply(file.path(out, files), perf_hash, character(1)),
    sep = ":", collapse = "\n"), algo = "sha256", serialize = FALSE)
}

rcv1_metadata_files <- function() c("config.rds", "config.json", "sources.csv",
  "raw_inputs.csv", "tasks.csv", "seeds.csv", "folds.rds", "software.json",
  "session_info.txt", "signature.sha256")

rcv1_write_metadata <- function(out, root, cfg, engine, data) {
  folds <- rcv1_folds(cfg, engine)
  rcv1_validate_design(data, folds, cfg)
  saveRDS(cfg, file.path(out, "config.rds"), version = 3)
  writeLines(perf_json(cfg), file.path(out, "config.json"))
  utils::write.csv(rcv1_sources(root, engine), file.path(out, "sources.csv"),
    row.names = FALSE)
  utils::write.csv(rcv1_raw_inputs(root, cfg), file.path(out, "raw_inputs.csv"),
    row.names = FALSE)
  utils::write.csv(rcv1_tasks(cfg), file.path(out, "tasks.csv"), row.names = FALSE)
  utils::write.csv(cfg$seeds, file.path(out, "seeds.csv"), row.names = FALSE)
  saveRDS(folds, file.path(out, "folds.rds"), version = 3)
  writeLines(perf_json(perf_software()), file.path(out, "software.json"))
  capture.output(sessionInfo(), .libPaths(), file = file.path(out, "session_info.txt"))
  writeLines(rcv1_signature(out), file.path(out, "signature.sha256"))
  files <- rcv1_metadata_files()
  utils::write.csv(data.frame(file = files,
    sha256 = unname(vapply(file.path(out, files), perf_hash, character(1)))),
    file.path(out, "metadata_manifest.csv"), row.names = FALSE)
  invisible(folds)
}

rcv1_verify_metadata <- function(out, root, dataset) {
  files <- rcv1_metadata_files()
  paths <- file.path(out, files)
  if (!all(file.exists(paths))) stop("Repeated-CV metadata missing")
  manifest <- utils::read.csv(file.path(out, "metadata_manifest.csv"),
    stringsAsFactors = FALSE)
  if (!identical(manifest$file, files) || anyDuplicated(manifest$file) ||
      !identical(manifest$sha256,
        unname(vapply(paths, perf_hash, character(1)))))
    stop("Repeated-CV metadata checksum mismatch")
  cfg <- readRDS(file.path(out, "config.rds"))
  engine <- rcv1_engine(root, dataset)
  if (!identical(cfg, rcv1_config(root, dataset, cfg$mode, engine)) ||
      !identical(paste(readLines(file.path(out, "config.json")), collapse = "\n"),
        perf_json(cfg))) stop("Repeated-CV configuration mismatch")
  sources <- utils::read.csv(file.path(out, "sources.csv"),
    stringsAsFactors = FALSE)
  if (!isTRUE(all.equal(sources, rcv1_sources(root, engine),
      check.attributes = FALSE)) ||
      !isTRUE(all.equal(utils::read.csv(file.path(out, "raw_inputs.csv")),
        rcv1_raw_inputs(root, cfg), check.attributes = FALSE)))
    stop("Repeated-CV source or raw-input signature mismatch")
  tasks <- utils::read.csv(file.path(out, "tasks.csv"))
  seeds <- utils::read.csv(file.path(out, "seeds.csv"))
  folds <- readRDS(file.path(out, "folds.rds"))
  if (!identical(tasks, rcv1_tasks(cfg)) || !identical(seeds, cfg$seeds) ||
      !identical(folds, rcv1_folds(cfg, engine)) ||
      !identical(paste(readLines(file.path(out, "software.json")), collapse = "\n"),
        perf_json(perf_software())) ||
      !identical(readLines(file.path(out, "signature.sha256")),
        rcv1_signature(out)))
    stop("Repeated-CV seeds, folds, software, or signature mismatch")
  data <- engine$qsar_v1_data(root, cfg$base)
  rcv1_validate_design(data, folds, cfg)
  list(cfg = cfg, engine = engine, data = data, tasks = tasks, folds = folds,
    signature = readLines(file.path(out, "signature.sha256")))
}

rcv1_unit_state <- function(state, repeat_id) {
  list(cfg = rcv1_rep_config(state$cfg, repeat_id), data = state$data,
    folds = state$folds[[repeat_id]], signature = state$signature)
}

rcv1_shard_path <- function(out, unit)
  file.path(out, "shards", sprintf("unit_%05d", unit))

rcv1_validate_shard <- function(shard, task, state) {
  if (!identical(as.integer(shard$unit), as.integer(task$unit)) ||
      !identical(as.integer(shard$repeat_id), as.integer(task$repeat_id)) ||
      !identical(as.integer(shard$outer_fold), as.integer(task$outer_fold)) ||
      !identical(shard$signature, state$signature))
    stop("Repeated-CV shard identity mismatch")
  state$engine$qsar_v1_validate_unit(shard$fit,
    as.integer(task$outer_fold), rcv1_unit_state(state, task$repeat_id))
  invisible(TRUE)
}

rcv1_read_shard <- function(out, unit, state) {
  path <- rcv1_shard_path(out, unit)
  file <- file.path(path, "result.rds")
  hash <- file.path(path, "result.sha256")
  if (!file.exists(file) || !file.exists(hash) ||
      !identical(readLines(hash), perf_hash(file)))
    stop("Missing or invalid repeated-CV shard: ", unit)
  shard <- readRDS(file)
  rcv1_validate_shard(shard, state$tasks[unit, , drop = FALSE], state)
  shard
}

rcv1_state <- function(out, root, dataset) {
  state <- rcv1_verify_metadata(out, root, dataset)
  paths <- list.dirs(file.path(out, "shards"), recursive = FALSE,
    full.names = FALSE)
  expected <- sprintf("unit_%05d", state$tasks$unit)
  if (!all(paths %in% expected)) stop("Unexpected or partial shard directory")
  failures <- list.files(file.path(out, "failures"), all.files = TRUE,
    no.. = TRUE)
  if (length(failures)) stop("Preserved failure evidence requires diagnosis")
  state$completed <- which(expected %in% paths)
  for (unit in state$completed) rcv1_read_shard(out, unit, state)
  state
}

rcv1_write_shard <- function(out, shard, state) {
  unit <- shard$unit
  rcv1_validate_shard(shard, state$tasks[unit, , drop = FALSE], state)
  path <- rcv1_shard_path(out, unit)
  if (file.exists(path)) stop("Refusing to overwrite completed shard")
  temp <- paste0(path, ".partial-", Sys.getpid())
  if (!dir.create(temp)) stop("Cannot stage repeated-CV shard")
  saveRDS(shard, file.path(temp, "result.rds"), version = 3,
    compress = "gzip")
  writeLines(perf_hash(file.path(temp, "result.rds")),
    file.path(temp, "result.sha256"))
  if (!file.rename(temp, path)) stop("Cannot commit repeated-CV shard")
  invisible(TRUE)
}

rcv1_progress <- function(out, state) {
  e <- new.env(parent = emptyenv())
  previous <- if (file.exists(file.path(out, "progress.json")))
    jsonlite::fromJSON(file.path(out, "progress.json")) else NULL
  e$previous_elapsed <- if (is.null(previous)) 0 else
    as.numeric(previous$elapsed_seconds)
  e$started <- proc.time()[[3L]]
  e$last <- -Inf
  e$phase <- "preflight"
  e$unit <- NA_integer_
  e$failed <- 0L
  e$completed <- state$completed
  e$before <- length(state$completed)
  e$durations <- numeric()
  e$emit <- function(detail = "", force = FALSE) {
    now <- proc.time()[[3L]]
    if (!force && now - e$last < 30) return(invisible(NULL))
    e$last <- now
    elapsed <- now - e$started
    total <- state$cfg$total_units
    done <- length(e$completed)
    running <- as.integer(!is.na(e$unit))
    eta <- if (done == total) 0 else if (length(e$durations))
      mean(e$durations) * (total - done) else NA_real_
    stamp <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    snapshot <- list(run_id = basename(out),
      slurm_job_id = Sys.getenv("SLURM_JOB_ID", unset = NA_character_),
      phase = e$phase,
      work_unit = "one outer repetition and fold, all inner CV and methods",
      total = total, validated_completed = done,
      validated_before_invocation = e$before,
      completed_this_invocation = done - e$before,
      running = running, failed = e$failed,
      pending = total - done - running - e$failed,
      percent = 100 * done / total,
      elapsed_seconds = e$previous_elapsed + elapsed,
      invocation_elapsed_seconds = elapsed,
      throughput_units_per_second = if (length(e$durations) && elapsed > 0)
        length(e$durations) / elapsed else NA_real_,
      estimated_remaining_seconds = eta,
      estimated_completion_utc = if (is.finite(eta))
        format(Sys.time() + eta, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC") else
        NA_character_,
      eta_scope = "computation only; aggregation and validation excluded",
      slurm_remaining_seconds = NA_real_,
      last_heartbeat_utc = stamp, current_work_unit = e$unit,
      most_recently_completed_unit = if (done) tail(e$completed, 1L) else
        NA_integer_, detail = detail)
    perf_atomic_json(snapshot, file.path(out, "progress.json"))
    history <- file.path(out, "progress.tsv")
    utils::write.table(as.data.frame(snapshot), history, sep = "\t",
      row.names = FALSE, col.names = !file.exists(history),
      append = file.exists(history), quote = TRUE, na = "NA")
    cat(sprintf("[%s] %s | %d/%d (%.1f%%) | unit %s | ETA %s | %s\n",
      stamp, e$phase, done, total, snapshot$percent, e$unit,
      if (is.finite(eta)) paste(round(eta), "s, computation only") else
        "estimating", detail))
    flush(stdout())
    invisible(snapshot)
  }
  e
}

rcv1_save_csv <- function(value, path) {
  temp <- paste0(path, ".tmp-", Sys.getpid())
  if (file.exists(temp)) stop("Existing aggregate staging file")
  utils::write.csv(value, temp, row.names = FALSE)
  if (file.exists(path)) {
    if (perf_hash(path) != perf_hash(temp))
      stop("Existing aggregate differs: ", path)
    if (!file.remove(temp)) stop("Cannot remove identical temporary file")
  } else if (!file.rename(temp, path)) stop("Cannot commit aggregate")
}

rcv1_aggregate <- function(out, state, write = FALSE) {
  if (length(state$completed) != state$cfg$total_units)
    stop("Repeated-CV run is incomplete")
  units <- lapply(state$tasks$unit, function(unit)
    rcv1_read_shard(out, unit, state))
  rows <- do.call(rbind, lapply(units, function(item) {
    x <- item$fit$rows
    x$unit <- item$unit
    x$repeat_id <- item$repeat_id
    x$outer_fold <- item$outer_fold
    x[, c("unit", "repeat_id", "outer_fold",
      setdiff(names(x), c("unit", "repeat_id", "outer_fold"))), drop = FALSE]
  }))
  predictions <- do.call(rbind, lapply(units, function(item) {
    x <- item$fit$predictions
    x$unit <- item$unit
    x$repeat_id <- item$repeat_id
    x$outer_fold <- item$outer_fold
    x[, c("unit", "repeat_id", "outer_fold",
      setdiff(names(x), c("unit", "repeat_id", "outer_fold"))), drop = FALSE]
  }))
  repeat_metrics <- do.call(rbind, lapply(seq_len(state$cfg$repeats),
    function(repeat_id) do.call(rbind, lapply(state$cfg$methods, function(method) {
      rr <- rows[rows$repeat_id == repeat_id & rows$method == method,
        , drop = FALSE]
      pp <- predictions[predictions$repeat_id == repeat_id &
        predictions$method == method, , drop = FALSE]
      expected <- which(state$folds[[repeat_id]]$outer %in%
        seq_len(state$cfg$active_folds))
      if (nrow(rr) != state$cfg$active_folds || nrow(pp) != length(expected) ||
          anyDuplicated(pp$sample_id) || !setequal(pp$sample_id, expected) ||
          !isTRUE(all.equal(pp$observed, state$data$y[pp$sample_id])))
        stop("Repeated-CV prediction coverage mismatch")
      complete <- state$cfg$mode == "production" && all(rr$ok) &&
        all(is.finite(pp$predicted)) && length(expected) == state$cfg$base$n
      error <- if (complete) abs(pp$observed - pp$predicted) else numeric()
      data.frame(repeat_id = repeat_id, method = method,
        expected_observations = state$cfg$base$n,
        available_predictions = sum(is.finite(pp$predicted)),
        complete_predictions = complete,
        failed_folds = sum(!rr$ok),
        partial_search_folds = sum(rr$partial_search),
        rmse = if (complete) sqrt(mean(error^2)) else NA_real_,
        mae = if (complete) mean(error) else NA_real_,
        medae = if (complete) stats::median(error) else NA_real_)
    }))))
  method_summary <- do.call(rbind, lapply(state$cfg$methods, function(method) {
    x <- repeat_metrics[repeat_metrics$method == method, , drop = FALSE]
    complete <- state$cfg$mode == "production" &&
      nrow(x) == state$cfg$repeats && all(x$complete_predictions)
    data.frame(method = method, repeats_expected = state$cfg$repeats,
      repeats_complete = sum(x$complete_predictions),
      all_repeats_complete = complete,
      failed_folds = sum(x$failed_folds),
      partial_search_folds = sum(x$partial_search_folds),
      mean_rmse = if (complete) mean(x$rmse) else NA_real_,
      median_rmse = if (complete) stats::median(x$rmse) else NA_real_,
      q25_rmse = if (complete) unname(stats::quantile(x$rmse, .25)) else NA_real_,
      q75_rmse = if (complete) unname(stats::quantile(x$rmse, .75)) else NA_real_,
      mean_mae = if (complete) mean(x$mae) else NA_real_,
      median_mae = if (complete) stats::median(x$mae) else NA_real_,
      mean_medae = if (complete) mean(x$medae) else NA_real_,
      median_medae = if (complete) stats::median(x$medae) else NA_real_)
  }))
  rownames(rows) <- rownames(predictions) <- rownames(repeat_metrics) <-
    rownames(method_summary) <- NULL
  value <- list(results = rows, predictions = predictions,
    repeat_metrics = repeat_metrics, method_summary = method_summary)
  if (write) for (name in names(value))
    rcv1_save_csv(value[[name]], file.path(out, paste0(name, ".csv")))
  value
}

rcv1_validate_complete <- function(out, root, dataset) {
  state <- rcv1_state(out, root, dataset)
  if (length(state$completed) != state$cfg$total_units ||
      !file.exists(file.path(out, "COMPLETED")) ||
      !identical(readLines(file.path(out, "COMPLETED")), state$signature))
    stop("Repeated-CV run is not complete")
  perf_verify_manifest(out)
  aggregate <- rcv1_aggregate(out, state)
  for (name in names(aggregate)) {
    staging <- tempfile(pattern = "rcv1_validate_", fileext = ".csv")
    utils::write.csv(aggregate[[name]], staging, row.names = FALSE)
    same <- identical(perf_hash(staging),
      perf_hash(file.path(out, paste0(name, ".csv"))))
    if (!file.remove(staging)) stop("Cannot remove validation staging CSV")
    if (!same)
      stop("Repeated-CV aggregate mismatch: ", name)
  }
  progress <- jsonlite::fromJSON(file.path(out, "progress.json"))
  if (progress$phase != "completed" ||
      progress$validated_completed != state$cfg$total_units ||
      progress$running != 0L || progress$failed != 0L ||
      progress$pending != 0L ||
      nrow(aggregate$results) != state$cfg$total_units * 7L ||
      nrow(aggregate$repeat_metrics) != state$cfg$repeats * 7L ||
      nrow(aggregate$method_summary) != 7L)
    stop("Repeated-CV final dimensions or progress mismatch")
  if (state$cfg$mode == "production" &&
      nrow(aggregate$predictions) != state$cfg$repeats *
        state$cfg$base$n * 7L)
    stop("Repeated-CV prediction count mismatch")
  cat(sprintf("COMPLETE: %s %s, %d validated units, %d method rows, %d predictions; raw/source/seed/shard/manifest gates verified.\n",
    dataset, state$cfg$mode, state$cfg$total_units,
    nrow(aggregate$results), nrow(aggregate$predictions)))
  invisible(state)
}

rcv1_main <- function(root, args) {
  opt <- list(dataset = NULL, mode = NULL, out = NULL, resume = FALSE,
    stop_after = Inf)
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key == "--resume") {
      opt$resume <- TRUE; i <- i + 1L; next
    }
    if (!key %in% c("--dataset", "--mode", "--out", "--stop-after") ||
        i == length(args) || startsWith(args[[i + 1L]], "--"))
      stop("Invalid repeated-CV CLI arguments")
    opt[[gsub("-", "_", substring(key, 3L))]] <- args[[i + 1L]]
    i <- i + 2L
  }
  if (is.null(opt$dataset) || !opt$dataset %in%
      c("diabetes", "qsar", "hitters") || is.null(opt$mode) ||
      !opt$mode %in% c("smoke", "production", "validate") || is.null(opt$out))
    stop("Use --dataset diabetes|qsar|hitters --mode smoke|production|validate --out DIRECTORY [--resume]")
  out <- if (startsWith(opt$out, "/")) opt$out else
    file.path(root, opt$out)
  if (opt$mode == "validate") {
    if (opt$resume || !is.infinite(opt$stop_after))
      stop("Validate accepts only dataset and output")
    rcv1_validate_complete(out, root, opt$dataset)
    return(invisible(NULL))
  }
  stop_after <- suppressWarnings(as.numeric(opt$stop_after))
  if (is.na(stop_after) || stop_after < 1 ||
      (!is.infinite(stop_after) && opt$mode != "smoke"))
    stop("--stop-after is available only for smoke testing")
  if (file.exists(out) && !opt$resume)
    stop("Output exists; refusing overwrite")
  if (!file.exists(out) && opt$resume) stop("Resume output is missing")
  perf_preflight(compile = TRUE, root = root)
  engine <- rcv1_engine(root, opt$dataset)
  cfg <- rcv1_config(root, opt$dataset, opt$mode, engine)
  data <- engine$qsar_v1_data(root, cfg$base)
  rcv1_validate_design(data, rcv1_folds(cfg, engine), cfg)
  if (!opt$resume) {
    if (!dir.create(out, recursive = TRUE) ||
        !dir.create(file.path(out, "shards")) ||
        !dir.create(file.path(out, "failures")))
      stop("Cannot create repeated-CV output")
  }
  lock <- perf_lock(out)
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  if (!opt$resume) rcv1_write_metadata(out, root, cfg, engine, data)
  state <- rcv1_state(out, root, opt$dataset)
  if (!identical(state$cfg$mode, opt$mode)) stop("Resume mode mismatch")
  if (file.exists(file.path(out, "COMPLETED"))) {
    rcv1_validate_complete(out, root, opt$dataset)
    return(invisible(NULL))
  }
  progress <- rcv1_progress(out, state)
  progress$phase <- "computation"
  progress$emit("validated run started", force = TRUE)
  for (unit in setdiff(state$tasks$unit, state$completed)) {
    task <- state$tasks[unit, , drop = FALSE]
    progress$unit <- unit
    progress$emit(sprintf("repeat %d outer fold %d started",
      task$repeat_id, task$outer_fold), force = TRUE)
    started <- proc.time()[[3L]]
    shard <- tryCatch({
      fit <- state$engine$qsar_v1_fit_unit(task$outer_fold,
        rcv1_unit_state(state, task$repeat_id),
        function(detail) progress$emit(detail))
      result <- list(unit = as.integer(unit),
        repeat_id = as.integer(task$repeat_id),
        outer_fold = as.integer(task$outer_fold),
        signature = state$signature, fit = fit)
      rcv1_write_shard(out, result, state)
      result
    }, error = function(e) {
      evidence <- list(unit = unit, repeat_id = task$repeat_id,
        outer_fold = task$outer_fold, signature = state$signature,
        message = conditionMessage(e), call = conditionCall(e))
      saveRDS(evidence, file.path(out, "failures",
        sprintf("unit_%05d_error.rds", unit)), version = 3)
      progress$unit <- NA_integer_
      progress$failed <- 1L
      progress$phase <- "failed"
      progress$emit("failure evidence preserved", force = TRUE)
      stop(e)
    })
    state$completed <- c(state$completed, unit)
    progress$completed <- state$completed
    progress$durations <- c(progress$durations,
      proc.time()[[3L]] - started)
    progress$unit <- NA_integer_
    progress$emit("outer unit validated", force = TRUE)
    if (length(progress$durations) >= stop_after) {
      progress$phase <- "interrupted_smoke"
      progress$emit("smoke checkpoint ready", force = TRUE)
      return(invisible(NULL))
    }
  }
  progress$phase <- "aggregation_and_validation"
  progress$emit("reading all validated shards", force = TRUE)
  rcv1_aggregate(out, state, write = TRUE)
  progress$phase <- "completed"
  progress$emit("all outputs aggregated", force = TRUE)
  perf_manifest(out)
  perf_verify_manifest(out)
  writeLines(state$signature, file.path(out, "COMPLETED"))
  rcv1_validate_complete(out, root, opt$dataset)
  invisible(NULL)
}
