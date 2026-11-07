# Versioned UCI QSAR aquatic-toxicity comparison; existing solvers are unchanged.
qsar_v1_methods <- function() c("ols", "lambert_cv", "lasso", "ridge",
  "elastic_net", "mcp", "scad")

qsar_v1_config <- function(root, mode) {
  stopifnot(mode %in% c("smoke", "production"))
  base <- perf_config(root, "confirmation")
  list(schema = "qsar_aquatic_realdata_v1", mode = mode, dataset_id = 505L,
    raw_dir = "data/qsar",
    archive_sha256 = "09e45d1b2b66c5bc871a64b87ad40e452623e1c2dc375e517e4769379a9db9cd",
    csv_sha256 = "cee55fadc4b021d3ff1567612a7ac32d7e212845a4ca9eb3ef983b78e0c17856",
    n = 546L, p = 8L, outer_folds = 5L,
    outer_units = if (mode == "smoke") 2L else 5L,
    outer_seed = 2026092301L, inner_seeds = 2026092302L + seq_len(5L),
    folds = 5L, lambda_points = if (mode == "smoke") 4L else 60L,
    elastic_net_alpha = base$elastic_net_alpha, mcp_gamma = base$mcp_gamma,
    scad_a = base$scad_a, lambert_shapes = base$lambert_shapes,
    methods = qsar_v1_methods()[-1L],
    max_coordinate_sweeps_per_start = base$max_coordinate_sweeps_per_start,
    normalized_kkt_tolerance = base$normalized_kkt_tolerance,
    support_reporting_tolerance = base$support_reporting_tolerance,
    implementation = base$implementation,
    response_transformation = "none", predictor_scaling = "training_fold_center_rms")
}

qsar_v1_source_paths <- function() c("R/penalties.R", "R/lambert_family.R",
  "R/performance_models.R", "R/performance_run.R",
  "R/qsar_aquatic_realdata_v1.R", "src/lambert_performance.cpp",
  "scripts/prepare_qsar_aquatic_v1.R", "scripts/qsar_aquatic_realdata_v1.R",
  "protocol/performance_v1.md", "protocol/qsar_aquatic_realdata_v1.md")

qsar_v1_sources <- function(root) {
  files <- qsar_v1_source_paths()
  paths <- file.path(root, files)
  data.frame(file = files, bytes = as.numeric(file.info(paths)$size),
    sha256 = unname(vapply(paths, perf_hash, character(1))))
}

qsar_v1_data <- function(root, cfg) {
  dir <- file.path(root, cfg$raw_dir)
  archive <- file.path(dir, "qsar_aquatic_toxicity.zip")
  csv <- file.path(dir, "qsar_aquatic_toxicity.csv")
  receipt_path <- file.path(dir, "ACQUISITION_RECEIPT.json")
  manifest_path <- file.path(dir, "manifest.csv")
  if (!all(file.exists(c(archive, csv, receipt_path, manifest_path))) ||
      perf_hash(archive) != cfg$archive_sha256 || perf_hash(csv) != cfg$csv_sha256)
    stop("Frozen UCI QSAR raw files are missing or changed")
  receipt <- jsonlite::fromJSON(receipt_path)
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  expected_files <- c(basename(archive), basename(csv))
  if (receipt$schema != "uci_qsar_aquatic_acquisition_v1" ||
      receipt$archive_sha256 != cfg$archive_sha256 ||
      receipt$csv_sha256 != cfg$csv_sha256 ||
      !identical(manifest$file, expected_files) ||
      !identical(manifest$sha256, c(cfg$archive_sha256, cfg$csv_sha256)) ||
      !identical(as.numeric(manifest$bytes), as.numeric(file.info(c(archive, csv))$size)))
    stop("UCI QSAR provenance mismatch")
  raw <- utils::read.csv(csv, header = FALSE, sep = ";",
    colClasses = rep("numeric", cfg$p + 1L))
  if (!identical(dim(raw), c(cfg$n, cfg$p + 1L)) ||
      any(!is.finite(as.matrix(raw)))) stop("UCI QSAR numeric schema changed")
  X <- as.matrix(raw[, seq_len(cfg$p), drop = FALSE]); y <- raw[[cfg$p + 1L]]
  if (any(apply(X, 2L, stats::sd) <= 0) ||
      qr(cbind(1, X))$rank != cfg$p + 1L) stop("UCI QSAR predictor rank/variance changed")
  colnames(X) <- paste0("x", seq_len(cfg$p))
  list(X = X, y = y, sample_id = seq_len(cfg$n))
}

qsar_v1_folds <- function(cfg) {
  RNGkind("Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  set.seed(cfg$outer_seed)
  outer <- sample(rep(seq_len(cfg$outer_folds), length.out = cfg$n))
  inner <- lapply(seq_len(cfg$outer_folds), function(k) {
    set.seed(cfg$inner_seeds[[k]])
    sample(rep(seq_len(cfg$folds), length.out = sum(outer != k)))
  })
  if (max(table(outer)) - min(table(outer)) > 1L ||
      any(vapply(inner, function(x) max(table(x)) - min(table(x)) > 1L, logical(1))))
    stop("Frozen CV partitions are not balanced")
  list(outer = outer, inner = inner)
}

qsar_v1_signature <- function(out) {
  files <- c("config.json", "sources.csv", "raw_inputs.csv", "folds.rds")
  digest::digest(paste(files, vapply(file.path(out, files), perf_hash, character(1)),
    sep = ":", collapse = "\n"), algo = "sha256", serialize = FALSE)
}

qsar_v1_metadata <- function(out, root, cfg, folds) {
  saveRDS(cfg, file.path(out, "config.rds"), version = 3)
  writeLines(perf_json(cfg), file.path(out, "config.json"))
  utils::write.csv(qsar_v1_sources(root), file.path(out, "sources.csv"), row.names = FALSE)
  raw_names <- c("qsar_aquatic_toxicity.zip", "qsar_aquatic_toxicity.csv",
    "ACQUISITION_RECEIPT.json", "manifest.csv")
  raw_paths <- file.path(root, cfg$raw_dir, raw_names)
  utils::write.csv(data.frame(file = raw_names,
    bytes = as.numeric(file.info(raw_paths)$size),
    sha256 = unname(vapply(raw_paths, perf_hash, character(1)))),
    file.path(out, "raw_inputs.csv"), row.names = FALSE)
  saveRDS(folds, file.path(out, "folds.rds"), version = 3)
  utils::write.csv(data.frame(sample_id = seq_len(cfg$n), outer_fold = folds$outer),
    file.path(out, "outer_folds.csv"), row.names = FALSE)
  writeLines(perf_json(perf_software()), file.path(out, "software.json"))
  capture.output(sessionInfo(), .libPaths(), file = file.path(out, "session_info.txt"))
  writeLines(qsar_v1_signature(out), file.path(out, "signature.sha256"))
  files <- c("config.rds", "config.json", "sources.csv", "raw_inputs.csv",
    "folds.rds", "outer_folds.csv", "software.json", "session_info.txt", "signature.sha256")
  utils::write.csv(data.frame(file = files,
    sha256 = unname(vapply(file.path(out, files), perf_hash, character(1)))),
    file.path(out, "metadata_manifest.csv"), row.names = FALSE)
}

qsar_v1_verify_metadata <- function(out, root) {
  meta <- utils::read.csv(file.path(out, "metadata_manifest.csv"), stringsAsFactors = FALSE)
  expected <- c("config.rds", "config.json", "sources.csv", "raw_inputs.csv",
    "folds.rds", "outer_folds.csv", "software.json", "session_info.txt", "signature.sha256")
  if (!identical(meta$file, expected) || anyDuplicated(meta$file) ||
      !identical(meta$sha256, unname(vapply(file.path(out, meta$file), perf_hash, character(1)))))
    stop("Metadata manifest mismatch")
  cfg <- readRDS(file.path(out, "config.rds"))
  if (!identical(cfg, qsar_v1_config(root, cfg$mode)) ||
      !identical(paste(readLines(file.path(out, "config.json")), collapse = "\n"), perf_json(cfg)))
    stop("Configuration mismatch")
  sources <- utils::read.csv(file.path(out, "sources.csv"), stringsAsFactors = FALSE)
  if (!isTRUE(all.equal(sources, qsar_v1_sources(root), check.attributes = FALSE)))
    stop("Source-code checksum mismatch")
  data <- qsar_v1_data(root, cfg)
  raw <- utils::read.csv(file.path(out, "raw_inputs.csv"), stringsAsFactors = FALSE)
  raw_paths <- file.path(root, cfg$raw_dir, raw$file)
  if (!identical(raw$file, c("qsar_aquatic_toxicity.zip", "qsar_aquatic_toxicity.csv",
      "ACQUISITION_RECEIPT.json", "manifest.csv")) ||
      !identical(as.numeric(raw$bytes), as.numeric(file.info(raw_paths)$size)) ||
      !identical(raw$sha256, unname(vapply(raw_paths, perf_hash, character(1)))))
    stop("Raw-input checksum mismatch")
  folds <- readRDS(file.path(out, "folds.rds"))
  outer_csv <- utils::read.csv(file.path(out, "outer_folds.csv"))
  if (!identical(folds, qsar_v1_folds(cfg)) ||
      !identical(outer_csv$sample_id, seq_len(cfg$n)) ||
      !identical(outer_csv$outer_fold, folds$outer)) stop("Frozen folds mismatch")
  if (!identical(paste(readLines(file.path(out, "software.json")), collapse = "\n"),
      perf_json(perf_software())) ||
      !identical(readLines(file.path(out, "signature.sha256")), qsar_v1_signature(out)))
    stop("Software or scientific signature mismatch")
  list(cfg = cfg, folds = folds, data = data,
    signature = readLines(file.path(out, "signature.sha256")))
}

qsar_v1_unit_path <- function(out, unit) file.path(out, "shards", sprintf("unit_%05d", unit))

qsar_v1_fit_unit <- function(unit, state, pulse) {
  cfg <- state$cfg; data <- state$data; folds <- state$folds
  train <- which(folds$outer != unit); test <- which(folds$outer == unit)
  started <- proc.time()[[3L]]
  fit <- perf_train_models(data$X[train, , drop = FALSE], data$y[train],
    folds$inner[[unit]], cfg, pulse)
  design <- cbind(intercept = 1, data$X[train, , drop = FALSE])
  ols <- stats::lm.fit(design, data$y[train])
  ols_ok <- ols$rank == ncol(design) && all(is.finite(ols$coefficients))
  rows <- predictions <- vector("list", length(qsar_v1_methods()))
  for (j in seq_along(qsar_v1_methods())) {
    method <- qsar_v1_methods()[[j]]
    model <- if (method == "ols") NULL else fit$models[[method]]
    ok <- if (method == "ols") ols_ok else isTRUE(model$ok)
    status <- if (method == "ols") if (ok) "full_rank" else "rank_deficient" else model$status
    pred <- rep(NA_real_, length(test)); beta <- rep(NA_real_, cfg$p)
    intercept <- NA_real_
    if (ok) {
      if (method == "ols") {
        intercept <- unname(ols$coefficients[[1L]])
        beta <- unname(ols$coefficients[-1L])
      } else {
        intercept <- model$intercept; beta <- model$beta
      }
      pred <- intercept + drop(data$X[test, , drop = FALSE] %*% beta)
      if (any(!is.finite(pred))) { ok <- FALSE; status <- "nonfinite_prediction"; pred[] <- NA_real_ }
    }
    selection <- if (method == "ols") NULL else model$selection
    denom <- sum((data$y[test] - mean(data$y[test]))^2)
    rows[[j]] <- data.frame(unit = unit, method = method, n_train = length(train),
      n_test = length(test), ok = ok, status = status,
      partial_search = if (method == "ols") FALSE else isTRUE(selection$partial),
      family = if (method == "ols") "ols" else if (isTRUE(selection$ok))
        as.character(model$setting$family) else NA_character_,
      shape = if (method != "ols" && isTRUE(selection$ok))
        model$setting$shape else NA_real_,
      lambda = if (method != "ols" && isTRUE(selection$ok))
        model$diagnostics$lambda else NA_real_,
      cv_mse = if (method != "ols" && isTRUE(selection$ok))
        selection$loss else NA_real_,
      kkt = if (method != "ols" && isTRUE(selection$ok))
        model$diagnostics$kkt else NA_real_,
      selected_size = if (ok) sum(abs(beta) > cfg$support_reporting_tolerance) else NA_integer_,
      rmse = if (ok) sqrt(mean((data$y[test] - pred)^2)) else NA_real_,
      mae = if (ok) mean(abs(data$y[test] - pred)) else NA_real_,
      r_squared = if (ok && denom > 0) 1 - sum((data$y[test] - pred)^2) / denom else NA_real_)
    predictions[[j]] <- data.frame(unit = unit, method = method,
      sample_id = test, observed = data$y[test], predicted = pred,
      status = status, stringsAsFactors = FALSE)
  }
  list(unit = unit, signature = state$signature, inner_folds = folds$inner[[unit]],
    test_ids = test, rows = do.call(rbind, rows),
    predictions = do.call(rbind, predictions), fit = fit,
    ols = list(rank = ols$rank, coefficients = ols$coefficients),
    elapsed_seconds = proc.time()[[3L]] - started)
}

qsar_v1_validate_unit <- function(result, unit, state) {
  cfg <- state$cfg; test <- which(state$folds$outer == unit)
  rows <- result$rows; pred <- result$predictions; methods <- qsar_v1_methods()
  if (!identical(result$unit, as.integer(unit)) ||
      !identical(result$signature, state$signature) ||
      !identical(result$inner_folds, state$folds$inner[[unit]]) ||
      !identical(result$test_ids, test) ||
      !identical(rows$method, methods) || any(rows$unit != unit) ||
      nrow(rows) != length(methods) || nrow(pred) != length(methods) * length(test) ||
      !identical(pred$method, rep(methods, each = length(test))) ||
      !identical(pred$sample_id, rep(test, times = length(methods))) ||
      !isTRUE(all.equal(pred$observed, rep(state$data$y[test], times = length(methods)))))
    stop("Shard identity, fold, or prediction coverage mismatch: ", unit)
  if (anyNA(rows$ok) || anyNA(rows$status) ||
      !identical(pred$status, rep(rows$status, each = length(test))) ||
      any(!is.finite(pred$predicted[rep(rows$ok, each = length(test))])) ||
      any(!is.na(pred$predicted[rep(!rows$ok, each = length(test))])) ||
      any(!is.finite(as.matrix(rows[rows$ok, c("rmse", "mae", "r_squared"), drop = FALSE]))) ||
      any(!is.na(rows$rmse[!rows$ok])) ||
      any(rows$kkt[rows$ok & rows$method != "ols"] > cfg$normalized_kkt_tolerance))
    stop("Shard numerical/status gate failed: ", unit)
  for (method in cfg$methods) {
    row <- rows[rows$method == method, , drop = FALSE]
    model <- result$fit$models[[method]]
    if (!identical(row$ok[[1L]], isTRUE(model$ok) && row$status[[1L]] != "nonfinite_prediction") ||
        !identical(row$partial_search[[1L]], isTRUE(model$selection$partial)))
      stop("Model-selection audit mismatch: ", unit, "/", method)
  }
  invisible(TRUE)
}

qsar_v1_read_unit <- function(out, unit, state) {
  path <- qsar_v1_unit_path(out, unit)
  file <- file.path(path, "result.rds")
  if (!file.exists(file) || !file.exists(file.path(path, "result.sha256")) ||
      !identical(readLines(file.path(path, "result.sha256")), perf_hash(file)))
    stop("Missing or invalid shard checksum: ", unit)
  result <- readRDS(file)
  qsar_v1_validate_unit(result, unit, state)
  result
}

qsar_v1_state <- function(out, root) {
  state <- qsar_v1_verify_metadata(out, root)
  paths <- list.dirs(file.path(out, "shards"), recursive = FALSE, full.names = FALSE)
  expected <- sprintf("unit_%05d", seq_len(state$cfg$outer_units))
  if (!all(paths %in% expected)) stop("Unexpected or partial shard directory")
  completed <- which(expected %in% paths)
  state$completed <- completed
  state$units <- setNames(lapply(completed, function(i) qsar_v1_read_unit(out, i, state)),
    as.character(completed))
  state
}

qsar_v1_write_unit <- function(out, result, state) {
  unit <- result$unit; qsar_v1_validate_unit(result, unit, state)
  path <- qsar_v1_unit_path(out, unit)
  if (file.exists(path)) stop("Refusing to overwrite committed shard")
  temp <- paste0(path, ".partial-", Sys.getpid())
  if (!dir.create(temp)) stop("Cannot stage shard")
  saveRDS(result, file.path(temp, "result.rds"), version = 3, compress = "gzip")
  writeLines(perf_hash(file.path(temp, "result.rds")), file.path(temp, "result.sha256"))
  if (!file.rename(temp, path)) stop("Cannot commit atomic shard")
}

qsar_v1_progress <- function(out, state) {
  e <- new.env(parent = emptyenv())
  e$started <- proc.time()[[3L]]; e$last <- -Inf; e$durations <- numeric()
  e$completed <- state$completed; e$unit <- NA_integer_; e$phase <- "preflight"
  e$emit <- function(detail = "", force = FALSE) {
    now <- proc.time()[[3L]]
    if (!force && now - e$last < 30) return(invisible(NULL))
    e$last <- now; n <- state$cfg$outer_units; done <- length(e$completed)
    eta <- if (done == n) 0 else if (length(e$durations))
      mean(e$durations) * (n - done) else NA_real_
    stamp <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    snapshot <- list(run_id = basename(out), slurm_job_id = NA_character_,
      phase = e$phase, work_unit = "one outer fold, all methods and inner CV",
      total = n, validated_completed = done,
      running = as.integer(!is.na(e$unit)), failed = 0L,
      pending = n - done - as.integer(!is.na(e$unit)),
      percent = 100 * done / n, elapsed_seconds = now - e$started,
      throughput_units_per_second = if (length(e$durations))
        length(e$durations) / (now - e$started) else NA_real_,
      estimated_remaining_seconds = eta,
      estimated_completion_utc = if (is.finite(eta))
        format(Sys.time() + eta, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC") else NA_character_,
      eta_scope = "fold computation only; final aggregation/validation excluded",
      slurm_remaining_seconds = NA_real_, last_heartbeat_utc = stamp,
      current_work_unit = e$unit,
      most_recently_completed_unit = if (done) tail(e$completed, 1L) else NA_integer_,
      detail = detail)
    perf_atomic_json(snapshot, file.path(out, "progress.json"))
    history <- file.path(out, "progress.tsv")
    utils::write.table(as.data.frame(snapshot), history, sep = "\t", row.names = FALSE,
      col.names = !file.exists(history), append = file.exists(history), quote = TRUE, na = "NA")
    cat(sprintf("[%s] %s | %d/%d (%.1f%%) | fold %s | ETA %s | %s\n",
      stamp, e$phase, done, n, snapshot$percent, e$unit,
      if (is.finite(eta)) paste(round(eta), "s, computation only") else "estimating",
      detail))
    flush(stdout())
  }
  e
}

qsar_v1_save_csv <- function(value, path) {
  temp <- paste0(path, ".tmp-", Sys.getpid())
  if (file.exists(temp)) stop("Existing aggregate staging file")
  utils::write.csv(value, temp, row.names = FALSE)
  if (file.exists(path)) {
    if (perf_hash(path) != perf_hash(temp)) stop("Existing aggregate differs: ", path)
    unlink(temp)
  } else if (!file.rename(temp, path)) stop("Cannot commit aggregate: ", path)
}

qsar_v1_aggregate <- function(out, state) {
  if (length(state$completed) != state$cfg$outer_units) stop("Incomplete folds")
  units <- lapply(seq_len(state$cfg$outer_units), function(i) qsar_v1_read_unit(out, i, state))
  rows <- do.call(rbind, lapply(units, `[[`, "rows"))
  pred <- do.call(rbind, lapply(units, `[[`, "predictions"))
  methods <- qsar_v1_methods()
  summary <- do.call(rbind, lapply(methods, function(method) {
    r <- rows[rows$method == method, , drop = FALSE]
    p <- pred[pred$method == method, , drop = FALSE]
    if (anyDuplicated(p$sample_id) ||
        !setequal(p$sample_id, which(state$folds$outer <= state$cfg$outer_units)))
      stop("Outer prediction coverage mismatch: ", method)
    valid <- all(r$ok) && all(is.finite(p$predicted))
    data.frame(method = method, folds = nrow(r), failed_folds = sum(!r$ok),
      partial_search_folds = sum(r$partial_search), n_test = nrow(p),
      pooled_rmse = if (valid) sqrt(mean((p$observed - p$predicted)^2)) else NA_real_,
      pooled_mae = if (valid) mean(abs(p$observed - p$predicted)) else NA_real_,
      pooled_r_squared = if (valid) 1 - sum((p$observed - p$predicted)^2) /
        sum((p$observed - mean(p$observed))^2) else NA_real_)
  }))
  qsar_v1_save_csv(rows, file.path(out, "results.csv"))
  qsar_v1_save_csv(pred, file.path(out, "predictions.csv"))
  qsar_v1_save_csv(summary, file.path(out, "summary.csv"))
  list(rows = rows, predictions = pred, summary = summary)
}

qsar_v1_validate_complete <- function(out, root) {
  state <- qsar_v1_state(out, root)
  if (length(state$completed) != state$cfg$outer_units ||
      !file.exists(file.path(out, "COMPLETED")) ||
      !identical(readLines(file.path(out, "COMPLETED")), state$signature))
    stop("Run is not complete")
  perf_verify_manifest(out)
  aggregate <- qsar_v1_aggregate(out, state)
  progress <- jsonlite::fromJSON(file.path(out, "progress.json"))
  if (progress$phase != "completed" || progress$validated_completed != state$cfg$outer_units ||
      progress$running != 0L || progress$failed != 0L || progress$pending != 0L)
    stop("Final progress mismatch")
  cat(sprintf("COMPLETE: %d outer folds, %d method rows, %d predictions; raw/source/seed/shard/manifest gates verified.\n",
    state$cfg$outer_units, nrow(aggregate$rows), nrow(aggregate$predictions)))
  print(aggregate$summary)
  invisible(state)
}

qsar_v1_main <- function(root, args) {
  opt <- list(mode = NULL, out = NULL, resume = FALSE, stop_after = Inf)
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key == "--resume") { opt$resume <- TRUE; i <- i + 1L; next }
    if (!key %in% c("--mode", "--out", "--stop-after") || i == length(args) ||
        startsWith(args[[i + 1L]], "--")) stop("Invalid CLI arguments")
    opt[[gsub("-", "_", substring(key, 3L))]] <- args[[i + 1L]]
    i <- i + 2L
  }
  if (is.null(opt$out) || is.null(opt$mode) ||
      !opt$mode %in% c("smoke", "production", "validate"))
    stop("Use --mode smoke|production|validate --out DIRECTORY [--resume]")
  out <- if (startsWith(opt$out, "/")) opt$out else file.path(root, opt$out)
  if (opt$mode == "validate") {
    if (opt$resume || !is.infinite(opt$stop_after)) stop("Validate takes only --out")
    qsar_v1_validate_complete(out, root); return(invisible(NULL))
  }
  stop_after <- suppressWarnings(as.numeric(opt$stop_after))
  if (is.na(stop_after) || stop_after < 1 ||
      (!is.infinite(stop_after) && opt$mode != "smoke"))
    stop("--stop-after is available only for smoke interruption testing")
  if (file.exists(out) && !opt$resume) stop("Output exists; refusing overwrite")
  if (!file.exists(out) && opt$resume) stop("Resume output is missing")
  perf_preflight(compile = TRUE, root = root)
  cfg <- qsar_v1_config(root, opt$mode)
  data <- qsar_v1_data(root, cfg)
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  if (!opt$resume) {
    if (!dir.create(out) || !dir.create(file.path(out, "shards")))
      stop("Cannot create fresh output")
  }
  lock <- perf_lock(out)
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  if (opt$resume) {
    state <- qsar_v1_state(out, root)
    if (state$cfg$mode != opt$mode) stop("Resume mode mismatch")
    if (file.exists(file.path(out, "COMPLETED"))) {
      qsar_v1_validate_complete(out, root); return(invisible(NULL))
    }
  } else {
    qsar_v1_metadata(out, root, cfg, qsar_v1_folds(cfg))
    state <- qsar_v1_state(out, root)
  }
  progress <- qsar_v1_progress(out, state)
  progress$phase <- "computation"; progress$emit("validated run started", force = TRUE)
  for (unit in setdiff(seq_len(cfg$outer_units), state$completed)) {
    progress$unit <- unit; progress$emit("outer fold started", force = TRUE)
    result <- qsar_v1_fit_unit(unit, state, function(detail) progress$emit(detail))
    qsar_v1_write_unit(out, result, state)
    state$completed <- c(state$completed, unit)
    progress$completed <- state$completed; progress$durations <- c(progress$durations, result$elapsed_seconds)
    progress$unit <- NA_integer_; progress$emit("outer fold validated", force = TRUE)
    if (length(progress$durations) >= stop_after) {
      progress$phase <- "interrupted_smoke"; progress$emit("smoke checkpoint ready", force = TRUE)
      return(invisible(NULL))
    }
  }
  progress$phase <- "aggregation_and_validation"; progress$emit("reading validated shards", force = TRUE)
  qsar_v1_aggregate(out, state)
  progress$phase <- "completed"; progress$emit("all outputs aggregated", force = TRUE)
  perf_manifest(out); perf_verify_manifest(out)
  writeLines(state$signature, file.path(out, "COMPLETED"))
  qsar_v1_validate_complete(out, root)
}
