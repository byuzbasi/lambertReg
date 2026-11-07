#!/usr/bin/env Rscript
root <- dirname(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(), value = TRUE)[1L]))))
setwd(root)

files <- c("penalties", "lambert_family", "performance_models", "performance_run",
  "qsar_aquatic_realdata_v1", "lambert_fixed_c1_realdata_v1",
  "lambert_fixed_c1_realdata_v2", "lambert_ic_pilot_v1",
  "lambert_bic_comparison_v1", "lambert_hbic_ebic_reanalysis_v1",
  "riboflavin_fixed_hd_v2", "spectral_ic_rescue_solver_v1",
  "scheetz_cv_ic_v1", "scheetz_cv_ic_run_v1")
for (f in file.path("R", paste0(files, ".R"))) source(f)

qsar_ic_unit <- function(i, e, state, baseline) {
  task <- state$tasks[i, , drop = FALSE]
  partition <- state$folds[[task$repeat_id]]$outer
  train <- which(partition != task$outer_fold)
  test <- which(partition == task$outer_fold)
  X <- state$data$X
  y <- state$data$y
  stopifnot(identical(test, baseline$sample_id[baseline$unit == i]),
    identical(y[test], baseline$observed[baseline$unit == i]),
    ncol(X) <= floor(length(train) / log(length(train))))
  prep <- perf_prepare(X[train, , drop = FALSE], y[train], e$cfg)
  paths <- lapply(seq_len(nrow(e$settings)), function(j) {
    setting <- e$settings[j, , drop = FALSE]
    original <- perf_path(prep, setting, e$cfg)
    path <- sch_rescue(original, prep, setting, e)$path
    candidates <- sch_candidates(path, prep, setting, e)
    valid <- path$diagnostics$ok & rfh_curvature_ok(candidates$curvature_status)
    list(path = path, candidates = candidates, valid = valid,
      rescued = sum(path$diagnostics$ok & !original$diagnostics$ok))
  })
  names(paths) <- paste0("0_", seq_along(paths))

  predictions <- vector("list", nrow(e$arms))
  rows <- lapply(seq_len(nrow(e$arms)), function(a) {
    arm <- e$arms[a, , drop = FALSE]
    choice <- sch_select(paths, arm$method, arm$criterion, e)
    ok <- !is.na(choice$index)
    lambda <- shape <- q <- size <- kkt <- NA_real_
    endpoint <- NA
    pred <- rep(NA_real_, length(test))
    if (ok) {
      z <- paths[[choice$setting]]
      k <- choice$index
      candidate <- z$candidates[k, , drop = FALSE]
      stopifnot(z$valid[k])
      lambda <- candidate$lambda
      shape <- e$settings$shape[choice$setting]
      q <- candidate$q
      size <- candidate$d
      kkt <- candidate$kkt
      endpoint <- k %in% c(1L, e$cfg$lambda_points)
      pred <- perf_predict(z$path$beta[, k], prep, X[test, , drop = FALSE])
      stopifnot(all(is.finite(pred)))
    }
    predictions[[a]] <<- data.frame(task, arm, sample_id = test,
      observed = y[test], predicted = pred, row.names = NULL)
    data.frame(task, arm, ok = ok, partial_search = choice$partial_search,
      setting = choice$setting, index = choice$index, shape = shape,
      lambda = lambda, score = choice$score, q = q, selected_size = size,
      kkt = kkt, path_endpoint = endpoint, n_train = length(train),
      n_test = length(test), rmse = sqrt(mean((y[test] - pred)^2)),
      mae = mean(abs(y[test] - pred)))
  })
  diagnostics <- lapply(seq_along(paths), function(j) {
    z <- paths[[j]]
    data.frame(task, setting = j, e$settings[j, , drop = FALSE],
      candidates = nrow(z$candidates), rejected = sum(!z$valid),
      rescued = z$rescued)
  })
  list(rows = sch_bind(rows), predictions = sch_bind(predictions),
    diagnostics = sch_bind(diagnostics))
}

qsar_ic_main <- function(args) {
  stopifnot(length(args) == 1L, args %in% c("smoke", "run", "resume"))
  threads <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")
  stopifnot(all(Sys.getenv(threads) == "1"))
  source_run <- "runs/cv"
  env <- lf2_environment(root)
  if (args == "smoke") {
    engine <- env$rcv1_engine(root, "qsar")
    cfg <- env$rcv1_config(root, "qsar", "production", engine)
    state <- list(cfg = cfg, data = engine$qsar_v1_data(root, cfg$base),
      tasks = env$rcv1_tasks(cfg), folds = env$rcv1_folds(cfg, engine))
    baseline <- do.call(rbind, lapply(c(1L, 500L), function(i) {
      task <- state$tasks[i, ]
      test <- which(state$folds[[task$repeat_id]]$outer == task$outer_fold)
      data.frame(unit = i, sample_id = test, observed = state$data$y[test])
    }))
  } else {
    state <- env$rcv1_verify_metadata(source_run, root, "qsar")
    stopifnot(state$cfg$total_units == 500L, state$cfg$base$p == 8L,
      identical(readLines(file.path(source_run, "COMPLETED")), state$signature))
    perf_verify_manifest(source_run)
    baseline <- read.csv(file.path(source_run, "predictions.csv"))
    baseline <- baseline[baseline$method == "ols", , drop = FALSE]
    cv <- read.csv(file.path(source_run, "results.csv"))
    stopifnot(nrow(baseline) == 54600L, nrow(cv) == 3500L, all(cv$ok))
  }

  e <- new.env(parent = globalenv())
  e$cfg <- state$cfg$base
  e$settings <- perf_settings(e$cfg)[c(8L, 5L, 1:4, 6:7), , drop = FALSE]
  stopifnot(identical(e$settings$shape, c(1, 1, 0, .25, .5, .75, 3, 3.7)),
    identical(as.character(e$settings$family),
      c("lambert", rep("enet", 5L), "mcp", "scad")))
  e$rescue <- list(chunk_sweeps = 250L, total_sweeps = 10000L,
    newton_steps = 20L, backtracks = 20L, armijo = 1e-4)
  e$arms <- sch_bind(lapply(e$cfg$methods, function(method) data.frame(method,
    criterion = if (method == "ridge") "BIC" else c("BIC", "EBIC_0.5", "EBIC_1", "HBIC"))))
  stopifnot(nrow(e$arms) == 21L)
  Rcpp::sourceCpp("src/lambert_performance.cpp", env = globalenv())
  numerical <- c("perf_capture", "perf_diagnose", "perf_hessian", "perf_penalty",
    "lambert_cd_cpp", "lic_choose", "lic_curvature", "lhx_scores")
  for (nm in numerical) e[[nm]] <- get(nm, envir = globalenv(), inherits = FALSE)

  if (args == "smoke") {
    for (i in c(1L, 500L)) {
      z <- qsar_ic_unit(i, e, state, baseline)
      stopifnot(nrow(z$rows) == 21L, all(z$rows$ok),
        all(z$rows$kkt <= e$cfg$normalized_kkt_tolerance))
    }
    cat("QSAR_IC_SMOKE_PASS: 2 splits, full 60-point grid, 21 selections per split.\n")
    return(invisible(NULL))
  }

  out <- "runs/information-criteria"
  metadata <- list(source_run = source_run, cv_signature = state$signature,
    cv_manifest_sha256 = perf_hash(file.path(source_run, "manifest.csv")),
    sources = perf_hash("scripts/qsar_information_criteria_v1.R"),
    calculation_sources = setNames(vapply(c(file.path("R", paste0(files, ".R")),
      "src/lambert_performance.cpp"), perf_hash, ""), c(files, "cpp")),
    settings = e$settings, rescue = e$rescue, workers = 1L, software = perf_software())
  if (args == "run") {
    if (!dir.create(out)) stop("Output exists. Use resume for a matching incomplete run.")
  }
  lock <- perf_lock(out)
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  if (args == "run") perf_atomic(file.path(out, "metadata.rds"), function(p) saveRDS(metadata, p))
  else stopifnot(!file.exists(file.path(out, "COMPLETED")),
    identical(metadata, readRDS(file.path(out, "metadata.rds"))))
  chunks <- vector("list", 10L)
  started <- proc.time()[[3]]
  for (j in seq_along(chunks)) {
    checkpoint <- file.path(out, sprintf("part_%02d", j))
    pack <- file.path(checkpoint, "records.rds")
    sha <- file.path(checkpoint, "records.sha256")
    units <- seq.int(50L * (j - 1L) + 1L, 50L * j)
    if (dir.exists(checkpoint)) {
      stopifnot(file.exists(pack), file.exists(sha),
        identical(readLines(sha), perf_hash(pack)))
      chunks[[j]] <- readRDS(pack)
    } else {
      records <- lapply(units, qsar_ic_unit, e = e, state = state, baseline = baseline)
      chunks[[j]] <- lapply(names(records[[1L]]), function(k)
        sch_bind(lapply(records, `[[`, k)))
      names(chunks[[j]]) <- names(records[[1L]])
      staging <- paste0(checkpoint, ".partial-", Sys.getpid())
      stopifnot(dir.create(staging))
      payload <- file.path(staging, "records.rds")
      saveRDS(chunks[[j]], payload)
      writeLines(perf_hash(payload), file.path(staging, "records.sha256"))
      stopifnot(file.rename(staging, checkpoint))
    }
    stopifnot(identical(unique(chunks[[j]]$rows$unit), units))
    cat(sprintf("QSAR IC: %d/500 splits committed, %.1f s this invocation.\n",
      tail(units, 1), proc.time()[[3]] - started))
    flush(stdout())
  }
  rows <- sch_bind(lapply(chunks, `[[`, "rows"))
  predictions <- sch_bind(lapply(chunks, `[[`, "predictions"))
  diagnostics <- sch_bind(lapply(chunks, `[[`, "diagnostics"))
  stopifnot(nrow(rows) == 10500L, nrow(predictions) == 1146600L,
    nrow(diagnostics) == 4000L, !anyDuplicated(rows[c("unit", "method", "criterion")]),
    all(rows$kkt[rows$ok] <= e$cfg$normalized_kkt_tolerance))

  cv$criterion <- ifelse(cv$method == "ols", "OLS", "CV")
  fields <- c("repeat_id", "outer_fold", "method", "criterion", "ok",
    "partial_search", "rmse", "mae", "selected_size")
  both <- rbind(cv[fields], rows[fields])
  metrics <- c("rmse", "mae", "selected_size")
  repeats <- sch_bind(lapply(split(both, both[c("method", "criterion", "repeat_id")],
    drop = TRUE), function(z) {
    stopifnot(identical(sort(z$outer_fold), 1:5))
    data.frame(z[1L, c("method", "criterion", "repeat_id")], ok = all(z$ok),
      as.list(colMeans(z[metrics])))
  }))
  summary <- sch_bind(lapply(split(repeats, repeats[c("method", "criterion")],
    drop = TRUE), function(z) {
    stopifnot(nrow(z) == 100L)
    data.frame(z[1L, c("method", "criterion")], complete_repeats = sum(z$ok),
      as.list(setNames(colMeans(z[metrics]), paste0("mean_", metrics))),
      as.list(setNames(vapply(z[metrics], sd, 0.0), paste0("sd_", metrics))))
  }))
  exports <- list(results = rows, predictions = predictions, path_diagnostics = diagnostics,
    repeat_metrics = repeats, summary = summary)
  for (nm in names(exports)) {
    file <- file.path(out, paste0(nm, ".csv"))
    if (file.exists(file)) {
      saved <- read.csv(file)
      stopifnot(isTRUE(all.equal(saved, exports[[nm]], check.attributes = FALSE,
        tolerance = 1e-12)))
    } else perf_atomic(file, function(p) write.csv(exports[[nm]], p, row.names = FALSE))
  }
  perf_manifest(out)
  perf_verify_manifest(out)
  perf_atomic(file.path(out, "COMPLETED"), function(p)
    writeLines(perf_hash(file.path(out, "manifest.csv")), p))
  print(summary, row.names = FALSE)
  cat("QSAR_IC_COMPLETE:", sum(!rows$ok), "failed selections;",
    sum(rows$partial_search), "partial searches.\n")
}

qsar_ic_main(commandArgs(trailingOnly = TRUE))
