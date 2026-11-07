#!/usr/bin/env Rscript
# Same split seeds, training-role preprocessing, path kernels, and selector as
# the manuscript. Each completed split is written as one checkpoint.
run <- function(args) {
  if (!length(args) || !args[1] %in% c("design", "smoke", "start", "resume"))
    stop("Usage: Rscript --vanilla run.R design | smoke/start/resume OUTPUT")
  mode <- args[1]
  root <- normalizePath(".")
  helpers <- c("checkpoint_packs_v1", "byrm_data_v1", "riboflavin_fixed_hd_v2",
    "lambert_bic_comparison_v1", "spectral_ic_rescue_solver_v1", "scheetz_cv_ic_v1",
    "scheetz_cv_ic_run_v1", "scheetz_cv_ic_cache_v1_a02", "scheetz_stream_collect_v1",
    "scheetz_cv_ic_main_v2", "magic_maize_preprocess_v1")
  prepare <- sch_prep
  for (f in helpers) sys.source(file.path("R", paste0(f, ".R")), envir = globalenv())
  assign("sch_prep", prepare, envir = globalenv())
  data <- mmp_source(root)
  stopifnot(identical(dim(data$X), c(529L, 54234L)), all(is.finite(data$y)))
  if (mode == "design") {
    cat("MAGIC maize: 529 lines, 54234 input markers, 50 splits of 423/106;",
        "5 inner folds and 8 path settings per training role; model fits: 0.\n")
    return(invisible(NULL))
  }
  if (length(args) != 2L) stop("Specify the output directory")
  out <- args[2]
  if (mode != "resume") {
    dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
    if (!dir.create(out)) stop("Output exists")
  } else if (!dir.exists(out)) stop("Missing output to resume")
  lock <- file.path(out, ".lock")
  if (!dir.create(lock)) stop("Output is locked")
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)

  e <- new.env(parent = globalenv())
  engines <- c("penalties", "lambert_family", "performance_models", "performance_acceleration_v1",
               "lambert_ic_pilot_v1", "lambert_hbic_ebic_reanalysis_v1")
  for (f in engines) sys.source(file.path("R", paste0(f, ".R")), envir = e)
  e$data <- data
  smoke <- mode == "smoke"
  e$view <- if (smoke) list(rows = seq_len(40L),
    columns = unique(round(seq(1, ncol(data$X), length.out = 128L)))) else
    list(rows = seq_len(nrow(data$X)), columns = seq_len(ncol(data$X)))
  n <- length(e$view$rows)
  e$repeats <- if (smoke) 1L else seq_len(50L)
  e$active_outer <- 1L
  e$folds <- lapply(e$repeats, function(r) {
    train <- sort(head(sch_sample(2026101801L + 1000L*r, seq_len(n)), floor(.8*n)))
    list(outer = ifelse(seq_len(n) %in% train, 0L, 1L),
         inner = list(sch_sample(2026101802L + 1000L*r, rep(1:5, length.out = length(train)))))
  })
  e$cfg <- e$perf_config(root, "confirmation")
  e$cfg$protocol <- "magic_maize_holdout_v1"
  e$cfg$methods <- c("lambert_c1", "lasso", "ridge", "elastic_net", "mcp", "scad")
  e$cfg$elastic_net_alpha <- c(0,.25,.5,.75,1)
  e$cfg$lambert_shapes <- 1
  e$cfg$mcp_gamma <- 3
  e$cfg$scad_a <- 3.7
  e$cfg$lambda_points <- if (smoke) 4L else 60L
  e$cfg$early_stop <- FALSE
  e$cfg$outer_repeats <- length(e$repeats)
  e$cfg$n_train <- floor(.8*n)
  e$cfg$n_test <- n - e$cfg$n_train
  e$cfg$p <- length(e$view$columns)
  e$rescue <- list(chunk_sweeps = 250L, total_sweeps = 10000L,
                   newton_steps = 20L, backtracks = 20L, armijo = 1e-4)
  e$settings <- data.frame(family = c("lambert", rep("enet",5), "mcp","scad"),
                           shape = c(1,1,0,.25,.5,.75,3,3.7))
  tasks <- sch_bind(lapply(e$repeats, function(r) data.frame(repeat_id = as.integer(r),
    outer_fold = 1L, role = rep(0:5, each=8L), setting = rep(1:8,6L))))
  e$tasks <- cbind(unit = seq_len(nrow(tasks)), tasks)
  predict <- e$perf_predict
  e$perf_predict <- function(beta, prep, X)
    predict(beta, prep, mmp_transform(prep$maize_context, X, FALSE))
  e$reference_path <- e$perf_path
  e$perf_path <- function(prep, setting, cfg, pulse = function(...) NULL, last = cfg$lambda_points) {
    if (setting$family == "enet" && setting$shape == 0)
      rfh_ridge_path(prep, setting, cfg, e, last)
    else e$pacc_path(prep, setting, cfg, pulse, last, e)
  }
  rfh_compile(e)
  paths <- c("run.R", list.files(c("R","protocol","data"), recursive=TRUE, full.names=TRUE))
  bind <- list(sources = vapply(paths, cpk_hash, ""),
    software = rfh_software()$packages, platform = R.version.string,
    cfg = e$cfg, tasks = e$tasks, folds = e$folds)
  binding <- file.path(out, "design.rds")
  if (mode == "resume") stopifnot(identical(readRDS(binding), bind))
  else saveRDS(bind, binding)
  write_atomic <- function(path, writer) {
    tmp <- tempfile(tmpdir = out)
    writer(tmp)
    if (!file.rename(tmp, path)) stop("Could not commit ", path)
  }
  pieces <- character()
  for (r in e$repeats) {
    path <- file.path(out, sprintf("split-%02d.rds", r))
    pieces <- c(pieces, path)
    if (file.exists(path)) {
      saved <- readRDS(path)
      stopifnot(saved$split == r, identical(saved$binding, cpk_object_hash(bind)),
                nrow(saved$results) == 27L, all(saved$results$repeat_id == r))
      next
    }
    e$memo <- new.env(parent = emptyenv())
    ids <- which(e$tasks$repeat_id == r)
    records <- setNames(lapply(ids, function(i) {
      fit <- sch_unit(i, e)
      sch_check_unit(fit, i, e)
      fit
    }), as.character(ids))
    part <- list2env(as.list(e), parent = parent.env(e))
    part$repeats <- r
    aggregate <- schm_collect_parts$group(records, part)
    stopifnot(nrow(aggregate$results) == 27L)
    aggregate$split <- r
    aggregate$binding <- cpk_object_hash(bind)
    write_atomic(path, function(p) saveRDS(aggregate, p))
    cat(sprintf("Saved %d/%d splits, 48 paths per split.\n", r, length(e$repeats)))
    flush(stdout())
  }
  completed <- lapply(pieces, readRDS)
  for (name in c("results", "predictions", "path_summary", "candidate_diagnostics"))
    write_atomic(file.path(out, paste0(name,".csv")), function(p)
      write.csv(sch_bind(lapply(completed, `[[`, name)), p, row.names = FALSE))
  rows <- sch_bind(lapply(completed, function(x) x$results))
  summary <- sch_bind(lapply(split(rows, rows[c("method", "criterion")], drop = TRUE), function(z) {
    stopifnot(nrow(z) == length(e$repeats))
    ans <- data.frame(method = z$method[1], criterion = z$criterion[1],
      complete = sum(z$ok), total = nrow(z), partial_searches = sum(z$partial_search))
    for (metric in c("rmse", "mae", "selected_size")) {
      ans[[paste0("mean_",metric)]] <- if (all(z$ok)) mean(z[[metric]]) else NA_real_
      ans[[paste0("sd_",metric)]] <- if (all(z$ok)) sd(z[[metric]]) else NA_real_
    }
    ans
  }))
  write_atomic(file.path(out, "summary.csv"), function(p) write.csv(summary, p, row.names = FALSE))
  manifest <- data.frame(file = list.files(out), stringsAsFactors = FALSE)
  manifest <- manifest[!manifest$file %in% c("COMPLETED", "manifest.csv"), , drop = FALSE]
  manifest$sha256 <- vapply(file.path(out, manifest$file), cpk_hash, "")
  write_atomic(file.path(out, "manifest.csv"), function(p) write.csv(manifest, p, row.names = FALSE))
  write_atomic(file.path(out, "COMPLETED"), function(p) writeLines(cpk_hash(binding), p))
}
sch_prep <- function(task, e) {
  key <- as.character(task$role)
  if (exists(key, e$memo, inherits = FALSE)) return(get(key, e$memo))
  s <- sch_split(task, e)
  X <- e$data$X[s$train, s$columns, drop = FALSE]
  context <- mmp_fit(X, e$data$annotation)
  prep <- e$perf_prepare(mmp_transform(context, X, FALSE), e$data$y[s$train], e$cfg)
  stopifnot(max(abs(prep$mx - context$center)) < 1e-12,
            max(abs(prep$sx - context$scale)) < 1e-12)
  prep$maize_context <- context
  assign(key, prep, e$memo)
  prep
}
threads <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")
do.call(Sys.setenv, setNames(as.list(rep("1", length(threads))), threads))
run(commandArgs(TRUE))
