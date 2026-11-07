# Isolated adapter: existing data, solvers and persistence stay unchanged.
lf1_labels <- function() c(ols = "OLS", lambert_c1 = "Lambert", lasso = "LASSO",
  ridge = "Ridge", elastic_net = "Elastic net", mcp = "MCP", scad = "SCAD")

lf1_summaries <- function(rows, cfg) {
  stopifnot(!anyDuplicated(rows[c("repeat_id", "outer_fold", "method")]),
    nrow(rows) == cfg$total_units * length(cfg$methods), !anyNA(rows$ok),
    setequal(rows$method, cfg$methods))
  repeats <- do.call(rbind, lapply(seq_len(cfg$repeats), function(r)
    do.call(rbind, lapply(cfg$methods, function(m) {
      z <- rows[rows$repeat_id == r & rows$method == m, ]
      stopifnot(nrow(z) == cfg$active_folds, setequal(z$outer_fold, seq_len(cfg$active_folds)))
      complete <- cfg$mode == "production" && all(z$ok) &&
        all(is.finite(as.matrix(z[c("rmse", "mae", "selected_size")])) )
      data.frame(repeat_id = r, method = m, method_label = unname(lf1_labels()[m]),
        complete_cv = complete, failed_folds = sum(!z$ok),
        partial_search_folds = sum(z$partial_search),
        cv_fold_mean_rmse = if (complete) mean(z$rmse) else NA_real_,
        cv_fold_mean_mae = if (complete) mean(z$mae) else NA_real_,
        cv_fold_mean_size = if (complete) mean(z$selected_size) else NA_real_)
    }))))
  summary <- do.call(rbind, lapply(cfg$methods, function(m) {
    z <- repeats[repeats$method == m, ]; complete <- all(z$complete_cv)
    row <- data.frame(method = m, method_label = unname(lf1_labels()[m]),
      repeats_expected = cfg$repeats, repeats_complete = sum(z$complete_cv),
      all_repeats_complete = complete, failed_folds = sum(z$failed_folds),
      partial_search_folds = sum(z$partial_search_folds))
    for (metric in c("rmse", "mae", "size")) {
      v <- z[[paste0("cv_fold_mean_", metric)]]
      row[[paste0("mean_", metric)]] <- if (complete) mean(v) else NA_real_
      row[[paste0("median_", metric)]] <- if (complete) median(v) else NA_real_
    }
    row
  }))
  rownames(repeats) <- rownames(summary) <- NULL
  list(repeat_metrics = repeats, method_summary = summary)
}

lf1_environment <- function(root) {
  env <- new.env(parent = globalenv())
  sys.source(file.path(root, "R/repeated_lowdim_realdata_v1.R"), envir = env)
  seed_tools <- new.env(parent = env)
  sys.source(file.path(root, "R/repeated_lowdim_realdata_v2.R"), envir = seed_tools)
  old_engine <- env$rcv1_engine; old_config <- env$rcv1_config
  old_sources <- env$rcv1_source_paths; old_aggregate <- env$rcv1_aggregate
  env$rcv1_engine <- function(root, dataset) {
    engine <- old_engine(root, dataset); base_config <- engine$qsar_v1_config
    engine$qsar_v1_methods <- function() names(lf1_labels())
    engine$qsar_v1_config <- function(root, mode) {
      cfg <- base_config(root, mode)
      cfg$schema <- paste0(cfg$schema, "_fixed_c1_v1")
      cfg$methods <- names(lf1_labels())[-1L]; cfg$lambert_shapes <- 1
      cfg$main_estimator <- list(label = "Lambert", method = "lambert_c1",
        shape = 1, tuned_parameter = "lambda", post_threshold = FALSE)
      cfg
    }
    original_validate <- engine$qsar_v1_validate_unit
    engine$qsar_v1_validate_unit <- function(result, unit, state) {
      original_validate(result, unit, state)
      row <- result$rows[result$rows$method == "lambert_c1", ]
      cv <- result$fit$cv
      stopifnot(nrow(row) == 1L, state$cfg$lambert_shapes == 1,
        all(cv$shape[cv$family == "lambert"] == 1),
        !any(c("lambert_cv", "lambert_c1_post_min") %in% result$rows$method))
      if (row$ok) stopifnot(row$shape == 1,
        identical(result$fit$models$lambert_c1$setting$family, "lambert"))
      invisible(TRUE)
    }
    engine
  }
  env$rcv1_config <- function(root, dataset, mode, engine = env$rcv1_engine(root, dataset)) {
    cfg <- old_config(root, dataset, mode, engine)
    cfg$schema <- "lambert_fixed_c1_realdata_v1"
    cfg$repeats <- if (mode == "smoke") 2L else 100L
    cfg$total_units <- as.integer(cfg$repeats * cfg$active_folds)
    seeds <- seed_tools$rcv2_seed_table(cfg$base, seq_len(cfg$repeats))
    cfg$seeds <- seeds[c("repeat_id", "outer_seed", paste0("inner_seed_", 1:5))]
    cfg$primary_summary <- "arithmetic_mean_of_five_fold_rmse_then_mean_and_median_over_100_repetitions"
    cfg$supportive_summaries <- c("fold_mean_mae", "fold_mean_selected_size")
    cfg
  }
  env$rcv1_source_paths <- function(engine) unique(c(old_sources(engine),
    "R/repeated_lowdim_realdata_v2.R", "R/lambert_fixed_c1_realdata_v1.R",
    "scripts/lambert_fixed_c1_realdata_v1.R", "tests/lambert_fixed_c1_v1_checks.R",
    "protocol/lambert_fixed_c1_v1.md"))
  env$rcv1_aggregate <- function(out, state, write = FALSE) {
    value <- old_aggregate(out, state, write = FALSE)
    summary <- lf1_summaries(value$results, state$cfg)
    value$repeat_metrics <- summary$repeat_metrics
    value$method_summary <- summary$method_summary
    if (write) for (name in names(value))
      env$rcv1_save_csv(value[[name]], file.path(out, paste0(name, ".csv")))
    value
  }
  env
}
