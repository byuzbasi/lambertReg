# Version-2 tail of the repeated low-dimensional comparison. This module is
# sourced after v1 in a private environment so its two bindings affect only
# the new run; the completed v1 output retains its original validator.

rcv2_prior_path <- function(root, dataset) {
  if (!dataset %in% c("diabetes", "qsar", "hitters"))
    stop("Unknown repeated-CV dataset")
  file.path(root, "outputs/real_data",
    paste0(dataset, "_repeated_cv_v1_a01"))
}

rcv2_seed_table <- function(base, global_ids) {
  global_ids <- as.integer(global_ids)
  if (!length(global_ids) || anyNA(global_ids) ||
      anyDuplicated(global_ids) || any(global_ids < 1L | global_ids > 100L))
    stop("Invalid global repetition IDs")
  offsets <- 1000 * (as.numeric(global_ids) - 1)
  if (any(base$outer_seed + offsets > .Machine$integer.max) ||
      any(max(base$inner_seeds) + offsets > .Machine$integer.max))
    stop("Repeated-CV seeds exceed the integer range")
  seeds <- data.frame(repeat_id = seq_along(global_ids),
    global_repeat_id = global_ids,
    outer_seed = as.integer(base$outer_seed + offsets))
  for (fold in seq_len(5L))
    seeds[[paste0("inner_seed_", fold)]] <-
      as.integer(base$inner_seeds[[fold]] + offsets)
  seed_values <- unlist(seeds[c("outer_seed", paste0("inner_seed_", 1:5))],
    use.names = FALSE)
  if (anyDuplicated(seed_values) || any(seed_values <= 0L))
    stop("Repeated-CV seed collision")
  seeds
}

rcv2_base_config <- rcv1_config
rcv1_config <- function(root, dataset, mode,
    engine = rcv1_engine(root, dataset)) {
  cfg <- rcv2_base_config(root, dataset, mode, engine)
  global_ids <- if (mode == "smoke") 11:12 else 11:100
  prior <- rcv2_prior_path(root, dataset)
  required <- file.path(prior, c("COMPLETED", "manifest.csv",
    "signature.sha256"))
  if (!all(file.exists(required)))
    stop("Validated first-ten production output is required")
  cfg$schema <- "repeated_lowdim_tail90_v2"
  cfg$repeats <- as.integer(length(global_ids))
  cfg$total_units <- as.integer(cfg$repeats * cfg$active_folds)
  cfg$global_repeat_ids <- as.integer(global_ids)
  cfg$seeds <- rcv2_seed_table(cfg$base, global_ids)
  cfg$prior <- list(path = file.path("outputs/real_data",
      paste0(dataset, "_repeated_cv_v1_a01")),
    signature = readLines(file.path(prior, "signature.sha256")),
    manifest_sha256 = perf_hash(file.path(prior, "manifest.csv")))
  cfg$primary_summary <- "five_outer_fold_rmse_arithmetic_mean; combine_with_prior_for_100"
  cfg$supportive_summaries <- "v1_pooled_oof_columns_are_tail_diagnostics_only"
  cfg
}

rcv2_base_source_paths <- rcv1_source_paths
rcv1_source_paths <- function(engine) unique(c(rcv2_base_source_paths(engine),
  "R/repeated_lowdim_realdata_v2.R",
  "scripts/repeated_lowdim_realdata_v2.R",
  "protocol/repeated_lowdim_realdata_v2.md"))

rcv2_validate_prior_link <- function(root, dataset, prior_state, cfg) {
  prior <- rcv2_prior_path(root, dataset)
  if (!identical(prior_state$cfg$mode, "production") ||
      prior_state$cfg$repeats != 10L ||
      prior_state$cfg$active_folds != 5L ||
      length(prior_state$completed) != 50L ||
      !identical(prior_state$signature, cfg$prior$signature) ||
      !identical(perf_hash(file.path(prior, "manifest.csv")),
        cfg$prior$manifest_sha256) ||
      !identical(as.integer(cfg$global_repeat_ids),
        if (cfg$mode == "smoke") 11:12 else 11:100))
    stop("First-ten run and version-2 tail are not linked")
  expected_old <- rcv2_seed_table(prior_state$cfg$base, 1:10)
  old_columns <- c("repeat_id", "outer_seed", paste0("inner_seed_", 1:5))
  if (!identical(prior_state$cfg$seeds, expected_old[old_columns]) ||
      !identical(cfg$seeds,
        rcv2_seed_table(cfg$base, cfg$global_repeat_ids)))
    stop("First-ten or tail seed schedule changed")
  if (cfg$mode == "production" &&
      !identical(prior_state$cfg$base, cfg$base))
    stop("First-ten and tail model configurations differ")
  invisible(TRUE)
}

rcv2_cv_values <- function(rows, dataset, repeats = 100L) {
  methods <- qsar_v1_methods()
  repeats <- as.integer(repeats)
  expected_n <- repeats * 5L * length(methods)
  if (!identical(nrow(rows), expected_n) ||
      !all(c("repeat_id", "outer_fold", "method", "ok", "rmse",
        "partial_search", "source_run") %in% names(rows)) ||
      anyNA(rows[, c("repeat_id", "outer_fold", "method", "ok",
        "partial_search", "source_run")]) ||
      !setequal(rows$repeat_id, seq_len(repeats)) ||
      !setequal(rows$outer_fold, seq_len(5L)) ||
      !setequal(rows$method, methods) ||
      anyDuplicated(rows[, c("repeat_id", "outer_fold", "method")]))
    stop("One-hundred-CV fold result coverage changed")
  values <- do.call(rbind, lapply(seq_len(repeats), function(repeat_id)
    do.call(rbind, lapply(methods, function(method) {
      fold <- rows[rows$repeat_id == repeat_id & rows$method == method,
        , drop = FALSE]
      if (nrow(fold) != 5L ||
          !identical(sort(fold$outer_fold), seq_len(5L)) ||
          length(unique(fold$source_run)) != 1L)
        stop("One-CV fold result coverage changed")
      complete <- all(fold$ok) && all(is.finite(fold$rmse))
      data.frame(dataset = dataset, method = method, repeat_id = repeat_id,
        source_run = fold$source_run[[1L]], folds_expected = 5L,
        folds_valid = sum(fold$ok & is.finite(fold$rmse)),
        failed_folds = sum(!fold$ok),
        partial_search_folds = sum(fold$partial_search),
        complete_cv = complete,
        cv_fold_mean_rmse = if (complete) mean(fold$rmse) else NA_real_,
        stringsAsFactors = FALSE)
    }))))
  rownames(values) <- NULL
  values
}

rcv2_method_summary <- function(values, repeats = 100L) {
  methods <- qsar_v1_methods()
  do.call(rbind, lapply(methods, function(method) {
    x <- values[values$method == method, , drop = FALSE]
    if (nrow(x) != repeats || anyDuplicated(x$repeat_id))
      stop("One-hundred-CV method coverage changed")
    complete <- all(x$complete_cv)
    z <- if (complete) x$cv_fold_mean_rmse else numeric()
    data.frame(dataset = x$dataset[[1L]], method = method,
      repeats_expected = as.integer(repeats),
      repeats_complete = sum(x$complete_cv),
      all_repeats_complete = complete,
      failed_folds = sum(x$failed_folds),
      partial_search_folds = sum(x$partial_search_folds),
      mean_100_cv_rmse = if (complete) mean(z) else NA_real_,
      median_100_cv_rmse = if (complete) stats::median(z) else NA_real_,
      q25_100_cv_rmse = if (complete)
        unname(stats::quantile(z, .25)) else NA_real_,
      q75_100_cv_rmse = if (complete)
        unname(stats::quantile(z, .75)) else NA_real_,
      stringsAsFactors = FALSE)
  }))
}
