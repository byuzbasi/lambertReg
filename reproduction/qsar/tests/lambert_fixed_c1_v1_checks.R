#!/usr/bin/env Rscript
# Design, scalar identity, aggregation and optional completed-smoke checks only.
driver <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L]))
root <- dirname(dirname(driver))
for (path in c("R/penalties.R", "R/lambert_family.R", "R/performance_models.R",
    "R/performance_run.R", "R/qsar_aquatic_realdata_v1.R", "R/diabetes_realdata_v1.R",
    "R/hitters_realdata_v1.R", "R/lambert_fixed_c1_realdata_v1.R")) source(file.path(root, path))
e <- lf1_environment(root); checks <- list()
gate <- function(name, passed) {
  stopifnot(isTRUE(passed)); checks[[length(checks) + 1L]] <<- data.frame(check = name, passed = TRUE)
}
for (dataset in c("diabetes", "qsar", "hitters")) {
  engine <- e$rcv1_engine(root, dataset)
  cfg <- e$rcv1_config(root, dataset, "production", engine)
  gate(paste0(dataset, "_fixed_shape"), identical(cfg$base$lambert_shapes, 1) &&
    identical(cfg$methods, names(lf1_labels())) && cfg$total_units == 500L && cfg$base$lambda_points == 60L)
  old <- readRDS(file.path(root, "outputs/real_data", paste0(dataset, "_repeated_cv_v1_a01"), "config.rds"))
  gate(paste0(dataset, "_same_seed_schedule"), identical(cfg$seeds[1:10, ], old$seeds))
  fields <- c("n", "p", "raw_dir", "folds", "mcp_gamma", "scad_a", "elastic_net_alpha",
    "support_reporting_tolerance", "normalized_kkt_tolerance", "implementation")
  gate(paste0(dataset, "_unchanged_comparators"), identical(cfg$base[fields], old$base[fields]))
  data <- engine$qsar_v1_data(root, cfg$base); folds <- e$rcv1_folds(cfg, engine)
  e$rcv1_validate_design(data, folds, cfg)
  gate(paste0(dataset, "_complete_design"), length(folds) == 100L && nrow(e$rcv1_tasks(cfg)) == 500L &&
    all(file.exists(file.path(root, e$rcv1_source_paths(engine)))))
  configs <- perf_settings(cfg$base)
  cv <- lapply(seq_len(nrow(configs)), function(i) list(loss = c(4, 2, 3, rep(5, 57))))
  chosen <- perf_choose(cv, configs, "lambert_c1", cfg$base)
  gate(paste0(dataset, "_lambda_only_cv"), chosen$ok && chosen$index == 2L &&
    configs$shape[chosen$setting] == 1 && length(chosen$eligible) == 1L)
}
t <- seq(0, 4, length.out = 101)
gate("existing_figure_is_fixed_c1", max(abs(lambert_penalty(t, 1, 1) - penalty_candidate(t))) < 1e-12 &&
  max(abs(lambert_threshold(t, 1, 1) - shrink_candidate(t))) < 1e-12)
gate("no_initial_positive_curvature", all(lambert_curvature(c(.001, .1, 1, 2), 1, 1) < 0))
cfg <- e$rcv1_config(root, "qsar", "production")
fixture <- expand.grid(outer_fold = 1:5, repeat_id = 1:100, method = cfg$methods, stringsAsFactors = FALSE)
fixture$ok <- TRUE; fixture$partial_search <- FALSE
fixture$rmse <- fixture$outer_fold + fixture$repeat_id / 100
fixture$mae <- fixture$rmse / 2; fixture$selected_size <- 3
s <- lf1_summaries(fixture, cfg)
gate("five_fold_arithmetic_mean", abs(s$repeat_metrics$cv_fold_mean_rmse[1] - 3.01) < 1e-12)
gate("hundred_mean_median", all(abs(s$method_summary$mean_rmse - 3.505) < 1e-12) &&
  all(abs(s$method_summary$median_rmse - 3.505) < 1e-12))
i <- which(fixture$method == "lambert_c1" & fixture$repeat_id == 100 & fixture$outer_fold == 5)
fixture$ok[i] <- FALSE; fixture$rmse[i] <- NA_real_
s <- lf1_summaries(fixture, cfg)$method_summary
gate("incomplete_summary_withheld", is.na(s$mean_rmse[s$method == "lambert_c1"]) &&
  s$repeats_complete[s$method == "lambert_c1"] == 99L)
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) {
  stopifnot(length(args) == 2L, args[1] == "--smoke")
  state <- e$rcv1_validate_complete(file.path(root, args[2]), root, "qsar")
  first <- e$rcv1_read_shard(file.path(root, args[2]), 1L, state)
  wrong <- first; wrong$signature <- "wrong"
  rejected <- inherits(try(e$rcv1_validate_shard(wrong, state$tasks[1, ], state), silent = TRUE), "try-error")
  gate("wrong_signature_rejected", rejected)
  cv <- first$fit$fit$cv; model <- first$fit$fit$models$lambert_c1
  gate("smoke_lambert_shape_one_only", all(cv$shape[cv$family == "lambert"] == 1))
  losses <- vapply(1:4, function(k) mean(cv$loss[cv$family == "lambert" & cv$index == k]), numeric(1))
  gate("smoke_inner_mse_selects_lambda", model$selection$index == which.min(losses))
  history <- readLines(file.path(root, args[2], "progress.tsv"))
  gate("smoke_resume_recorded", any(grepl("interrupted_smoke", history)))
}
print(do.call(rbind, checks), row.names = FALSE)
cat("PASS:", length(checks), "fixed-c1 checks; production fits: 0.\n")
