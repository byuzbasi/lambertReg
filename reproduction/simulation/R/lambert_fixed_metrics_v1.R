# Evaluation only: no tuning, coefficient changes or model fitting.
lfb_add_metrics <- function(rows) {
  rows$tn <- NA_integer_; rows$mcc <- rows$excess_risk <- NA_real_
  rows$mcc_zero_denominator <- NA
  good <- which(rows$ok)
  if (!length(good)) return(rows)
  z <- rows[good, ]; tn <- z$p - z$tp - z$fp - z$fn
  counts <- cbind(z$tp, z$fp, tn, z$fn)
  stopifnot(all(is.finite(counts)), all(counts >= 0), all(counts == floor(counts)),
    all(z$tp + z$fp == z$selected_size), all(z$sigma2 > 0))
  den <- sqrt((z$tp + z$fp) * (z$tp + z$fn) * (tn + z$fp) * (tn + z$fn))
  rows$tn[good] <- as.integer(tn)
  rows$mcc_zero_denominator[good] <- den == 0
  rows$mcc[good] <- vapply(seq_along(good), function(i)
    mltools::mcc(TP = z$tp[i], FP = z$fp[i], TN = tn[i], FN = z$fn[i]), numeric(1))
  expected <- ifelse(den == 0, 0, (z$tp * tn - z$fp * z$fn) / den)
  er <- z$analytic_risk / z$sigma2 - 1
  stopifnot(all(is.finite(er)), all(er >= -1e-12),
    max(abs(rows$mcc[good] - expected)) < 1e-12)
  rows$excess_risk[good] <- pmax(er, 0)
  rows
}

lfb_summary <- function(v, expected = length(v), monte_carlo = FALSE) {
  z <- v[is.finite(v)]; n <- length(z); complete <- n == expected
  se <- if (complete && monte_carlo && n > 1L) sd(z) / sqrt(n) else NA_real_
  data.frame(expected = expected, available = n, complete = complete,
    mean = if (complete && n) mean(z) else NA_real_,
    median = if (complete && n) median(z) else NA_real_,
    q25 = if (complete && n) unname(quantile(z, .25)) else NA_real_,
    q75 = if (complete && n) unname(quantile(z, .75)) else NA_real_,
    mcse = se,
    conditional_mean = if (!complete && n) mean(z) else NA_real_)
}

lfb_stability <- function(Z) {
  stopifnot(is.matrix(Z), !anyNA(Z), all(Z %in% c(0, 1)), ncol(Z) > 0L)
  lfb_stability_counts(colSums(Z), nrow(Z))
}

lfb_stability_counts <- function(selected, M) {
  stopifnot(length(selected) > 0, M >= 0, all(selected >= 0 & selected <= M),
    all(selected == floor(selected)))
  q <- if (M) selected / M else rep(NA_real_, length(selected))
  qbar <- mean(q)
  status <- if (M < 2) "insufficient_fits" else if (qbar == 0 || qbar == 1)
    "all_empty_or_all_full_uninformative" else "estimated"
  value <- if (status == "estimated") 1 - mean(M * q * (1-q)/(M-1)) /
    (qbar * (1-qbar)) else NA_real_
  list(estimate = value, status = status, frequency = q, valid_fits = M)
}

lfb_numerical <- function(fit, method, row) {
  if (method == "ols") return(data.frame(method = method, ok = row$ok,
    status = row$status, partial_search = FALSE, cv_rows = 0L, valid_cv_rows = 0L,
    cv_coverage = NA_real_, candidates = 0L, valid_candidates = 0L,
    candidate_coverage = NA_real_, lower_endpoint = NA, upper_endpoint = NA,
    tuning_seconds = NA_real_, refit_seconds = NA_real_))
  cv <- fit$cv; eligible <- fit$models[[method]]$selection$eligible
  cv <- cv[cv$setting %in% eligible, ]
  groups <- split(seq_len(nrow(cv)), interaction(cv$setting, cv$index, drop = TRUE))
  valid <- cv$ok & is.finite(cv$loss)
  full <- vapply(groups, function(i) length(i) == 5L && all(valid[i]), logical(1))
  model <- fit$models[[method]]
  val <- function(name, fallback = NA_real_) {
    x <- row[[name]]; if (is.null(x)) x <- model[[name]]
    if (is.null(x)) fallback else x
  }
  index <- model$selection$index
  lower <- if (is.null(index)) NA else index == max(cv$index)
  upper <- if (is.null(index)) NA else index == 1L
  data.frame(method = method, ok = row$ok, status = row$status,
    partial_search = row$partial_search, cv_rows = nrow(cv), valid_cv_rows = sum(valid),
    cv_coverage = mean(valid), candidates = length(groups), valid_candidates = sum(full),
    candidate_coverage = mean(full), lower_endpoint = lower,
    upper_endpoint = upper, tuning_seconds = val("tuning_seconds"),
    refit_seconds = val("refit_seconds"))
}
