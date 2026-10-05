cv.lambert <- function(x, y, nfolds = 5L, foldid = NULL,
                       lambda_fraction = 10^seq(0, -3, length.out = 60L),
                       max_sweeps = 10000L, kkt_tol = 1e-7) {
  cl <- match.call(); d <- .lr_data(x, y); n <- nrow(d$x)
  lambda_fraction <- .lr_grid(lambda_fraction, "lambda_fraction", TRUE)
  cfg <- .lr_config(max_sweeps, kkt_tol, lambda_fraction = lambda_fraction)
  if (is.null(foldid)) {
    .lr_positive(nfolds, "nfolds", TRUE)
    if (nfolds < 2L || nfolds > n) stop("nfolds must be between 2 and nrow(x).")
    foldid <- sample(rep(seq_len(nfolds), length.out = n))
  } else {
    if (!is.numeric(foldid) || is.complex(foldid) || !is.null(dim(foldid)) ||
        length(foldid) != n || any(!is.finite(foldid)) ||
        any(foldid != floor(foldid)) || any(foldid < 1) || any(foldid > n))
      stop("foldid must contain one integer fold label per observation.")
    K <- max(foldid)
    if (K < 2 || !identical(sort(unique(as.integer(foldid))), seq_len(as.integer(K))))
      stop("foldid labels must be consecutive integers starting at 1, with at least 2 folds.")
    if (!missing(nfolds) && (!identical(as.double(nfolds), as.double(K))))
      stop("nfolds does not match foldid.")
    nfolds <- as.integer(K)
  }
  foldid <- as.integer(foldid)
  if (any(n - tabulate(foldid, nfolds) < 3L))
    stop("Every training fold must contain at least 3 observations.")
  setting <- list(family = "lambert", shape = 1)
  # Prepare every training fold first: input problems fail before path fitting.
  preps <- lapply(seq_len(nfolds), function(f) {
    tryCatch(.lr_prepare(d$x[foldid != f, , drop = FALSE], d$y[foldid != f], cfg),
      error = function(e) stop("Training fold ", f, ": ", conditionMessage(e), call. = FALSE))
  })
  fullprep <- .lr_prepare(d$x, d$y, cfg)
  paths <- vector("list", nfolds); records <- attempts <- vector("list", nfolds)
  losses <- matrix(NA_real_, length(lambda_fraction), nfolds)
  for (f in seq_len(nfolds)) {
    pp <- perf_path(preps[[f]], setting, cfg); paths[[f]] <- pp
    for (k in which(pp$diagnostics$ok)) {
      pred <- perf_predict(pp$beta[, k], preps[[f]], d$x[foldid == f, , drop = FALSE])
      if (all(is.finite(pred))) losses[k, f] <- mean((d$y[foldid == f] - pred)^2)
    }
    records[[f]] <- cbind(fold = f, pp$diagnostics, loss = losses[, f])
    attempts[[f]] <- cbind(fold = f, pp$attempts)
  }
  cvm <- rowMeans(losses); cvm[!is.finite(cvm)] <- Inf
  valid <- which(is.finite(cvm)); partial <- any(!is.finite(cvm))
  index <- NA_integer_; ok <- FALSE; status <- "no_valid_cv"
  # Compute the full training path once, in the original coordinate order.
  # Later positions cannot change the selected prefix or the CV scores.
  fullpath <- perf_path(fullprep, setting, cfg)
  refit <- .lr_object(fullpath, fullprep, cfg, colnames(x), cl)
  nzero <- rep(NA_integer_, length(lambda_fraction))
  verified <- which(refit$verified)
  nzero[verified] <- colSums(abs(refit$beta[, verified, drop = FALSE]) > 1e-8)
  selected_curvature <- NULL
  if (length(valid)) {
    best <- min(cvm[valid])
    index <- which(cvm <= best + cfg$implementation$tie_tol * (1 + best))[1L]
    selected_curvature <- lapply(seq_len(nfolds), function(f)
      perf_hessian(paths[[f]]$beta[, index], preps[[f]],
        paths[[f]]$diagnostics$lambda[index], setting, cfg))
    negative <- refit$diagnostics$ok[index] &&
      any(vapply(selected_curvature, function(h) h$status == "negative_curvature", logical(1)))
    ok <- refit$verified[index] && !negative
    status <- if (negative) "negative_curvature" else refit$status[index]
  }
  fit <- structure(list(call = cl, index = index, lambda_fraction = lambda_fraction,
    lambda = fullprep$lambda0 * lambda_fraction,
    lambda.min = if (is.na(index)) NA_real_ else fullprep$lambda0 * lambda_fraction[index],
    fraction.min = if (is.na(index)) NA_real_ else lambda_fraction[index],
    cvm = cvm, fold_loss = losses, foldid = foldid, nfolds = nfolds,
    cv_diagnostics = do.call(rbind, records), attempts = do.call(rbind, attempts),
    preprocessing = lapply(c(preps, list(fullprep)), function(p) p[c("mx", "my", "sx", "lambda0")]),
    selected_curvature = selected_curvature, lambert.fit = refit,
    nzero = nzero, support_tol = 1e-8,
    full_path_complete = all(refit$verified),
    ok = ok, status = status, partial_search = partial, control = cfg, shape = 1),
    class = "cv.lambert")
  if (partial) warning("CV search is incomplete; unavailable candidates were not eligible. See cv_diagnostics.", call. = FALSE)
  if (any(!refit$verified)) warning("Some full-training path fits are unverified; their variable counts are unavailable.", call. = FALSE)
  if (!ok) warning("Selected fit is unavailable or unverified: ", status, call. = FALSE)
  fit
}
