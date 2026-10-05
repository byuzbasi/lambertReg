# Public interface: fixed c = 1; training-only centering and RMS scaling.
.lr_data <- function(x, y) {
  if (!is.matrix(x) || !is.numeric(x) || is.complex(x) ||
      nrow(x) < 3L || ncol(x) < 1L || any(!is.finite(x)))
    stop("x must be a finite numeric matrix with at least 3 rows and 1 column.")
  if (!is.numeric(y) || is.complex(y) || !is.null(dim(y)) ||
      length(y) != nrow(x) || any(!is.finite(y)))
    stop("y must be a finite numeric vector with nrow(x) entries.")
  if (!is.null(colnames(x)) && (anyNA(colnames(x)) || any(!nzchar(colnames(x))) ||
      anyDuplicated(colnames(x)))) stop("Predictor names must be nonempty and unique.")
  storage.mode(x) <- "double"
  list(x = x, y = as.double(y))
}

.lr_positive <- function(x, name, integer = FALSE) {
  if (!is.numeric(x) || is.complex(x) || length(x) != 1L || !is.finite(x) ||
      x <= 0 || (integer && (x != floor(x) || x > .Machine$integer.max)))
    stop(name, " must be a finite positive ", if (integer) "integer." else "number.")
  invisible(TRUE)
}

.lr_grid <- function(x, name, relative = FALSE) {
  if (!is.numeric(x) || is.complex(x) || !is.null(dim(x)) || !length(x) ||
      any(!is.finite(x)) || any(x <= 0) || any(diff(x) >= 0) ||
      (relative && any(x > 1)))
    stop(name, " must be a strictly decreasing positive finite vector",
         if (relative) " with entries <= 1." else ".")
  if (any(!is.finite(x * exp(1))) || any(!is.finite(x^2 * expm1(2)/4)) || any(x^2 * expm1(2)/4 <= 0))
    stop(name, " exceeds the representable numerical range.")
  as.double(x)
}

.lr_config <- function(max_sweeps, kkt_tol, lambda = NULL,
                       lambda_fraction = 10^seq(0, -3, length.out = 60L)) {
  .lr_positive(max_sweeps, "max_sweeps", TRUE)
  .lr_positive(kkt_tol, "kkt_tol")
  list(lambda = lambda, lambda_fraction = lambda_fraction,
       lambda_points = if (is.null(lambda)) length(lambda_fraction) else length(lambda),
       max_coordinate_sweeps_per_start = as.integer(max_sweeps),
       normalized_kkt_tolerance = kkt_tol,
       implementation = list(tie_tol = 1e-10, standardization_tol = 1e-12,
         curvature_join_tol = 1e-8, curvature_negative_tol = 1e-8))
}

.lr_prepare <- function(x, y, cfg) {
  p <- perf_prepare(x, y, cfg)
  .lr_grid(perf_lambdas(p, NULL, cfg), "derived lambda")
  p
}

.lr_object <- function(path, prep, cfg, xnames, call) {
  setting <- list(family = "lambert", shape = 1)
  curvature <- lapply(seq_len(ncol(path$beta)), function(k) {
    if (!path$diagnostics$ok[k]) return(list(status = "unavailable",
      min_eigenvalue = NA_real_, condition = NA_real_))
    perf_hessian(path$beta[, k], prep, path$diagnostics$lambda[k], setting, cfg)
  })
  curvature_status <- vapply(curvature, `[[`, character(1), "status")
  verified <- path$diagnostics$ok & curvature_status != "negative_curvature"
  status <- path$diagnostics$status
  status[path$diagnostics$ok & !verified] <- "negative_curvature"
  beta <- sweep(path$beta, 1L, prep$sx, "/")
  intercept <- prep$my - colSums(sweep(beta, 1L, prep$mx, "*"))
  rownames(beta) <- if (is.null(xnames)) paste0("V", seq_len(nrow(beta))) else xnames
  colnames(beta) <- paste0("s", seq_len(ncol(beta)))
  names(intercept) <- colnames(beta)
  structure(list(call = call, beta = beta, intercept = intercept,
    lambda = path$diagnostics$lambda,
    lambda_fraction = path$diagnostics$lambda/prep$lambda0,
    beta_standardized = path$beta, diagnostics = path$diagnostics,
    attempts = path$attempts, curvature = curvature, verified = verified,
    status = status, center = prep$mx, scale = prep$sx,
    response_mean = prep$my, lambda_max = prep$lambda0,
    nobs = nrow(prep$Z), nvars = ncol(prep$Z), xnames = xnames,
    control = cfg, shape = 1, elapsed = path$elapsed), class = "lambert")
}

lambert <- function(x, y, lambda = NULL,
                    lambda_fraction = 10^seq(0, -3, length.out = 60L),
                    max_sweeps = 10000L, kkt_tol = 1e-7) {
  cl <- match.call(); d <- .lr_data(x, y)
  if (!is.null(lambda) && !missing(lambda_fraction))
    stop("Supply lambda or lambda_fraction, not both.")
  if (!is.null(lambda)) lambda <- .lr_grid(lambda, "lambda")
  lambda_fraction <- .lr_grid(lambda_fraction, "lambda_fraction", TRUE)
  cfg <- .lr_config(max_sweeps, kkt_tol, lambda, lambda_fraction)
  prep <- .lr_prepare(d$x, d$y, cfg)
  path <- perf_path(prep, list(family = "lambert", shape = 1), cfg)
  fit <- .lr_object(path, prep, cfg, colnames(x), cl)
  if (any(!fit$verified)) warning(sum(!fit$verified),
    " path fits are unverified; inspect diagnostics and status.", call. = FALSE)
  fit
}

lambert_penalty <- function(t, lambda) {
  .lr_positive(lambda, "lambda")
  .lr_grid(lambda, "lambda")
  if (!is.numeric(t) || is.complex(t) || !is.null(dim(t)) || any(!is.finite(t)) || any(t < 0))
    stop("t must be a finite nonnegative numeric vector.")
  lambert_values_cpp(as.double(t), lambda, 1)$penalty
}

lambert_threshold <- function(z, lambda) {
  .lr_positive(lambda, "lambda")
  .lr_grid(lambda, "lambda")
  if (!is.numeric(z) || is.complex(z) || !is.null(dim(z)) || any(!is.finite(z)))
    stop("z must be a finite numeric vector.")
  a <- abs(z); out <- numeric(length(z)); cut <- lambda * exp(1)
  mid <- a > lambda & a < cut
  out[mid] <- z[mid] * log1p((a[mid] - lambda)/lambda)
  tail <- a >= cut; out[tail] <- z[tail]
  out
}
