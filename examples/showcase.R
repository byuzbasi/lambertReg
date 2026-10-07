# Small synthetic demonstration; not a Monte Carlo study or a paper result.
# Source this file, then call lambert_showcase("lambert-example-v1").
lambert_showcase <- function(out) {
  if (length(out) != 1L || !is.character(out) || !nzchar(out))
    stop("Supply one output directory.")
  if (file.exists(out)) stop("Refusing existing output: ", out)
  if (!requireNamespace("lambertReg", quietly = TRUE)) stop("Install lambertReg first.")
  if (!requireNamespace("digest", quietly = TRUE)) stop("The example manifest requires digest.")
  dir.create(out, recursive = TRUE)
  set.seed(42)
  x <- matrix(rnorm(100 * 8), 100, 8)
  colnames(x) <- paste0("x", seq_len(ncol(x)))
  y <- 2 * x[, 1] - x[, 2] + 0.5 * x[, 3] + rnorm(100)
  train <- seq_len(80)
  foldid <- sample(rep(1:5, length.out = length(train)))
  fit <- lambertReg::cv.lambert(x[train, , drop = FALSE], y[train], foldid = foldid)
  stopifnot(fit$ok, !fit$partial_search, fit$full_path_complete)
  pred <- drop(stats::predict(fit, newx = x[-train, , drop = FALSE]))
  summary <- data.frame(lambda = fit$lambda.min, selected = fit$nzero[fit$index],
    cv_mse = fit$cvm[fit$index], test_rmse = sqrt(mean((y[-train] - pred)^2)),
    test_mae = mean(abs(y[-train] - pred)))
  write.csv(x, file.path(out, "predictors.csv"), row.names = FALSE)
  write.csv(data.frame(id = seq_len(100), y = y,
    sample = ifelse(seq_len(100) %in% train, "train", "test"),
    fold = c(foldid, rep(NA_integer_, 20))), file.path(out, "observations.csv"), row.names = FALSE)
  write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
  write.csv(data.frame(lambda = fit$lambda, fraction = fit$lambda_fraction,
    cv_mse = fit$cvm, selected = fit$nzero), file.path(out, "cv_path.csv"), row.names = FALSE)
  write.csv(data.frame(predictor = rownames(fit$lambert.fit$beta),
    fit$lambert.fit$beta, check.names = FALSE), file.path(out, "coefficients.csv"), row.names = FALSE)
  saveRDS(fit, file.path(out, "fit.rds"))
  saveRDS(list(seed = 42L, train = train, foldid = foldid,
    lambda_fraction = fit$lambda_fraction, controls = fit$control,
    package_version = as.character(utils::packageVersion("lambertReg"))), file.path(out, "configuration.rds"))
  draw <- function(name, code) {
    grDevices::png(file.path(out, name), width = 1600, height = 1000, res = 160)
    on.exit(grDevices::dev.off())
    force(code)
  }
  draw("cv.png", graphics::plot(fit))
  draw("coefficients.png", graphics::plot(fit, type = "coefficients", labels = TRUE))
  versions <- c(paste("R", getRversion()), paste("Platform", R.version$platform),
    paste("lambertReg", utils::packageVersion("lambertReg")),
    paste("Rcpp", utils::packageVersion("Rcpp")), paste("digest", utils::packageVersion("digest")))
  writeLines(versions, file.path(out, "versions.txt"))
  files <- list.files(out, full.names = TRUE)
  manifest <- data.frame(file = basename(files), bytes = file.info(files)$size,
    sha256 = vapply(files, function(p) digest::digest(file = p, algo = "sha256"), ""))
  write.csv(manifest, file.path(out, "manifest.csv"), row.names = FALSE)
  writeLines("Synthetic example verified; 80 training and 20 held-out observations; seed 42.", file.path(out, "COMPLETED"))
  invisible(list(fit = fit, summary = summary, directory = out))
}
