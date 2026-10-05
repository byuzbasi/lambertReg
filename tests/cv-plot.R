library(lambertReg)
set.seed(221)
x <- matrix(rnorm(60 * 7), 60, 7)
y <- 3 + 2 * x[, 1] - x[, 3] + rnorm(60)
foldid <- rep(1:5, each = 12)
cv <- cv.lambert(x, y, foldid = foldid, lambda_fraction = 10^seq(0, -3, length.out = 7))
stopifnot(length(cv$nzero) == 7L, ncol(cv$lambert.fit$beta) == 7L,
  all(cv$nzero == colSums(abs(cv$lambert.fit$beta) > 1e-8)),
  cv$nzero[1] == 0L, cv$support_tol == 1e-8, cv$full_path_complete)
pd <- getFromNamespace(".lr_cv_plot_data", "lambertReg")
d <- pd(cv, "fraction"); a <- pd(cv, "lambda")
stopifnot(identical(d$nzero, cv$nzero), identical(d$index, seq_len(7L)),
  identical(d$error, cv$cvm), isTRUE(all.equal(d$x, log(cv$lambda_fraction))),
  isTRUE(all.equal(a$x, log(cv$lambda))), identical(d$label, as.character(cv$nzero)))
# Descriptive SE uses sample SD across folds and requires complete losses.
expected_se <- apply(cv$fold_loss, 1L, stats::sd)/sqrt(ncol(cv$fold_loss))
stopifnot(isTRUE(all.equal(d$se, expected_se)),
  isTRUE(all.equal(d$lower, cv$cvm - expected_se)),
  isTRUE(all.equal(d$upper, cv$cvm + expected_se)))
known <- cv
known$fold_loss[1, ] <- 1:5; known$cvm[1] <- 3
known$fold_loss[2, ] <- 2; known$cvm[2] <- 2
known$fold_loss[3, 1] <- NA_real_
k <- pd(known, "fraction")
stopifnot(abs(k$se[1] - sqrt(.5)) < 1e-14, k$se[2] == 0,
  is.na(k$se[3]), is.na(k$lower[3]), is.na(k$upper[3]))
# Plotting uses stored values only: no RNG consumption or changes to the fit.
before <- serialize(cv, NULL); rng <- .Random.seed
path <- tempfile(fileext = ".pdf"); grDevices::pdf(path, width = 9, height = 6)
par_before <- graphics::par(c("mar", "mgp", "las", "col.axis"))
stopifnot(identical(plot(cv), cv))
stopifnot(identical(graphics::par(c("mar", "mgp", "las", "col.axis")), par_before))
plot(cv, xaxis = "lambda", main = "Custom CV", col = "navy", label_cex = .7)
plot(cv, show_nzero = FALSE)
plot(cv, error_bars = FALSE)
# Legacy objects and unverified positions cannot be reported as zero counts.
legacy <- cv; legacy$nzero <- NULL
warned <- FALSE
withCallingHandlers(plot(legacy), warning = function(w) {
  warned <<- grepl("older CV object", conditionMessage(w)); invokeRestart("muffleWarning")
})
stopifnot(warned, all(pd(legacy, "fraction")$label == "--"))
unverified <- cv; unverified$nzero[2] <- NA_integer_; unverified$cvm[2] <- Inf
stopifnot(pd(unverified, "fraction")$label[2] == "--",
  is.na(pd(unverified, "fraction")$error[2]))
plot(unverified)
unverified_selected <- cv; unverified_selected$ok <- FALSE
status_warned <- FALSE
withCallingHandlers(plot(unverified_selected), warning = function(w) {
  status_warned <<- grepl("full-training fit is unverified", conditionMessage(w))
  invokeRestart("muffleWarning")
})
stopifnot(status_warned)
no_losses <- cv; no_losses$fold_loss <- NULL
se_warned <- FALSE
withCallingHandlers(plot(no_losses), warning = function(w) {
  se_warned <<- grepl("error bars are unavailable", conditionMessage(w))
  invokeRestart("muffleWarning")
})
stopifnot(se_warned, all(is.na(pd(no_losses, "fraction")$se)))
# A one-position path is supported without assuming an interval of positions.
single <- cv.lambert(x, y, foldid = foldid, lambda_fraction = 1)
plot(single)
grDevices::dev.off()
stopifnot(identical(before, serialize(cv, NULL)), identical(rng, .Random.seed),
  file.info(path)$size > 0)
failed <- suppressWarnings(cv.lambert(x, y, foldid = foldid,
  lambda_fraction = c(1, .1, .01), max_sweeps = 1))
stopifnot(all(is.na(failed$nzero[!failed$lambert.fit$verified])),
  all(failed$nzero[failed$lambert.fit$verified] >= 0))
cat("PASS: full-path counts, axis alignment, legacy objects, failure gaps, descriptive SE and plot state.\n")

# Two display styles and coefficient dispatch must preserve fits and RNG.
before <- serialize(cv, NULL); rng <- .Random.seed
path_data <- getFromNamespace(".lr_path_plot_data", "lambertReg")
b <- path_data(cv$lambert.fit, "lambda")
stopifnot(identical(b$beta, cv$lambert.fit$beta),
  identical(b$x, log(cv$lambda)), nrow(b$beta) == ncol(x))
# A deliberately invalid middle position is a gap, not a segment across it.
gapped <- cv$lambert.fit; gapped$verified[3] <- FALSE
g <- path_data(gapped, "fraction")
stopifnot(all(is.na(g$beta[,3])), identical(g$beta[,-3],gapped$beta[,-3]))
# No accidental switch to standardized coefficients on a nonstandard-scale fit.
nonstandard <- cv$lambert.fit
nonstandard$beta <- sweep(nonstandard$beta, 1, seq_len(nrow(nonstandard$beta)), "/")
stopifnot(identical(path_data(nonstandard,"lambda")$beta, nonstandard$beta),
  !identical(nonstandard$beta, nonstandard$beta_standardized))
file <- tempfile(fileext = ".pdf"); grDevices::pdf(file,width=10,height=7)
prior <- graphics::par(no.readonly=TRUE)
stopifnot(identical(plot(cv,style="detailed"),cv),
  identical(plot(cv,style="paper"),cv),
  identical(plot(cv,type="coefficients",labels=TRUE),cv),
  identical(plot(cv,type="coefficients",style="paper"),cv),
  identical(plot(cv$lambert.fit,labels=TRUE,col="navy",lwd=2),cv$lambert.fit))
plot(gapped,style="paper")
plot(single,type="coefficients",labels=TRUE)
plot(cv,type="coefficients",main="Custom paths",xlab="Penalty",ylab="Beta")
stopifnot(isTRUE(all.equal(prior,graphics::par(no.readonly=TRUE))))
warned_path <- FALSE
withCallingHandlers(plot(unverified_selected,type="coefficients"),warning=function(w) {
  warned_path <<- grepl("unverified",conditionMessage(w));invokeRestart("muffleWarning")
})
stopifnot(warned_path)
for (expr in list(quote(plot(cv,style="unknown")),quote(plot(cv,type="unknown")),
  quote(plot(cv,type="coefficients",labels=NA)),quote(plot(cv$lambert.fit,selected=0)))) {
  stopifnot(inherits(try(eval(expr),silent=TRUE),"try-error"))
}
grDevices::dev.off()
stopifnot(identical(before,serialize(cv,NULL)),identical(rng,.Random.seed))
cat("PASS: detailed/paper styles, coefficient dispatch, original-scale values, gaps and plot state.\n")
