.lr_indices <- function(object, index) {
  if (is.null(index)) index <- seq_along(object$lambda)
  if (!is.numeric(index) || is.complex(index) || !length(index) ||
      any(!is.finite(index)) || any(index != floor(index)) ||
      any(index < 1 | index > length(object$lambda))) stop("Invalid path index.")
  if (any(!object$verified[index]))
    stop("Requested path includes unverified fits; inspect status and diagnostics.")
  as.integer(index)
}

.lr_dots <- function(...) {
  if (length(list(...))) stop("Unused arguments in ... .")
}

coef.lambert <- function(object, index = NULL, ...) {
  .lr_dots(...); index <- .lr_indices(object, index)
  rbind("(Intercept)" = object$intercept[index], object$beta[, index, drop = FALSE])
}

predict.lambert <- function(object, newx, index = NULL, ...) {
  .lr_dots(...); index <- .lr_indices(object, index)
  if (!is.matrix(newx) || !is.numeric(newx) || is.complex(newx) ||
      ncol(newx) != object$nvars || any(!is.finite(newx)))
    stop("newx must be a finite numeric matrix with the training predictor count.")
  if (!is.null(object$xnames) && !identical(colnames(newx), object$xnames))
    stop("newx column names and order must match training x.")
  sweep(newx %*% object$beta[, index, drop = FALSE], 2L, object$intercept[index], "+")
}

coef.cv.lambert <- function(object, ...) {
  .lr_dots(...)
  if (!isTRUE(object$ok)) stop("Selected fit is unverified: ", object$status)
  coef(object$lambert.fit, index = object$index)
}

predict.cv.lambert <- function(object, newx, ...) {
  .lr_dots(...)
  if (!isTRUE(object$ok)) stop("Selected fit is unverified: ", object$status)
  predict(object$lambert.fit, newx = newx, index = object$index)
}

print.lambert <- function(x, ...) {
  .lr_dots(...)
  cat("Lambert regression (fixed shape c = 1)\n",
      x$nobs, " observations; ", x$nvars, " predictors; ", length(x$lambda),
      " penalty levels\n", sum(x$verified), "/", length(x$lambda),
      " numerically verified fits (not global-optimum certificates)\n", sep = "")
  invisible(x)
}

print.cv.lambert <- function(x, ...) {
  .lr_dots(...)
  cat("Cross-validated Lambert regression (", x$nfolds, " folds)\n",
    "Selected index: ", x$index, "; full-training lambda: ", format(x$lambda.min),
    "\nStatus: ", x$status, "; partial search: ", x$partial_search, "\n", sep = "")
  invisible(x)
}

# Plotting reads the original-scale coefficient matrix, never standardized betas.
.lr_path_plot_data <- function(x, xaxis) {
  position <- if (xaxis == "fraction") x$lambda_fraction else x$lambda
  beta <- x$beta
  if (!is.matrix(beta) || !is.numeric(beta) || ncol(beta) != length(position) ||
      length(x$verified) != length(position) || anyNA(x$verified) ||
      any(!is.finite(position)) || any(position <= 0)) stop("Invalid coefficient path.")
  if (!any(x$verified)) stop("No verified path fits to plot.")
  if (any(!is.finite(beta[, x$verified, drop = FALSE]))) stop("Nonfinite verified coefficients.")
  beta[, !x$verified] <- NA_real_
  list(x = log(position), beta = beta,
    labels = if (is.null(rownames(beta))) paste0("V", seq_len(nrow(beta))) else rownames(beta))
}

# Compute limits and ticks together so rounding is applied only once.
.lr_axis_spec <- function(x, limits = NULL, reverse = FALSE) {
  if (is.null(limits)) {
    limits <- range(x, finite = TRUE)
    if (diff(limits) == 0) limits <- limits + c(-.5, .5)
    ticks <- pretty(limits, n = 5)
    limits <- range(ticks)
    if (reverse) limits <- rev(limits)
  } else {
    ticks <- pretty(range(limits), n = 5)
    ticks <- ticks[ticks >= min(limits) & ticks <= max(limits)]
  }
  list(limits = limits, ticks = ticks)
}

.lr_plot_axes <- function(horizontal, vertical) {
  graphics::axis(1, at = horizontal$ticks)
  graphics::axis(2, at = vertical$ticks)
  graphics::box(bty = "l")
}

plot.lambert <- function(x, ..., xaxis = c("lambda", "fraction"),
                         style = c("detailed", "paper"), labels = FALSE,
                         col = NULL, selected = NULL, selected_col = "black",
                         main = NULL, xlab = NULL, ylab = "Original-scale coefficient") {
  xaxis <- match.arg(xaxis); style <- match.arg(style)
  if (!is.logical(labels) || length(labels) != 1L || is.na(labels)) stop("labels must be TRUE or FALSE.")
  d <- .lr_path_plot_data(x, xaxis)
  if (!is.null(selected) && (length(selected) != 1L || !is.numeric(selected) ||
      is.na(selected) || !selected %in% seq_along(d$x))) stop("Invalid selected path index.")
  dots <- list(...)
  if (length(dots) && (is.null(names(dots)) || any(!nzchar(names(dots)))))
    stop("Additional graphical arguments must be named.")
  if (any(names(dots) %in% c("x", "y", "axes", "ann")) ||
      (!is.null(dots$log) && !identical(dots$log, ""))) stop("Path coordinates are managed internally and already log-transformed.")
  if (is.null(col)) col <- grDevices::hcl.colors(nrow(d$beta), "Dark 3")
  horizontal <- .lr_axis_spec(d$x, dots$xlim, reverse = TRUE)
  vertical <- .lr_axis_spec(d$beta, dots$ylim)
  dots$xlim <- horizontal$limits; dots$ylim <- vertical$limits
  if (is.null(dots$xaxs)) dots$xaxs <- "i"
  if (is.null(dots$yaxs)) dots$yaxs <- "i"
  if (is.null(dots$type)) dots$type <- if (length(d$x) == 1L) "p" else "l"
  if (is.null(dots$lty)) dots$lty <- 1
  if (is.null(dots$lwd)) dots$lwd <- 1.5
  if (is.null(xlab)) xlab <- if (xaxis == "fraction")
    expression(log(lambda/lambda[max])) else expression(log(lambda))
  if (missing(main) && style == "detailed") main <- "Lambert coefficient paths"
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)
  right <- if (labels) max(5, max(nchar(d$labels))*.6 + 2) else 1.5
  graphics::par(mar = c(if (style == "detailed") 5.5 else 4.5, 5, 3, right),
    mgp = c(2.8, .7, 0), las = 1)
  do.call(graphics::matplot, c(list(x = d$x, y = t(d$beta), col = col,
    axes = FALSE, ann = FALSE), dots))
  .lr_plot_axes(horizontal, vertical)
  usr <- graphics::par("usr")
  graphics::abline(h = 0, col = "gray75", lty = 3)
  if (!is.null(selected)) graphics::abline(v = d$x[selected], col = selected_col, lty = 2)
  if (labels) graphics::legend(usr[2] + .03*diff(usr[1:2]), usr[4], xjust = 0,
    legend = d$labels, col = rep(col, length.out = nrow(d$beta)),
    lty = rep(dots$lty, length.out = nrow(d$beta)), lwd = dots$lwd,
    cex = .75, bty = "n", xpd = NA)
  graphics::title(xlab = xlab, ylab = ylab, main = main)
  if (style == "detailed") {
    note <- "Original-scale coefficients; intercept excluded"
    if (!is.null(selected)) note <- paste0(note, " | CV-selected lambda = ",
      format(signif(x$lambda[selected], 4)))
    if (any(!x$verified)) note <- paste0(note, " | Gaps: unverified fits")
    graphics::mtext(note, side = 1, line = 4, cex = .7)
  }
  invisible(x)
}

# A missing count is never interpreted as an empty selected model.
.lr_cv_plot_data <- function(x, xaxis) {
  position <- if (xaxis == "fraction") x$lambda_fraction else x$lambda
  if (length(position) != length(x$cvm) || any(!is.finite(position)) || any(position <= 0))
    stop("CV object has an invalid penalty path.")
  nzero <- x$nzero
  if (is.null(nzero)) nzero <- rep(NA_integer_, length(position))
  if (length(nzero) != length(position)) stop("CV object has incompatible variable counts.")
  se <- rep(NA_real_, length(position))
  losses <- x$fold_loss
  if (!is.null(losses)) {
    if (!is.matrix(losses) || !is.numeric(losses) ||
        nrow(losses) != length(position) || ncol(losses) < 2L)
      stop("CV object has incompatible fold losses.")
    complete <- apply(is.finite(losses), 1L, all) & is.finite(x$cvm)
    if (any(complete)) se[complete] <- apply(losses[complete, , drop = FALSE],
      1L, stats::sd)/sqrt(ncol(losses))
  }
  error <- replace(x$cvm, !is.finite(x$cvm), NA_real_)
  data.frame(index = seq_along(position), x = log(position),
    error = error, se = se, lower = error - se, upper = error + se, nzero = nzero,
    label = ifelse(is.na(nzero), "--", as.character(nzero)))
}

plot.cv.lambert <- function(x, xaxis = c("fraction", "lambda"),
                            show_nzero = TRUE, label_cex = NULL,
                            col = "#315C76", selected_col = "#C65D32",
                            main = NULL, xlab = NULL,
                            ylab = "Cross-validation MSE", error_bars = TRUE,
                            bar_col = "#B5C4CE", style = c("detailed", "paper"),
                            type = c("cv", "coefficients"), labels = FALSE, ...) {
  xaxis <- match.arg(xaxis); style <- match.arg(style); type <- match.arg(type)
  if (type == "coefficients") {
    if (is.null(x$lambert.fit)) stop("This CV object has no stored coefficient path.")
    if (!isTRUE(x$ok)) warning("The CV-selected full-training fit is unverified; inspect status before prediction.", call. = FALSE)
    path_args <- list(x = x$lambert.fit, xaxis = xaxis, style = style,
      labels = labels, selected = x$index, ...)
    if (!missing(selected_col)) path_args$selected_col <- selected_col
    if (!missing(col)) path_args$col <- col
    if (!missing(main)) path_args$main <- main
    if (!missing(xlab)) path_args$xlab <- xlab
    if (!missing(ylab)) path_args$ylab <- ylab
    do.call(plot.lambert, path_args)
    return(invisible(x))
  }
  if (missing(main) && style == "detailed") main <- "Lambert cross-validation"
  if (!is.logical(show_nzero) || length(show_nzero) != 1L || is.na(show_nzero))
    stop("show_nzero must be TRUE or FALSE.")
  if (!is.logical(error_bars) || length(error_bars) != 1L || is.na(error_bars))
    stop("error_bars must be TRUE or FALSE.")
  if (!is.null(label_cex)) .lr_positive(label_cex, "label_cex")
  if (!any(is.finite(x$cvm))) stop("No finite CV scores to plot.")
  d <- .lr_cv_plot_data(x, xaxis)
  if (show_nzero && is.null(x$nzero))
    warning("Variable counts are absent from this older CV object; shown as --. No models were fitted.", call. = FALSE)
  if (error_bars && any(is.finite(d$error) & !is.finite(d$se)))
    warning("Some fold-based error bars are unavailable; they are omitted, not replaced by zero.", call. = FALSE)
  dots <- list(...)
  if (length(dots) && (is.null(names(dots)) || any(!nzchar(names(dots)))))
    stop("Additional graphical arguments must be named.")
  if (!is.null(dots$log) && !identical(dots$log, ""))
    stop("Do not log-transform these axes again; use xaxis to select the penalty scale.")
  if (any(names(dots) %in% c("x", "y", "type", "axes", "ann")))
    stop("The CV plot manages x, y, type, axes and ann internally.")
  horizontal <- .lr_axis_spec(d$x, dots$xlim)
  vertical <- .lr_axis_spec(
    if (error_bars) c(d$error, d$lower, d$upper) else d$error, dots$ylim)
  dots$xlim <- horizontal$limits; dots$ylim <- vertical$limits
  if (is.null(dots$xaxs)) dots$xaxs <- "i"
  if (is.null(dots$yaxs)) dots$yaxs <- "i"
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)
  graphics::par(mar = c(if (style == "detailed") 5.5 else 4.5, 4.5,
    if (style == "detailed") { if (show_nzero) 4.7 else 3 } else { if (show_nzero) 3.5 else 2 }, 1.5),
    mgp = c(2.8, .7, 0), tcl = -.25, las = 1,
    col.axis = "black", col.lab = "black", fg = "#BDC9D1")
  do.call(graphics::plot, c(list(x = d$x, y = d$error, type = "n",
    axes = FALSE, ann = FALSE), dots))
  graphics::abline(h = vertical$ticks, col = "#EAF0F3", lwd = .7)
  .lr_plot_axes(horizontal, vertical)
  if (error_bars) {
    keep <- is.finite(d$lower) & is.finite(d$upper)
    cap <- abs(diff(graphics::par("usr")[1:2])) * .003
    graphics::segments(d$x[keep], d$lower[keep], d$x[keep], d$upper[keep],
      col = bar_col, lwd = 1)
    graphics::segments(d$x[keep] - cap, d$lower[keep], d$x[keep] + cap, d$lower[keep],
      col = bar_col, lwd = 1)
    graphics::segments(d$x[keep] - cap, d$upper[keep], d$x[keep] + cap, d$upper[keep],
      col = bar_col, lwd = 1)
  }
  ord <- order(d$x)
  graphics::lines(d$x[ord], d$error[ord], col = col, lwd = 1.4)
  graphics::points(d$x, d$error, pch = 21, bg = "white", col = col, cex = .65, xpd = NA)
  chosen <- x$index
  selected <- length(chosen) == 1L && !is.na(chosen) && chosen %in% d$index &&
    is.finite(d$error[chosen])
  if (selected) {
    graphics::abline(v = d$x[chosen], col = selected_col, lty = 2, lwd = 1.1)
    graphics::points(d$x[chosen], d$error[chosen], pch = 21, bg = selected_col,
      col = "white", cex = 1.15, xpd = NA)
    if (!isTRUE(x$ok))
      warning("The CV-selected full-training fit is unverified; inspect status before prediction.", call. = FALSE)
  }
  if (show_nzero) {
    if (is.null(label_cex)) {
      label_cex <- .8
      if (nrow(d) > 1L) {
        ord <- order(d$x)
        widths <- abs(graphics::strwidth(d$label[ord], cex = 1))
        spacing <- diff(d$x[ord])
        needed <- (widths[-length(widths)] + widths[-1L])/2
        label_cex <- min(.8, .85 * min(spacing/needed))
      }
    }
    # Draw every label (including repeats); axis() can silently thin labels.
    usr <- graphics::par("usr")
    top <- usr[4] + .8 * graphics::strheight("M", cex = label_cex)
    graphics::text(d$x, top, labels = d$label, srt = 0, adj = c(.5, 0),
      cex = label_cex, col = "black", xpd = NA)
    graphics::mtext("Selected variables", side = 3, line = 1.2, cex = .85, col = "black")
  }
  if (is.null(xlab)) xlab <- if (xaxis == "fraction")
    expression(log(lambda/lambda[max])) else expression(log(lambda))
  graphics::title(xlab = xlab, ylab = ylab)
  if (!is.null(main)) graphics::mtext(main, side = 3,
    line = if (style == "detailed") { if (show_nzero) 3 else 1.4 } else {
      if (show_nzero) 2.6 else 1 }, font = 2, adj = 0, col = "black")
  if (style == "detailed") {
    note <- if (selected) paste0("Selected lambda = ", format(signif(x$lambda.min, 4)),
      " | CV MSE = ", format(signif(d$error[chosen], 4)),
      " | Selected = ", d$label[chosen]) else "No CV-selected candidate"
    if (!isTRUE(x$ok)) note <- paste0(note, " | Unverified selected fit")
    if (show_nzero && anyNA(d$nzero)) note <- paste0(note, " | -- unavailable")
    graphics::mtext(note, side = 1, line = 4, cex = .72, adj = 0, col = "black")
  }
  invisible(x)
}
