# Score eligible candidates from the same paths used for constrained BIC.
lhx_methods <- function() c("lambert_c1","lasso","mcp","scad")
lhx_scores <- function(rss,d,n,p) {
  stopifnot(n>exp(1),p>=1,length(rss)==length(d))
  good <- is.finite(rss) & rss>0 & is.finite(d) & d>=0 & d<=p & d==floor(d)
  z <- matrix(Inf,length(rss),4,dimnames=list(NULL,c("BIC","HBIC","EBIC_0.5","EBIC_1")))
  base <- n*log(rss[good]/n); q <- d[good]
  z[good,] <- cbind(base+q*log(n),base+q*log(p)*log(log(n)),
    base+q*log(n)+lchoose(p,q),base+q*log(n)+2*lchoose(p,q))
  z
}

score_ic <- function(task, cfg, paths, candidates, bic, prep, test, truth, e) {
  fields <- e$lic_metric_fields()
  rows <- lapply(lhx_methods(), function(method) {
    j <- e$lbc_settings_for(e$perf_settings(cfg), method)
    stopifnot(length(j) == 1L)
    cc <- candidates[[j]]
    scores <- lhx_scores(cc$rss, cc$d, task$n, task$p)
    z <- bic[bic$method == method & bic$criterion == "BIC", , drop = FALSE]
    pick <- function(criterion)
      e$lic_choose(scores[, criterion], cc$eligible, cfg$implementation$tie_tol)
    stopifnot(identical(pick("BIC"), z$index))
    do.call(rbind, lapply(c("HBIC", "EBIC_0.5", "EBIC_1"), function(criterion) {
      k <- pick(criterion)
      row <- z
      row$criterion <- criterion
      row$index <- k
      row$ok <- !is.na(k)
      row$criterion_score <- if (row$ok) scores[k, criterion] else NA_real_
      if (row$ok) {
        row$lambda <- cc$lambda[k]
        row$q <- cc$q[k]
        row$bic <- cc$BIC[k]
        row$kkt <- cc$kkt[k]
        row$curvature_status <- cc$curvature_status[k]
        row$cap_boundary <- cc$d[k] == floor(task$n/log(task$n))
        row$path_endpoint <- k %in% c(1L,cfg$lambda_points)
        row[fields] <- e$lic_metrics(paths[[j]]$beta[,k], prep, test, truth, cfg, e)[fields]
      } else {
        row[fields] <- NA
        row$status <- "no_valid_candidate"
      }
      row
    }))
  })
  do.call(rbind, rows)
}
