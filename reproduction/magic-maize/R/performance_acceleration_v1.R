# Opt-in execution layer. Frozen methods, folds, grids and acceptance tolerances
# are supplied by the caller; existing runners do not load this file.
pacc_environment <- function(root, accelerated = TRUE, early_stop = FALSE, compile = TRUE) {
  stopifnot(is.logical(early_stop), length(early_stop)==1L, !is.na(early_stop),
    is.logical(accelerated), length(accelerated)==1L, !is.na(accelerated),
    accelerated || !early_stop)
  e <- new.env(parent=globalenv())
  for (f in c("R/penalties.R", "R/lambert_family.R", "R/performance_models.R",
      "R/performance_run.R", "R/lambert_fixed_metrics_v1.R", "R/lambert_fixed_simulation_v1.R"))
    sys.source(file.path(root,f),e)
  fixed <- e$lfb_sim_environment(root); parent.env(fixed) <- e
  e$perf_config <- fixed$perf_config
  if (compile) {
    e$perf_preflight(FALSE,root)
    Rcpp::sourceCpp(file.path(root,"src/lambert_performance.cpp"),env=e,
      cacheDir=file.path(tempdir(),"pacc_rcpp"),showOutput=FALSE)
    if (accelerated) Rcpp::sourceCpp(file.path(root,"src/lambert_performance_fast_v1.cpp"),
      env=e,cacheDir=file.path(tempdir(),"pacc_rcpp"),showOutput=FALSE)
  }
  if (accelerated) {
    e$reference_path <- e$perf_path
    # The ENET branch is delegated unchanged and does not recurse through perf_path.
    e$lambert_cd_cpp <- e$lambert_cd_fast_v1
    e$lambert_values_cpp <- e$lambert_values_fast_v1
    e$perf_path <- function(prep,setting,cfg,pulse=function(...)NULL,last=cfg$lambda_points)
      pacc_path(prep,setting,cfg,pulse,last,e)
    e$perf_train_models <- function(X,y,folds,cfg,pulse=function(...)NULL)
      pacc_train(X,y,folds,cfg,pulse,e,early_stop)
  }
  e
}

# A nonnegative partial fold sum / total folds bounds the final mean from below.
# A full-path incumbent and a positive roundoff buffer are required. Missing
# losses, negative values and an absent incumbent can never certify pruning.
pacc_tail_certificate <- function(partial, k, incumbent, folds, tie_tol) {
  stopifnot(folds>=2L, k>=1L, k<=length(partial), is.finite(tie_tol),tie_tol>=0)
  lb <- partial[seq.int(k,length(partial))]/folds
  guard <- 64*.Machine$double.eps*pmax(1,abs(lb),abs(incumbent))
  accepted <- is.finite(incumbent) && incumbent>=0 && all(is.finite(lb)) &&
    all(lb>=0) && all(lb>incumbent+tie_tol*(1+incumbent)+guard)
  list(accepted=accepted, lower_bound=lb, incumbent=incumbent, roundoff_buffer=guard)
}

pacc_path <- function(prep,setting,cfg,pulse,last,engine,before=NULL,after=NULL) {
  began <- proc.time()[[3]]
  if (setting$family=="enet") {
    z <- engine$reference_path(prep,setting,cfg,pulse,last)
    z$attempts$executed <- TRUE; z$attempts$reused <- FALSE
    z$acceleration <- list(solver_seconds=z$elapsed,diagnose_seconds=NA_real_,
      duplicate_starts=0L,starts_executed=1L,pruned_lambda=0L)
    return(z)
  }
  lambda <- engine$perf_lambdas(prep,setting,cfg)[seq_len(last)]
  beta <- matrix(NA_real_,ncol(prep$Z),last); rows <- vector("list",last)
  attempts <- list(); warm <- numeric(ncol(prep$Z))
  solver_seconds <- diagnose_seconds <- 0; duplicate_starts <- starts_executed <- 0L
  certificate <- NULL
  for (k in seq_len(last)) {
    certificate <- if (is.null(before)) NULL else before(k)
    if (!is.null(certificate) && certificate$accepted) {
      for (l in seq.int(k,last)) rows[[l]] <- data.frame(index=l,lambda=lambda[l],ok=FALSE,
        status="certified_cv_pruned",kkt=NA_real_,solver_lambda=lambda[l],solver_shape=setting$shape,
        objective=NA_real_,iterations=NA_real_,selected_start="not_computed",elapsed=0,
        max_increase=NA_real_,warning="",error="")
      certificate$first_index <- k
      break
    }
    candidates <- list()
    for (start_id in c("warm","zero")) {
      if (start_id=="zero" && all(warm==0) && candidates$warm$row$ok) {
        candidate <- candidates$warm; candidate$row$selected_start <- "zero"
        candidate$row$elapsed <- 0
        duplicate_starts <- duplicate_starts+1L
        audit <- cbind(candidate$row,executed=FALSE,reused=TRUE)
      } else {
        start <- if (start_id=="warm") warm+0 else numeric(length(warm))
        ts <- proc.time()[[3]]
        cap <- engine$perf_capture(function() {
          if (setting$family=="lambert") engine$lambert_cd_cpp(prep$Z,prep$y,start,lambda[k],
            setting$shape,cfg$max_coordinate_sweeps_per_start,cfg$normalized_kkt_tolerance/2,pulse)
          else ncvreg::ncvfit(prep$Z,prep$y,init=start,penalty=toupper(setting$family),
            gamma=setting$shape,lambda=lambda[k],eps=cfg$implementation$ncvreg_eps,
            max.iter=cfg$implementation$ncvreg_max_iter,warn=TRUE)
        })
        solver_seconds <- solver_seconds+proc.time()[[3]]-ts
        starts_executed <- starts_executed+1L
        fit <- cap$value; failed <- inherits(fit,"perf_error")
        b <- if (failed) rep(NA_real_,length(warm)) else as.numeric(fit$beta)
        td <- proc.time()[[3]]
        dd <- engine$perf_diagnose(b,prep,lambda[k],setting,cfg)
        diagnose_seconds <- diagnose_seconds+proc.time()[[3]]-td
        iter <- if (failed) NA_real_ else if (setting$family=="lambert") fit$iterations else fit$iter
        solver_ok <- !failed && if (setting$family=="lambert") fit$converged else iter<cfg$implementation$ncvreg_max_iter
        increased <- !failed && setting$family=="lambert" && fit$increase_failed
        ok <- solver_ok && dd$ok && !increased
        status <- if (failed) "solver_error" else if (increased) "objective_increased" else
          if (!solver_ok) "budget_exhausted" else if (!dd$ok) "kkt_not_met" else "verified_stationary"
        row <- data.frame(index=k,lambda=lambda[k],ok=ok,status=status,kkt=dd$kkt,
          solver_lambda=lambda[k],solver_shape=setting$shape,objective=dd$objective,
          iterations=iter,selected_start=start_id,elapsed=proc.time()[[3]]-ts,
          max_increase=if (!failed && setting$family=="lambert") fit$max_increase else NA_real_,
          warning=cap$warnings,error=if(failed) fit$message else "")
        candidate <- list(beta=b,row=row); audit <- cbind(row,executed=TRUE,reused=FALSE)
      }
      candidates[[start_id]] <- candidate
      attempts[[length(attempts)+1L]] <- audit
      pulse(paste(setting$family,setting$shape,"lambda",k,start_id))
    }
    good <- which(vapply(candidates,function(x)x$row$ok,logical(1)))
    if (length(good)) {
      q <- vapply(candidates[good],function(x)x$row$objective,numeric(1)); best <- min(q)
      selected <- good[which(q<=best+cfg$implementation$tie_tol*(1+abs(best)))[1]]
      beta[,k] <- candidates[[selected]]$beta; rows[[k]] <- candidates[[selected]]$row; warm <- beta[,k]
    } else { rows[[k]] <- candidates[[1]]$row; rows[[k]]$selected_start <- "none_verified"; warm[] <- 0 }
    if (!is.null(after)) after(k,beta[,k],rows[[k]])
  }
  list(beta=beta,diagnostics=do.call(rbind,rows),attempts=do.call(rbind,attempts),
    elapsed=proc.time()[[3]]-began,certificate=certificate,
    acceleration=list(solver_seconds=solver_seconds,diagnose_seconds=diagnose_seconds,
      duplicate_starts=duplicate_starts,starts_executed=starts_executed,
      pruned_lambda=sum(vapply(rows,function(x)x$status=="certified_cv_pruned",logical(1)))))
}

pacc_train <- function(X,y,folds,cfg,pulse,engine,early_stop=FALSE) {
  settings <- engine$perf_settings(cfg); K <- cfg$folds; J <- nrow(settings); L <- cfg$lambda_points
  t0 <- proc.time()[[3]]
  preps <- lapply(seq_len(K),function(f)engine$perf_prepare(X[folds!=f,,drop=FALSE],y[folds!=f],cfg))
  fullprep <- engine$perf_prepare(X,y,cfg)
  validation_X <- lapply(seq_len(K),function(f)
    sweep(sweep(X[folds==f,,drop=FALSE],2,preps[[f]]$mx,"-"),2,preps[[f]]$sx,"/"))
  preprocessing_seconds <- proc.time()[[3]]-t0; prediction_seconds <- 0
  paths <- lapply(seq_len(J),function(j)vector("list",K)); cv <- vector("list",J)
  losses <- lapply(seq_len(J),function(j)matrix(NA_real_,L,K))
  records <- attempts <- timings <- certificates <- list(); times <- matrix(0,J,K)
  for (f in seq_len(K)) for (j in seq_len(J)) {
    loss <- rep(NA_real_,L); incumbent <- Inf
    partial <- if (f==K) rowSums(losses[[j]][,seq_len(K-1L),drop=FALSE]) else rep(NA_real_,L)
    after <- function(k,b,row) {
      if (row$ok) {
        tp <- proc.time()[[3]]
        pred <- preps[[f]]$my+drop(validation_X[[f]]%*%b)
        if (all(is.finite(pred))) loss[k] <<- mean((y[folds==f]-pred)^2)
        prediction_seconds <<- prediction_seconds+proc.time()[[3]]-tp
      }
      if (f==K && is.finite(partial[k]) && is.finite(loss[k]))
        incumbent <<- min(incumbent,(partial[k]+loss[k])/K)
    }
    before <- if (early_stop && f==K) function(k)
      pacc_tail_certificate(partial,k,incumbent,K,cfg$implementation$tie_tol) else NULL
    prediction_before <- prediction_seconds
    pp <- pacc_path(preps[[f]],settings[j,],cfg,pulse,L,engine,before,after)
    path_seconds <- pp$elapsed
    if (settings$family[j]!="enet") path_seconds <- max(0,path_seconds-(prediction_seconds-prediction_before))
    if (settings$family[j]=="enet") for (k in seq_len(L)) after(k,pp$beta[,k],pp$diagnostics[k,])
    paths[[j]][[f]] <- pp; times[j,f] <- path_seconds; losses[[j]][,f] <- loss
    records[[length(records)+1L]] <- cbind(fold=f,setting=j,family=settings$family[j],
      shape=settings$shape[j],pp$diagnostics,loss=loss)
    attempts[[length(attempts)+1L]] <- cbind(phase="cv",fold=f,setting=j,
      family=settings$family[j],shape=settings$shape[j],pp$attempts)
    timings[[length(timings)+1L]] <- cbind(phase="cv",fold=f,setting=j,family=settings$family[j],as.data.frame(pp$acceleration))
    if (!is.null(pp$certificate) && isTRUE(pp$certificate$accepted)) certificates[[length(certificates)+1L]] <-
      data.frame(fold=f,setting=j,index=seq.int(pp$certificate$first_index,L),
        lower_bound=pp$certificate$lower_bound,incumbent=pp$certificate$incumbent,
        roundoff_buffer=pp$certificate$roundoff_buffer)
  }
  cv_records <- do.call(rbind,records)
  for (j in seq_len(J)) { loss <- rowMeans(losses[[j]]); loss[!is.finite(loss)] <- Inf; cv[[j]] <- list(loss=loss) }
  selections <- setNames(lapply(cfg$methods,function(m)engine$perf_choose(cv,settings,m,cfg)),cfg$methods)
  # Keep the legacy partial-search indicator; add its cause, never silently
  # relabel uncomputed points as converged or as numerical failures.
  for (m in names(selections)) {
    z <- cv_records[cv_records$setting %in% selections[[m]]$eligible,]
    selections[[m]]$certified_pruned <- any(z$status=="certified_cv_pruned")
    selections[[m]]$numerically_partial <- any(!is.finite(z$loss) & z$status!="certified_cv_pruned")
  }
  fullpaths <- vector("list",J)
  for (j in seq_len(J)) {
    needed <- Filter(function(s)isTRUE(s$ok) && s$setting==j,selections)
    if (length(needed)) {
      last <- max(vapply(needed,function(s)s$index,numeric(1)))
      pp <- pacc_path(fullprep,settings[j,],cfg,pulse,last,engine); fullpaths[[j]] <- pp
      attempts[[length(attempts)+1L]] <- cbind(phase="full",fold=0,setting=j,
        family=settings$family[j],shape=settings$shape[j],pp$attempts)
      timings[[length(timings)+1L]] <- cbind(phase="full",fold=0,setting=j,family=settings$family[j],as.data.frame(pp$acceleration))
    }
  }
  td <- proc.time()[[3]]; models <- list()
  for (method in cfg$methods) {
    ss <- selections[[method]]
    model <- list(selection=ss,status="no_valid_cv",ok=FALSE,tuning_seconds=sum(times[ss$eligible,,drop=FALSE]))
    if (ss$ok) {
      j <- ss$setting; k <- ss$index; pp <- fullpaths[[j]]
      model$setting <- settings[j,]; model$diagnostics <- pp$diagnostics[k,]; model$refit_seconds <- pp$elapsed
      model$cv_curvature <- lapply(seq_len(K),function(f)engine$perf_hessian(paths[[j]][[f]]$beta[,k],
        preps[[f]],paths[[j]][[f]]$diagnostics$lambda[k],settings[j,],cfg))
      if (pp$diagnostics$ok[k]) {
        model$curvature <- engine$perf_hessian(pp$beta[,k],fullprep,pp$diagnostics$lambda[k],settings[j,],cfg)
        negative <- any(vapply(c(model$cv_curvature,list(model$curvature)),function(h)h$status=="negative_curvature",logical(1)))
        model$ok <- !negative; model$status <- if(negative) "negative_curvature" else "verified_stationary"
        model$beta <- pp$beta[,k]/fullprep$sx; model$intercept <- fullprep$my-sum(fullprep$mx*model$beta)
      } else model$status <- pp$diagnostics$status[k]
    }
    models[[method]] <- model
  }
  selected_folds <- lapply(selections,function(ss)if(!ss$ok) NULL else lapply(seq_len(K),function(f)
    list(beta=paths[[ss$setting]][[f]]$beta[,ss$index],mx=preps[[f]]$mx,sx=preps[[f]]$sx,my=preps[[f]]$my)))
  list(models=models,cv=cv_records,attempts=do.call(rbind,attempts),selected_folds=selected_folds,folds=folds,
    preprocessing=lapply(c(preps,list(fullprep)),function(x)x[c("mx","my","sx","lambda0","spectral")]),
    timings=list(cv_by_setting_fold=times,full_by_setting=vapply(fullpaths,function(x)if(is.null(x))0 else x$elapsed,numeric(1))),
    shared_computation_note="Method times include full reusable paths; do not sum across methods.",
    acceleration=list(early_stop=early_stop,paths=do.call(rbind,timings),certificates=do.call(rbind,certificates),
      preprocessing_seconds=preprocessing_seconds,prediction_seconds=prediction_seconds,
      selected_diagnostics_seconds=proc.time()[[3]]-td,total_seconds=proc.time()[[3]]-t0))
}
