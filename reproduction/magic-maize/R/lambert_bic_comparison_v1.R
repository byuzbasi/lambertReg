# Versioned adapter: historical estimators, CV outputs and pilot remain frozen.
lbc_spec <- function(root) {
  x <- readLines(file.path(root,"protocol/lambert_bic_comparison_v1.md"))
  a <- which(x=="```json"); b <- which(seq_along(x)>a & x=="```")[1]
  jsonlite::fromJSON(paste(x[(a+1):(b-1)],collapse="\n"))
}
lbc_sources <- function() c("protocol/lambert_bic_comparison_v1.md",
  "protocol/lambert_simulation_main_v2.md","R/penalties.R","R/lambert_family.R",
  "R/performance_models.R","R/performance_run.R","R/performance_truba_v1.R",
  "R/lambert_fixed_metrics_v1.R","R/performance_acceleration_v1.R",
  "R/lambert_simulation_main_v2.R","R/lambert_ic_pilot_v1.R",
  "src/lambert_performance_fast_v1.cpp","R/lambert_bic_comparison_v1.R",
  "R/lambert_bic_execution_v1.R","scripts/lambert_bic_comparison_v1.R",
  "scripts/analyze_lambert_bic_comparison_v1.R","scripts/run_lambert_bic_comparison_v1.sh",
  "tests/lambert_bic_comparison_v1_checks.R")

lbc_df <- function(Z,active,lambda,alpha,cache=NULL) {
  stopifnot(is.finite(lambda),lambda>0,alpha>=0,alpha<1)
  if(!length(active))return(0)
  key <- paste(active,collapse=",")
  ev <- if(!is.null(cache))get0(key,envir=cache,inherits=FALSE) else NULL
  if(is.null(ev)) {
    A <- Z[,active,drop=FALSE]
    G <- if(ncol(A)<=nrow(A))crossprod(A)/nrow(A) else tcrossprod(A)/nrow(A)
    ev <- eigen(G,symmetric=TRUE,only.values=TRUE)$values
    stopifnot(min(ev)>= -1e-10*max(1,max(abs(ev))))
    ev <- pmax(ev,0)
    if(!is.null(cache))assign(key,ev,envir=cache)
  }
  sum(ev/(ev+lambda*(1-alpha)))
}
lbc_candidates <- function(path,prep,setting,cfg,e,cache) {
  z <- path$diagnostics
  z$d <- z$q <- z$rss <- z$BIC <- NA_real_
  z$curvature_status <- "not_assessed"; z$cap_excluded <- z$eligible <- FALSE
  z$complexity <- if(setting$family=="enet" && setting$shape<1)"effective_df" else "active_count"
  for(k in which(z$ok)) {
    b <- path$beta[,k]; dd <- e$perf_diagnose(b,prep,z$lambda[k],setting,cfg)
    stopifnot(dd$ok,abs(dd$objective-z$objective[k])<=1e-10*(1+abs(dd$objective)))
    z$d[k] <- sum(abs(b/prep$sx)>cfg$support_reporting_tolerance)
    z$q[k] <- if(z$complexity[k]=="effective_df")
      lbc_df(prep$Z,if(setting$shape==0)seq_along(b) else which(b!=0),
        z$lambda[k],setting$shape,cache) else z$d[k]
    z$rss[k] <- sum((prep$y-drop(prep$Z%*%b))^2)
    z$cap_excluded[k] <- z$q[k]>floor(nrow(prep$Z)/log(nrow(prep$Z)))
    if(is.finite(z$rss[k]) && z$rss[k]>0)
      z$BIC[k] <- nrow(prep$Z)*log(z$rss[k]/nrow(prep$Z))+z$q[k]*log(nrow(prep$Z))
    if(!z$cap_excluded[k]) {
      z$curvature_status[k] <- if(setting$family=="enet")"convex_psd" else
        e$lic_curvature(b,prep,z$lambda[k],setting,cfg,e)$status
      z$eligible[k] <- z$curvature_status[k] %in%
        c("convex_psd","empty_active","nonnegative_local","nonnegative_join_lower") && is.finite(z$BIC[k])
    }
  }
  z
}
lbc_settings_for <- function(settings,method) switch(method,
  lambert_c1=which(settings$family=="lambert"),mcp=which(settings$family=="mcp"),
  scad=which(settings$family=="scad"),lasso=which(settings$family=="enet" & settings$shape==1),
  ridge=which(settings$family=="enet" & settings$shape==0),
  elastic_net=which(settings$family=="enet"))
lbc_select <- function(candidates,settings,method,cfg,e) {
  ids <- lbc_settings_for(settings,method)
  ids <- ids[order(-settings$shape[ids])]
  z <- do.call(rbind,lapply(ids,function(j)cbind(setting=j,candidates[[j]])))
  k <- e$lic_choose(z$BIC,z$eligible,cfg$implementation$tie_tol)
  list(setting=if(is.na(k))NA_integer_ else z$setting[k],
    index=if(is.na(k))NA_integer_ else z$index[k],
    partial=any(!z$ok | (z$ok & !z$cap_excluded & !z$eligible)),
    unavailable=sum(!z$ok | (z$ok & !z$cap_excluded & !z$eligible)),
    excluded=sum(z$cap_excluded),valid=sum(z$eligible))
}
lbc_oracle <- function(prep,cfg) {
  A <- cfg$support_indices; Z <- prep$Z[,A,drop=FALSE]; qrZ <- qr(Z)
  b <- rep(0,ncol(prep$Z)); ok <- qrZ$rank==length(A)
  if(ok)b[A] <- qr.coef(qrZ,prep$y) else b[] <- NA_real_
  kkt <- if(ok)max(abs(crossprod(Z,drop(Z%*%b[A])-prep$y)/nrow(Z)))/(1+max(abs(prep$score))) else Inf
  list(beta=b,rank=qrZ$rank,kkt=kkt,ok=ok && all(is.finite(b)) && kkt<=cfg$normalized_kkt_tolerance)
}
lbc_oracle_check <- function(oracle,prep,cfg) {
  A <- cfg$support_indices; Z <- prep$Z[,A,drop=FALSE]; rank <- qr(Z)$rank
  stopifnot(oracle$rank==rank,length(oracle$beta)==ncol(prep$Z))
  if(rank<length(A)) { stopifnot(!oracle$ok,all(is.na(oracle$beta))); return(invisible(TRUE)) }
  b <- oracle$beta
  kkt <- max(abs(crossprod(Z,drop(Z%*%b[A])-prep$y)/nrow(Z)))/(1+max(abs(prep$score)))
  stopifnot(all(is.finite(b)),all(b[-A]==0),isTRUE(all.equal(kkt,oracle$kkt,tolerance=1e-12)),
    identical(oracle$ok,kkt<=cfg$normalized_kkt_tolerance))
  invisible(TRUE)
}
lbc_rows <- function(task,cfg,paths,candidates,old,oracle,prep,test,truth,e) {
  settings <- e$perf_settings(cfg); rows <- list()
  make <- function(method,criterion) {
    r <- task[c("unit","source_unit","scenario","repetition","n","p","rho","snr","profile")]
    r$method <- method; r$criterion <- criterion; r$ok <- FALSE; r$status <- "no_valid_bic"
    r$partial_search <- FALSE; r$reference_support <- method=="oracle"
    r$setting <- r$index <- NA_integer_; r$family <- r$complexity <- r$curvature_status <- NA_character_
    r$shape <- r$lambda <- r$q <- r$bic <- r$kkt <- r$cv_loss <- NA_real_
    r$unavailable_numerical <- r$excluded_cap <- r$valid_candidates <- 0L
    r$cap_boundary <- r$path_endpoint <- FALSE
    for(f in e$lic_metric_fields())r[[f]] <- NA_real_
    r
  }
  for(m in cfg$methods)for(criterion in c("CV","BIC")) {
    r <- make(m,criterion); cv <- old[old$method==m,,drop=FALSE]; iscv <- criterion=="CV"
    ss <- lbc_select(candidates,settings,m,cfg,e)
    j <- if(iscv && cv$ok)which(settings$family==cv$family & settings$shape==cv$shape) else ss$setting
    k <- if(iscv)cv$index else ss$index
    r$ok <- if(iscv)cv$ok else !is.na(k)
    r$partial_search <- if(iscv)cv$partial_search else ss$partial
    if(!iscv) {
      r$unavailable_numerical <- ss$unavailable; r$excluded_cap <- ss$excluded; r$valid_candidates <- ss$valid
    }
    if(r$ok) {
      stopifnot(length(j)==1L); cc <- candidates[[j]][k,]
      r$setting <- j; r$index <- k; r$family <- settings$family[j]; r$shape <- settings$shape[j]
      r$lambda <- cc$lambda; r$q <- cc$q; r$complexity <- cc$complexity
      r$kkt <- cc$kkt; r$status <- "verified_stationary"; r$curvature_status <- cc$curvature_status
      r$path_endpoint <- k %in% c(1L,cfg$lambda_points)
      if(!iscv) {
        r$bic <- cc$BIC; r$cap_boundary <- abs(cc$q-floor(task$n/log(task$n)))<1e-8
        r[e$lic_metric_fields()] <- e$lic_metrics(paths[[j]]$beta[,k],prep,test,truth,cfg,e)[e$lic_metric_fields()]
      }
    }
    if(iscv) {
      r[e$lic_metric_fields()] <- cv[e$lic_metric_fields()]
      r$status <- cv$status; r$lambda <- cv$lambda; r$kkt <- cv$kkt; r$cv_loss <- cv$cv_loss
      r$curvature_status <- cv$curvature_status
    }
    rows[[length(rows)+1L]] <- r
  }
  r <- make("oracle","Oracle"); r$ok <- oracle$ok; r$q <- cfg$support_size
  r$complexity <- "known_support"; r$kkt <- oracle$kkt
  r$status <- if(oracle$ok)"verified_oracle" else "oracle_rank_or_score_failure"
  if(oracle$ok) {
    z <- e$lic_metrics(oracle$beta,prep,test,truth,cfg,e)
    z$selected_size <- z$tp <- cfg$support_size; z$fp <- z$fn <- 0L
    z$tpr <- 1; z$fdr <- 0; z$exact_support <- TRUE
    z <- e$lfb_add_metrics(z); r[e$lic_metric_fields()] <- z[e$lic_metric_fields()]
  }
  rows[[length(rows)+1L]] <- r
  ans <- do.call(rbind,rows); rownames(ans) <- NULL; ans
}

lbc_environment <- function(root,stage="production") {
  stage <- match.arg(stage,c("production","smoke","replay"))
  e <- new.env(parent=environment(lbc_environment))
  for(f in c("R/penalties.R","R/lambert_family.R","R/performance_models.R","R/performance_run.R",
    "R/performance_truba_v1.R","R/lambert_fixed_metrics_v1.R","R/performance_acceleration_v1.R",
    "R/lambert_simulation_main_v2.R","R/lambert_ic_pilot_v1.R"))sys.source(file.path(root,f),e)
  e$lambert_values_cpp <- function(t,lambda,shape)list(penalty=e$lambert_penalty(t,lambda,shape),
    slope=e$lambert_slope(t,lambda,shape),curvature=e$lambert_curvature(t,lambda,shape))
  spec <- lbc_spec(root); input <- file.path(root,spec$baseline)
  stopifnot(e$perf_hash(file.path(input,"manifest.csv"))==spec$baseline_manifest_sha256)
  e$perf_verify_manifest(input)
  stopifnot(readLines(file.path(input,"COMPLETED"))==spec$baseline_manifest_sha256)
  archived <- readRDS(file.path(input,"identity.rds"))$sources
  used <- intersect(lbc_sources(),archived$file)
  stopifnot(all(vapply(file.path(root,used),e$perf_hash,"")==archived$sha256[match(used,archived$file)]))
  cfg <- e$lms_config(root,"main","fast"); all_tasks <- e$lms_tasks(cfg)
  stopifnot(isTRUE(all.equal(all_tasks,read.csv(file.path(input,"tasks.csv")),check.attributes=FALSE)))
  tasks <- all_tasks[all_tasks$profile==spec$profile & all_tasks$p %in% spec$p &
    all_tasks$rho %in% spec$rho & all_tasks$snr %in% spec$snr & all_tasks$repetition<=spec$repetitions,,drop=FALSE]
  tasks$source_unit <- tasks$unit; tasks$unit <- seq_len(nrow(tasks)); rownames(tasks) <- NULL
  stopifnot(nrow(tasks)==4800L,length(unique(tasks$scenario))==16L,all(table(tasks$scenario)==300L))
  old <- read.csv(file.path(input,"evaluation.csv"))
  old <- old[old$unit %in% tasks$source_unit & old$method %in% spec$methods,,drop=FALSE]
  stopifnot(nrow(old)==28800L,all(old$ok),!anyDuplicated(paste(old$unit,old$method)))
  cfg$stage <- stage; cfg$bic <- spec; cfg$ic <- list(reproduction_tolerance=spec$reproduction_tolerance)
  if(stage=="replay")tasks <- tasks[c(1L,4501L),,drop=FALSE]
  if(stage=="smoke") {
    tasks <- tasks[rep(c(1L,4501L),length.out=spec$smoke_tasks),,drop=FALSE]
    tasks$n <- cfg$n_train <- spec$smoke_n; tasks$p <- rep(spec$smoke_p,length.out=nrow(tasks))
    tasks$stream_index <- spec$smoke_stream_offset+seq_len(nrow(tasks)); tasks$source_unit <- 0L
    cfg$n_test <- spec$smoke_n_test; cfg$lambda_points <- spec$smoke_lambda_points
  }
  tasks$unit <- seq_len(nrow(tasks)); rownames(tasks) <- NULL
  files <- c("manifest.csv","COMPLETED","evaluation.csv","tasks.csv","identity.rds")
  cfg$inputs <- data.frame(file=file.path(spec$baseline,files),bytes=as.numeric(file.info(file.path(input,files))$size),
    sha256=unname(vapply(file.path(input,files),e$perf_hash,"")))
  e$perf_config <- function(root,stage) { stopifnot(stage==cfg$stage); cfg }
  e$perf_tasks <- function(cfg)tasks
  e$perf_source_paths <- lbc_sources
  original_software <- e$perf_software
  e$perf_software <- function() {
    z <- original_software()
    for(p in c("mltools","kableExtra"))z$packages[[p]] <- as.character(utils::packageVersion(p))
    z
  }
  original_preflight <- e$perf_preflight
  e$perf_preflight <- function(compile=TRUE,root) {
    original_preflight(FALSE,root)
    for(p in names(wanted <- c(glmnet="5.0",ncvreg="3.16.0",mltools="0.3.5")))
      if(!requireNamespace(p,quietly=TRUE)||as.character(utils::packageVersion(p))!=wanted[[p]])stop("Required version: ",p," ",wanted[[p]])
    if(!requireNamespace("kableExtra",quietly=TRUE))stop("Required installed package: kableExtra")
    if(compile) {
      Rcpp::sourceCpp(file.path(root,"src/lambert_performance_fast_v1.cpp"),env=e,
        cacheDir=file.path(tempdir(),"lbc_v1_rcpp"),showOutput=FALSE)
      e$lambert_cd_cpp <- e$lambert_cd_fast_v1; e$lambert_values_cpp <- e$lambert_values_fast_v1
    }
    invisible(TRUE)
  }
  e$reference_path <- e$perf_path
  e$perf_path <- function(prep,setting,cfg,pulse=function(...)NULL,last=cfg$lambda_points)
    e$pacc_path(prep,setting,cfg,pulse,last,e)
  e$perf_train_models <- function(X,y,folds,cfg,pulse=function(...)NULL)e$pacc_train(X,y,folds,cfg,pulse,e,FALSE)
  smoke_cv <- e$perf_run_task
  e$perf_run_task <- function(task,cfg,seeds,pulse=function(...)NULL) {
    began <- proc.time()[[3]]; data <- e$perf_training(task,cfg,seeds); prep <- e$perf_prepare(data$X,data$y,cfg)
    cv <- if(stage=="smoke")e$lfb_add_metrics(smoke_cv(task,cfg,seeds,pulse)$rows) else
      old[old$unit==task$source_unit,,drop=FALSE]
    settings <- e$perf_settings(cfg); paths <- candidates <- vector("list",nrow(settings))
    cache <- new.env(parent=emptyenv())
    for(j in seq_len(nrow(settings))) {
      paths[[j]] <- e$perf_path(prep,settings[j,],cfg,pulse)
      candidates[[j]] <- lbc_candidates(paths[[j]],prep,settings[j,],cfg,e,cache)
    }
    oracle <- lbc_oracle(prep,cfg)
    test <- e$perf_test_data(task,cfg,seeds,data$truth)
    checks <- lbc_reproduction(paths,cv,prep,test,data$truth,cfg,e)
    list(task=task,seeds=seeds,paths=paths,candidates=candidates,oracle=oracle,archived_cv=cv,
      reproduction=checks,data_hashes=list(training=data$hashes,test=test$hashes),
      rows=lbc_rows(task,cfg,paths,candidates,cv,oracle,prep,test,data$truth,e),elapsed=proc.time()[[3]]-began)
  }
  e$perf_validate_unit <- function(result,task,cfg,seeds,signature) {
    stopifnot(identical(result$signature,signature),identical(result$task,task),identical(result$seeds,seeds))
    data <- e$perf_training(task,cfg,seeds); prep <- e$perf_prepare(data$X,data$y,cfg)
    test <- e$perf_test_data(task,cfg,seeds,data$truth)
    stopifnot(identical(result$data_hashes,list(training=data$hashes,test=test$hashes)))
    cv <- if(stage=="smoke")result$archived_cv else old[old$unit==task$source_unit,,drop=FALSE]
    stopifnot(identical(cv,result$archived_cv),nrow(cv)==6L)
    settings <- e$perf_settings(cfg); cc <- vector("list",nrow(settings)); cache <- new.env(parent=emptyenv())
    stopifnot(length(result$paths)==nrow(settings))
    for(j in seq_len(nrow(settings))) {
      pp <- result$paths[[j]]
      stopifnot(identical(dim(pp$beta),as.integer(c(task$p,cfg$lambda_points))),
        identical(pp$diagnostics$index,seq_len(cfg$lambda_points)),
        isTRUE(all.equal(pp$diagnostics$lambda,e$perf_lambdas(prep,settings[j,],cfg),tolerance=1e-14)))
      cc[[j]] <- lbc_candidates(pp,prep,settings[j,],cfg,e,cache)
      stopifnot(isTRUE(all.equal(cc[[j]],result$candidates[[j]],tolerance=1e-12)))
    }
    checks <- lbc_reproduction(result$paths,cv,prep,test,data$truth,cfg,e)
    stopifnot(all(checks$passed),isTRUE(all.equal(checks,result$reproduction,tolerance=1e-12)))
    lbc_oracle_check(result$oracle,prep,cfg)
    expected <- lbc_rows(task,cfg,result$paths,cc,cv,result$oracle,prep,test,data$truth,e)
    stopifnot(nrow(result$rows)==13L,isTRUE(all.equal(result$rows,expected,tolerance=1e-12)),
      all(expected$q[expected$ok & expected$criterion=="BIC"]<=floor(task$n/log(task$n))))
    invisible(TRUE)
  }
  # Existing dispatcher, atomic shards and explicit failure preservation.
  e$perf_truba_write_unit <- function(out,result,state) {
    e$accept_started <- as.numeric(Sys.time()); e$accept_unit <- result$task$unit
    e$perf_write_unit(out,result,state)
  }
  e$perf_truba_read_unit <- function(out,i,state) {
    r <- e$perf_read_unit(out,i,state)
    # Scheduling cost includes parent validation/write/readback; stored scientific
    # payload and the worker's measured fit time are never modified.
    if(identical(i,e$accept_unit)) {
      r$elapsed <- r$elapsed+max(0,as.numeric(Sys.time())-e$accept_started)
      e$accept_unit <- NULL
    }
    r
  }
  e$perf_truba_progress_emit <- lbc_progress_emit
  environment(e$perf_truba_progress_emit) <- e
  e
}
lbc_reproduction <- function(paths,cv,prep,test,truth,cfg,e) {
  settings <- e$perf_settings(cfg)
  do.call(rbind,lapply(cfg$methods,function(m) {
    old <- cv[cv$method==m,,drop=FALSE]
    if(!old$ok)return(data.frame(method=m,field="archived_failure_retained",actual=NA_real_,expected=NA_real_,relative_error=0,passed=TRUE))
    j <- which(settings$family==old$family & settings$shape==old$shape)
    stopifnot(length(j)==1L)
    if(!paths[[j]]$diagnostics$ok[old$index])return(data.frame(method=m,field="CV_path_fit_unavailable",
      actual=NA_real_,expected=NA_real_,relative_error=Inf,passed=FALSE))
    cbind(method=m,e$lic_cv_check(paths[[j]],old,prep,test,truth,cfg,e))
  }))
}
