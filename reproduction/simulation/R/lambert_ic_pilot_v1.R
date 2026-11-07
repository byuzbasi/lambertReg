# Isolated IC reanalysis. Existing engines, CV outputs and manuscript are read-only.
lic_read_config <- function(root) {
  x <- readLines(file.path(root,"protocol/lambert_ic_pilot_v1.md"))
  a <- which(x=="```json"); b <- which(seq_along(x)>a & x=="```")[1]
  jsonlite::fromJSON(paste(x[(a+1):(b-1)],collapse="\n"))
}
lic_sources <- function() c("protocol/lambert_ic_pilot_v1.md",
  "protocol/lambert_simulation_main_v2.md","R/penalties.R","R/lambert_family.R",
  "R/performance_models.R","R/performance_run.R","R/lambert_fixed_metrics_v1.R",
  "R/performance_acceleration_v1.R","R/lambert_simulation_main_v2.R",
  "R/lambert_ic_pilot_v1.R","src/lambert_performance_fast_v1.cpp",
  "scripts/lambert_ic_pilot_v1.R","scripts/run_lambert_ic_pilot_v1.sh",
  "scripts/analyze_lambert_ic_pilot_v1.R",
  "tests/lambert_ic_pilot_v1_checks.R")

# Unknown Gaussian variance is profiled at RSS/n. No true variance or test data.
lic_scores <- function(rss,d,n,p,gamma=c(0,.5,1)) {
  stopifnot(length(rss)==length(d),n>1,p>=1,all(d>=0 & d<=p & d==floor(d)))
  score <- outer(n*log(rss/n)+d*log(n),rep(1,length(gamma))) +
    outer(2*lchoose(p,d),gamma)
  score[!is.finite(rss)|rss<=0,] <- Inf
  colnames(score) <- c("BIC","EBIC_0.5","EBIC_1")
  score
}
lic_choose <- function(score,eligible,tol) {
  good <- which(eligible & is.finite(score))
  if(!length(good))return(NA_integer_)
  best <- min(score[good]); good[which(score[good]<=best+tol*(1+abs(best)))[1]]
}

# At joins use a conservative lower limiting active Hessian, never ignore a join.
lic_curvature <- function(b,prep,lambda,setting,cfg,e) {
  h <- e$perf_hessian(b,prep,lambda,setting,cfg)
  if(h$status!="join_unavailable")return(h)
  active <- which(b!=0); t <- abs(b[active]); family <- setting$family
  cuts <- switch(family,lambert=exp(setting$shape)*lambda,
    mcp=setting$shape*lambda,scad=c(lambda,setting$shape*lambda))
  near <- apply(abs(outer(t,cuts,"-")),1,min)<=
    cfg$implementation$curvature_join_tol*max(1,lambda)
  curvature <- e$perf_penalty(t,lambda,setting)$curvature
  curvature[near] <- switch(family,lambert=-.5,mcp=-1/setting$shape,
    scad=-1/(setting$shape-1))
  H <- crossprod(prep$Z[,active,drop=FALSE])/nrow(prep$Z)
  diag(H) <- diag(H)+curvature
  ev <- eigen(H,symmetric=TRUE,only.values=TRUE)$values
  ok <- min(ev)>= -cfg$implementation$curvature_negative_tol*max(1,max(abs(ev)))
  list(status=if(ok)"nonnegative_join_lower" else "join_unresolved",
    min_eigenvalue=min(ev),condition=if(min(ev)>0)max(ev)/min(ev) else Inf)
}
lic_candidates <- function(path,prep,setting,cfg,e) {
  z <- path$diagnostics; L <- nrow(z)
  z$d <- z$rss <- NA_real_; z$curvature_status <- "not_assessed"
  z$min_eigenvalue <- NA_real_; z$size_excluded <- FALSE; z$eligible <- FALSE
  for(k in which(z$ok)) {
    b <- path$beta[,k]
    dd <- e$perf_diagnose(b,prep,z$lambda[k],setting,cfg)
    stopifnot(dd$ok,abs(dd$objective-z$objective[k])<=1e-10*(1+abs(dd$objective)))
    z$d[k] <- sum(abs(b/prep$sx)>cfg$support_reporting_tolerance)
    z$rss[k] <- sum((prep$y-drop(prep$Z%*%b))^2)
    z$size_excluded[k] <- z$d[k]>floor(nrow(prep$Z)/log(nrow(prep$Z)))
    if(!z$size_excluded[k]) {
      h <- lic_curvature(b,prep,z$lambda[k],setting,cfg,e)
      z$curvature_status[k] <- h$status; z$min_eigenvalue[k] <- h$min_eigenvalue
      z$eligible[k] <- h$status %in% c("empty_active","nonnegative_local","nonnegative_join_lower") &&
        is.finite(z$rss[k]) && z$rss[k]>0
    }
  }
  scores <- matrix(Inf,L,3,dimnames=list(NULL,c("BIC","EBIC_0.5","EBIC_1")))
  ids <- which(z$ok)
  if(length(ids))scores[ids,] <- lic_scores(z$rss[ids],z$d[ids],nrow(prep$Z),ncol(prep$Z),cfg$ic$gamma)
  cbind(z,as.data.frame(scores))
}
lic_metrics <- function(b,prep,test,truth,cfg,e) {
  beta <- b/prep$sx; intercept <- prep$my-sum(prep$mx*beta)
  pred <- intercept+drop(test$X%*%beta); active <- which(abs(beta)>cfg$support_reporting_tolerance)
  delta <- beta-truth$beta; tp <- sum(active %in% cfg$support_indices); fp <- length(active)-tp
  fn <- cfg$support_size-tp; mse <- mean((test$y-pred)^2)
  z <- data.frame(ok=TRUE,p=length(beta),tp=tp,fp=fp,fn=fn,
    selected_size=length(active),sigma2=truth$sigma2,
    mse=mse,normalized_mse=mse/truth$sigma2,coefficient_error=sum(delta^2),
    analytic_risk=truth$sigma2+intercept^2+drop(crossprod(delta,truth$Sigma%*%delta)),
    tpr=tp/cfg$support_size,fdr=fp/max(1,length(active)),
    exact_support=setequal(active,cfg$support_indices))
  e$lfb_add_metrics(z)
}
lic_metric_fields <- function() c("mse","normalized_mse","coefficient_error","excess_risk",
  "selected_size","tp","fp","fn","tn","tpr","fdr","mcc","exact_support","sigma2","mcc_zero_denominator")
lic_cv_check <- function(path,cv,prep,test,truth,cfg,e) {
  k <- cv$index
  stopifnot(nrow(cv)==1L,cv$ok,length(k)==1L,k>=1,k<=nrow(path$diagnostics),path$diagnostics$ok[k])
  met <- lic_metrics(path$beta[,k],prep,test,truth,cfg,e)
  actual <- c(lambda=path$diagnostics$lambda[k],objective=path$diagnostics$objective[k],
    unlist(met[lic_metric_fields()]))
  expected <- c(lambda=cv$lambda,objective=cv$objective,unlist(cv[lic_metric_fields()]))
  exact <- names(actual) %in% c("selected_size","tp","fp","fn","tn","exact_support","mcc_zero_denominator")
  error <- abs(actual-expected)/(1+abs(expected))
  passed <- is.finite(error) & ifelse(exact,actual==expected,error<=cfg$ic$reproduction_tolerance)
  data.frame(field=names(actual),actual=as.numeric(actual),expected=as.numeric(expected),
    relative_error=as.numeric(error),passed=as.logical(passed),row.names=NULL)
}
lic_rows <- function(task,cfg,paths,candidates,cv,prep,test,truth,e) {
  rows <- list()
  for(m in cfg$methods) {
    old <- cv[cv$method==m,,drop=FALSE]; path <- paths[[m]]; cc <- candidates[[m]]
    for(criterion in cfg$ic$criteria) {
      iscv <- criterion=="CV"
      k <- if(iscv)old$index else lic_choose(cc[[criterion]],cc$eligible,cfg$implementation$tie_tol)
      numerical <- sum(!cc$ok | (cc$ok & !cc$size_excluded & !cc$eligible))
      row <- data.frame(unit=task$unit,source_unit=task$source_unit,scenario=task$scenario,
        repetition=task$repetition,n=task$n,p=task$p,rho=task$rho,snr=task$snr,profile=task$profile,
        method=m,criterion=criterion,ok=!is.na(k),status=if(is.na(k))"no_valid_ic" else "verified_stationary",
        partial_search=if(iscv)old$partial_search else numerical>0,
        unavailable_numerical=if(iscv)NA_integer_ else numerical,
        excluded_size=if(iscv)0L else sum(cc$size_excluded),
        valid_candidates=if(iscv)NA_integer_ else sum(cc$eligible),
        index=k,lambda=if(is.na(k))NA_real_ else path$diagnostics$lambda[k],
        criterion_score=if(iscv||is.na(k))NA_real_ else cc[[criterion]][k],
        cv_loss=if(iscv)old$cv_loss else NA_real_,
        kkt=if(is.na(k))NA_real_ else path$diagnostics$kkt[k],
        curvature_status=if(iscv)old$curvature_status else if(is.na(k))"unavailable" else cc$curvature_status[k])
      for(f in lic_metric_fields())row[[f]] <- NA_real_
      if(iscv) {
        # Preserve the saved CV metrics verbatim; reproduction is an audit only.
        row[lic_metric_fields()] <- old[lic_metric_fields()]
        row$lambda <- old$lambda; row$kkt <- old$kkt
      } else if(!is.na(k)) row[lic_metric_fields()] <-
        lic_metrics(path$beta[,k],prep,test,truth,cfg,e)[lic_metric_fields()]
      rows[[length(rows)+1L]] <- row
    }
  }
  ans <- do.call(rbind,rows); rownames(ans) <- NULL; ans
}

lic_environment <- function(root,stage="pilot") {
  stage <- match.arg(stage,c("smoke","pilot")); e <- new.env(parent=environment(lic_environment))
  for(f in c("R/penalties.R","R/lambert_family.R","R/performance_models.R","R/performance_run.R",
    "R/lambert_fixed_metrics_v1.R","R/performance_acceleration_v1.R","R/lambert_simulation_main_v2.R"))
    sys.source(file.path(root,f),e)
  # Validation uses the existing scalar R formulas and never fits a model.
  e$lambert_values_cpp <- function(t,lambda,shape)list(
    penalty=e$lambert_penalty(t,lambda,shape),slope=e$lambert_slope(t,lambda,shape),
    curvature=e$lambert_curvature(t,lambda,shape))
  ic <- lic_read_config(root); input <- file.path(root,ic$baseline)
  stopifnot(e$perf_hash(file.path(input,"manifest.csv"))==ic$baseline_manifest_sha256)
  e$perf_verify_manifest(input)
  stopifnot(readLines(file.path(input,"COMPLETED"))==ic$baseline_manifest_sha256)
  archived_sources <- readRDS(file.path(input,"identity.rds"))$sources
  used <- intersect(lic_sources(),archived_sources$file)
  expected_sources <- archived_sources[match(used,archived_sources$file),]
  stopifnot(all(vapply(file.path(root,used),e$perf_hash,"")==expected_sources$sha256))
  cfg <- e$lms_config(root,"main","fast"); all_tasks <- e$lms_tasks(cfg)
  archived_tasks <- read.csv(file.path(input,"tasks.csv"))
  stopifnot(isTRUE(all.equal(all_tasks,archived_tasks,check.attributes=FALSE)))
  tasks <- all_tasks[all_tasks$profile==ic$profile & all_tasks$snr %in% ic$snr &
    all_tasks$p %in% ic$p & all_tasks$rho %in% ic$rho & all_tasks$repetition<=ic$repetitions,,drop=FALSE]
  tasks$source_unit <- tasks$unit; tasks$unit <- seq_len(nrow(tasks)); rownames(tasks) <- NULL
  stopifnot(nrow(tasks)==160L,length(unique(tasks$scenario))==16L,all(table(tasks$scenario)==10))
  cv <- read.csv(file.path(input,"evaluation.csv"))
  cv <- cv[cv$unit %in% tasks$source_unit & cv$method %in% ic$methods,,drop=FALSE]
  stopifnot(nrow(cv)==480L,all(cv$ok),!anyDuplicated(paste(cv$unit,cv$method)))
  cfg$stage <- stage; cfg$ic <- ic; cfg$methods <- ic$methods
  if(stage=="smoke") {
    tasks <- tasks[c(1,nrow(tasks)),,drop=FALSE]; tasks$unit <- 1:2
    tasks$n <- ic$smoke_n; tasks$p <- ic$smoke_p; tasks$stream_index <- ic$smoke_stream_offset+1:2
    tasks$source_unit <- 0L; rownames(tasks) <- NULL # zero denotes synthetic smoke only
    cfg$n_train <- ic$smoke_n; cfg$n_test <- ic$smoke_n_test; cfg$lambda_points <- ic$smoke_lambda_points
  }
  files <- c("manifest.csv","COMPLETED","evaluation.csv","tasks.csv","identity.rds")
  cfg$inputs <- data.frame(file=file.path(ic$baseline,files),
    bytes=as.numeric(file.info(file.path(input,files))$size),
    sha256=unname(vapply(file.path(input,files),e$perf_hash,"")))
  e$perf_config <- function(root,stage) { stopifnot(stage==cfg$stage); cfg }
  e$perf_tasks <- function(cfg) tasks
  e$perf_source_paths <- lic_sources
  original_software <- e$perf_software
  e$perf_software <- function() {
    z <- original_software(); z$packages$mltools <- as.character(utils::packageVersion("mltools")); z
  }
  original_preflight <- e$perf_preflight
  e$perf_preflight <- function(compile=TRUE,root) {
    original_preflight(FALSE,root)
    for(pkg in names(wanted <- c(glmnet="5.0",ncvreg="3.16.0",mltools="0.3.5")))
      if(!requireNamespace(pkg,quietly=TRUE)||as.character(utils::packageVersion(pkg))!=wanted[[pkg]])
        stop("Required package version: ",pkg," ",wanted[[pkg]])
    if(compile) {
      Rcpp::sourceCpp(file.path(root,"src/lambert_performance_fast_v1.cpp"),env=e,
        cacheDir=file.path(tempdir(),"lic_v1_rcpp"),showOutput=FALSE)
      e$lambert_cd_cpp <- e$lambert_cd_fast_v1; e$lambert_values_cpp <- e$lambert_values_fast_v1
    }
    invisible(TRUE)
  }
  e$reference_path <- e$perf_path
  e$perf_path <- function(prep,setting,cfg,pulse=function(...)NULL,last=cfg$lambda_points)
    e$pacc_path(prep,setting,cfg,pulse,last,e)
  e$perf_settings <- function(cfg) data.frame(family=c("lambert","mcp","scad"),shape=c(1,3,3.7))
  e$perf_train_models <- function(X,y,folds,cfg,pulse=function(...)NULL)
    e$pacc_train(X,y,folds,cfg,pulse,e,FALSE)
  smoke_cv <- e$perf_run_task
  e$perf_run_task <- function(task,cfg,seeds,pulse=function(...)NULL) {
    began <- proc.time()[[3]]; data <- e$perf_training(task,cfg,seeds)
    prep <- e$perf_prepare(data$X,data$y,cfg)
    old <- if(stage=="pilot")cv[cv$unit==task$source_unit,,drop=FALSE] else
      e$lfb_add_metrics(smoke_cv(task,cfg,seeds,pulse)$rows)
    settings <- e$perf_settings(cfg); paths <- candidates <- list()
    for(j in seq_along(cfg$methods)) {
      m <- cfg$methods[j]; pulse(paste("full-training path",m))
      paths[[m]] <- e$perf_path(prep,settings[j,],cfg,pulse)
      candidates[[m]] <- lic_candidates(paths[[m]],prep,settings[j,],cfg,e)
    }
    # Test outcomes become available only after all candidate scores are frozen.
    test <- e$perf_test_data(task,cfg,seeds,data$truth)
    checks <- do.call(rbind,lapply(cfg$methods,function(m)cbind(method=m,
      lic_cv_check(paths[[m]],old[old$method==m,,drop=FALSE],prep,test,data$truth,cfg,e))))
    if(!all(checks$passed)) {
      failure <- structure(list(message="Archived CV reproduction failed; no IC comparison accepted",
        call=NULL,evidence=list(task=task,paths=paths,archived_cv=old,reproduction=checks)),
        class=c("lic_reproduction_error","error","condition"))
      stop(failure)
    }
    list(task=task,seeds=seeds,paths=paths,candidates=candidates,archived_cv=old,
      reproduction=checks,data_hashes=list(training=data$hashes,test=test$hashes),
      rows=lic_rows(task,cfg,paths,candidates,old,prep,test,data$truth,e),
      elapsed=proc.time()[[3]]-began)
  }
  e$perf_validate_unit <- function(result,task,cfg,seeds,signature) {
    stopifnot(identical(result$signature,signature),identical(result$task,task),
      identical(result$seeds,seeds),identical(names(result$paths),cfg$methods))
    data <- e$perf_training(task,cfg,seeds); prep <- e$perf_prepare(data$X,data$y,cfg)
    test <- e$perf_test_data(task,cfg,seeds,data$truth)
    stopifnot(identical(result$data_hashes,list(training=data$hashes,test=test$hashes)))
    old <- if(stage=="pilot")cv[cv$unit==task$source_unit,,drop=FALSE] else result$archived_cv
    stopifnot(identical(old,result$archived_cv))
    cc <- list(); settings <- e$perf_settings(cfg)
    for(j in seq_along(cfg$methods)) {
      m <- cfg$methods[j]; pp <- result$paths[[m]]
      stopifnot(identical(dim(pp$beta),as.integer(c(task$p,cfg$lambda_points))),
        identical(pp$diagnostics$index,seq_len(cfg$lambda_points)),
        isTRUE(all.equal(pp$diagnostics$lambda,e$perf_lambdas(prep,settings[j,],cfg),tolerance=1e-14)))
      cc[[m]] <- lic_candidates(pp,prep,settings[j,],cfg,e)
      stopifnot(isTRUE(all.equal(cc[[m]],result$candidates[[m]],tolerance=1e-12)))
      check <- lic_cv_check(pp,old[old$method==m,,drop=FALSE],prep,test,data$truth,cfg,e)
      stopifnot(all(check$passed),isTRUE(all.equal(check,
        result$reproduction[result$reproduction$method==m,setdiff(names(result$reproduction),"method")],
        tolerance=1e-12,check.attributes=FALSE)))
    }
    expected <- lic_rows(task,cfg,result$paths,cc,old,prep,test,data$truth,e)
    stopifnot(nrow(result$rows)==12L,isTRUE(all.equal(result$rows,expected,tolerance=1e-12)),
      !anyDuplicated(paste(result$rows$method,result$rows$criterion)),
      all(result$rows$selected_size[result$rows$ok & result$rows$criterion!="CV"]<=floor(task$n/log(task$n))))
    invisible(TRUE)
  }
  base_validate <- e$perf_validate
  e$perf_validate <- function(out,root,complete=TRUE,quiet=FALSE) {
    z <- base_validate(out,root,complete,TRUE)
    if(complete && !quiet)cat(sprintf("IC_PILOT_COMPLETE: %d validated datasets, %d method/criterion rows; %d failed selections; CV preserved.\n",
      nrow(z$tasks),nrow(z$rows),sum(!z$rows$ok)))
    invisible(z)
  }
  e$perf_progress <- function(out,state)lic_progress(out,state,e)
  e
}

lic_progress <- function(out,state,e) {
  z <- new.env(parent=emptyenv()); z$completed <- state$completed; z$phase <- "preflight"
  z$unit <- NA_integer_; z$running <- z$failed <- 0L; z$began <- proc.time()[[3]]; z$last <- -Inf
  z$durations <- numeric(); z$cost_class <- character()
  z$last_completed <- if(length(z$completed))tail(z$completed,1) else NA_integer_
  before <- length(z$completed)
  z$emit <- function(detail="",force=FALSE) {
    now <- proc.time()[[3]]; if(!force && now-z$last<state$cfg$ic$progress_seconds)return(invisible(NULL))
    z$last <- now; elapsed <- now-z$began; pending <- setdiff(state$tasks$unit,z$completed)
    cls <- with(state$tasks,paste(p,rho,snr,profile,sep="/")); means <- tapply(z$durations,z$cost_class,mean)
    eta <- if(!length(pending))0 else if(length(means)&&all(cls[pending]%in%names(means)))sum(means[cls[pending]]) else NA_real_
    # Sequential ETA counts the active dataset conservatively as a whole task.
    stamp <- format(Sys.time(),"%Y-%m-%dT%H:%M:%SZ",tz="UTC")
    snapshot <- list(run_id=basename(out),slurm_job_id=NA_character_,phase=z$phase,
      work_unit="one regenerated dataset: three full training paths, IC scores and archived CV audit",
      total=nrow(state$tasks),validated_completed=length(z$completed),validated_before_invocation=before,
      completed_this_invocation=length(z$completed)-before,running=z$running,failed=z$failed,
      pending=nrow(state$tasks)-length(z$completed)-z$running-z$failed,
      percent_completed=100*length(z$completed)/nrow(state$tasks),elapsed_seconds=elapsed,
      throughput_units_per_second=if(elapsed>0)(length(z$completed)-before)/elapsed else NA_real_,
      estimated_remaining_seconds=eta,eta_status=if(is.finite(eta))"estimated" else "estimating_or_unknown",
      estimated_completion_utc=if(is.finite(eta))format(Sys.time()+eta,"%Y-%m-%dT%H:%M:%SZ",tz="UTC") else NA_character_,
      eta_scope="current-invocation scenario timings; path computation and scoring only; per-unit and final validation excluded",
      slurm_remaining_seconds=NA_real_,last_heartbeat_utc=stamp,current_work_unit=z$unit,
      most_recently_completed_unit=z$last_completed,detail=detail)
    e$perf_atomic_json(snapshot,file.path(out,"progress.json"))
    history <- file.path(out,"progress.tsv")
    write.table(as.data.frame(snapshot),history,sep="\t",row.names=FALSE,col.names=!file.exists(history),append=file.exists(history),quote=TRUE)
    cat(sprintf("[%s] %s | %d/%d | running %d failed %d | ETA %s | %s\n",stamp,z$phase,
      length(z$completed),nrow(state$tasks),z$running,z$failed,
      if(is.finite(eta))paste(round(eta),"s") else "unknown",detail)); flush(stdout())
    invisible(snapshot)
  }
  z
}

lic_main <- function(root,args) {
  opt <- list(stage="pilot",action="design",out=NULL,stop_after=Inf)
  for(arg in args) {
    if(!grepl("^--[a-z-]+=.+$",arg))stop("Use --key=value arguments")
    key <- gsub("-","_",sub("=.*$","",substring(arg,3)))
    if(!key %in% names(opt))stop("Unknown option: ",key)
    opt[[key]] <- sub("^[^=]*=","",arg)
  }
  stopifnot(!is.null(opt$out),opt$stage %in% c("pilot","smoke"),
    opt$action %in% c("start","resume","validate","design"))
  if(nzchar(Sys.getenv("SLURM_JOB_ID")))stop("This approved pilot is local only")
  out <- if(startsWith(opt$out,"/"))opt$out else file.path(root,opt$out)
  if(opt$action %in% c("start","design") && file.exists(out))stop("Refusing overwrite: ",out)
  e <- lic_environment(root,opt$stage)
  if(opt$stage=="pilot" && opt$action %in% c("start","resume")) {
    gate_dir <- file.path(root,"outputs/validation/lambert_ic_pilot_v1_a01")
    if(!file.exists(file.path(gate_dir,"COMPLETED")))stop("Local smoke/replay validation gate missing")
    e$perf_verify_manifest(gate_dir)
    stopifnot(readLines(file.path(gate_dir,"COMPLETED"))==e$perf_hash(file.path(gate_dir,"manifest.csv")))
    gate <- readRDS(file.path(gate_dir,"gate.rds"))
    stopifnot(isTRUE(gate$passed),identical(gate$sources,e$perf_sources(root)),
      identical(gate$software,e$perf_software()),identical(gate$inputs,e$perf_config(root,"pilot")$inputs))
  }
  mode <- if(opt$action %in% c("design","validate"))opt$action else opt$stage
  argv <- c("--mode",mode,"--stage",opt$stage,"--out",out,"--stop-after",as.character(opt$stop_after))
  if(opt$action=="resume")argv <- c(argv,"--resume")
  e$commandArgs <- function(trailingOnly=FALSE)argv
  tryCatch(e$perf_main(root),lic_reproduction_error=function(err) {
    if(dir.exists(out))saveRDS(err$evidence,file.path(out,paste0("reproduction_failure_",Sys.getpid(),".rds")))
    stop(err)
  })
}

# Analysis only: paired summaries never substitute a complete-case mean for an
# unconditional result. Exclusions and original CV partial searches stay visible.
lic_analyze <- function(root,input,out,stage="pilot") {
  stopifnot(!file.exists(out)); e <- lic_environment(root,stage)
  e$perf_preflight(FALSE,root); state <- e$perf_validate(input,root,TRUE)
  rows <- state$rows; metrics <- c("normalized_mse","mcc","tpr","fdr",
    "exact_support","selected_size","coefficient_error","excess_risk")
  summarize <- function(v,expected) {
    v <- as.numeric(v); good <- is.finite(v); n <- sum(good); complete <- n==expected
    data.frame(expected=expected,available=n,complete=complete,
      mean=if(complete)mean(v) else NA_real_,mcse=if(complete&&n>1)sd(v)/sqrt(n) else NA_real_,
      conditional_mean=if(!complete&&n)mean(v[good]) else NA_real_)
  }
  summaries <- pairs <- list()
  for(s in unique(rows$scenario)) {
    z <- rows[rows$scenario==s,,drop=FALSE]; reps <- sort(unique(z$repetition))
    for(m in state$cfg$methods)for(criterion in state$cfg$ic$criteria) {
      r <- z[z$method==m & z$criterion==criterion,,drop=FALSE]
      for(metric in metrics)summaries[[length(summaries)+1L]] <- cbind(
        data.frame(scenario=s,method=m,criterion=criterion,metric=metric),
        summarize(r[[metric]],length(reps)))
    }
    add_pair <- function(a,b,label) {
      a <- a[order(a$repetition),]; b <- b[order(b$repetition),]
      stopifnot(identical(a$repetition,b$repetition))
      for(scope in c("all_pairs","complete_search_pairs"))for(metric in metrics) {
        valid <- a$ok & b$ok
        if(scope=="complete_search_pairs")valid <- valid & !a$partial_search & !b$partial_search
        delta <- as.numeric(a[[metric]])-as.numeric(b[[metric]])
        delta[!valid] <- NA_real_
        pairs[[length(pairs)+1L]] <<- cbind(data.frame(scenario=s,contrast=label,
          metric=metric,scope=scope,direction="first_minus_second"),summarize(delta,length(reps)))
      }
    }
    for(m in state$cfg$methods)for(criterion in state$cfg$ic$criteria[-1])
      add_pair(z[z$method==m & z$criterion==criterion,],z[z$method==m & z$criterion=="CV",],
        paste(m,criterion,"minus",m,"CV"))
    for(criterion in state$cfg$ic$criteria)for(other in c("mcp","scad"))
      add_pair(z[z$method=="lambert_c1" & z$criterion==criterion,],z[z$method==other & z$criterion==criterion,],
        paste("lambert_c1",criterion,"minus",other,criterion))
  }
  summary <- do.call(rbind,summaries); paired <- do.call(rbind,pairs)
  stopifnot(dir.create(out,recursive=TRUE))
  write.csv(summary,file.path(out,"summary.csv"),row.names=FALSE)
  write.csv(paired,file.path(out,"paired.csv"),row.names=FALSE)
  write.csv(rows,file.path(out,"selections.csv"),row.names=FALSE)
  cells <- unique(rows[c("scenario","n","p","rho","snr","profile")])
  write.csv(cells,file.path(out,"scenarios.csv"),row.names=FALSE)
  numerical <- aggregate(cbind(failed=as.integer(!rows$ok),partial=as.integer(rows$partial_search)),
    rows[c("method","criterion")],sum)
  write.csv(numerical,file.path(out,"numerical_summary.csv"),row.names=FALSE)
  for(metric in metrics) {
    z <- merge(cells,summary[summary$metric==metric,],by="scenario",sort=FALSE)
    z <- z[order(z$p,z$rho,z$snr,z$method,z$criterion),]
    tab <- data.frame(p=z$p,rho=z$rho,SNR=z$snr,Method=z$method,Criterion=z$criterion,
      Estimate=ifelse(z$complete,sprintf("%.4f (%.4f)",z$mean,z$mcse),"Incomplete"),
      Available=paste0(z$available,"/",z$expected))
    latex <- knitr::kable(tab,format="latex",booktabs=TRUE,longtable=TRUE,
      caption=paste(metric,"in the exploratory IC pilot; arithmetic averages (MCSE)."))
    latex <- kableExtra::column_spec(latex,4,width="2cm")
    latex <- kableExtra::collapse_rows(latex,columns=1:4,latex_hline="major",valign="top")
    writeLines(as.character(latex),file.path(out,paste0(metric,".tex")))
  }
  writeLines(c("Exploratory IC pilot; original CV retained; no new fits in analysis.",
    "BIC/EBIC use penalized RSS and a shared size cap. Pilot results do not establish oracle behavior.",
    "Paired differences are first minus second. Negative NMSE/FDR and positive MCC/TPR favor the first entry.",
    "Incomplete unconditional summaries are NA; conditional means and counts are separate.",
    capture.output(knitr::kable(numerical))),file.path(out,"report.txt"))
  e$perf_atomic_json(list(input_manifest_sha256=e$perf_hash(file.path(input,"manifest.csv")),
    new_model_fits=0L,source_hashes=e$perf_sources(root)),file.path(out,"provenance.json"))
  e$perf_manifest(out); e$perf_verify_manifest(out)
  writeLines(e$perf_hash(file.path(out,"manifest.csv")),file.path(out,"COMPLETED"))
  cat("IC_ANALYSIS_COMPLETE: means, MCSE, paired contrasts and numerical coverage; no new fits.\n")
  invisible(summary)
}
