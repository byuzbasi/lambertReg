# Isolated adapter; historical engines, datasets and results are read-only.
rfh_sources <- function() c("protocol/riboflavin_fixed_hd_v2.md",
  "protocol/riboflavin_realdata_model_v1.md", "R/penalties.R", "R/lambert_family.R",
  "R/performance_models.R", "R/performance_acceleration_v1.R",
  "R/lambert_ic_pilot_v1.R", "R/riboflavin_data_v1.R",
  "R/riboflavin_realdata_model_v1.R", "R/checkpoint_packs_v1.R",
  "src/lambert_performance_fast_v1.cpp", "R/riboflavin_fixed_hd_v2.R",
  "R/riboflavin_fixed_hd_execution_v2.R", "scripts/riboflavin_fixed_hd_v2.R",
  "scripts/run_riboflavin_fixed_hd_v2.sh", "tests/riboflavin_fixed_hd_v2_checks.R")
rfh_software <- function() list(R=R.version.string, platform=R.version$platform,
  packages=setNames(vapply(c("Rcpp","glmnet","ncvreg","digest","jsonlite"),
    function(p)as.character(utils::packageVersion(p)),""),
    c("Rcpp","glmnet","ncvreg","digest","jsonlite")))
rfh_files <- function(root,files) data.frame(file=files,
  bytes=as.numeric(file.info(file.path(root,files))$size),
  sha256=unname(vapply(file.path(root,files),cpk_hash,"")))
rfh_environment <- function(root,stage="production",compile=FALSE) {
  stopifnot(stage %in% c("production","smoke"))
  e <- new.env(parent=environment(rfh_environment)); e$root <- root; e$stage <- stage
  for(f in c("penalties.R","lambert_family.R","performance_models.R",
    "performance_acceleration_v1.R","lambert_ic_pilot_v1.R","riboflavin_data_v1.R",
    "riboflavin_realdata_model_v1.R"))sys.source(file.path(root,"R",f),e)
  legacy <- e$ribo_config(root,stage); data <- e$ribo_validate_prepared(root)
  e$view <- e$ribo_stage_view(data,legacy); e$tasks <- e$ribo_tasks(legacy)
  e$seeds <- e$ribo_unit_seeds(legacy,e$tasks)
  keep <- c("stage","folds","raw_relpath","processed_relpath","sample_count","feature_count",
    "outer_repeats","outer_folds","inner_folds","lambda_points","elastic_net_alpha",
    "max_coordinate_sweeps_per_start","normalized_kkt_tolerance","support_reporting_tolerance",
    "fold_master_seed","unit_master_seed","implementation","active_repeats",
    "active_outer_folds","active_feature_count")
  cfg <- legacy[keep]; cfg$protocol <- "riboflavin_fixed_hd_v2"
  cfg$fold_protocol <- "riboflavin_realdata_model_v1"
  cfg$expected_production_units <- 50L; cfg$expected_production_rows <- 500L
  cfg$expected_production_predictions <- 7100L
  cfg$methods <- c("lambert_c1","lasso","ridge","elastic_net","mcp","scad")
  cfg$hbic_methods <- c("lambert_c1","lasso","mcp","scad")
  cfg$lambert_shapes <- 1; cfg$mcp_gamma <- 3; cfg$scad_a <- 3.7
  cfg$workers <- if(stage=="production")5L else 2L
  cfg$early_stop <- FALSE; cfg$hbic_multiplier <- "log(log(n))"
  e$cfg <- cfg; e$settings <- e$perf_settings(cfg)
  e$mapping <- cpk_mapping(e$tasks$unit)
  e$sources <- rfh_files(root,rfh_sources())
  prepared <- file.path(cfg$processed_relpath,c("riboflavin.rds","outer_folds.csv",
    "inner_folds.csv","fold_seed_hashes.csv","manifest.csv","PREPARED"))
  e$inputs <- rfh_files(root,c(cfg$raw_relpath,prepared))
  e$software <- rfh_software()
  e$scientific_signature <- cpk_object_hash(list(cfg=cfg,tasks=e$tasks,seeds=e$seeds,
    data=e$ribo_data_signatures(e$view),inputs=e$inputs))
  e$execution_signature <- cpk_object_hash(list(sources=e$sources,software=e$software))
  e$reference_path <- e$perf_path
  e$lambert_values_cpp <- function(t,lambda,shape)list(penalty=e$lambert_penalty(t,lambda,shape),
    slope=e$lambert_slope(t,lambda,shape),curvature=e$lambert_curvature(t,lambda,shape))
  e$perf_path <- function(prep,setting,cfg,pulse=function(...)NULL,last=cfg$lambda_points) {
    if(setting$family=="enet" && setting$shape==0)
      return(rfh_ridge_path(prep,setting,cfg,e,last))
    e$pacc_path(prep,setting,cfg,pulse,last,e)
  }
  if(compile)rfh_compile(e)
  e
}
rfh_compile <- function(e) {
  Rcpp::sourceCpp(file.path(e$root,"src/lambert_performance_fast_v1.cpp"),env=e,
    cacheDir=file.path(tempdir(),"rfh_v2_rcpp"),showOutput=FALSE)
  e$lambert_cd_cpp <- e$lambert_cd_fast_v1; e$lambert_values_cpp <- e$lambert_values_fast_v1
  invisible(e)
}

# Exact same Gaussian objective as the existing path; n x n dual eigensystem.
rfh_ridge_path <- function(prep,setting,cfg,e,last=cfg$lambda_points) {
  started <- proc.time()[[3]]; n <- nrow(prep$Z)
  lambda <- e$perf_lambdas(prep,setting,cfg)[seq_len(last)]
  eg <- eigen(tcrossprod(prep$Z)/n,symmetric=TRUE)
  stopifnot(min(eg$values)>= -1e-10*max(1,max(eg$values)))
  d <- pmax(eg$values,0); uy <- drop(crossprod(eg$vectors,prep$y))
  beta <- crossprod(prep$Z,eg$vectors %*% (uy/(outer(d,lambda,"+"))))/n
  rows <- lapply(seq_along(lambda),function(k) {
    z <- e$perf_diagnose(beta[,k],prep,lambda[k],setting,cfg)
    data.frame(index=k,lambda=lambda[k],ok=z$ok,
      status=if(z$ok)"verified_stationary" else "kkt_not_met",kkt=z$kkt,
      solver_lambda=lambda[k],solver_shape=0,objective=z$objective,iterations=0,
      selected_start="dual_ridge",elapsed=NA_real_,max_increase=NA_real_,warning="",error="")
  })
  rows <- do.call(rbind,rows)
  list(beta=beta,diagnostics=rows,attempts=cbind(rows,executed=TRUE,reused=FALSE),
    elapsed=proc.time()[[3]]-started)
}
rfh_curvature <- function(b,prep,lambda,setting,e) {
  if(setting$family=="enet")return(list(status="convex_psd",min_eigenvalue=NA_real_))
  e$lic_curvature(b,prep,lambda,setting,e$cfg,e)
}
rfh_curvature_ok <- function(x)x %in% c("convex_psd","empty_active",
  "nonnegative_local","nonnegative_join_lower")
rfh_hbic <- function(path,prep,setting,e) {
  z <- path$diagnostics; n <- nrow(prep$Z); p <- ncol(prep$Z)
  z$d <- z$rss <- z$HBIC <- NA_real_; z$cap <- floor(n/log(n))
  z$cap_excluded <- z$eligible <- FALSE; z$curvature_status <- "not_assessed"
  for(k in which(z$ok)) {
    b <- path$beta[,k]; dd <- e$perf_diagnose(b,prep,z$lambda[k],setting,e$cfg)
    stopifnot(dd$ok,abs(dd$objective-z$objective[k])<1e-10*(1+abs(dd$objective)))
    z$d[k] <- sum(abs(b/prep$sx)>e$cfg$support_reporting_tolerance)
    z$rss[k] <- sum((prep$y-drop(prep$Z%*%b))^2)
    if(is.finite(z$rss[k]) && z$rss[k]>0)
      z$HBIC[k] <- n*log(z$rss[k]/n)+z$d[k]*log(p)*log(log(n))
    z$cap_excluded[k] <- z$d[k]>z$cap[k]
    if(!z$cap_excluded[k]) {
      z$curvature_status[k] <- rfh_curvature(b,prep,z$lambda[k],setting,e)$status
      z$eligible[k] <- rfh_curvature_ok(z$curvature_status[k]) && is.finite(z$HBIC[k])
    }
  }
  z
}
rfh_preps <- function(X,y,folds,e) c(lapply(1:5,function(f)
  e$perf_prepare(X[folds!=f,,drop=FALSE],y[folds!=f],e$cfg)),
  list(e$perf_prepare(X,y,e$cfg)))
rfh_cv_choices <- function(cv,e) {
  loss <- lapply(seq_len(nrow(e$settings)),function(j) {
    ll <- vapply(1:5,function(f)cv$loss[cv$setting==j & cv$fold==f],
      numeric(e$cfg$lambda_points))
    x <- rowMeans(ll); x[!is.finite(x)] <- Inf; list(loss=x)
  })
  setNames(lapply(e$cfg$methods,function(m)e$perf_choose(loss,e$settings,m,e$cfg)),e$cfg$methods)
}
rfh_hbic_setting <- function(m,e) switch(m,
  lambert_c1=which(e$settings$family=="lambert"),lasso=which(e$settings$family=="enet" & e$settings$shape==1),
  mcp=which(e$settings$family=="mcp"),scad=which(e$settings$family=="scad"))
rfh_run_unit <- function(i,e)e$perf_with_seed(e$seeds[[i]],function()rfh_fit_unit(i,e))
rfh_fit_unit <- function(i,e) {
  began <- proc.time()[[3]]; task <- e$tasks[i,,drop=FALSE]
  split <- e$ribo_task_split(e$view,task,e$cfg)
  X <- e$view$x[split$train,,drop=FALSE]; y <- e$view$y[split$train]; folds <- split$inner_folds
  preps <- rfh_preps(X,y,folds,e); J <- nrow(e$settings); L <- e$cfg$lambda_points
  paths <- vector("list",5*J); cv <- validation_predictions <- attempts <- list()
  for(f in 1:5)for(j in seq_len(J)) {
    a <- (f-1)*J+j; pp <- e$perf_path(preps[[f]],e$settings[j,],e$cfg)
    Ztest <- sweep(sweep(X[folds==f,,drop=FALSE],2,preps[[f]]$mx,"-"),2,preps[[f]]$sx,"/")
    pred <- preps[[f]]$my+Ztest%*%pp$beta; pred[,!pp$diagnostics$ok] <- NA_real_
    loss <- colMeans((y[folds==f]-pred)^2)
    cv[[a]] <- cbind(fold=f,setting=j,pp$diagnostics,loss=loss)
    validation_predictions[[a]] <- pred; paths[[a]] <- pp$beta
    attempts[[a]] <- cbind(phase="inner",fold=f,setting=j,pp$attempts)
  }
  cv <- do.call(rbind,cv); choices <- rfh_cv_choices(cv,e)
  selected_inner <- lapply(choices,function(ss)if(!ss$ok)NULL else lapply(1:5,function(f)
    paths[[(f-1)*J+ss$setting]][,ss$index]))
  paths <- NULL
  full <- lapply(seq_len(J),function(j)e$perf_path(preps[[6]],e$settings[j,],e$cfg))
  hbic <- setNames(lapply(e$cfg$hbic_methods,function(m) {
    j <- rfh_hbic_setting(m,e); rfh_hbic(full[[j]],preps[[6]],e$settings[j,],e)
  }),e$cfg$hbic_methods)
  z <- list(task=task,seed=e$seeds[[i]],split=split,cv=cv,
    validation_predictions=validation_predictions,selected_inner=selected_inner,
    cv_attempts=do.call(rbind,attempts),full=full,hbic=hbic,
    preprocessing=lapply(preps,function(p)p[c("mx","my","sx","lambda0","spectral")]))
  z$output <- rfh_outputs(z,preps,e); z$elapsed <- proc.time()[[3]]-began
  z
}

rfh_outputs <- function(z,preps,e) {
  choices <- rfh_cv_choices(z$cv,e); pp <- preps[[6]]
  arms <- rbind(data.frame(criterion="CV",method=e$cfg$methods),
    data.frame(criterion="HBIC",method=e$cfg$hbic_methods))
  rows <- predictions <- coefficients <- vector("list",nrow(arms))
  for(a in seq_len(nrow(arms))) {
    m <- arms$method[a]; cr <- arms$criterion[a]; ss <- choices[[m]]
    if(cr=="HBIC") {
      cand <- z$hbic[[m]]; k <- e$lic_choose(cand$HBIC,cand$eligible,e$cfg$implementation$tie_tol)
      ss <- list(ok=!is.na(k),setting=rfh_hbic_setting(m,e),index=k,
        partial=any(!cand$ok | (cand$ok & !cand$cap_excluded & !cand$eligible)))
    }
    r <- cbind(z$task,arms[a,,drop=FALSE],ok=FALSE,status="no_eligible_candidate",
      partial_search=ss$partial,n_train=length(z$split$train_index),n_test=length(z$split$test_index),
      setting=NA_integer_,index=NA_integer_,lambda=NA_real_,shape=NA_real_,
      cv_loss=NA_real_,hbic=NA_real_,kkt=NA_real_,curvature="not_assessed",
      selected_size=NA_integer_,rmse=NA_real_,mae=NA_real_,path_endpoint=FALSE)
    b <- rep(NA_real_,ncol(pp$Z)); intercept <- NA_real_
    if(ss$ok) {
      j <- ss$setting; k <- ss$index; path <- z$full[[j]]; diag <- path$diagnostics[k,]
      r$setting <- j; r$index <- k; r$lambda <- diag$lambda; r$shape <- e$settings$shape[j]
      r$kkt <- diag$kkt; r$status <- diag$status; r$path_endpoint <- k %in% c(1,e$cfg$lambda_points)
      if(cr=="CV")r$cv_loss <- ss$loss else r$hbic <- z$hbic[[m]]$HBIC[k]
      if(diag$ok) {
        h <- rfh_curvature(path$beta[,k],pp,diag$lambda,e$settings[j,],e)$status
        hs <- if(cr=="CV")vapply(1:5,function(f)rfh_curvature(z$selected_inner[[m]][[f]],preps[[f]],
          z$cv$lambda[z$cv$fold==f & z$cv$setting==j & z$cv$index==k],e$settings[j,],e)$status,"") else character()
        r$curvature <- h; r$ok <- all(rfh_curvature_ok(c(h,hs)))
        r$status <- if(r$ok)"verified_stationary" else "curvature_not_accepted"
        if(r$ok) { b <- path$beta[,k]/pp$sx; intercept <- pp$my-sum(pp$mx*b) }
      }
    }
    pred <- intercept+drop(e$view$x[z$split$test,,drop=FALSE]%*%b)
    if(r$ok) {
      met <- e$ribo_fold_metrics(e$view$y[z$split$test],pred)
      r$rmse <- met[["rmse"]]; r$mae <- met[["mae"]]
      r$selected_size <- sum(abs(b)>e$cfg$support_reporting_tolerance)
    }
    rows[[a]] <- r
    predictions[[a]] <- data.frame(unit=z$task$unit,repeat_id=z$task$repeat_id,
      outer_fold=z$task$outer_fold,criterion=cr,method=m,sample_index=z$split$test_index,
      sample_id=e$view$sample_id[z$split$test],observed=e$view$y[z$split$test],predicted=pred,ok=r$ok)
    coefficients[[a]] <- list(criterion=cr,method=m,beta=b,intercept=intercept)
  }
  list(rows=do.call(rbind,rows),predictions=do.call(rbind,predictions),coefficients=coefficients)
}
