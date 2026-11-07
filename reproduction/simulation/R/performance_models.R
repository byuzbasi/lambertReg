# Statistical implementation of protocol/performance_v1.md.
perf_config <- function(root,stage) {
  stopifnot(stage %in% c("smoke","pilot","confirmation"))
  lines <- readLines(file.path(root,"protocol/performance_v1.md"),warn=FALSE)
  fence <- strrep(intToUtf8(96),3)
  start <- which(lines==paste0(fence,"json")); end <- which(lines==fence)
  cfg <- jsonlite::fromJSON(paste(lines[seq(start+1,min(end[end>start])-1)],collapse="\n"))
  cfg$stage <- stage
  cfg$implementation <- list(ncvreg_eps=1e-12,ncvreg_max_iter=10000L,
    glmnet_thresh=1e-16,glmnet_maxit=1000000L,progress_seconds=30,
    tie_tol=1e-10,objective_increase_tol=1e-9,curvature_join_tol=1e-8,
    curvature_negative_tol=1e-8,standardization_tol=1e-12,
    execution_order="replicate_then_scenario",version="performance_v1")
  if(stage=="smoke") {
    cfg$n_train<-40L; cfg$p<-c(20L,60L); cfg$rho<-c(0,0.9)
    cfg$snr<-1; cfg$profiles<-"heterogeneous"; cfg$n_test<-40L
    cfg$lambda_points<-4L; cfg$repetitions<-1L
  } else cfg$repetitions<-cfg[[paste0(stage,"_repetitions")]]
  cfg
}
perf_tasks <- function(cfg) {
  rows<-list(); s<-0L; u<-0L
  for(p in cfg$p) for(rho in cfg$rho) for(snr in cfg$snr) for(profile in cfg$profiles) {
    s<-s+1L
    for(r in seq_len(cfg$repetitions)) {
      u<-u+1L
      stream<-if(cfg$stage=="smoke")cfg$smoke_stream_indices[u] else cfg[[paste0(cfg$stage,"_stream_offset")]]+u
      rows[[u]]<-data.frame(unit=u,scenario=s,repetition=r,n=cfg$n_train,p=p,
        rho=rho,snr=snr,profile=profile,stream_index=stream)
    }
  }
  do.call(rbind,rows)
}
perf_seeds <- function(cfg,tasks) {
  RNGkind("L'Ecuyer-CMRG",normal.kind="Inversion",sample.kind="Rejection")
  set.seed(cfg$master_seed); state<-.Random.seed
  states<-vector("list",nrow(tasks)); j<-1L
  for(i in seq_len(max(tasks$stream_index))) {
    if(j<=nrow(tasks)&&i==tasks$stream_index[j]) {
      sub<-vector("list",5L); sub[[1]]<-state
      for(k in 2:5) sub[[k]]<-parallel::nextRNGSubStream(sub[[k-1]])
      states[[j]]<-sub; j<-j+1L
    }
    state<-parallel::nextRNGStream(state)
  }
  states
}
perf_with_seed <- function(seed,fun) {
  old<-get0(".Random.seed",envir=.GlobalEnv)
  on.exit(if(is.null(old))rm(".Random.seed",envir=.GlobalEnv) else assign(".Random.seed",old,envir=.GlobalEnv))
  assign(".Random.seed",seed,envir=.GlobalEnv); fun()
}
perf_x <- function(n,p,rho,seed) perf_with_seed(seed,function() {
  X<-matrix(rnorm(n*p),n,p)
  if(p>1) for(j in 2:p) X[,j]<-rho*X[,j-1]+sqrt(1-rho^2)*X[,j]
  X
})
perf_truth <- function(task,cfg) {
  beta<-numeric(task$p)
  beta[cfg$support_indices]<-cfg$signs*cfg[[paste0("magnitudes_",task$profile)]]
  beta<-beta/sqrt(sum(beta^2))
  Sigma<-task$rho^abs(outer(seq_len(task$p),seq_len(task$p),"-"))
  list(beta=beta,Sigma=Sigma,sigma2=drop(crossprod(beta,Sigma%*%beta))/task$snr)
}
perf_numeric_hash <- function(x) digest::digest(writeBin(as.double(x),raw(),size=8,endian="little"),
  algo="sha256",serialize=FALSE)
perf_training <- function(task,cfg,seeds) {
  truth<-perf_truth(task,cfg); X<-perf_x(task$n,task$p,task$rho,seeds[[1]])
  y<-drop(X%*%truth$beta)+perf_with_seed(seeds[[2]],function()rnorm(task$n,sd=sqrt(truth$sigma2)))
  folds<-perf_with_seed(seeds[[3]],function()sample(rep(seq_len(cfg$folds),each=task$n/cfg$folds)))
  list(X=X,y=y,folds=folds,truth=truth,
    hashes=c(X=perf_numeric_hash(X),y=perf_numeric_hash(y),folds=perf_numeric_hash(folds)))
}
perf_test_data <- function(task,cfg,seeds,truth) {
  X<-perf_x(cfg$n_test,task$p,task$rho,seeds[[4]])
  y<-drop(X%*%truth$beta)+perf_with_seed(seeds[[5]],function()rnorm(cfg$n_test,sd=sqrt(truth$sigma2)))
  list(X=X,y=y,hashes=c(X=perf_numeric_hash(X),y=perf_numeric_hash(y)))
}
perf_prepare <- function(X,y,cfg) {
  mx<-colMeans(X); my<-mean(y); Z<-sweep(X,2,mx,"-"); sx<-sqrt(colMeans(Z^2))
  if(any(!is.finite(sx))||any(sx<=cfg$implementation$standardization_tol)||any(!is.finite(y)))
    stop("data_validation_failed: nonfinite or constant training inputs")
  Z<-sweep(Z,2,sx,"/"); yc<-y-my; score<-drop(crossprod(Z,yc)/nrow(Z)); lambda0<-max(abs(score))
  if(!is.finite(lambda0)||lambda0<=0)stop("data_validation_failed: lambda0")
  list(Z=Z,y=yc,mx=mx,my=my,sx=sx,score=score,lambda0=lambda0,
    spectral=max(svd(Z,nu=0,nv=0)$d)^2/nrow(Z))
}
perf_settings <- function(cfg) rbind(data.frame(family="enet",shape=cfg$elastic_net_alpha),
  data.frame(family="mcp",shape=cfg$mcp_gamma),data.frame(family="scad",shape=cfg$scad_a),
  data.frame(family="lambert",shape=cfg$lambert_shapes))
perf_lambdas <- function(prep,setting,cfg) {
  if(setting$family=="enet"&&setting$shape==0)
    prep$spectral*10^seq(3,-4,length.out=cfg$lambda_points)
  else prep$lambda0*10^seq(0,-3,length.out=cfg$lambda_points)/
    if(setting$family=="enet")setting$shape else 1
}
perf_penalty <- function(t,lambda,setting) {
  if(setting$family=="lambert")return(lambert_values_cpp(t,lambda,setting$shape))
  if(setting$family=="enet")return(list(penalty=lambda*(setting$shape*t+(1-setting$shape)*t^2/2),
    slope=lambda*(setting$shape+(1-setting$shape)*t),curvature=rep(lambda*(1-setting$shape),length(t))))
  list(penalty=penalty_value(t,lambda,setting$family,setting$shape),
    slope=penalty_slope(t,lambda,setting$family,setting$shape),
    curvature=penalty_curvature(t,lambda,setting$family,setting$shape))
}
perf_diagnose <- function(beta,prep,lambda,setting,cfg) {
  if(length(beta)!=ncol(prep$Z)||any(!is.finite(beta)))return(list(ok=FALSE,kkt=Inf,objective=NA_real_))
  r<-drop(prep$Z%*%beta)-prep$y; gradient<-drop(crossprod(prep$Z,r)/nrow(prep$Z))
  pp<-perf_penalty(abs(beta),lambda,setting)
  entry<-if(setting$family=="enet")lambda*setting$shape else lambda
  residual<-pmax(abs(gradient)-entry,0); active<-beta!=0
  residual[active]<-abs(gradient[active]+sign(beta[active])*pp$slope[active])
  kkt<-max(residual)/(1+max(abs(prep$score))); q<-mean(r^2)/2+sum(pp$penalty)
  list(ok=is.finite(q)&&is.finite(kkt)&&kkt<=cfg$normalized_kkt_tolerance,kkt=kkt,objective=q)
}
perf_hessian <- function(beta,prep,lambda,setting,cfg) {
  active<-which(beta!=0)
  if(!length(active))return(list(status="empty_active",min_eigenvalue=NA_real_,condition=NA_real_))
  t<-abs(beta[active])
  cuts<-switch(setting$family,lambert=lambda*exp(setting$shape),mcp=lambda*setting$shape,
    scad=c(lambda,lambda*setting$shape),enet=numeric())
  if(length(cuts)&&any(abs(outer(t,cuts,"-"))<=cfg$implementation$curvature_join_tol*max(1,lambda)))
    return(list(status="join_unavailable",min_eigenvalue=NA_real_,condition=NA_real_))
  H<-crossprod(prep$Z[,active,drop=FALSE])/nrow(prep$Z); diag(H)<-diag(H)+perf_penalty(t,lambda,setting)$curvature
  ee<-eigen(H,symmetric=TRUE,only.values=TRUE)$values
  list(status=if(min(ee)< -cfg$implementation$curvature_negative_tol*max(1,max(abs(ee))))
    "negative_curvature" else "nonnegative_local",min_eigenvalue=min(ee),
    condition=if(min(ee)>0)max(ee)/min(ee) else Inf)
}
perf_predict <- function(beta,prep,X)prep$my+drop(sweep(sweep(X,2,prep$mx,"-"),2,prep$sx,"/")%*%beta)
perf_capture <- function(fun) {
  warnings<-character()
  value<-tryCatch(withCallingHandlers(fun(),warning=function(w) {
    warnings<<-c(warnings,conditionMessage(w)); invokeRestart("muffleWarning")
  }),error=function(e)structure(list(message=conditionMessage(e)),class="perf_error"))
  list(value=value,warnings=paste(unique(warnings),collapse=" | "))
}

perf_path <- function(prep,setting,cfg,pulse=function(...)NULL,last=cfg$lambda_points) {
  lambda<-perf_lambdas(prep,setting,cfg)[seq_len(last)]
  betas<-matrix(NA_real_,ncol(prep$Z),last)
  diagrows<-vector("list",last); attempts<-list(); warm<-numeric(ncol(prep$Z))
  total_start<-proc.time()[[3]]
  if(setting$family=="enet") {
    old<-glmnet::glmnet.control()
    on.exit(do.call(glmnet::glmnet.control,old),add=TRUE)
    glmnet::glmnet.control(fdev=0,devmax=1)
    # The fast Gaussian backend rescales y internally. Its original-scale
    # quadratic penalty is lambda_api*(1-alpha_api)/(2*response_rms).
    # Map API parameters so the approved objective is exactly preserved.
    response_rms<-sqrt(mean(prep$y^2))
    api_factor<-setting$shape+response_rms*(1-setting$shape)
    api_alpha<-setting$shape/api_factor
    api_lambda<-lambda*api_factor
    cap<-perf_capture(function()glmnet::glmnet(prep$Z,prep$y,family="gaussian",
      alpha=api_alpha,lambda=api_lambda,standardize=FALSE,intercept=FALSE,
      control=list(thresh=cfg$implementation$glmnet_thresh,maxit=cfg$implementation$glmnet_maxit,
        dfmax=ncol(prep$Z)+1,pmax=ncol(prep$Z),fdev=0,devmax=1)))
    fit<-cap$value; failed<-inherits(fit,"perf_error")
    for(k in seq_along(lambda)) {
      found<-if(failed)integer() else which(abs(fit$lambda-api_lambda[k])<=1e-12*api_lambda[k])
      beta<-if(length(found))as.numeric(fit$beta[,found[1]]) else rep(NA_real_,ncol(prep$Z))
      dd<-perf_diagnose(beta,prep,lambda[k],setting,cfg)
      ok<-dd$ok&&length(found)>0&&!failed&&fit$jerr<=0
      status<-if(failed)"solver_error" else if(!length(found))"budget_exhausted" else
        if(!ok)"kkt_not_met" else "verified_stationary"
      betas[,k]<-beta
      diagrows[[k]]<-data.frame(index=k,lambda=lambda[k],ok=ok,status=status,kkt=dd$kkt,
        solver_lambda=api_lambda[k],solver_shape=api_alpha,
        objective=dd$objective,iterations=if(failed)NA else fit$npasses,
        selected_start="glmnet_path",elapsed=NA_real_,max_increase=NA_real_,
        warning=cap$warnings,error=if(failed)fit$message else "")
      pulse(paste("glmnet",setting$shape,"lambda",k))
    }
    attempts<-diagrows
  } else {
    for(k in seq_along(lambda)) {
      candidates<-list()
      for(start_id in c("warm","zero")) {
        start<-if(start_id=="warm")warm+0 else numeric(length(warm))
        began<-proc.time()[[3]]
        cap<-perf_capture(function() {
          if(setting$family=="lambert")
            lambert_cd_cpp(prep$Z,prep$y,start,lambda[k],setting$shape,
              cfg$max_coordinate_sweeps_per_start,cfg$normalized_kkt_tolerance/2,pulse)
          else ncvreg::ncvfit(prep$Z,prep$y,init=start,penalty=toupper(setting$family),
            gamma=setting$shape,lambda=lambda[k],eps=cfg$implementation$ncvreg_eps,
            max.iter=cfg$implementation$ncvreg_max_iter,warn=TRUE)
        })
        fit<-cap$value; failed<-inherits(fit,"perf_error")
        beta<-if(failed)rep(NA_real_,length(warm)) else as.numeric(fit$beta)
        dd<-perf_diagnose(beta,prep,lambda[k],setting,cfg)
        iter<-if(failed)NA_real_ else if(setting$family=="lambert")fit$iterations else fit$iter
        solver_ok<-!failed&&(if(setting$family=="lambert")fit$converged else iter<cfg$implementation$ncvreg_max_iter)
        increased<-!failed&&setting$family=="lambert"&&fit$increase_failed
        ok<-solver_ok&&dd$ok&&!increased
        status<-if(failed)"solver_error" else if(increased)"objective_increased" else
          if(!solver_ok)"budget_exhausted" else if(!dd$ok)"kkt_not_met" else "verified_stationary"
        row<-data.frame(index=k,lambda=lambda[k],ok=ok,status=status,kkt=dd$kkt,
          solver_lambda=lambda[k],solver_shape=setting$shape,
          objective=dd$objective,iterations=iter,selected_start=start_id,
          elapsed=proc.time()[[3]]-began,
          max_increase=if(!failed&&setting$family=="lambert")fit$max_increase else NA_real_,
          warning=cap$warnings,error=if(failed)fit$message else "")
        candidates[[start_id]]<-list(beta=beta,row=row)
        attempts[[length(attempts)+1L]]<-row
        pulse(paste(setting$family,setting$shape,"lambda",k,start_id))
      }
      good<-which(vapply(candidates,function(x)x$row$ok,logical(1)))
      if(length(good)) {
        q<-vapply(candidates[good],function(x)x$row$objective,numeric(1)); best<-min(q)
        selected<-good[which(q<=best+cfg$implementation$tie_tol*(1+abs(best)))[1]]
        betas[,k]<-candidates[[selected]]$beta
        diagrows[[k]]<-candidates[[selected]]$row; warm<-betas[,k]
      } else {
        diagrows[[k]]<-candidates[[1]]$row
        diagrows[[k]]$selected_start<-"none_verified"; warm[]<-0
      }
    }
  }
  list(beta=betas,diagnostics=do.call(rbind,diagrows),attempts=do.call(rbind,attempts),
    elapsed=proc.time()[[3]]-total_start)
}
perf_choose <- function(cv,settings,method,cfg) {
  eligible<-switch(method,
    lasso=which(settings$family=="enet"&settings$shape==1),
    ridge=which(settings$family=="enet"&settings$shape==0),
    elastic_net=which(settings$family=="enet"),mcp=which(settings$family=="mcp"),
    scad=which(settings$family=="scad"),lambert_cv=which(settings$family=="lambert"),
    lambert_c1=which(settings$family=="lambert"&settings$shape==1),
    lambert_c1p5=which(settings$family=="lambert"&settings$shape==1.5),
    lambert_c2=which(settings$family=="lambert"&settings$shape==2))
  best_each<-lapply(eligible,function(j) {
    loss<-cv[[j]]$loss; valid<-which(is.finite(loss))
    if(!length(valid))return(NULL)
    best<-min(loss[valid]); k<-which(loss<=best+cfg$implementation$tie_tol*(1+best))[1]
    list(setting=j,index=k,loss=loss[k])
  })
  best_each<-Filter(Negate(is.null),best_each)
  partial<-any(vapply(cv[eligible],function(x)any(!is.finite(x$loss)),logical(1)))
  if(!length(best_each))return(list(ok=FALSE,partial=partial,eligible=eligible))
  values<-vapply(best_each,function(x)x$loss,numeric(1)); best<-min(values)
  tied<-best_each[which(values<=best+cfg$implementation$tie_tol*(1+best))]
  ord<-if(startsWith(method,"lambert"))order(vapply(tied,function(x)settings$shape[x$setting],numeric(1))) else
    order(-vapply(tied,function(x)settings$shape[x$setting],numeric(1)))
  c(list(ok=TRUE,partial=partial,eligible=eligible),tied[[ord[1]]])
}

# Only training observations and folds are accessible to model selection.
perf_train_models <- function(X,y,folds,cfg,pulse=function(...)NULL) {
  settings<-perf_settings(cfg); K<-cfg$folds; J<-nrow(settings)
  preps<-lapply(seq_len(K),function(f)perf_prepare(X[folds!=f,,drop=FALSE],y[folds!=f],cfg))
  paths<-lapply(seq_len(J),function(j)vector("list",K)); cv<-vector("list",J)
  records<-attempts<-list(); times<-matrix(0,J,K)
  for(f in seq_len(K))for(j in seq_len(J)) {
    pp<-perf_path(preps[[f]],settings[j,],cfg,pulse)
    paths[[j]][[f]]<-pp; times[j,f]<-pp$elapsed
    loss<-rep(NA_real_,cfg$lambda_points)
    for(k in which(pp$diagnostics$ok)) {
      pred<-perf_predict(pp$beta[,k],preps[[f]],X[folds==f,,drop=FALSE])
      if(all(is.finite(pred)))loss[k]<-mean((y[folds==f]-pred)^2)
    }
    records[[length(records)+1L]]<-cbind(fold=f,setting=j,
      family=settings$family[j],shape=settings$shape[j],pp$diagnostics,loss=loss)
    attempts[[length(attempts)+1L]]<-cbind(phase="cv",fold=f,setting=j,
      family=settings$family[j],shape=settings$shape[j],pp$attempts)
  }
  cv_records<-do.call(rbind,records)
  for(j in seq_len(J)) {
    ll<-vapply(seq_len(K),function(f)cv_records$loss[cv_records$setting==j&cv_records$fold==f],
      numeric(cfg$lambda_points))
    loss<-rowMeans(ll); loss[!is.finite(loss)]<-Inf; cv[[j]]<-list(loss=loss)
  }
  selections<-setNames(lapply(cfg$methods,function(m)perf_choose(cv,settings,m,cfg)),cfg$methods)
  fullprep<-perf_prepare(X,y,cfg); fullpaths<-vector("list",J)
  for(j in seq_len(J)) {
    needed<-Filter(function(s)isTRUE(s$ok)&&s$setting==j,selections)
    if(length(needed)) {
      last<-max(vapply(needed,function(s)s$index,numeric(1)))
      pp<-perf_path(fullprep,settings[j,],cfg,pulse,last=last); fullpaths[[j]]<-pp
      attempts[[length(attempts)+1L]]<-cbind(phase="full",fold=0,setting=j,
        family=settings$family[j],shape=settings$shape[j],pp$attempts)
    }
  }
  models<-list()
  for(method in cfg$methods) {
    ss<-selections[[method]]
    model<-list(selection=ss,status="no_valid_cv",ok=FALSE,tuning_seconds=sum(times[ss$eligible,,drop=FALSE]))
    if(ss$ok) {
      j<-ss$setting; k<-ss$index; pp<-fullpaths[[j]]
      model$setting<-settings[j,]; model$diagnostics<-pp$diagnostics[k,]; model$refit_seconds<-pp$elapsed
      model$cv_curvature<-lapply(seq_len(K),function(f)
        perf_hessian(paths[[j]][[f]]$beta[,k],preps[[f]],paths[[j]][[f]]$diagnostics$lambda[k],settings[j,],cfg))
      if(pp$diagnostics$ok[k]) {
        model$curvature<-perf_hessian(pp$beta[,k],fullprep,pp$diagnostics$lambda[k],settings[j,],cfg)
        negative<-any(vapply(c(model$cv_curvature,list(model$curvature)),function(h)h$status=="negative_curvature",logical(1)))
        model$ok<-!negative; model$status<-if(negative)"negative_curvature" else "verified_stationary"
        model$beta<-pp$beta[,k]/fullprep$sx
        model$intercept<-fullprep$my-sum(fullprep$mx*model$beta)
      } else model$status<-pp$diagnostics$status[k]
    }
    models[[method]]<-model
  }
  selected_folds<-lapply(selections,function(ss)if(!ss$ok)NULL else
    lapply(seq_len(K),function(f)list(beta=paths[[ss$setting]][[f]]$beta[,ss$index],
      mx=preps[[f]]$mx,sx=preps[[f]]$sx,my=preps[[f]]$my)))
  list(models=models,cv=cv_records,attempts=do.call(rbind,attempts),
    selected_folds=selected_folds,folds=folds,
    preprocessing=lapply(c(preps,list(fullprep)),function(x)x[c("mx","my","sx","lambda0","spectral")]),
    timings=list(cv_by_setting_fold=times,full_by_setting=vapply(fullpaths,function(x)if(is.null(x))0 else x$elapsed,numeric(1))),
    shared_computation_note="Method times include full reusable paths; do not sum across methods.")
}
perf_run_task <- function(task,cfg,seeds,pulse=function(...)NULL) {
  began<-proc.time()[[3]]; data<-perf_training(task,cfg,seeds)
  fitted<-perf_train_models(data$X,data$y,data$folds,cfg,pulse)
  test<-perf_test_data(task,cfg,seeds,data$truth)
  rows<-lapply(cfg$methods,function(method) {
    model<-fitted$models[[method]]; ss<-model$selection
    row<-data.frame(unit=task$unit,scenario=task$scenario,repetition=task$repetition,
      method=method,n=task$n,p=task$p,rho=task$rho,snr=task$snr,profile=task$profile,
      ok=model$ok,status=model$status,partial_search=ss$partial,
      family=if(ss$ok)model$setting$family else NA_character_,shape=if(ss$ok)model$setting$shape else NA_real_,
      index=if(ss$ok)ss$index else NA_integer_,lambda=if(ss$ok)model$diagnostics$lambda else NA_real_,
      cv_loss=if(ss$ok)ss$loss else NA_real_,lower_endpoint=if(ss$ok)ss$index==cfg$lambda_points else NA,
      upper_endpoint=if(ss$ok)ss$index==1 else NA,kkt=if(ss$ok)model$diagnostics$kkt else NA_real_,
      objective=if(ss$ok)model$diagnostics$objective else NA_real_,
      min_eigenvalue=if(is.null(model$curvature))NA_real_ else model$curvature$min_eigenvalue,
      condition=if(is.null(model$curvature))NA_real_ else model$curvature$condition,
      curvature_status=if(is.null(model$curvature))NA_character_ else model$curvature$status,
      tuning_seconds=model$tuning_seconds,refit_seconds=if(is.null(model$refit_seconds))NA_real_ else model$refit_seconds,
      mse=NA_real_,normalized_mse=NA_real_,coefficient_error=NA_real_,analytic_risk=NA_real_,
      tp=NA_integer_,fp=NA_integer_,fn=NA_integer_,tpr=NA_real_,fdr=NA_real_,
      selected_size=NA_integer_,exact_support=NA,sigma2=data$truth$sigma2)
    if(model$ok) {
      pred<-model$intercept+drop(test$X%*%model$beta)
      if(any(!is.finite(pred))) { row$ok<-FALSE; row$status<-"nonfinite_prediction"; return(row) }
      difference<-model$beta-data$truth$beta; active<-which(abs(model$beta)>cfg$support_reporting_tolerance)
      row$mse<-mean((test$y-pred)^2); row$normalized_mse<-row$mse/data$truth$sigma2
      row$coefficient_error<-sum(difference^2)
      row$analytic_risk<-data$truth$sigma2+model$intercept^2+drop(crossprod(difference,data$truth$Sigma%*%difference))
      row$tp<-sum(active%in%cfg$support_indices); row$fp<-length(active)-row$tp; row$fn<-cfg$support_size-row$tp
      row$tpr<-row$tp/cfg$support_size; row$fdr<-row$fp/max(1,length(active))
      row$selected_size<-length(active); row$exact_support<-setequal(active,cfg$support_indices)
    }
    row
  })
  list(task=task,rows=do.call(rbind,rows),fit=fitted,seeds=seeds,
    data_hashes=list(training=data$hashes,test=test$hashes),truth=list(beta=data$truth$beta,sigma2=data$truth$sigma2),
    elapsed=proc.time()[[3]]-began)
}
