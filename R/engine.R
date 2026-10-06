# Lambert-only engine with equivalent computational shortcuts; see inst/PROVENANCE.json.
perf_prepare <- function(X,y,cfg) {
  mx<-colMeans(X); my<-mean(y); Z<-sweep(X,2,mx,"-"); sx<-sqrt(colMeans(Z^2))
  if(any(!is.finite(sx))||any(sx<=cfg$implementation$standardization_tol)||any(!is.finite(y)))
    stop("data_validation_failed: nonfinite or constant training inputs")
  Z<-sweep(Z,2,sx,"/"); yc<-y-my; score<-drop(crossprod(Z,yc)/nrow(Z)); lambda0<-max(abs(score))
  if(!is.finite(lambda0)||lambda0<=0)stop("data_validation_failed: lambda0")
  list(Z=Z,y=yc,mx=mx,my=my,sx=sx,score=score,lambda0=lambda0)
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
  reused<-logical(2L*last)
  total_start<-proc.time()[[3]]
    for(k in seq_along(lambda)) {
      candidates<-list()
      for(start_id in c("warm","zero")) {
        # Reuse only a verified result from an exactly identical initial vector.
        if(start_id=="zero" && all(warm==0) && isTRUE(candidates$warm$row$ok)) {
          candidate<-candidates$warm
          candidate$row$selected_start<-"zero"; candidate$row$elapsed<-0
          candidates[[start_id]]<-candidate
          attempts[[length(attempts)+1L]]<-candidate$row
          reused[length(attempts)]<-TRUE
          pulse(paste(setting$family,setting$shape,"lambda",k,start_id))
          next
        }
        start<-if(start_id=="warm")warm+0 else numeric(length(warm))
        began<-proc.time()[[3]]
        cap<-perf_capture(function() {
          lambert_cd_cpp(prep$Z,prep$y,start,lambda[k],1,
            cfg$max_coordinate_sweeps_per_start,cfg$normalized_kkt_tolerance/2,pulse)
        })
        fit<-cap$value; failed<-inherits(fit,"perf_error")
        beta<-if(failed)rep(NA_real_,length(warm)) else as.numeric(fit$beta)
        dd<-perf_diagnose(beta,prep,lambda[k],setting,cfg)
        iter<-if(failed)NA_real_ else fit$iterations
        solver_ok<-!failed&&fit$converged
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
  attempt_rows<-do.call(rbind,attempts)
  attempt_rows$executed<-!reused; attempt_rows$reused<-reused
  list(beta=betas,diagnostics=do.call(rbind,diagrows),attempts=attempt_rows,
    elapsed=proc.time()[[3]]-total_start)
}

perf_penalty <- function(t,lambda,setting) lambert_values_cpp(t,lambda,1)
perf_lambdas <- function(prep,setting,cfg) {
  if (!is.null(cfg$lambda)) cfg$lambda else prep$lambda0 * cfg$lambda_fraction
}

