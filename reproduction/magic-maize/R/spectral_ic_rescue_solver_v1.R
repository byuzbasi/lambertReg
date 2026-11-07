# Pilot-only hybrid, following the safeguarded active-set correction used in
# lambert_active_set_pilot_v1.R. Historical solvers and acceptance tolerances stay frozen.
sir_newton <- function(beta,prep,lambda,setting,e) {
  active <- which(beta!=0); steps <- 0L; status <- "empty_active"
  if(!length(active))return(list(beta=beta,steps=steps,status=status))
  signs <- sign(beta[active]); gram <- crossprod(prep$Z[,active,drop=FALSE])/nrow(prep$Z)
  for(step in seq_len(e$rescue$newton_steps)) {
    dd <- e$perf_diagnose(beta,prep,lambda,setting,e$cfg)
    if(dd$ok) { status<-"kkt_met";break }
    pp <- e$perf_penalty(abs(beta[active]),lambda,setting)
    if(any(!is.finite(pp$curvature))) { status<-"curvature_unavailable";break }
    H <- gram; diag(H)<-diag(H)+pp$curvature
    # A strictly positive Hessian is a proposal safeguard, not a relaxed final gate.
    if(min(eigen(H,symmetric=TRUE,only.values=TRUE)$values)<=0) {
      status<-"nonpositive_hessian";break
    }
    ch <- tryCatch(chol(H),error=function(err)NULL)
    if(is.null(ch)) { status<-"factorization_failed";break }
    residual <- drop(prep$Z%*%beta)-prep$y
    gradient <- drop(crossprod(prep$Z[,active,drop=FALSE],residual))/nrow(prep$Z)+signs*pp$slope
    direction <- -backsolve(ch,forwardsolve(t(ch),gradient)); descent<-sum(gradient*direction)
    if(!is.finite(descent)||descent>=0) { status<-"non_descent_direction";break }
    accepted <- FALSE
    for(backtrack in 0:e$rescue$backtracks) {
      alpha<-2^(-backtrack); candidate<-beta; candidate[active]<-beta[active]+alpha*direction
      if(any(sign(candidate[active])!=signs))next
      next_dd<-e$perf_diagnose(candidate,prep,lambda,setting,e$cfg)
      if(is.finite(next_dd$objective) && next_dd$objective<=dd$objective+e$rescue$armijo*alpha*descent) {
        accepted<-TRUE;break
      }
    }
    if(!accepted) { status<-"line_search_failed";break }
    beta<-candidate;steps<-steps+1L;status<-"step_budget_exhausted"
  }
  list(beta=beta,steps=steps,status=status)
}

sir_fit <- function(init,prep,lambda,setting,e) {
  started<-proc.time()[[3]]; beta<-as.numeric(init); sweeps<-newton_steps<-0L
  warnings<-character(); error<-""; status<-"rescue_budget_exhausted"; trace<-list()
  diagnostic<-e$perf_diagnose(beta,prep,lambda,setting,e$cfg)
  record<-function(stage,detail) {
    trace[[length(trace)+1L]]<<-data.frame(stage=stage,sweeps=sweeps,newton_steps=newton_steps,
      objective=diagnostic$objective,kkt=diagnostic$kkt,detail=detail)
  }
  record("initial","initialized")
  for(block in seq_len(e$rescue$total_sweeps %/% e$rescue$chunk_sweeps)) {
    before<-diagnostic$objective
    fit<-e$perf_capture(function() {
      if(setting$family=="lambert")e$lambert_cd_cpp(prep$Z,prep$y,beta,lambda,setting$shape,
        e$rescue$chunk_sweeps,e$cfg$normalized_kkt_tolerance/2,function(...)NULL)
      else ncvreg::ncvfit(prep$Z,prep$y,init=beta,lambda=lambda,
        penalty=if(setting$family=="enet")"lasso" else toupper(setting$family),
        gamma=if(setting$family=="enet")3 else setting$shape,alpha=1,
        eps=e$cfg$implementation$ncvreg_eps,max.iter=e$rescue$chunk_sweeps,warn=TRUE)
    })
    warnings<-unique(c(warnings,fit$warnings)); value<-fit$value
    if(inherits(value,"perf_error")) { status<-"solver_error";error<-value$message;break }
    candidate<-as.numeric(value$beta)
    sweeps<-sweeps+as.integer(if(setting$family=="lambert")value$iterations else value$iter)
    dd<-e$perf_diagnose(candidate,prep,lambda,setting,e$cfg)
    if(!is.finite(dd$objective) || isTRUE(value$increase_failed) ||
      dd$objective>before+e$cfg$implementation$objective_increase_tol*(1+abs(before))) {
      status<-"coordinate_objective_gate_failed";break
    }
    beta<-candidate;diagnostic<-dd;record("coordinate","updated")
    correction<-sir_newton(beta,prep,lambda,setting,e)
    dd<-e$perf_diagnose(correction$beta,prep,lambda,setting,e$cfg)
    if(!is.finite(dd$objective) || dd$objective>diagnostic$objective+
      e$cfg$implementation$objective_increase_tol*(1+abs(diagnostic$objective))) {
      status<-"newton_objective_gate_failed";break
    }
    beta<-correction$beta;diagnostic<-dd;newton_steps<-newton_steps+correction$steps
    record("newton",correction$status)
    if(dd$ok && rfh_curvature_ok(rfh_curvature(beta,prep,lambda,setting,e)$status)) {
      status<-"verified_stationary";break
    }
  }
  list(ok=identical(status,"verified_stationary"),status=status,beta=beta,
    kkt=diagnostic$kkt,objective=diagnostic$objective,sweeps=sweeps,newton_steps=newton_steps,
    seconds=proc.time()[[3]]-started,warning=paste(warnings[nzchar(warnings)],collapse=" | "),
    error=error,trace=do.call(rbind,trace))
}

sir_unit <- function(i,e) {
  started<-proc.time()[[3]]; old<-e$baseline[[i]]; z<-old; task<-z$task
  prep<-sip_prep(task$dataset,e); setting<-e$settings[task$setting,,drop=FALSE]
  attempts<-list(); warm<-numeric(ncol(prep$Z)); targets<-which(!old$path$diagnostics$ok)
  for(k in seq_len(e$cfg$lambda_points)) {
    if(old$path$diagnostics$ok[k]) { warm<-old$path$beta[,k];next }
    lambda<-old$path$diagnostics$lambda[k]
    starts<-list(warm=warm+0,zero=numeric(length(warm)))
    fits<-lapply(starts,function(b)sir_fit(b,prep,lambda,setting,e))
    good<-which(vapply(fits,function(f)f$ok,FALSE))
    chosen<-NA_integer_
    if(length(good)) {
      q<-vapply(fits[good],function(f)f$objective,1.0);best<-min(q)
      chosen<-good[which(q<=best+e$cfg$implementation$tie_tol*(1+abs(best)))[1]]
      fit<-fits[[chosen]];z$path$beta[,k]<-fit$beta;warm<-fit$beta
      d<-z$path$diagnostics;d$ok[k]<-TRUE;d$status[k]<-"verified_stationary"
      d$kkt[k]<-fit$kkt;d$objective[k]<-fit$objective;d$iterations[k]<-fit$sweeps
      d$selected_start[k]<-paste0("rescue_",names(fits)[chosen]);d$elapsed[k]<-sum(vapply(fits,function(f)f$seconds,1.0))
      d$max_increase[k]<-max(c(0,diff(fit$trace$objective)))
      d$warning[k]<-fit$warning;d$error[k]<-fit$error;z$path$diagnostics<-d
    }
    # Failed endpoints are retained as evidence, never propagated as warm starts.
    attempts[[as.character(k)]]<-list(index=k,initial=starts,fits=fits,chosen=chosen)
  }
  z$path_seconds<-proc.time()[[3]]-started;z$path$elapsed<-z$path_seconds
  scoring<-proc.time()[[3]];z$candidates<-sip_score(z$path,prep,setting,e)
  z$selected<-sip_choices(z$candidates,z$path,prep,e);z$scoring_seconds<-proc.time()[[3]]-scoring
  z$elapsed<-proc.time()[[3]]-started
  z$rescue<-list(baseline_hash=cpk_object_hash(old),targets=targets,attempts=attempts)
  z
}
