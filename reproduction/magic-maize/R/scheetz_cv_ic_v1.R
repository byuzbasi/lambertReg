# New dataset adapter; existing numerical sources and previous results are immutable.
sch_raw <- function()"data/raw/scheetz2006_iowabiostat_20260921/Scheetz2006.rds"
sch_prepared <- function()"data/processed/scheetz_cv_ic_v1_a01"
sch_validation <- function()"outputs/validation/scheetz_cv_ic_v1_a02"
sch_sources <- function()c("protocol/scheetz_cv_ic_pilot_v1.md","protocol/performance_v1.md",
  "R/penalties.R","R/lambert_family.R","R/performance_models.R","R/performance_acceleration_v1.R",
  "R/lambert_ic_pilot_v1.R","R/lambert_bic_comparison_v1.R","R/lambert_hbic_ebic_reanalysis_v1.R",
  "R/riboflavin_fixed_hd_v2.R","R/riboflavin_fixed_hd_execution_v2.R","R/spectral_ic_rescue_solver_v1.R",
  "R/spectral_ic_rescue_v1.R","R/checkpoint_packs_v1.R","src/lambert_performance_fast_v1.cpp",
  "R/lambert_bic_release_v1.R","R/scheetz_cv_ic_v1.R","R/scheetz_cv_ic_run_v1.R","R/scheetz_cv_ic_release_v1.R",
  "scripts/scheetz_cv_ic_v1.R","tests/scheetz_cv_ic_v1_checks.R","tests/scheetz_cv_ic_linux_v1.R",
  "truba/run_scheetz_cv_ic_v1.slurm","truba/preflight_scheetz_cv_ic_v1.sh",
  "truba/build_scheetz_cv_ic_v1_bundle.R","truba/README_TRUBA_SCHEETZ_CV_IC_V1.md")
sch_close <- function(a,b)isTRUE(all.equal(a,b,tolerance=1e-11,check.attributes=TRUE))
sch_number <- function(a,b)all(is.finite(a)==is.finite(b)) &&
  all(abs(a[is.finite(a)]-b[is.finite(b)])<=1e-11*(1+abs(b[is.finite(b)])))
sch_bind <- function(x) {z<-do.call(rbind,x);rownames(z)<-NULL;z}
sch_text_hash <- function(x)digest::digest(paste(x,collapse="\n"),algo="sha256",serialize=FALSE)
sch_load <- function(root) {
  stopifnot(cpk_hash(file.path(root,sch_raw()))=="1b3450f6ed59afce27dd162e4563aa27c6512e41f65735b5fa448e0643d5d6b4")
  d<-readRDS(file.path(root,sch_raw()))
  stopifnot(identical(dim(d$X),c(120L,18975L)),length(d$y)==120L,
    all(is.finite(d$X)),all(is.finite(d$y)),!anyDuplicated(colnames(d$X)))
  d
}
sch_sample <- function(seed,labels) {
  kind<-RNGkind();old<-get0(".Random.seed",envir=.GlobalEnv)
  on.exit({do.call(RNGkind,as.list(kind));if(is.null(old)) {
    if(exists(".Random.seed",.GlobalEnv,inherits=FALSE))rm(".Random.seed",envir=.GlobalEnv)
  } else assign(".Random.seed",old,.GlobalEnv)},add=TRUE)
  RNGkind("Mersenne-Twister","Inversion","Rejection");set.seed(seed);sample(labels)
}
sch_folds <- function(n=120L,repeats=2L)lapply(seq_len(repeats),function(r) {
  outer<-sch_sample(2026101801L+1000L*(r-1L),rep(1:5,length.out=n))
  list(outer=outer,inner=lapply(1:5,function(k)
    sch_sample(2026101802L+1000L*(r-1L)+k,rep(1:5,length.out=sum(outer!=k)))))
})
sch_prepare <- function(root) {
  sch_load(root);out<-file.path(root,sch_prepared());stopifnot(!file.exists(out));dir.create(out,recursive=TRUE)
  design<-list(schema="scheetz_cv_ic_design_v1",sample_index=1:120,folds=sch_folds(),
    seeds=data.frame(repeat_id=1:2,outer_seed=2026101801L+1000L*(0:1),inner_seed_base=2026101802L+1000L*(0:1)),
    raw_sha256=cpk_hash(file.path(root,sch_raw())))
  saveRDS(design,file.path(out,"design.rds"),version=3)
  mf<-rfh_files(root,file.path(sch_prepared(),"design.rds"));mf$file<-basename(mf$file)
  write.csv(mf,file.path(out,"manifest.csv"),row.names=FALSE)
  invisible(out)
}
sch_environment <- function(root,stage="pilot",compile=FALSE) {
  stopifnot(stage %in% c("pilot","smoke"))
  e<-new.env(parent=environment(sch_environment));e$root<-root;e$stage<-stage
  for(f in c("penalties.R","lambert_family.R","performance_models.R","performance_acceleration_v1.R",
    "lambert_ic_pilot_v1.R","lambert_hbic_ebic_reanalysis_v1.R"))sys.source(file.path(root,"R",f),e)
  e$data<-sch_load(root);e$view<-list(rows=1:120,columns=1:18975)
  prepared<-file.path(root,sch_prepared());mf<-read.csv(file.path(prepared,"manifest.csv"))
  stopifnot(identical(mf$file,"design.rds"),mf$bytes==file.info(file.path(prepared,mf$file))$size,
    mf$sha256==cpk_hash(file.path(prepared,mf$file)))
  design<-readRDS(file.path(prepared,"design.rds"));stopifnot(identical(design$folds,sch_folds()),
    design$raw_sha256==cpk_hash(file.path(root,sch_raw())))
  e$folds<-design$folds;e$repeats<-1:2;e$active_outer<-1:5;e$workers<-56L
  e$cfg<-e$perf_config(root,"confirmation")
  e$cfg$methods<-c("lambert_c1","lasso","ridge","elastic_net","mcp","scad")
  e$cfg$elastic_net_alpha<-c(0,.25,.5,.75,1);e$cfg$lambert_shapes<-1;e$cfg$mcp_gamma<-3;e$cfg$scad_a<-3.7
  e$cfg$lambda_points<-60L;e$cfg$early_stop<-FALSE
  e$rescue<-list(chunk_sweeps=250L,total_sweeps=10000L,newton_steps=20L,backtracks=20L,armijo=1e-4)
  if(stage=="smoke") {
    e$view<-list(rows=1:40,columns=1:24);e$folds<-sch_folds(40L,1L);e$repeats<-1L;e$active_outer<-1:2;e$workers<-2L
    e$cfg$lambda_points<-4L;e$cfg$max_coordinate_sweeps_per_start<-200L;e$cfg$implementation$ncvreg_max_iter<-200L
    e$rescue$chunk_sweeps<-50L;e$rescue$total_sweeps<-200L;e$rescue$newton_steps<-5L
  }
  e$settings<-data.frame(family=c("lambert",rep("enet",5),"mcp","scad"),shape=c(1,1,0,.25,.5,.75,3,3.7))
  e$tasks<-sch_bind(lapply(e$repeats,function(r)sch_bind(lapply(e$active_outer,function(k)
    data.frame(repeat_id=r,outer_fold=k,role=rep(0:5,each=8),setting=rep(1:8,6))))))
  e$tasks<-cbind(unit=seq_len(nrow(e$tasks)),e$tasks);e$mapping<-cpk_mapping(e$tasks$unit)
  e$sources<-rfh_files(root,sch_sources());e$inputs<-rfh_files(root,c(sch_raw(),file.path(sch_prepared(),c("design.rds","manifest.csv"))))
  e$software<-rfh_software()
  e$scientific_signature<-sch_text_hash(c("scheetz_cv_ic_v1",stage,e$sources$sha256,e$inputs$sha256))
  e$execution_signature<-sch_text_hash(c(e$scientific_signature,as.character(e$workers),"two_adjacent_packs_v1"))
  e$reference_path<-e$perf_path
  e$lambert_values_cpp<-function(t,lambda,shape)list(penalty=e$lambert_penalty(t,lambda,shape),
    slope=e$lambert_slope(t,lambda,shape),curvature=e$lambert_curvature(t,lambda,shape))
  e$perf_path<-function(prep,setting,cfg,pulse=function(...)NULL,last=cfg$lambda_points)
    if(setting$family=="enet" && setting$shape==0)rfh_ridge_path(prep,setting,cfg,e,last)
    else e$pacc_path(prep,setting,cfg,pulse,last,e)
  if(compile)rfh_compile(e)
  e
}
sch_split <- function(task,e) {
  outer<-e$folds[[task$repeat_id]]$outer;train<-which(outer!=task$outer_fold);test<-which(outer==task$outer_fold)
  if(task$role>0L) {
    inner<-e$folds[[task$repeat_id]]$inner[[task$outer_fold]]
    test<-train[inner==task$role];train<-train[inner!=task$role]
  }
  list(train=e$view$rows[train],test=e$view$rows[test],columns=e$view$columns)
}
sch_prep <- function(task,e) {
  s<-sch_split(task,e);e$perf_prepare(e$data$X[s$train,s$columns,drop=FALSE],e$data$y[s$train],e$cfg)
}
sch_newton <- function(beta,prep,lambda,setting,e) {
  active<-which(beta!=0);n<-nrow(prep$Z)
  if(setting$family!="enet" || setting$shape>=1 || length(active)<=n)
    return(sir_newton(beta,prep,lambda,setting,e))
  # Exact Woodbury solve for the same positive ENET Hessian, with an n x n system.
  Z<-prep$Z[,active,drop=FALSE];tau<-lambda*(1-setting$shape)
  G<-tcrossprod(Z)/n;diag(G)<-diag(G)+tau;ch<-chol(G)
  signs<-sign(beta[active]);steps<-0L;status<-"step_budget_exhausted"
  for(step in seq_len(e$rescue$newton_steps)) {
    dd<-e$perf_diagnose(beta,prep,lambda,setting,e$cfg)
    if(dd$ok) {status<-"kkt_met";break}
    gradient<-drop(crossprod(Z,drop(prep$Z%*%beta)-prep$y))/n+
      lambda*(setting$shape*signs+(1-setting$shape)*beta[active])
    dual<-backsolve(ch,forwardsolve(t(ch),drop(Z%*%gradient)))
    direction<--(gradient-drop(crossprod(Z,dual))/n)/tau;descent<-sum(gradient*direction)
    if(!is.finite(descent)||descent>=0) {status<-"non_descent_direction";break}
    accepted<-FALSE
    for(backtrack in 0:e$rescue$backtracks) {
      a<-2^(-backtrack);candidate<-beta;candidate[active]<-beta[active]+a*direction
      if(any(sign(candidate[active])!=signs))next
      next_dd<-e$perf_diagnose(candidate,prep,lambda,setting,e$cfg)
      if(is.finite(next_dd$objective) && next_dd$objective<=dd$objective+e$rescue$armijo*a*descent) {
        accepted<-TRUE;break
      }
    }
    if(!accepted) {status<-"line_search_failed";break}
    beta<-candidate;steps<-steps+1L
  }
  list(beta=beta,steps=steps,status=status)
}
sch_hybrid <- function(init,prep,lambda,setting,e) {
  # Reuse the frozen hybrid for Lambert/MCP/SCAD and LASSO.
  if(setting$family!="enet" || setting$shape==1)return(sir_fit(init,prep,lambda,setting,e))
  stopifnot(setting$shape>0,setting$shape<1)
  # Only the restartable coordinate block differs: same ENET objective, supplied alpha.
  f<-sir_fit;changes<-0L;private<-new.env(parent=environment(sir_fit));private$sir_newton<-sch_newton
  environment(f)<-private
  replace_alpha<-function(x) {
    if(!is.call(x))return(x)
    if(identical(x[[1]],quote(`::`)))return(x)
    if(identical(x[[1]],quote(ncvreg::ncvfit))) {x$alpha<-quote(setting$shape);changes<<-changes+1L}
    for(i in seq_along(x)[-1])x[i]<-list(replace_alpha(x[[i]]))
    x
  }
  body(f)<-replace_alpha(body(f));stopifnot(changes==1L);f(init,prep,lambda,setting,e)
}
sch_rescue <- function(old,prep,setting,e) {
  path<-old;attempts<-list();warm<-numeric(ncol(prep$Z))
  for(k in seq_len(e$cfg$lambda_points)) {
    if(old$diagnostics$ok[k]) {warm<-old$beta[,k];next}
    if(setting$family=="enet" && setting$shape==0)next
    init<-list(warm=warm+0,zero=numeric(length(warm)));lambda<-old$diagnostics$lambda[k]
    fits<-lapply(init,function(b)sch_hybrid(b,prep,lambda,setting,e))
    good<-which(vapply(fits,function(f)f$ok,FALSE));chosen<-NA_integer_
    if(length(good)) {
      q<-vapply(fits[good],function(f)f$objective,1.0);best<-min(q)
      chosen<-good[which(q<=best+e$cfg$implementation$tie_tol*(1+abs(best)))[1]];fit<-fits[[chosen]]
      path$beta[,k]<-fit$beta;warm<-fit$beta;d<-path$diagnostics
      d$ok[k]<-TRUE;d$status[k]<-"verified_stationary";d$kkt[k]<-fit$kkt;d$objective[k]<-fit$objective
      d$iterations[k]<-fit$sweeps;d$selected_start[k]<-paste0("rescue_",names(fits)[chosen]);path$diagnostics<-d
    }
    attempts[[as.character(k)]]<-list(initial=init,fits=fits,chosen=chosen)
  }
  list(path=path,attempts=attempts)
}
sch_candidates <- function(path,prep,setting,e) {
  z<-lbc_candidates(path,prep,setting,e$cfg,e,new.env(parent=emptyenv()))
  z$eligible_BIC<-z$eligible;z$eligible_sparse<-FALSE
  z$EBIC_0.5<-z$EBIC_1<-z$HBIC<-NA_real_
  if(!(setting$family=="enet" && setting$shape==0)) {
    scores<-e$lhx_scores(z$rss,z$d,nrow(prep$Z),ncol(prep$Z))
    for(nm in c("EBIC_0.5","EBIC_1","HBIC"))z[[nm]]<-scores[,nm]
    z$eligible_sparse<-z$ok & is.finite(z$rss) & z$rss>0 & z$d<=floor(nrow(prep$Z)/log(nrow(prep$Z))) &
      rfh_curvature_ok(z$curvature_status)
  }
  z
}
sch_unit <- function(i,e) {
  began<-proc.time()[[3]];task<-e$tasks[i,,drop=FALSE];split<-sch_split(task,e);prep<-sch_prep(task,e)
  setting<-e$settings[task$setting,,drop=FALSE];baseline<-e$perf_path(prep,setting,e$cfg)
  rescued<-sch_rescue(baseline,prep,setting,e);path<-rescued$path
  curvature<-rep("not_assessed",e$cfg$lambda_points)
  for(k in which(path$diagnostics$ok))curvature[k]<-rfh_curvature(path$beta[,k],prep,path$diagnostics$lambda[k],setting,e)$status
  valid<-path$diagnostics$ok & rfh_curvature_ok(curvature)
  predictions<-matrix(NA_real_,length(split$test),e$cfg$lambda_points)
  if(task$role>0L)for(k in which(valid))predictions[,k]<-e$perf_predict(path$beta[,k],prep,
    e$data$X[split$test,split$columns,drop=FALSE])
  list(task=task,split=split,preprocessing=prep[c("mx","sx","my","lambda0","spectral")],
    baseline=baseline,path=path,rescue=rescued$attempts,curvature=curvature,valid=valid,
    candidates=if(task$role==0L)sch_candidates(path,prep,setting,e) else NULL,
    predictions=if(task$role>0L)predictions else NULL,
    loss=if(task$role>0L)unname(colMeans((e$data$y[split$test]-predictions)^2)) else NULL,
    elapsed=proc.time()[[3]]-began)
}
