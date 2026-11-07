# Independent analysis adapter. Historical numerical functions stay unchanged.
lhx_input <- function() paste0("incoming/lambert_bic_truba_v1_a01_results_20261017/",
  "extracted/outputs/performance/lambert_bic_simulation_truba_v1_a01")
lhx_sources <- function() c("protocol/lambert_hbic_ebic_reanalysis_v1.md",
  "R/lambert_hbic_ebic_reanalysis_v1.R","scripts/analyze_lambert_hbic_ebic_v1.R",
  "tests/lambert_hbic_ebic_v1_checks.R")
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
lhx_manifest_check <- function(path) {
  mf <- read.csv(file.path(path,"manifest.csv"),stringsAsFactors=FALSE)
  stopifnot(!anyDuplicated(mf$file))
  for(i in seq_len(nrow(mf))) {
    f <- file.path(path,mf$file[i]);stopifnot(file.exists(f),cpk_hash(f)==mf$sha256[i])
    if("bytes" %in% names(mf))stopifnot(file.info(f)$size==mf$bytes[i])
  }
  invisible(mf)
}
lhx_environment <- function(root) {
  e <- new.env(parent=globalenv())
  for(f in c("R/penalties.R","R/lambert_family.R","R/performance_models.R","R/performance_run.R",
    "R/lambert_fixed_metrics_v1.R","R/lambert_ic_pilot_v1.R","R/lambert_bic_comparison_v1.R"))sys.source(file.path(root,f),e)
  e$lambert_values_cpp <- function(t,lambda,shape)list(penalty=e$lambert_penalty(t,lambda,shape),
    slope=e$lambert_slope(t,lambda,shape),curvature=e$lambert_curvature(t,lambda,shape))
  for(nm in ls(e)[grepl("^perf_(path|train_models|run_task|fit)|^lbc_oracle$",ls(e))])
    assign(nm,function(...)stop("Model fitting is forbidden in archived-path reanalysis"),e)
  e$input <- file.path(root,lhx_input());e$meta <- readRDS(file.path(e$input,"metadata.rds"))
  e$cfg <- e$meta$cfg;e$settings <- e$perf_settings(e$cfg)
  proof <- file.path(root,"outputs/validation/lambert_bic_import_v1_a01")
  lhx_manifest_check(proof);receipt <- jsonlite::fromJSON(file.path(proof,"receipt.json"))
  stopifnot(receipt$accepted,receipt$datasets==4800,receipt$rows==62400,
    receipt$input_manifest_sha256==cpk_hash(file.path(e$input,"manifest.csv")),
    receipt$scientific_signature==e$meta$scientific_signature)
  e$mf <- read.csv(file.path(e$input,"manifest.csv"));pm <- read.delim(file.path(e$input,"manifest.tsv"))
  stopifnot(!anyDuplicated(e$mf$file),all(pm$validated),sum(pm$kind=="pack")==10L,
    all(pm$records[pm$kind=="pack"]==480L),identical(pm$sha256,e$mf$sha256[match(pm$file,e$mf$file)]))
  # Large packs are hashed once immediately before reading.
  for(i in which(!grepl("^checkpoints/",e$mf$file)))
    stopifnot(cpk_hash(file.path(e$input,e$mf$file[i]))==e$mf$sha256[i])
  port <- file.path(root,"outputs/diagnostics/lambert_bic_receipt_portability_v1_a01")
  hf <- file.path(port,"record_hashes_r43.csv")
  stopifnot(cpk_hash(hf)==receipt$record_receipts$evidence_sha256,
    cpk_hash(file.path(port,"verify_r43_receipts.R"))==receipt$record_receipts$verifier_sha256,
    cpk_hash(file.path(port,"r43_receipts_runtime.txt"))==receipt$record_receipts$runtime_receipt_sha256)
  e$rh <- read.csv(hf);stopifnot(nrow(e$rh)==4800,!anyDuplicated(e$rh$unit),
    all(e$rh$accepted),identical(e$rh$actual_sha256,e$rh$expected_sha256))
  for(i in seq_len(nrow(e$meta$sources)))stopifnot(cpk_hash(file.path(root,e$meta$sources$file[i]))==e$meta$sources$sha256[i])
  stopifnot(identical(e$perf_seeds(e$cfg,e$meta$tasks),e$meta$seeds),
    nrow(e$meta$tasks)==4800L,all(table(e$meta$tasks$scenario)==300L),
    readLines(file.path(e$input,"COMPLETED"))==e$meta$scientific_signature,
    !file.exists(file.path(e$input,"failure_evidence.rds")))
  e$base <- readRDS(file.path(e$input,"results.rds"));e$base_by_unit <- split(e$base,e$base$unit)
  stopifnot(nrow(e$base)==62400L,all(e$base$ok),!anyDuplicated(e$base[c("unit","method","criterion")]))
  e$source_hashes <- setNames(vapply(file.path(root,lhx_sources()),cpk_hash,""),lhx_sources())
  e$contract <- list(schema="archived_hbic_ebic_v1",input_manifest=receipt$input_manifest_sha256,
    record_receipt=receipt$record_receipts$evidence_sha256,sources=e$source_hashes,
    methods=lhx_methods(),criteria=c("HBIC","EBIC_0.5","EBIC_1"),model_fits=0L)
  e
}
lhx_pack <- function(e,batch) {
  file <- cpk_pack_file(batch);w <- match(file,e$mf$file);path <- file.path(e$input,file)
  stopifnot(file.info(path)$size==e$mf$bytes[w],cpk_hash(path)==e$mf$sha256[w])
  p <- readRDS(path);ids <- e$meta$mapping$unit[e$meta$mapping$batch==batch];h <- e$rh[e$rh$batch==batch,]
  stopifnot(identical(p$schema,"checkpoint_pack_v1"),identical(p$batch,as.integer(batch)),
    identical(p$scientific_signature,e$meta$scientific_signature),
    identical(p$execution_signature,e$meta$execution_signature),identical(p$units,ids),
    identical(names(p$records),as.character(ids)),identical(p$receipts$unit,ids),all(p$receipts$accepted),
    identical(h$unit,ids),identical(h$actual_sha256,p$receipts$record_sha256),all(h$pack_sha256==e$mf$sha256[w]))
  p
}
lhx_compare <- function(a,b,tol=1e-6) {
  stopifnot(identical(names(a),names(b)))
  error <- abs(as.numeric(a[1,])-as.numeric(b[1,]))/(1+abs(as.numeric(b[1,])))
  stopifnot(all(is.finite(error)),all(error<=tol));max(error)
}
lhx_unit <- function(z,e) {
  task <- z$task;cfg <- e$cfg;i <- task$unit;old <- e$base_by_unit[[as.character(i)]]
  stopifnot(identical(task,e$meta$tasks[i,,drop=FALSE]),identical(z$seeds,e$meta$seeds[[i]]),
    identical(z$signature,e$meta$scientific_signature),identical(z$execution_signature,e$meta$execution_signature),
    all(z$cv_audit$passed),all(z$rows$ok),all(z$rows$kkt<=cfg$normalized_kkt_tolerance))
  zr <- z$rows[order(z$rows$method,z$rows$criterion),];rownames(zr)<-NULL
  oo <- old[order(old$method,old$criterion),];rownames(oo)<-NULL;stopifnot(identical(zr,oo))
  data <- e$perf_training(task,cfg,z$seeds);prep <- e$perf_prepare(data$X,data$y,cfg)
  test <- e$perf_test_data(task,cfg,z$seeds,data$truth);e$lbc_oracle_check(z$oracle,prep,cfg)
  fields <- e$lic_metric_fields();new <- paths <- audits <- list()
  for(method in lhx_methods()) {
    j <- e$lbc_settings_for(e$settings,method);stopifnot(length(j)==1)
    pp <- z$paths[[j]];cc <- z$candidates[[j]];setting <- e$settings[j,]
    stopifnot(identical(dim(pp$beta),c(task$p,cfg$lambda_points)),identical(cc$index,seq_len(cfg$lambda_points)),
      identical(cc$lambda,pp$diagnostics$lambda),all(diff(cc$lambda)<0),
      all(cc$complexity=="active_count"),!anyNA(cc$eligible),all(cc$ok[cc$eligible]))
    g <- which(cc$ok);beta <- sweep(pp$beta[,g,drop=FALSE],1,prep$sx,"/")
    active <- abs(beta)>cfg$support_reporting_tolerance;d <- colSums(active)
    rss <- colSums((prep$y-prep$Z%*%pp$beta[,g,drop=FALSE])^2)
    rss_err <- max(abs(rss-cc$rss[g])/(1+abs(cc$rss[g])))
    stopifnot(all(d==cc$d[g]),all(cc$q[g]==d),rss_err<=cfg$bic$reproduction_tolerance,
      identical(cc$cap_excluded[g],d>floor(task$n/log(task$n))),
      all(cc$curvature_status[cc$eligible] %in% c("convex_psd","empty_active","nonnegative_local","nonnegative_join_lower")))
    scores <- lhx_scores(cc$rss,cc$d,task$n,task$p)
    stopifnot(max(abs(scores[g,"BIC"]-cc$BIC[g])/(1+abs(cc$BIC[g])))<1e-12)
    kk <- vapply(colnames(scores),function(nm)e$lic_choose(scores[,nm],cc$eligible,cfg$implementation$tie_tol),integer(1))
    bic <- old[old$method==method & old$criterion=="BIC",]
    stopifnot(length(kk)==4L,!anyNA(kk),kk["BIC"]==bic$index,j==bic$setting)
    metrics <- list();diffs <- numeric();curvature <- character()
    for(k in unique(kk)) {
      b <- pp$beta[,k];dd <- e$perf_diagnose(b,prep,cc$lambda[k],setting,cfg)
      h <- if(setting$family=="enet")list(status="convex_psd") else e$lic_curvature(b,prep,cc$lambda[k],setting,cfg,e)
      stopifnot(dd$ok,abs(dd$objective-cc$objective[k])<=1e-10*(1+abs(dd$objective)),
        h$status %in% c("convex_psd","empty_active","nonnegative_local","nonnegative_join_lower"))
      metrics[[as.character(k)]] <- e$lic_metrics(b,prep,test,data$truth,cfg,e)[fields]
      diffs[as.character(k)] <- max(abs((b-z$oracle$beta)/prep$sx))/(1+max(abs(z$oracle$beta/prep$sx)))
      curvature[as.character(k)] <- h$status
    }
    met_err <- lhx_compare(metrics[[as.character(kk["BIC"])]],bic[fields],cfg$bic$reproduction_tolerance)
    target <- seq_len(task$p) %in% cfg$support_indices
    exact <- colSums(active!=target)==0;path_exact <- any(exact & cc$eligible[g])
    near <- vapply(seq_along(g),function(l)max(abs(beta[,l]-z$oracle$beta/prep$sx))/
      (1+max(abs(z$oracle$beta/prep$sx)))<=1e-6,logical(1))
    paths[[method]] <- cbind(task,method=method,eligible=sum(cc$eligible),
      support_on_path=path_exact,oracle_near_on_path=any(near & cc$eligible[g]),
      exact_support_points=sum(exact & cc$eligible[g]),partial_search=bic$partial_search)
    for(criterion in setdiff(colnames(scores),"BIC")) {
      k <- kk[criterion];r <- bic;r$criterion <- criterion;r$index <- unname(k)
      r$lambda <- cc$lambda[k];r$q <- cc$q[k];r$bic <- cc$BIC[k];r$kkt <- cc$kkt[k]
      r$curvature_status <- curvature[as.character(k)];r$cap_boundary <- cc$d[k]==floor(task$n/log(task$n))
      r$path_endpoint <- k %in% c(1L,cfg$lambda_points);r[fields] <- metrics[[as.character(k)]]
      r$criterion_score <- scores[k,criterion];r$support_on_path <- path_exact
      r$oracle_relative_max_difference <- diffs[as.character(k)]
      r$oracle_numerical_agreement <- diffs[as.character(k)]<=1e-6
      stopifnot(!r$exact_support || path_exact);new[[paste(method,criterion)]] <- r
    }
    audits[[method]] <- data.frame(unit=i,method=method,rss_relative_error=rss_err,
      bic_metric_relative_error=met_err,bic_choice_unchanged=TRUE,selected_kkt_curvature_passed=TRUE)
  }
  list(rows=do.call(rbind,new),paths=do.call(rbind,paths),audit=do.call(rbind,audits))
}
lhx_summarize <- function(rows) {
  metrics <- c("normalized_mse","mcc","fdr","tpr","coefficient_error","excess_risk","exact_support","selected_size")
  groups <- split(rows,interaction(rows$scenario,rows$method,rows$criterion,drop=TRUE))
  ans <- lapply(groups,function(z)do.call(rbind,lapply(metrics,function(m) {
    v <- as.numeric(z[[m]]);v[!z$ok]<-NA_real_;s<-lfb_summary(v,300L,TRUE)
    cbind(z[1,c("scenario","n","p","rho","snr","profile","method","criterion")],metric=m,
      s[c("expected","available","complete","mean","mcse","conditional_mean")])
  })))
  x<-do.call(rbind,ans);rownames(x)<-NULL;x
}
lhx_paired <- function(rows) {
  pairs <- list();id<-0L
  for(s in unique(rows$scenario))for(metric in c("normalized_mse","mcc","fdr","tpr","exact_support")) {
    dat<-rows[rows$scenario==s,]
    compare<-function(a,b,kind,scope) {
      z<-merge(a,b,by="unit",suffixes=c(".a",".b"));good<-z$ok.a & z$ok.b
      if(scope=="complete_search_pairs")good<-good & !z$partial_search.a & !z$partial_search.b
      d<-as.numeric(z[[paste0(metric,".a")]][good])-as.numeric(z[[paste0(metric,".b")]][good])
      id<<-id+1L;pairs[[id]]<<-data.frame(scenario=s,p=dat$p[1],rho=dat$rho[1],snr=dat$snr[1],
        metric=metric,comparison=kind,scope=scope,method_a=a$method[1],criterion_a=a$criterion[1],
        method_b=b$method[1],criterion_b=b$criterion[1],expected=300L,available=length(d),
        mean_difference=if(length(d))mean(d) else NA_real_,mcse=if(length(d)>1)sd(d)/sqrt(length(d)) else NA_real_)
    }
    for(scope in c("all_pairs","complete_search_pairs")) {
      for(m in lhx_methods())for(c in c("HBIC","EBIC_0.5","EBIC_1"))for(b in c("CV","BIC"))
        compare(dat[dat$method==m & dat$criterion==c,],dat[dat$method==m & dat$criterion==b,],"criterion",scope)
      for(c in c("CV","BIC","HBIC","EBIC_0.5","EBIC_1"))for(m in c("lasso","mcp","scad"))
        compare(dat[dat$method=="lambert_c1" & dat$criterion==c,],dat[dat$method==m & dat$criterion==c,],"method",scope)
    }
  }
  do.call(rbind,pairs)
}
