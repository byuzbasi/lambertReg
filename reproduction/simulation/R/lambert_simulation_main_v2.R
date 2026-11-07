# Versioned integration only: established models and persistence stay untouched.
lms_sources <- function() c(
  "protocol/lambert_simulation_main_v2.md", "protocol/performance_acceleration_v1.md",
  "R/penalties.R", "R/lambert_family.R", "R/performance_models.R", "R/performance_run.R",
  "R/performance_truba_v1.R", "R/lambert_fixed_truba_v1.R", "R/lambert_fixed_metrics_v1.R",
  "R/performance_acceleration_v1.R", "R/lambert_simulation_main_v2.R",
  "src/lambert_performance.cpp", "src/lambert_performance_fast_v1.cpp",
  "scripts/pilot_performance_acceleration_v1.R", "scripts/lambert_simulation_main_v2.R",
  "scripts/run_lambert_main_v2_pilot.sh", "scripts/validate_lambert_simulation_main_v2.R",
  "scripts/analyze_lambert_simulation_main_v2.R", "scripts/analyze_lambert_fixed_benchmark_v1.R",
  "tests/lambert_simulation_main_v2_checks.R")

lms_config <- function(root, study_stage, backend) {
  lines <- readLines(file.path(root, "protocol/lambert_simulation_main_v2.md"))
  a <- which(lines == "```json"); b <- which(seq_along(lines) > a & lines == "```")[1]
  cfg <- jsonlite::fromJSON(paste(lines[(a+1):(b-1)], collapse="\n"))
  cfg$study_stage <- match.arg(study_stage, c("smoke", "pilot", "main"))
  cfg$stage <- if (study_stage == "main") "confirmation" else "smoke"
  cfg$execution <- list(backend=backend, early_stop=backend=="fast_early",
    paired_variants=if (backend=="paired") c("reference","fast","fast_early") else character())
  if (study_stage != "main") cfg$repetitions <- 1L
  if (study_stage == "smoke") {
    cfg$n_train <- 40L; cfg$p <- c(20L,60L); cfg$n_test <- 40L; cfg$lambda_points <- 4L
  }
  cfg
}

lms_tasks <- function(cfg) {
  if (cfg$study_stage == "main") {
    cells <- list(); k <- 0L
    for (p in cfg$p) for (rho in cfg$rho) for (snr in cfg$snr) for (profile in cfg$profiles) {
      k <- k+1L; cells[[k]] <- data.frame(p=p,rho=rho,snr=snr,profile=profile)
    }
    for (p in cfg$stress$p) for (rho in cfg$stress$rho) {
      k <- k+1L; cells[[k]] <- data.frame(p=p,rho=rho,snr=cfg$stress$snr,profile=cfg$stress$profile)
    }
    cells <- do.call(rbind,cells)
  } else {
    cells <- cfg$pilot_cases
    if (cfg$study_stage == "smoke") cells$p <- ifelse(cells$p==125,20L,60L)
  }
  scenario <- rep(seq_len(nrow(cells)),each=cfg$repetitions)
  z <- cells[scenario,,drop=FALSE]; u <- seq_len(nrow(z))
  ans <- data.frame(unit=u,scenario=scenario,repetition=rep(seq_len(cfg$repetitions),nrow(cells)),
    n=cfg$n_train,p=z$p,rho=z$rho,snr=z$snr,profile=z$profile,
    stream_index=cfg$stream_offsets[[cfg$study_stage]]+u)
  rownames(ans) <- NULL; ans
}

lms_certificates <- function(result, cfg, enabled) {
  cv <- result$fit$cv; pruned <- cv$status=="certified_cv_pruned"
  cert <- result$fit$acceleration$certificates
  if (!any(pruned)) { stopifnot(is.null(cert) || nrow(cert)==0L); return(invisible(TRUE)) }
  stopifnot(enabled, !is.null(cert), nrow(cert)==sum(pruned),all(cert$fold==cfg$folds),
    all(!cv$ok[pruned]),all(is.na(cv$loss[pruned])),
    all(cv$family[pruned] %in% c("lambert","mcp","scad")))
  key <- function(x) paste(x$fold,x$setting,x$index,sep="/")
  stopifnot(!anyDuplicated(key(cert)),setequal(key(cert),key(cv[pruned,])))
  for (j in unique(cert$setting)) {
    cc <- cert[cert$setting==j,,drop=FALSE]; cc <- cc[order(cc$index),]
    first <- min(cc$index)
    stopifnot(first>1L,identical(as.integer(cc$index),seq.int(first,cfg$lambda_points)))
    z <- cv[cv$setting==j,,drop=FALSE]
    values <- vapply(seq_len(first-1L),function(i) {
      loss <- z$loss[z$index==i]
      if (length(loss)==cfg$folds && all(is.finite(loss))) mean(loss) else Inf
    },numeric(1))
    incumbent <- min(values)
    lb <- vapply(cc$index,function(i) {
      loss <- z$loss[z$index==i & z$fold<cfg$folds]
      stopifnot(length(loss)==cfg$folds-1L,all(is.finite(loss)),all(loss>=0))
      sum(loss)/cfg$folds
    },numeric(1))
    buffer <- 64*.Machine$double.eps*pmax(1,abs(lb),abs(incumbent))
    # The path uses (sum(first K-1 losses)+last)/K; R's mean can differ by an ulp.
    # Check that recomputation agrees at roundoff scale, retaining the strict bound.
    stopifnot(is.finite(incumbent),length(unique(cc$incumbent))==1L,
      all(abs(cc$incumbent-incumbent)<=64*.Machine$double.eps*max(1,abs(incumbent))),
      isTRUE(all.equal(cc$lower_bound,lb,tolerance=1e-13)),
      isTRUE(all.equal(cc$roundoff_buffer,buffer,tolerance=1e-13)),
      all(lb>pmax(incumbent,cc$incumbent)+cfg$implementation$tie_tol*(1+incumbent)+buffer))
  }
  invisible(TRUE)
}

lms_environment <- function(root, study_stage="smoke", backend="reference", pilot=NULL) {
  stopifnot(study_stage %in% c("smoke","pilot","main"),
    backend %in% c("reference","fast","fast_early","paired"),
    study_stage != "pilot" || backend=="paired", study_stage != "main" || backend!="paired")
  e <- new.env(parent=environment(lms_environment))
  for (f in c("R/penalties.R","R/lambert_family.R","R/performance_models.R","R/performance_run.R",
      "R/performance_truba_v1.R","R/lambert_fixed_truba_v1.R","R/lambert_fixed_metrics_v1.R",
      "R/performance_acceleration_v1.R","scripts/pilot_performance_acceleration_v1.R"))
    sys.source(file.path(root,f),e)
  e$perf_config <- function(root,stage) {
    cfg <- lms_config(root,study_stage,backend); stopifnot(stage==cfg$stage); cfg
  }
  e$perf_tasks <- lms_tasks
  e$perf_source_paths <- e$perf_truba_source_paths <- lms_sources
  e$perf_truba_release_id <- function() "lambert_simulation_main_v2_preparation"
  e$perf_truba_scan <- e$lft_scan
  software <- e$perf_software
  e$perf_software <- function() {
    x <- software(); x$packages$mltools <- as.character(utils::packageVersion("mltools")); x
  }
  base_preflight <- e$perf_preflight; compiled <- FALSE; engines <- NULL
  e$perf_preflight <- function(compile=TRUE,root) {
    base_preflight(FALSE,root)
    for (pkg in names(wanted <- c(glmnet="5.0",ncvreg="3.16.0",mltools="0.3.5")))
      if (!requireNamespace(pkg,quietly=TRUE) || as.character(utils::packageVersion(pkg))!=wanted[[pkg]])
        stop("Required installed package version: ",pkg," ",wanted[[pkg]])
    if (compile && !compiled) {
      Rcpp::sourceCpp(file.path(root,"src/lambert_performance.cpp"),env=e,
        cacheDir=file.path(tempdir(),"lms_v2_rcpp"),showOutput=FALSE)
      if (backend %in% c("fast","fast_early")) {
        Rcpp::sourceCpp(file.path(root,"src/lambert_performance_fast_v1.cpp"),env=e,
          cacheDir=file.path(tempdir(),"lms_v2_rcpp"),showOutput=FALSE)
        e$lambert_cd_cpp <- e$lambert_cd_fast_v1; e$lambert_values_cpp <- e$lambert_values_fast_v1
      }
      if (backend=="paired") {
        # Variant environments use the same scientific configuration passed at fit time.
        engines <<- setNames(lapply(c("reference","fast","fast_early"),function(v)
          lms_environment(root,if(study_stage=="pilot") "main" else "smoke",v)),
          c("reference","fast","fast_early"))
        for (engine in engines) engine$perf_preflight(TRUE,root)
      }
      compiled <<- TRUE
    }
    invisible(TRUE)
  }
  e$perf_truba_runtime_preflight <- function(root,stage,require_release=FALSE,require_truba_runtime=FALSE) {
    e$perf_preflight(FALSE,root)
    if (require_release) e$perf_truba_validate_release(root)
    invisible(e$perf_software())
  }
  e$perf_truba_validate_release <- function(root,quiet=FALSE) {
    if (study_stage!="main") return(invisible(TRUE))
    if (is.null(pilot)) stop("Main run requires --pilot=<completed source-bound timing pilot>")
    gate <- lms_pilot_report(root,pilot,write=FALSE)
    if (!isTRUE(gate$accepted)) stop("Main run closed: pilot equivalence or selected-fit numerical gate failed")
    invisible(gate)
  }
  e$perf_truba_operational_guard <- function(stage,workers,max_seconds,task_limit,interval,...) {
    stopifnot(.Platform$OS.type=="unix",workers==floor(workers),is.finite(workers),workers>=1,
      workers<=if(study_stage=="smoke")2L else 56L,
      is.finite(max_seconds),max_seconds>0,max_seconds<=if(study_stage=="smoke")300 else 252000,
      length(task_limit)==1L,!is.na(task_limit),task_limit>=1,
      is.infinite(task_limit) || task_limit==floor(task_limit),
      is.finite(interval),interval>0,interval<=60)
    if (nzchar(Sys.getenv("SLURM_JOB_ID")))
      stop("TRUBA resource/bundle release follows the validated local timing pilot; do not submit this preparation copy")
  }
  if (backend %in% c("fast","fast_early")) {
    e$reference_path <- e$perf_path
    e$perf_path <- function(prep,setting,cfg,pulse=function(...)NULL,last=cfg$lambda_points)
      e$pacc_path(prep,setting,cfg,pulse,last,e)
    e$perf_train_models <- function(X,y,folds,cfg,pulse=function(...)NULL)
      e$pacc_train(X,y,folds,cfg,pulse,e,backend=="fast_early")
  }
  base_run <- e$perf_run_task; base_validate <- e$perf_validate_unit
  e$perf_run_task <- function(task,cfg,seeds,pulse=function(...)NULL) {
    if (backend!="paired") {
      z <- base_run(task,cfg,seeds,pulse); z$rows <- e$lfb_add_metrics(z$rows); return(z)
    }
    stopifnot(compiled,!is.null(engines))
    order <- names(engines)[((0:2+task$unit-1L)%%3L)+1L]
    results <- list(); timing <- list(); started <- proc.time()[[3]]
    for (v in order) {
      t0 <- proc.time()[[3]]
      results[[v]] <- engines[[v]]$perf_run_task(task,cfg,seeds,pulse)
      a <- results[[v]]$fit$acceleration
      timing[[v]] <- data.frame(unit=task$unit,variant=v,elapsed_seconds=proc.time()[[3]]-t0,
        failed_methods=sum(!results[[v]]$rows$ok),
        partial_searches=sum(results[[v]]$rows$partial_search),
        reused_starts=if(is.null(a))0L else sum(a$paths$duplicate_starts),
        pruned_lambda=if(is.null(a))0L else sum(a$paths$pruned_lambda))
    }
    z <- results$reference
    z$paired <- list(results=results,order=order,timing=do.call(rbind,timing))
    rownames(z$paired$timing) <- NULL; z$elapsed <- proc.time()[[3]]-started; z
  }
  validate_single <- function(z,task,cfg,seeds,signature,early) {
    base_validate(z,task,cfg,seeds,signature)
    fields <- c("tn","mcc","mcc_zero_denominator","excess_risk")
    stopifnot(isTRUE(all.equal(z$rows[fields],e$lfb_add_metrics(z$rows)[fields],tolerance=1e-13)))
    settings <- e$perf_settings(cfg); cv <- z$fit$cv
    stopifnot(all(cv$shape==settings$shape[cv$setting]),all(cv$family==settings$family[cv$setting]))
    for (m in cfg$methods) if (z$fit$models[[m]]$ok) {
      support <- abs(z$fit$models[[m]]$beta)>cfg$support_reporting_tolerance
      truth <- z$truth$beta!=0; r <- z$rows[z$rows$method==m,]
      stopifnot(r$tp==sum(support & truth),r$fp==sum(support & !truth),
        r$fn==sum(!support & truth),r$tn==sum(!support & !truth))
    }
    lms_certificates(z,cfg,early)
  }
  e$perf_validate_unit <- function(z,task,cfg,seeds,signature) {
    validate_single(z,task,cfg,seeds,signature,backend=="fast_early")
    if (backend=="paired") {
      rr <- z$paired$results; stopifnot(setequal(names(rr),c("reference","fast","fast_early")))
      data <- e$perf_training(task,cfg,seeds)
      test <- e$perf_test_data(task,cfg,seeds,data$truth)
      stopifnot(identical(z$rows,rr$reference$rows),identical(z$fit,rr$reference$fit))
      for (v in names(rr)) {
        rr[[v]]$signature <- signature
        validate_single(rr[[v]],task,cfg,seeds,signature,v=="fast_early")
        e$pacc_compare_result(rr$reference,rr[[v]],cfg)
        for (m in cfg$methods) if (rr$reference$fit$models[[m]]$ok) {
          a <- rr$reference$fit$models[[m]]; b <- rr[[v]]$fit$models[[m]]
          stopifnot(isTRUE(all.equal(a$intercept+drop(test$X%*%a$beta),
            b$intercept+drop(test$X%*%b$beta),tolerance=1e-11)))
        }
      }
      tm <- z$paired$timing
      stopifnot(nrow(tm)==3L,setequal(tm$variant,names(rr)),all(tm$elapsed_seconds>=0),
        all(is.finite(tm$elapsed_seconds)),all(tm$unit==task$unit))
    }
    invisible(TRUE)
  }
  e
}

lms_pilot_report <- function(root,input,write=TRUE) {
  cfg <- readRDS(file.path(input,"config.rds"))
  stopifnot(cfg$study_stage %in% c("smoke","pilot"),cfg$execution$backend=="paired")
  e <- lms_environment(root,cfg$study_stage,"paired")
  state <- e$perf_truba_scan(input,root,complete=TRUE,portable=TRUE,quiet=TRUE)
  timing <- do.call(rbind,lapply(state$tasks$unit,function(i)
    e$perf_truba_read_unit(input,i,state)$paired$timing))
  totals <- aggregate(elapsed_seconds~variant,timing,sum)
  totals$speedup <- totals$elapsed_seconds[totals$variant=="reference"]/totals$elapsed_seconds
  accepted <- cfg$study_stage=="pilot" && all(timing$failed_methods==0)
  receipt <- list(schema="lambert_main_v2_pilot_gate",accepted=accepted,
    study_stage=cfg$study_stage,equivalence_verified=TRUE,selected_fits_verified=all(timing$failed_methods==0),
    sources=state$cfg$protocol,source_inventory=e$perf_sources(root),
    input_manifest_sha256=e$perf_hash(file.path(input,"manifest.csv")),
    main_task_count=15600L,early_stop_default=FALSE,timing=timing,totals=totals)
  if (write) { print(totals); cat("PAIRED_EQUIVALENCE=TRUE; MAIN_PILOT_GATE=",accepted,"\n",sep="") }
  invisible(receipt)
}

lms_main <- function(root,args) {
  if (any(!grepl("^--[a-z-]+=.+$",args))) stop("Use --key=value options")
  keys <- sub("^--([^=]+)=.*$","\\1",args)
  x <- list(stage="main",action="design",backend="reference",out=NULL,pilot=NULL,
    workers="1",`max-seconds`="252000",`task-limit`="Inf",interval="30")
  stopifnot(!anyDuplicated(keys),all(keys %in% names(x)))
  x[keys] <- as.list(sub("^--[^=]+=","",args))
  if (x$stage %in% c("smoke","pilot") && !"backend" %in% keys) x$backend <- "paired"
  if (x$stage=="smoke" && !"max-seconds" %in% keys) x$`max-seconds` <- "300"
  stopifnot(x$action %in% c("design","preflight","start","resume","verify","portable-verify"))
  e <- lms_environment(root,x$stage,x$backend,x$pilot)
  stage <- if(x$stage=="main") "confirmation" else "smoke"
  if (x$action=="design") {
    cfg <- e$perf_config(root,stage); tasks <- e$perf_tasks(cfg)
    print(unique(tasks[c("scenario","n","p","rho","snr","profile")]))
    cat("DESIGN:",nrow(tasks),"units;",nrow(tasks)*length(cfg$methods),
      "reference-method rows; stage:",x$stage,"; backend:",x$backend,"; model fits: 0\n")
    return(invisible(0L))
  }
  if (x$action=="preflight") { e$perf_preflight(FALSE,root); print(e$perf_software()); return(invisible(0L)) }
  stopifnot(!is.null(x$out))
  out <- if(startsWith(x$out,"/"))x$out else file.path(root,x$out)
  if (x$action %in% c("verify","portable-verify")) {
    e$perf_truba_scan(out,root,complete=TRUE,portable=x$action=="portable-verify")
    if (x$backend=="paired") lms_pilot_report(root,out)
    return(invisible(0L))
  }
  workers <- as.numeric(x$workers); seconds <- as.numeric(x$`max-seconds`)
  limit <- as.numeric(x$`task-limit`); interval <- as.numeric(x$interval)
  e$perf_truba_operational_guard(stage,workers,seconds,limit,interval)
  if (x$action=="start" && file.exists(out)) stop("Fresh start refuses existing output")
  if (x$action=="resume" && !dir.exists(out)) stop("Resume requires existing output")
  if (x$stage=="main") e$perf_truba_validate_release(root)
  e$perf_preflight(TRUE,root)
  z <- e$perf_truba_run(root,out,stage,workers,seconds,limit,interval,resume=x$action=="resume")
  if (z$complete && x$backend=="paired") lms_pilot_report(root,out)
  invisible(if(isTRUE(z$complete))0L else 3L)
}
