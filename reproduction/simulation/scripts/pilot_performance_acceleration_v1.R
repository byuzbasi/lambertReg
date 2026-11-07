#!/usr/bin/env Rscript
# Bounded eight-case executable smoke, never a production simulation runner.
pacc_pilot_sources <- c("R/performance_acceleration_v1.R","src/lambert_performance_fast_v1.cpp",
  "scripts/pilot_performance_acceleration_v1.R","tests/performance_acceleration_v1_checks.R",
  "protocol/performance_acceleration_v1.md","R/penalties.R","R/lambert_family.R",
  "R/performance_models.R","R/performance_run.R","src/lambert_performance.cpp",
  "R/lambert_fixed_simulation_v1.R","R/lambert_fixed_metrics_v1.R","protocol/lambert_superiority_v2.md")

pacc_compare_result <- function(a,b,cfg) {
  eq <- function(x,y) isTRUE(all.equal(x,y,tolerance=1e-11,check.attributes=FALSE))
  stopifnot(identical(a$data_hashes,b$data_hashes),identical(a$seeds,b$seeds),
    identical(a$fit$preprocessing,b$fit$preprocessing))
  fields <- setdiff(names(a$rows),c("tuning_seconds","refit_seconds","partial_search"))
  stopifnot(eq(a$rows[fields],b$rows[fields]))
  for (m in cfg$methods) {
    x <- a$fit$models[[m]]; y <- b$fit$models[[m]]
    stopifnot(eq(x$selection[c("ok","eligible","setting","index","loss")],
      y$selection[c("ok","eligible","setting","index","loss")]))
    if (x$ok) stopifnot(eq(x$beta,y$beta),eq(x$intercept,y$intercept),
      identical(abs(x$beta)>cfg$support_reporting_tolerance,abs(y$beta)>cfg$support_reporting_tolerance))
  }
  aa <- a$fit$cv; bb <- b$fit$cv
  stopifnot(identical(aa[c("fold","setting","index")],bb[c("fold","setting","index")]))
  computed <- bb$status!="certified_cv_pruned"
  fields <- setdiff(names(aa),"elapsed")
  stopifnot(eq(aa[computed,fields],bb[computed,fields]))
  cert <- b$fit$acceleration$certificates
  if (any(!computed)) {
    stopifnot(!is.null(cert),nrow(cert)==sum(!computed),all(cert$fold==cfg$folds))
    for (i in seq_len(nrow(cert))) {
      losses <- aa$loss[aa$setting==cert$setting[i] & aa$index==cert$index[i]]
      # A failed reference candidate is also ineligible. Every computable
      # pruned candidate must be worse than the recorded full-fold incumbent.
      if (all(is.finite(losses))) stopifnot(mean(losses)>cert$incumbent[i]+
        cfg$implementation$tie_tol*(1+cert$incumbent[i]))
    }
  }
  invisible(TRUE)
}

pacc_pilot_main <- function(root,args) {
  source(file.path(root,"R/performance_acceleration_v1.R"))
  opt <- list(mode="design",out=NULL,resume="false",`stop-after`="Inf")
  stopifnot(length(args)%%2L==0)
  if (length(args)) {
    keys <- sub("^--","",args[seq(1,length(args),2)])
    stopifnot(!anyDuplicated(keys),all(keys %in% names(opt)))
    opt[keys] <- as.list(args[seq(2,length(args),2)])
  }
  stopifnot(opt$mode %in% c("design","smoke","validate"),opt$resume %in% c("true","false"))
  stop_after <- as.numeric(opt$`stop-after`)
  stopifnot(length(stop_after)==1L,!is.na(stop_after),stop_after>=1,
    is.infinite(stop_after) || (stop_after==floor(stop_after) && stop_after<=8L && opt$mode=="smoke"))
  base <- pacc_environment(root,FALSE,compile=FALSE)
  cfg <- base$perf_config(root,"smoke"); original <- base$perf_tasks(cfg)
  old_seeds <- base$perf_seeds(cfg,original)
  tasks <- list(); seeds <- list(); i <- 0L
  # Preserve the old smoke streams; use common innovations across rho within p.
  # These are execution checks, not independent Monte Carlo replications.
  for (p in cfg$p) for (rho in c(0,.3,.6,.9)) {
    i <- i+1L; k <- which(original$p==p)[1]
    task <- original[k,,drop=FALSE]; task$unit <- i; task$scenario <- i; task$rho <- rho
    tasks[[i]] <- task; seeds[[i]] <- old_seeds[[k]]
  }
  tasks <- do.call(rbind,tasks); rownames(tasks) <- NULL
  stopifnot(nrow(tasks)==8L,cfg$n_train==40,cfg$lambda_points==4,cfg$folds==5)
  if (opt$mode=="design") { print(tasks); cat("SMOKE ONLY: 8 paired cases, 3 execution variants; production runs unsupported.\n"); return(invisible(NULL)) }
  stopifnot(!is.null(opt$out))
  out <- if (startsWith(opt$out,"/")) opt$out else file.path(root,opt$out)
  sources <- data.frame(file=pacc_pilot_sources,
    bytes=file.info(file.path(root,pacc_pilot_sources))$size,
    sha256=unname(vapply(file.path(root,pacc_pilot_sources),base$perf_hash,character(1))))
  meta <- list(schema="performance_acceleration_v1_smoke",cfg=cfg,tasks=tasks,seeds=seeds,
    sources=sources,software=base$perf_software(),variants=c("reference","fast","fast_early"),
    scientific_scope="small execution equivalence checks; not simulation evidence")
  signature <- digest::digest(serialize(meta,NULL,version=3),algo="sha256",serialize=FALSE)
  read_unit <- function(i) {
    path <- file.path(out,"shards",sprintf("unit_%05d",i)); file <- file.path(path,"result.rds")
    stopifnot(file.exists(file),identical(readLines(file.path(path,"result.sha256")),base$perf_hash(file)))
    z <- readRDS(file)
    stopifnot(identical(z$signature,signature),identical(z$unit,as.integer(i)))
    pacc_compare_result(z$results$reference,z$results$fast,cfg)
    pacc_compare_result(z$results$reference,z$results$fast_early,cfg)
    z
  }
  validate <- function() {
    stopifnot(identical(readLines(file.path(out,"COMPLETED")),signature))
    base$perf_verify_manifest(out)
    z <- lapply(tasks$unit,read_unit)
    expected <- do.call(rbind,lapply(z,function(x)x$summary))
    stopifnot(isTRUE(all.equal(read.csv(file.path(out,"summary.csv")),expected,check.attributes=FALSE,tolerance=1e-12)))
    progress <- jsonlite::fromJSON(file.path(out,"progress.json"))
    stopifnot(progress$phase=="completed",progress$validated_completed==8,progress$running==0,
      progress$failed==0,progress$pending==0)
    cat("ACCELERATION_SMOKE_VALID: 8 paired cases, 24 execution results; model/selection/prediction/KKT/manifest checks passed.\n")
  }
  if (opt$mode=="validate") { validate(); return(invisible(NULL)) }
  if (file.exists(out) && opt$resume!="true") stop("Refusing existing output; use --resume true")
  if (!file.exists(out) && opt$resume=="true") stop("Cannot resume missing output")
  if (file.exists(file.path(out,"COMPLETED"))) { validate(); return(invisible(NULL)) }
  if (file.exists(paste0(out,".run_lock"))) stop("Existing run lock; inspect before recovery")
  dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE)
  stopifnot(dir.create(paste0(out,".run_lock")))
  on.exit(unlink(paste0(out,".run_lock"),recursive=TRUE),add=TRUE)
  writeLines(as.character(Sys.getpid()),file.path(paste0(out,".run_lock"),"pid"))
  if (opt$resume=="false") {
    dir.create(out); dir.create(file.path(out,"shards"))
    saveRDS(meta,file.path(out,"metadata.rds"),version=3)
    writeLines(signature,file.path(out,"signature.sha256"))
    write.csv(tasks,file.path(out,"tasks.csv"),row.names=FALSE)
    saveRDS(seeds,file.path(out,"seeds.rds"),version=3)
    write.csv(sources,file.path(out,"sources.csv"),row.names=FALSE)
  } else stopifnot(identical(readRDS(file.path(out,"metadata.rds")),meta),
    identical(readLines(file.path(out,"signature.sha256")),signature))
  dirs <- list.files(file.path(out,"shards"),all.files=TRUE,no..=TRUE)
  stopifnot(all(dirs %in% sprintf("unit_%05d",tasks$unit)))
  done <- tasks$unit[sprintf("unit_%05d",tasks$unit) %in% dirs]
  invisible(lapply(done,read_unit))
  state <- list(cfg=cfg,tasks=tasks,completed=done)
  progress <- base$perf_progress(out,state)
  progress$phase <- "smoke_preflight"; progress$emit("compiling reference and acceleration backends",TRUE)
  engines <- list(reference=pacc_environment(root,FALSE),fast=pacc_environment(root,TRUE),
    fast_early=pacc_environment(root,TRUE,TRUE))
  invocation_count <- 0L
  for (i in setdiff(tasks$unit,done)) {
    t0 <- proc.time()[[3]]; progress$unit <- i; progress$running <- 1L; progress$phase <- "smoke_computation"
    progress$emit("paired reference / fast / optional early-stop CV",TRUE)
    results <- list(); elapsed <- numeric(3); names(elapsed) <- names(engines)
    order <- names(engines)[((0:2+i-1L)%%3L)+1L]
    for (mode in order) {
      t1 <- proc.time()[[3]]
      r <- engines[[mode]]$perf_run_task(tasks[i,,drop=FALSE],cfg,seeds[[i]],progress$emit)
      r$rows <- base$lfb_add_metrics(r$rows)
      elapsed[mode] <- proc.time()[[3]]-t1; results[[mode]] <- r
    }
    pacc_compare_result(results$reference,results$fast,cfg)
    pacc_compare_result(results$reference,results$fast_early,cfg)
    summary <- do.call(rbind,lapply(names(engines),function(mode) {
      a <- results[[mode]]$fit$acceleration
      data.frame(unit=i,p=tasks$p[i],rho=tasks$rho[i],mode=mode,elapsed_seconds=elapsed[mode],
        speedup_vs_reference=elapsed["reference"]/elapsed[mode],
        failed_method_fits=sum(!results[[mode]]$rows$ok),
        reused_starts=if(is.null(a))0L else sum(a$paths$duplicate_starts),
        pruned_lambda=if(is.null(a))0L else sum(a$paths$pruned_lambda),
        preprocessing_seconds=if(is.null(a))NA_real_ else a$preprocessing_seconds,
        prediction_seconds=if(is.null(a))NA_real_ else a$prediction_seconds,
        solver_seconds=if(is.null(a))NA_real_ else sum(a$paths$solver_seconds),
        diagnostic_seconds=if(is.null(a))NA_real_ else sum(a$paths$diagnose_seconds,na.rm=TRUE),
        selected_diagnostics_seconds=if(is.null(a))NA_real_ else a$selected_diagnostics_seconds)
    }))
    rownames(summary) <- NULL
    z <- list(signature=signature,unit=as.integer(i),results=results,summary=summary,order=order)
    stage <- file.path(out,"shards",sprintf("unit_%05d.partial",i)); stopifnot(dir.create(stage))
    saveRDS(z,file.path(stage,"result.rds"),version=3)
    writeLines(base$perf_hash(file.path(stage,"result.rds")),file.path(stage,"result.sha256"))
    stopifnot(file.rename(stage,file.path(out,"shards",sprintf("unit_%05d",i))))
    read_unit(i)
    progress$completed <- c(progress$completed,i); progress$running <- 0L; progress$last_completed <- i
    progress$durations <- c(progress$durations,proc.time()[[3]]-t0)
    progress$cost_class <- c(progress$cost_class,with(tasks[i,],paste(p,rho,snr,profile,sep="/")))
    progress$emit("paired unit validated",TRUE); invocation_count <- invocation_count+1L
    if (invocation_count>=stop_after && length(progress$completed)<8L) {
      progress$phase <- "checkpointed"; progress$emit("intentional smoke checkpoint",TRUE)
      return(invisible(NULL))
    }
  }
  progress$phase <- "aggregation_and_validation"; progress$emit("validating all paired outputs",TRUE)
  z <- lapply(tasks$unit,read_unit); summary <- do.call(rbind,lapply(z,function(x)x$summary))
  write.csv(summary,file.path(out,"summary.csv"),row.names=FALSE)
  text <- c("Small execution-equivalence smoke; not a production benchmark.",
    "Same frozen smoke streams, methods and 5-fold CV; 4 lambda values; rho=0,.3,.6,.9; n=40; p=20,60.",
    "One timing per case/variant with rotated execution order; no general speedup or risk guarantee.",
    "Wall timings include data, all solvers, diagnostics and predictions; compile time excluded.",
    "Solver timings for ENET include its existing path diagnostics. Sub-timers are not an exhaustive partition of wall time.",
    "Reference and accelerated failed fits are retained. Certified pruning is distinct from convergence failure.",
    "early_stop defaults to FALSE; it can only discard a certified last-fold suffix for Lambert/MCP/SCAD.",
    "Selected fits retain reference KKT, objective, curvature and support gates.","",
    capture.output(knitr::kable(summary,digits=4,format="pipe")))
  writeLines(text,file.path(out,"report.txt"))
  progress$phase <- "completed"; progress$unit <- NA_integer_; progress$emit("all paired smoke cases verified",TRUE)
  base$perf_manifest(out); writeLines(signature,file.path(out,"COMPLETED")); validate()
}

if (sys.nframe()==0L) {
  driver <- normalizePath(sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1]))
  pacc_pilot_main(dirname(dirname(driver)),commandArgs(TRUE))
}
