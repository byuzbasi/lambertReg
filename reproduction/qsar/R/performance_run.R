# Versioned local persistence. Existing reliability-v1 sources stay frozen.
perf_hash <- function(path)digest::digest(file=path,algo="sha256")
perf_json <- function(x)as.character(jsonlite::toJSON(x,auto_unbox=TRUE,digits=17,
  pretty=TRUE,na="null",null="null"))
perf_atomic <- function(path,writer) {
  tmp<-paste0(path,".tmp-",Sys.getpid())
  if(file.exists(tmp))stop("Existing temporary file: ",tmp)
  writer(tmp)
  if(!file.rename(tmp,path))stop("Atomic rename failed: ",path)
}
perf_atomic_json <- function(x,path)perf_atomic(path,function(p)writeLines(perf_json(x),p))
perf_source_paths <- function()c("R/penalties.R","R/lambert_family.R","R/performance_models.R",
  "R/performance_run.R","src/lambert_performance.cpp","scripts/performance_v1.R",
  "tests/performance_v1_checks.R","protocol/performance_v1.md")
perf_software <- function() {
  pkgs<-c("glmnet","ncvreg","digest","jsonlite","knitr","Rcpp")
  vars<-c("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","BLIS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS","NUMEXPR_NUM_THREADS")
  list(R=R.version.string,platform=R.version$platform,
    packages=setNames(lapply(pkgs,function(p)as.character(utils::packageVersion(p))),pkgs),
    external=as.list(extSoftVersion()),libraries=.libPaths(),threads=as.list(Sys.getenv(vars)),
    compiler=system2(file.path(R.home("bin"),"R"),c("CMD","config","CXX"),stdout=TRUE))
}
perf_sources <- function(root) {
  files<-perf_source_paths()
  data.frame(file=files,bytes=as.numeric(file.info(file.path(root,files))$size),
    sha256=unname(vapply(file.path(root,files),perf_hash,character(1))))
}
perf_signature <- function(out) {
  files<-c("config.json","sources.csv","software.json","tasks.csv","seeds.rds")
  text<-paste(files,vapply(file.path(out,files),perf_hash,character(1)),sep=":")
  digest::digest(paste(text,collapse="\n"),algo="sha256",serialize=FALSE)
}
perf_metadata <- function(out,root,cfg,tasks,seeds) {
  saveRDS(cfg,file.path(out,"config.rds"),version=3)
  writeLines(perf_json(cfg),file.path(out,"config.json"))
  write.csv(perf_sources(root),file.path(out,"sources.csv"),row.names=FALSE)
  writeLines(perf_json(perf_software()),file.path(out,"software.json"))
  capture.output(sessionInfo(),.libPaths(),file=file.path(out,"session_info.txt"))
  write.csv(tasks,file.path(out,"tasks.csv"),row.names=FALSE)
  saveRDS(seeds,file.path(out,"seeds.rds"),version=3)
  writeLines(perf_signature(out),file.path(out,"signature.sha256"))
  files<-c("config.json","config.rds","sources.csv","software.json","tasks.csv",
    "seeds.rds","signature.sha256","session_info.txt")
  write.csv(data.frame(file=files,sha256=unname(vapply(file.path(out,files),perf_hash,character(1)))),
    file.path(out,"metadata_manifest.csv"),row.names=FALSE)
}
perf_verify_metadata <- function(out,root) {
  mf<-read.csv(file.path(out,"metadata_manifest.csv"),stringsAsFactors=FALSE)
  expected<-c("config.json","config.rds","sources.csv","software.json","tasks.csv",
    "seeds.rds","signature.sha256","session_info.txt")
  if(!identical(mf$file,expected)||anyDuplicated(mf$file))stop("Metadata file set mismatch")
  if(!identical(unname(vapply(file.path(out,mf$file),perf_hash,character(1))),mf$sha256))
    stop("Metadata checksum mismatch")
  cfg<-readRDS(file.path(out,"config.rds"))
  if(!identical(perf_json(cfg),perf_json(perf_config(root,cfg$stage)))||
     !identical(paste(readLines(file.path(out,"config.json")),collapse="\n"),perf_json(cfg)))
    stop("Configuration mismatch")
  sources<-read.csv(file.path(out,"sources.csv"),stringsAsFactors=FALSE); sources$bytes<-as.numeric(sources$bytes)
  if(!identical(sources,perf_sources(root)))stop("Source signature mismatch")
  if(!identical(paste(readLines(file.path(out,"software.json")),collapse="\n"),perf_json(perf_software())))
    stop("Software/thread environment mismatch")
  signature<-readLines(file.path(out,"signature.sha256"))
  if(!identical(signature,perf_signature(out)))stop("Scientific signature mismatch")
  tasks<-perf_tasks(cfg); stored<-read.csv(file.path(out,"tasks.csv"),stringsAsFactors=FALSE)
  if(!isTRUE(all.equal(tasks,stored,check.attributes=FALSE)))stop("Task grid mismatch")
  seeds<-readRDS(file.path(out,"seeds.rds"))
  if(!identical(seeds,perf_seeds(cfg,tasks)))stop("Frozen RNG states mismatch")
  list(cfg=cfg,tasks=tasks,seeds=seeds,signature=signature)
}
perf_unit_path <- function(out,i)file.path(out,"shards",sprintf("unit_%05d",i))
perf_validate_unit <- function(result,task,cfg,seeds,signature) {
  if(!identical(result$signature,signature)||!identical(result$task,task)||
     !identical(result$seeds,seeds))stop("Shard signature/task/seed mismatch")
  r<-result$rows; cv<-result$fit$cv
  if(nrow(r)!=length(cfg$methods)||!identical(r$method,cfg$methods)||
     any(r$unit!=task$unit)||any(r$scenario!=task$scenario)||any(r$repetition!=task$repetition))
    stop("Result dimensions or task keys mismatch")
  allowed<-c("verified_stationary","no_valid_cv","negative_curvature","solver_error",
    "budget_exhausted","kkt_not_met","objective_increased","nonfinite_prediction")
  if(anyNA(r$ok)||!all(r$status%in%allowed)||!all(r$ok==(r$status=="verified_stationary")))
    stop("Invalid scientific status")
  fields<-c("mse","normalized_mse","coefficient_error","analytic_risk","tpr","fdr")
  if(any(!is.finite(as.matrix(r[r$ok,fields,drop=FALSE])))||
     any(as.matrix(r[r$ok,fields,drop=FALSE])<0)||any(r$kkt[r$ok]>cfg$normalized_kkt_tolerance)||
     any(!is.na(r$mse[!r$ok])))stop("Invalid successful/failed fit metrics")
  J<-nrow(perf_settings(cfg)); L<-cfg$lambda_points
  keys<-paste(cv$fold,cv$setting,cv$index)
  expected<-expand.grid(index=seq_len(L),setting=seq_len(J),fold=seq_len(cfg$folds))
  if(nrow(cv)!=cfg$folds*J*L||anyDuplicated(keys)||
     !setequal(keys,paste(expected$fold,expected$setting,expected$index)))stop("CV coverage mismatch")
  if(any(cv$kkt[cv$ok]>cfg$normalized_kkt_tolerance)||
     any(!is.na(cv$loss[!cv$ok])))stop("CV failure was silently scored")
  if(length(result$fit$folds)!=task$n||any(table(result$fit$folds)!=task$n/cfg$folds))
    stop("Invalid CV folds")
  if(any(!grepl("^[a-f0-9]{64}$",unlist(result$data_hashes))))stop("Invalid input checksums")
  for(m in cfg$methods) {
    row<-r[r$method==m,]; model<-result$fit$models[[m]]
    if(!identical(row$partial_search,model$selection$partial))stop("Selection audit mismatch")
    if(row$ok) {
      if(!model$ok||length(model$beta)!=task$p||any(!is.finite(model$beta))||
         !is.finite(model$intercept)||model$curvature$status=="negative_curvature"||
         any(vapply(model$cv_curvature,function(h)h$status=="negative_curvature",logical(1))))
        stop("Selected-model numerical gate failed")
      ss<-model$selection
      selcv<-cv[cv$setting==ss$setting&cv$index==ss$index,]
      if(nrow(selcv)!=cfg$folds||!all(selcv$ok)||any(!is.finite(selcv$loss))||
         abs(mean(selcv$loss)-row$cv_loss)>1e-12*(1+abs(row$cv_loss)))
        stop("Selected CV score mismatch")
    }
  }
  invisible(TRUE)
}
perf_read_unit <- function(out,i,state) {
  path<-perf_unit_path(out,i); file<-file.path(path,"result.rds")
  if(!file.exists(file)||!file.exists(file.path(path,"result.sha256")))stop("Incomplete shard: ",path)
  if(!identical(readLines(file.path(path,"result.sha256")),perf_hash(file)))stop("Shard checksum mismatch: ",path)
  result<-readRDS(file)
  perf_validate_unit(result,state$tasks[i,,drop=FALSE],state$cfg,state$seeds[[i]],state$signature)
  result
}
perf_write_unit <- function(out,result,state) {
  i<-result$task$unit
  perf_validate_unit(result,state$tasks[i,,drop=FALSE],state$cfg,state$seeds[[i]],state$signature)
  path<-perf_unit_path(out,i)
  if(file.exists(path))stop("Refusing to overwrite a committed unit")
  tmp<-paste0(path,".partial-",Sys.getpid())
  if(!dir.create(tmp))stop("Cannot create atomic shard staging directory")
  saveRDS(result,file.path(tmp,"result.rds"),version=3,compress="gzip")
  writeLines(perf_hash(file.path(tmp,"result.rds")),file.path(tmp,"result.sha256"))
  if(!file.rename(tmp,path))stop("Cannot commit shard")
}
perf_files <- function(out) {
  paths<-list.files(out,recursive=TRUE,all.files=TRUE,no..=TRUE)
  paths[!dir.exists(file.path(out,paths)) & !paths%in%c("manifest.csv","COMPLETED") &
    !grepl("\\.partial-|\\.tmp-",paths)]
}
perf_manifest <- function(out) {
  files<-perf_files(out); paths<-file.path(out,files)
  mf<-data.frame(file=files,bytes=file.info(paths)$size,sha256=unname(vapply(paths,perf_hash,character(1))))
  perf_atomic(file.path(out,"manifest.csv"),function(p)write.csv(mf,p,row.names=FALSE))
}
perf_verify_manifest <- function(out) {
  mf<-read.csv(file.path(out,"manifest.csv"),stringsAsFactors=FALSE)
  if(anyDuplicated(mf$file)||!setequal(mf$file,perf_files(out))||
     any(grepl("(^/|(^|/)\\.\\.(/|$))",mf$file)))stop("Manifest file set mismatch")
  paths<-file.path(out,mf$file)
  if(any(!file.exists(paths))||any(file.info(paths)$size!=mf$bytes)||
     !identical(unname(vapply(paths,perf_hash,character(1))),mf$sha256))stop("Output checksum mismatch")
}
perf_validate <- function(out,root,complete=TRUE,quiet=FALSE) {
  state<-perf_verify_metadata(out,root)
  present<-vapply(state$tasks$unit,function(i)dir.exists(perf_unit_path(out,i)),logical(1))
  actual<-list.dirs(file.path(out,"shards"),recursive=FALSE,full.names=FALSE)
  committed<-actual[!grepl("\\.partial-",actual)]
  if(!setequal(committed,sprintf("unit_%05d",which(present))))stop("Unexpected shard directories")
  rows<-list()
  for(i in which(present))rows[[length(rows)+1L]]<-perf_read_unit(out,i,state)$rows
  state$completed<-which(present); state$rows<-if(length(rows))do.call(rbind,rows) else NULL
  if(complete) {
    if(!all(present)||!file.exists(file.path(out,"COMPLETED")))stop("Run incomplete")
    perf_verify_manifest(out)
    aggregate<-readRDS(file.path(out,"results.rds"))
    if(!identical(aggregate,state$rows))stop("Aggregate differs from validated shards")
    csv<-read.csv(file.path(out,"results.csv"),stringsAsFactors=FALSE)
    if(!isTRUE(all.equal(csv,aggregate,tolerance=1e-12,check.attributes=FALSE)))stop("CSV differs from shard results")
    status<-jsonlite::fromJSON(file.path(out,"progress.json"))
    if(status$phase!="completed"||status$validated_completed!=nrow(state$tasks)||
       status$running!=0||status$failed!=0||status$pending!=0)stop("Final progress inconsistent")
    if(!identical(readLines(file.path(out,"COMPLETED")),state$signature))stop("Completion signature mismatch")
    if(!quiet)cat(sprintf("COMPLETE: %d work units, %d method rows, CV/seed/source/manifest gates verified.\nFailed method fits retained: %d; partial searches: %d.\n",
      nrow(state$tasks),nrow(aggregate),sum(!aggregate$ok),sum(aggregate$partial_search)))
  }
  invisible(state)
}
perf_lock <- function(out) {
  path<-paste0(out,".lock")
  if(!dir.create(path,showWarnings=FALSE)) {
    owner<-tryCatch(jsonlite::fromJSON(file.path(path,"owner.json")),error=function(e)NULL)
    if(is.null(owner)||!identical(owner$host,unname(Sys.info()["nodename"]))||
       isTRUE(tools::pskill(as.integer(owner$pid),signal=0L)))stop("Active or unidentifiable run lock: ",path)
    stale<-paste0(path,".stale-",format(Sys.time(),"%Y%m%d%H%M%S"),"-",Sys.getpid())
    if(!file.rename(path,stale)||!dir.create(path))stop("Cannot preserve and recover stale lock")
  }
  writeLines(perf_json(list(host=unname(Sys.info()["nodename"]),pid=Sys.getpid())),
    file.path(path,"owner.json"))
  path
}
perf_progress <- function(out,state) {
  e<-new.env(parent=emptyenv()); e$completed<-state$completed; e$phase<-"preflight"
  e$unit<-NA_integer_; e$running<-0L; e$failed<-0L; e$began<-proc.time()[[3]]; e$last<--Inf
  e$durations<-numeric(); e$cost_class<-character(); e$last_completed<-if(length(e$completed))tail(e$completed,1) else NA
  e$emit<-function(detail="",force=FALSE) {
    now<-proc.time()[[3]]
    if(!force&&now-e$last<state$cfg$implementation$progress_seconds)return(invisible(NULL))
    e$last<-now; elapsed<-now-e$began
    pending<-setdiff(state$tasks$unit,e$completed)
    # Separate p/rho/profile/SNR classes: no invented ETA for unmeasured costs.
    cls<-with(state$tasks,paste(p,rho,snr,profile,sep="/"))
    means<-tapply(e$durations,e$cost_class,mean)
    eta<-if(!length(pending))0 else if(length(means)&&all(cls[pending]%in%names(means)))
      sum(means[cls[pending]]) else NA_real_
    timestamp<-format(Sys.time(),"%Y-%m-%dT%H:%M:%SZ",tz="UTC")
    snapshot<-list(run_id=basename(out),slurm_job_id=Sys.getenv("SLURM_JOB_ID",unset=NA_character_),
      phase=e$phase,work_unit="one scenario/replication including all methods and CV",
      total=nrow(state$tasks),validated_completed=length(e$completed),running=e$running,failed=e$failed,
      pending=nrow(state$tasks)-length(e$completed)-e$running-e$failed,
      percent=100*length(e$completed)/nrow(state$tasks),elapsed_seconds=elapsed,
      throughput_units_per_second=if(length(e$durations)&&elapsed>0)length(e$durations)/elapsed else NA_real_,
      estimated_remaining_seconds=eta,
      estimated_completion_utc=if(is.finite(eta))format(Sys.time()+eta,"%Y-%m-%dT%H:%M:%SZ",tz="UTC") else NA_character_,
      eta_scope="current computation phase; aggregation and validation excluded",
      slurm_remaining_seconds=NA_real_,last_heartbeat_utc=timestamp,
      current_work_unit=e$unit,most_recently_completed_unit=e$last_completed,detail=detail)
    perf_atomic_json(snapshot,file.path(out,"progress.json"))
    record<-as.data.frame(snapshot,stringsAsFactors=FALSE)
    history<-file.path(out,"progress.tsv")
    write.table(record,history,sep="\t",row.names=FALSE,col.names=!file.exists(history),
      append=file.exists(history),quote=TRUE,na="NA")
    cat(sprintf("[%s] %s | %d/%d (%.1f%%) | unit %s | ETA %s | %s\n",timestamp,e$phase,
      length(e$completed),nrow(state$tasks),snapshot$percent,e$unit,
      if(is.finite(eta))paste(round(eta),"s, computation only") else "estimating",detail))
    flush(stdout()); invisible(snapshot)
  }
  e
}
perf_preflight <- function(compile=TRUE,root) {
  for(p in c("glmnet","ncvreg","digest","jsonlite","knitr","Rcpp"))
    if(!requireNamespace(p,quietly=TRUE))stop("Required installed package missing: ",p)
  vars<-c("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","BLIS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS","NUMEXPR_NUM_THREADS")
  if(any(Sys.getenv(vars)!="1"))stop("Set all six BLAS/OpenMP thread variables to 1 before starting R.")
  if(compile)Rcpp::sourceCpp(file.path(root,"src/lambert_performance.cpp"),
    env=.GlobalEnv,cacheDir=file.path(tempdir(),"performance_rcpp"),showOutput=FALSE)
}
perf_main <- function(root) {
  args<-commandArgs(trailingOnly=TRUE)
  flags<-c("--resume"); values<-c("--mode","--stage","--out","--stop-after")
  opt<-list(mode="design",stage="pilot",resume=FALSE,stop_after=Inf)
  i<-1L
  while(i<=length(args)) {
    key<-args[i]
    if(key%in%flags) { opt$resume<-TRUE; i<-i+1L; next }
    if(!key%in%values||i==length(args)||startsWith(args[i+1],"--"))stop("Invalid CLI arguments")
    opt[[gsub("-","_",substring(key,3))]]<-args[i+1]; i<-i+2L
  }
  if(is.null(opt$out)||!opt$mode%in%c("design","check","smoke","pilot","confirmation","validate"))
    stop("Use --mode design|check|smoke|pilot|confirmation|validate --out DIRECTORY")
  out<-if(startsWith(opt$out,"/"))opt$out else file.path(root,opt$out)
  stop_after<-as.numeric(opt$stop_after)
  if(is.na(stop_after)||stop_after<1||(!is.infinite(stop_after)&&opt$mode!="smoke"))
    stop("--stop-after is for smoke interruption testing only")
  if(opt$mode!="validate"&&file.exists(out)&&!opt$resume)
    stop("Output exists; refusing overwrite. Use a new version or --resume.")
  perf_preflight(compile=!opt$mode%in%c("design","validate"),root=root)
  if(opt$mode=="validate") { perf_validate(out,root); return(invisible(NULL)) }
  stage<-if(opt$mode%in%c("design","check"))opt$stage else opt$mode
  cfg<-perf_config(root,stage)
  if(opt$resume&&opt$mode%in%c("design","check"))stop("Design/check outputs are not resumable")
  dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE)
  lock<-perf_lock(out)
  on.exit(unlink(lock,recursive=TRUE),add=TRUE) # Only this driver's ephemeral lock.
  if(opt$resume) {
    if(!dir.exists(out))stop("Cannot resume missing run")
    state<-perf_validate(out,root,complete=FALSE)
    if(!identical(perf_json(state$cfg),perf_json(cfg)))stop("Resume stage mismatch")
    if(file.exists(file.path(out,"COMPLETED"))) { perf_validate(out,root); return(invisible(NULL)) }
  } else {
    if(!dir.create(out)||!dir.create(file.path(out,"shards")))stop("Cannot create fresh run")
    tasks<-perf_tasks(cfg); seeds<-perf_seeds(cfg,tasks)
    perf_metadata(out,root,cfg,tasks,seeds)
    state<-perf_verify_metadata(out,root); state$completed<-integer()
  }
  if(opt$mode%in%c("design","check")) {
    if(opt$mode=="check") {
      source(file.path(root,"tests/performance_v1_checks.R"))
      checks<-perf_checks(root)
      write.csv(checks,file.path(out,"checks.csv"),row.names=FALSE)
    }
    perf_manifest(out); perf_verify_manifest(out)
    cat(sprintf("DESIGN VERIFIED: %d %s work units; no production fits.\n",nrow(state$tasks),stage))
    return(invisible(NULL))
  }
  progress<-perf_progress(out,state); progress$emit(force=TRUE); new_units<-0L
  tryCatch({
    queue<-state$tasks$unit[order(state$tasks$repetition,state$tasks$scenario)]
    for(i in setdiff(queue,state$completed)) {
      progress$unit<-i; progress$phase<-"computation"; progress$running<-1L
      progress$emit("starting work unit",force=TRUE)
      result<-perf_run_task(state$tasks[i,,drop=FALSE],cfg,state$seeds[[i]],progress$emit)
      result$signature<-state$signature
      perf_write_unit(out,result,state); perf_read_unit(out,i,state)
      progress$completed<-c(progress$completed,i); progress$running<-0L; progress$last_completed<-i
      progress$durations<-c(progress$durations,result$elapsed)
      t<-state$tasks[i,]; progress$cost_class<-c(progress$cost_class,paste(t$p,t$rho,t$snr,t$profile,sep="/"))
      progress$emit("work unit validated",force=TRUE); new_units<-new_units+1L
      if(new_units>=stop_after&&length(progress$completed)<nrow(state$tasks)) {
        progress$phase<-"checkpointed"; progress$emit("intentional smoke stop; resume supported",force=TRUE)
        return(invisible(NULL))
      }
    }
    progress$phase<-"aggregation"; progress$unit<-NA_integer_; progress$emit("validating shards",force=TRUE)
    checked<-perf_validate(out,root,complete=FALSE)
    perf_atomic(file.path(out,"results.rds"),function(p)saveRDS(checked$rows,p,version=3))
    perf_atomic(file.path(out,"results.csv"),function(p)write.csv(checked$rows,p,row.names=FALSE))
    progress$phase<-"completed"; progress$emit("all result/status rows validated",force=TRUE)
    perf_manifest(out); perf_verify_manifest(out)
    writeLines(state$signature,file.path(out,"COMPLETED"))
    perf_validate(out,root)
  },error=function(e) {
    progress$phase<-"failed"; progress$running<-0L
    progress$failed<-as.integer(length(progress$completed)<nrow(state$tasks))
    progress$emit(conditionMessage(e),force=TRUE)
    filename<-paste0("failure_",format(Sys.time(),"%Y%m%d%H%M%S"),"_",Sys.getpid(),".json")
    perf_atomic_json(list(unit=progress$unit,error=conditionMessage(e)),file.path(out,filename))
    stop(e)
  })
}
