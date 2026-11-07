# One writer over the shared checkpoint helper; two adjacent packs may compute.
sch_check_unit <- function(z,i,e) {
  stopifnot(identical(z$task,e$tasks[i,,drop=FALSE]),identical(z$split,sch_split(z$task,e)),
    is.finite(z$elapsed),z$elapsed>=0)
  prep<-sch_prep(z$task,e);setting<-e$settings[z$task$setting,,drop=FALSE];L<-e$cfg$lambda_points
  stopifnot(sch_close(z$preprocessing,prep[c("mx","sx","my","lambda0","spectral")]),
    identical(dim(z$path$beta),c(ncol(prep$Z),L)),
    sch_number(z$path$diagnostics$lambda,e$perf_lambdas(prep,setting,e$cfg)),
    identical(z$path$diagnostics$index,seq_len(L)),all(diff(z$path$diagnostics$lambda)<0))
  for(nm in c("baseline","path")) {
    p<-z[[nm]]
    for(k in which(p$diagnostics$ok)) {
      d<-e$perf_diagnose(p$beta[,k],prep,p$diagnostics$lambda[k],setting,e$cfg)
      stopifnot(d$ok,sch_number(d$kkt,p$diagnostics$kkt[k]),sch_number(d$objective,p$diagnostics$objective[k]))
    }
  }
  keep<-which(z$baseline$diagnostics$ok);targets<-which(!z$baseline$diagnostics$ok)
  stopifnot(identical(z$path$beta[,keep,drop=FALSE],z$baseline$beta[,keep,drop=FALSE]),
    identical(z$path$diagnostics[keep,],z$baseline$diagnostics[keep,]),
    identical(z$path$attempts,z$baseline$attempts))
  warm<-numeric(ncol(prep$Z))
  for(k in seq_len(L)) {
    if(k %in% keep) {warm<-z$baseline$beta[,k];next}
    if(setting$family=="enet" && setting$shape==0)next
    a<-z$rescue[[as.character(k)]];stopifnot(identical(a$initial,list(warm=warm+0,zero=numeric(length(warm)))))
    for(nm in c("warm","zero")) {
      f<-a$fits[[nm]];dd<-e$perf_diagnose(f$beta,prep,z$path$diagnostics$lambda[k],setting,e$cfg)
      tr<-f$trace
      stopifnot(sch_number(dd$kkt,f$kkt),sch_number(dd$objective,f$objective),
        all(diff(tr$objective)<=e$cfg$implementation$objective_increase_tol*(1+abs(head(tr$objective,-1)))),
        f$sweeps<=e$rescue$total_sweeps,
        f$newton_steps<=e$rescue$newton_steps*(e$rescue$total_sweeps %/% e$rescue$chunk_sweeps),
        f$ok==identical(f$status,"verified_stationary"))
      if(f$ok)stopifnot(dd$ok,rfh_curvature_ok(rfh_curvature(f$beta,prep,z$path$diagnostics$lambda[k],setting,e)$status))
    }
    good<-which(vapply(a$fits,function(f)f$ok,FALSE));chosen<-NA_integer_
    if(length(good)) {
      q<-vapply(a$fits[good],function(f)f$objective,1.0);best<-min(q)
      chosen<-good[which(q<=best+e$cfg$implementation$tie_tol*(1+abs(best)))[1]]
    }
    stopifnot(identical(chosen,a$chosen),z$path$diagnostics$ok[k]==!is.na(chosen))
    if(!is.na(chosen)) {stopifnot(identical(z$path$beta[,k],a$fits[[chosen]]$beta));warm<-a$fits[[chosen]]$beta}
  }
  curvature<-rep("not_assessed",L)
  for(k in which(z$path$diagnostics$ok))curvature[k]<-rfh_curvature(z$path$beta[,k],prep,z$path$diagnostics$lambda[k],setting,e)$status
  stopifnot(identical(curvature,z$curvature),identical(z$valid,z$path$diagnostics$ok & rfh_curvature_ok(curvature)))
  if(z$task$role==0L)stopifnot(sch_close(z$candidates,sch_candidates(z$path,prep,setting,e))) else {
    pred<-matrix(NA_real_,length(z$split$test),L)
    for(k in which(z$valid))pred[,k]<-e$perf_predict(z$path$beta[,k],prep,
      e$data$X[z$split$test,z$split$columns,drop=FALSE])
    stopifnot(sch_close(z$predictions,pred),sch_close(z$loss,unname(colMeans((e$data$y[z$split$test]-pred)^2))))
  }
  invisible(TRUE)
}
sch_metadata <- function(e)list(schema="scheetz_cv_ic_run_v1",scientific_signature=e$scientific_signature,
  execution_signature=e$execution_signature,mapping=e$mapping,stage=e$stage,cfg=e$cfg,rescue=e$rescue,
  tasks=e$tasks,settings=e$settings,view=e$view,folds=e$folds,sources=e$sources,inputs=e$inputs,software=e$software)
sch_unchanged <- function(e) {
  stopifnot(identical(e$sources,rfh_files(e$root,sch_sources())),identical(e$inputs,rfh_files(e$root,e$inputs$file)))
  invisible(TRUE)
}
sch_scan <- function(out,e,recover=FALSE,portable=FALSE,pulse=function(...)NULL) {
  m<-cpk_metadata(out,e$scientific_signature,e$execution_signature);expected<-sch_metadata(e)
  same<-identical(m$software,e$software);expected$software<-m$software
  stopifnot(identical(m,expected));if(!portable && !same)stop("Runtime differs; use portable verification, not resume")
  verify<-function(z,i) {sch_check_unit(z,i,e);pulse(i)}
  if(!portable || same)return(cpk_scan(out,e$scientific_signature,e$execution_signature,verify,recover=recover))
  stopifnot(!recover)
  mf<-cpk_manifest(out);packs<-mf[mf$kind=="pack",];packs<-packs[order(packs$batch),]
  stopifnot(identical(as.integer(packs$batch),seq_len(nrow(packs))),
    setequal(list.files(file.path(out,"checkpoints")),basename(packs$file)))
  completed<-integer()
  for(k in seq_len(nrow(packs))) {
    b<-packs$batch[k];p<-readRDS(file.path(out,packs$file[k]));ids<-m$mapping$unit[m$mapping$batch==b]
    stopifnot(p$schema=="checkpoint_pack_v1",identical(p$units,ids),identical(names(p$records),as.character(ids)),
      p$scientific_signature==m$scientific_signature,p$execution_signature==m$execution_signature,
      identical(p$receipts$unit,ids),all(p$receipts$accepted),length(p$records)==length(ids))
    for(i in ids)verify(p$records[[as.character(i)]],i)
    completed<-c(completed,ids)
  }
  list(metadata=m,manifest=mf,completed=completed,batches=as.integer(packs$batch))
}
sch_ids <- function(method,criterion,e) {
  s<-e$settings
  ids<-switch(method,lambert_c1=1L,lasso=2L,ridge=3L,elastic_net=which(s$family=="enet"),mcp=7L,scad=8L)
  if(method=="elastic_net" && !criterion %in% c("CV","BIC"))ids<-ids[s$shape[ids]>0]
  ids[order(-s$shape[ids])]
}
sch_select <- function(paths,method,criterion,e) {
  ids<-sch_ids(method,criterion,e);rows<-list()
  for(j in ids) {
    outer<-paths[[paste(0,j,sep="_")]]
    if(criterion=="CV") {
      inner<-lapply(1:5,function(f)paths[[paste(f,j,sep="_")]])
      loss<-rowMeans(vapply(inner,function(z)z$loss,numeric(e$cfg$lambda_points)))
      eligible<-Reduce(`&`,lapply(inner,function(z)z$valid)) & is.finite(loss)
      partial<-any(!eligible)
    } else {
      cc<-outer$candidates;loss<-cc[[criterion]]
      eligible<-if(criterion=="BIC")cc$eligible_BIC else cc$eligible_sparse
      partial<-any(!outer$valid)
    }
    rows[[length(rows)+1L]]<-data.frame(setting=j,index=seq_len(e$cfg$lambda_points),score=loss,
      eligible=eligible,partial_search=partial)
  }
  all<-sch_bind(rows);k<-e$lic_choose(all$score,all$eligible,e$cfg$implementation$tie_tol)
  list(setting=if(is.na(k))NA_integer_ else all$setting[k],index=if(is.na(k))NA_integer_ else all$index[k],
    score=if(is.na(k))NA_real_ else all$score[k],partial_search=any(all$partial_search))
}
sch_collect <- function(out,e,portable=FALSE,pulse=function(...)NULL) {
  state<-sch_scan(out,e,portable=portable,pulse=pulse);stopifnot(identical(state$completed,e$tasks$unit))
  records<-unlist(lapply(state$batches,function(b)readRDS(file.path(out,cpk_pack_file(b)))$records),recursive=FALSE)
  rows<-predictions<-diagnostics<-candidate_rows<-list()
  for(z in records) {
    diagnostics[[length(diagnostics)+1L]]<-cbind(z$task,elapsed=z$elapsed,baseline_verified=sum(z$baseline$diagnostics$ok),
      rescued=sum(z$path$diagnostics$ok)-sum(z$baseline$diagnostics$ok),failed=sum(!z$path$diagnostics$ok),
      curvature_rejected=sum(z$path$diagnostics$ok & !z$valid),max_kkt=if(any(z$valid))max(z$path$diagnostics$kkt[z$valid]) else NA_real_)
    candidate_rows[[length(candidate_rows)+1L]]<-cbind(z$task[rep(1,e$cfg$lambda_points),],z$path$diagnostics,
      valid=z$valid,curvature=z$curvature,loss=if(z$task$role>0L)z$loss else NA_real_)
  }
  for(r in e$repeats)for(f in e$active_outer) {
    ids<-which(e$tasks$repeat_id==r & e$tasks$outer_fold==f)
    paths<-setNames(records[as.character(ids)],paste(e$tasks$role[ids],e$tasks$setting[ids],sep="_"))
    for(method in e$cfg$methods)for(cr in if(method=="ridge")c("CV","BIC") else c("CV","BIC","EBIC_0.5","EBIC_1","HBIC")) {
      ss<-sch_select(paths,method,cr,e);ok<-!is.na(ss$index);n_test<-length(paths[[1]]$split$test)
      base<-paths[[1]];split<-base$split;lambda<-q<-d<-kkt<-NA_real_;alpha<-NA_real_;curvature<-"not_available"
      pred<-rep(NA_real_,n_test);endpoint<-cap_boundary<-NA
      if(ok) {
        z<-paths[[paste(0,ss$setting,sep="_")]];ok<-z$valid[ss$index];split<-z$split
        j<-ss$setting;k<-ss$index;lambda<-z$path$diagnostics$lambda[k];alpha<-e$settings$shape[j]
        cc<-z$candidates;q<-cc$q[k];d<-cc$d[k];kkt<-z$path$diagnostics$kkt[k];curvature<-z$curvature[k]
        endpoint<-k %in% c(1L,e$cfg$lambda_points)
        cap_boundary<-if(cr=="CV")FALSE else abs((if(cr=="BIC")q else d)-floor(length(split$train)/log(length(split$train))))<1e-10
        if(ok) {prep<-sch_prep(z$task,e);pred<-e$perf_predict(z$path$beta[,k],prep,e$data$X[split$test,split$columns,drop=FALSE])}
      }
      y<-e$data$y[split$test]
      rows[[length(rows)+1L]]<-data.frame(repeat_id=r,outer_fold=f,method=method,criterion=cr,ok=ok,
        partial_search=ss$partial_search,index=ss$index,setting=ss$setting,shape=alpha,lambda=lambda,score=ss$score,
        q=q,selected_size=d,kkt=kkt,curvature=curvature,path_endpoint=endpoint,cap_boundary=cap_boundary,
        n_train=length(split$train),n_test=length(y),rmse=if(ok)sqrt(mean((y-pred)^2)) else NA_real_,mae=if(ok)mean(abs(y-pred)) else NA_real_)
      predictions[[length(predictions)+1L]]<-data.frame(repeat_id=r,outer_fold=f,method=method,criterion=cr,
        sample_index=split$test,observed=y,predicted=pred,ok=ok)
    }
  }
  rows<-sch_bind(rows);predictions<-sch_bind(predictions);repeats<-summary<-list()
  arms<-unique(rows[c("method","criterion")])
  for(a in seq_len(nrow(arms))) {
    rr<-rows[rows$method==arms$method[a] & rows$criterion==arms$criterion[a],]
    rp<-sch_bind(lapply(e$repeats,function(r) {
      z<-rr[rr$repeat_id==r,];ok<-nrow(z)==5L && all(z$ok)
      data.frame(method=arms$method[a],criterion=arms$criterion[a],repeat_id=r,complete=ok,
        rmse=if(ok)mean(z$rmse) else NA_real_,mae=if(ok)mean(z$mae) else NA_real_,
        selected_size=if(ok)mean(z$selected_size) else NA_real_)
    }))
    su<-data.frame(method=arms$method[a],criterion=arms$criterion[a],repeats_expected=length(e$repeats),
      repeats_complete=sum(rp$complete),failed_folds=sum(!rr$ok),partial_search_folds=sum(rr$partial_search),
      endpoint_folds=sum(rr$path_endpoint,na.rm=TRUE),cap_boundary_folds=sum(rr$cap_boundary,na.rm=TRUE))
    for(nm in c("rmse","mae","selected_size")) {
      su[[paste0("mean_",nm)]]<-if(all(rp$complete))mean(rp[[nm]]) else NA_real_
      su[[paste0("sd_",nm)]]<-if(all(rp$complete))sd(rp[[nm]]) else NA_real_
    }
    repeats[[a]]<-rp;summary[[a]]<-su
  }
  stopifnot(nrow(rows)==27L*length(e$repeats)*length(e$active_outer),
    nrow(predictions)==27L*sum(vapply(e$repeats,function(r)sum(e$folds[[r]]$outer %in% e$active_outer),1L)),
    !anyDuplicated(rows[c("repeat_id","outer_fold","method","criterion")]),
    !anyDuplicated(predictions[c("repeat_id","method","criterion","sample_index")]))
  list(results=rows,predictions=predictions,path_summary=sch_bind(diagnostics),candidate_diagnostics=sch_bind(candidate_rows),
    repeat_metrics=sch_bind(repeats),summary=sch_bind(summary))
}
sch_check_exports <- function(out,z) {
  mf<-cpk_manifest(out);stopifnot(setequal(mf$file[mf$kind=="aggregate"],paste0(names(z),".csv")))
  for(nm in names(z))stopifnot(isTRUE(all.equal(rfh_read_csv(file.path(out,paste0(nm,".csv")),z[[nm]]),
    z[[nm]],tolerance=1e-11,check.attributes=FALSE)))
  invisible(TRUE)
}
sch_verify <- function(out,e,portable=FALSE) {
  sch_unchanged(e);stopifnot(readLines(file.path(out,"COMPLETED"))==cpk_hash(file.path(out,"manifest.tsv")))
  last<-as.numeric(Sys.time());z<-sch_collect(out,e,portable,pulse=function(i) {
    if(as.numeric(Sys.time())-last>=30) {cat("validation | path",i,"/",nrow(e$tasks),"; no fits\n");flush(stdout());last<<-as.numeric(Sys.time())}
  });sch_check_exports(out,z)
  cat(sprintf("SCHEETZ_CV_IC_COMPLETE: %d paths, %d result rows, %d predictions; failed selected fits=%d; partial searches=%d; model fits in verification=0.\n",
    nrow(e$tasks),nrow(z$results),nrow(z$predictions),sum(!z$results$ok),sum(z$results$partial_search)))
  print(z$summary,row.names=FALSE);invisible(z)
}
sch_emit <- function(p,phase="computation",detail="",force=FALSE) {
  now<-as.numeric(Sys.time());if(!force && now-p$last<30)return(invisible(NULL));p$last<-now
  e<-p$e;N<-nrow(e$tasks);done<-union(p$committed,p$accepted);todo<-setdiff(e$tasks$unit,done);eta<-NA_real_
  cls<-paste(e$tasks$setting,ifelse(e$tasks$role==0,"outer","inner"),sep="_")
  if(phase=="completed")eta<-0 else if(phase=="computation" && length(todo)) {
    cost<-tapply(p$durations,cls,function(x)if(sum(is.finite(x))>=2L)mean(x[is.finite(x)]) else NA_real_)
    cc<-cost[cls[todo]];if(all(is.finite(cc)))eta<-max(max(cc),sum(cc)/e$workers)
  }
  end<-suppressWarnings(as.numeric(Sys.getenv("SLURM_JOB_END_TIME")))
  z<-list(schema="scheetz_cv_ic_progress_v1",run_id=basename(p$out),slurm_job_id=if(nzchar(Sys.getenv("SLURM_JOB_ID")))Sys.getenv("SLURM_JOB_ID") else NA_character_,
    scientific_signature=e$scientific_signature,execution_signature=e$execution_signature,phase=phase,
    work_unit="one repeat/outer fold/inner role/penalty-setting path",workers=e$workers,total=N,
    validated_completed=length(done),durably_committed=length(p$committed),validated_before_invocation=p$before,
    completed_this_invocation=length(p$accepted),running=length(p$running),failed=length(p$failed),
    pending=N-length(done)-length(p$running)-length(p$failed),percent_completed=100*length(p$committed)/N,
    percent_validated=100*length(done)/N,elapsed_seconds=p$previous+now-p$started,invocation_elapsed_seconds=now-p$started,
    throughput_units_per_second=length(p$accepted)/max(.001,now-p$started),estimated_remaining_seconds=eta,
    estimated_completion_utc=if(is.finite(eta))format(as.POSIXct(now+eta,origin="1970-01-01",tz="UTC"),"%Y-%m-%dT%H:%M:%SZ") else NA_character_,
    eta_status=if(is.finite(eta))"estimated" else "estimating_or_unknown",eta_scope="computation and per-path checks; pack writes, final replay and archive excluded",
    slurm_remaining_seconds=if(is.finite(end))max(0,end-now) else NA_real_,last_heartbeat_utc=format(Sys.time(),"%Y-%m-%dT%H:%M:%SZ",tz="UTC"),
    peak_running=p$peak,current_work_units=paste(p$running,collapse=","),most_recently_completed_unit=p$latest,detail=detail)
  cpk_atomic(file.path(p$out,"progress.json"),function(f)jsonlite::write_json(z,f,auto_unbox=TRUE,pretty=TRUE,na="null",digits=17),replace=TRUE)
  f<-file.path(p$out,"progress.tsv");write.table(as.data.frame(z),f,sep="\t",row.names=FALSE,quote=TRUE,col.names=!file.exists(f),append=file.exists(f))
  cat(sprintf("[%s] %s | validated %d/%d saved %d/%d running %d failed %d | ETA %s | %s\n",z$last_heartbeat_utc,
    phase,length(done),N,length(p$committed),N,z$running,z$failed,if(is.finite(eta))paste(round(eta),"s, computation only") else "unknown",detail));flush(stdout())
}
sch_run <- function(out,e,action="start",stop_batches=Inf,max_seconds=Inf,require_release=TRUE) {
  stopifnot(action %in% c("start","resume"),stop_batches>=1,max_seconds>0)
  if(e$stage=="pilot")sch_preflight(e,arf=TRUE,storage=TRUE)
  if(require_release)sch_gate(e)
  if(action=="start" && file.exists(out))stop("Refusing overwrite: ",out)
  if(action=="resume" && !dir.exists(out))stop("Missing run")
  dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE);lock<-paste0(out,".run_lock")
  if(!dir.create(lock,showWarnings=FALSE))stop("Existing run lock; inspect owner before resume")
  saveRDS(list(pid=Sys.getpid(),job=Sys.getenv("SLURM_JOB_ID"),host=Sys.info()[["nodename"]],started=Sys.time()),file.path(lock,"owner.rds"))
  on.exit(unlink(lock,recursive=TRUE),add=TRUE)
  if(action=="start")cpk_create(out,sch_metadata(e))
  state<-sch_scan(out,e,recover=action=="resume");if(file.exists(file.path(out,"failure_evidence.rds")))stop("Failure evidence requires inspection")
  if(file.exists(file.path(out,"COMPLETED")))return(invisible(sch_verify(out,e)))
  p<-rfh_progress(out,e,state$completed);p$durations<-rep(NA_real_,nrow(e$tasks));p$peak<-0L;buffer<-receipts<-list()
  remaining<-setdiff(e$tasks$unit,state$completed);next_position<-1L
  next_batch<-if(length(state$batches))max(state$batches)+1L else 1L;initial_batch<-next_batch;paused<-FALSE
  commit_ready<-function() {
    while(next_batch<=max(e$mapping$batch)) {
      ids<-e$mapping$unit[e$mapping$batch==next_batch]
      if(!all(as.character(ids) %in% names(buffer)))break
      sch_emit(p,"checkpoint_commit",paste("saving pack",next_batch),TRUE)
      cpk_commit(out,state$metadata,next_batch,buffer[as.character(ids)],sch_bind(receipts[as.character(ids)]))
      p$committed<-union(p$committed,ids);buffer[as.character(ids)]<<-NULL;receipts[as.character(ids)]<<-NULL
      next_batch<<-next_batch+1L;sch_emit(p,detail="validated pack committed",force=TRUE)
    }
  }
  sch_emit(p,detail="frozen CV/IC pilot started",force=TRUE)
  tryCatch({
    accept<-function(z,i) {
      began<-proc.time()[[3]];sch_check_unit(z,i,e);buffer[[as.character(i)]]<<-z;receipts[[as.character(i)]]<<-cpk_receipt(z,i)
      p$accepted<-union(p$accepted,i);p$running<-setdiff(p$running,i);p$latest<-i
      p$durations[i]<-z$elapsed+proc.time()[[3]]-began;commit_ready()
    }
    pulse<-function(event,i,running) {
      if(event=="start") {p$running<-union(p$running,i);next_position<<-next_position+1L}
      if(event %in% c("accepted","failed"))p$running<-setdiff(p$running,i)
      if(event=="failed")p$failed<-union(p$failed,i)
      peak_event<-event=="start" && length(p$running)>p$peak && length(p$running)==e$workers
      p$peak<-max(p$peak,length(p$running))
      sch_emit(p,detail=paste(event,if(is.finite(i))i else ""),force=peak_event||event=="failed")
    }
    launch<-function() {
      if(as.numeric(Sys.time())-p$started>=max_seconds || next_batch-initial_batch>=stop_batches)return(FALSE)
      if(next_position>length(remaining))return(FALSE)
      b<-e$mapping$batch[match(remaining[next_position],e$mapping$unit)]
      b<=next_batch+1L && b<initial_batch+stop_batches
    }
    result<-cpk_dispatch(remaining,e$workers,function(i)sch_unit(i,e),accept,pulse,launch)
    if(length(result$failed)) {
      cpk_atomic(file.path(out,"failure_evidence.rds"),function(f)saveRDS(result$failed,f,version=3))
      stop("Worker/validation failure; evidence preserved")
    }
    commit_ready();sch_unchanged(e)
    if(length(p$committed)<nrow(e$tasks)) {
      p$running<-integer();sch_emit(p,"checkpointed","dispatch stopped; uncommitted work requires recomputation",TRUE);return(invisible(NULL))
    }
    sch_emit(p,"aggregation_and_validation","replaying all paths and selections; no new fits",TRUE)
    z<-sch_collect(out,e,pulse=function(i)sch_emit(p,"aggregation_and_validation",paste("replayed path",i)))
    for(nm in names(z)) {
      file<-paste0(nm,".csv");mf<-cpk_manifest(out);f<-file.path(out,file)
      if(!file %in% mf$file) {
        cpk_atomic(f,function(p)write.csv(z[[nm]],p,row.names=FALSE,na="NA"));cpk_register(out,cpk_row(out,file,"aggregate",records=nrow(z[[nm]])))
      } else stopifnot(isTRUE(all.equal(rfh_read_csv(f,z[[nm]]),z[[nm]],tolerance=1e-11,check.attributes=FALSE)))
    }
    sch_unchanged(e);sch_check_exports(out,z)
    cpk_atomic(file.path(out,"COMPLETED"),function(f)writeLines(cpk_hash(file.path(out,"manifest.tsv")),f))
    sch_emit(p,"completed","paths, choices, held-out predictions and manifests verified",TRUE)
    cat("SCHEETZ_CV_IC_COMPLETE:",nrow(e$tasks),"paths;",nrow(z$results),"rows;",nrow(z$predictions),"predictions\n")
    print(z$summary,row.names=FALSE);invisible(z)
  },interrupt=function(err) {p$running<-integer();sch_emit(p,"checkpointed","interrupted; verified packs preserved",TRUE);stop("Interrupted")},
  error=function(err) {p$running<-integer();sch_emit(p,"failed",conditionMessage(err),TRUE);stop(err)})
}
