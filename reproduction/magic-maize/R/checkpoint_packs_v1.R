# Reusable single-coordinator, bounded ten-pack persistence. No model code.
cpk_hash <- function(path) digest::digest(file=path,algo="sha256")
cpk_object_hash <- function(x) digest::digest(x,algo="sha256",serializeVersion=3)
cpk_atomic <- function(path,writer,replace=FALSE) {
  if(!replace && file.exists(path))stop("Refusing overwrite: ",path)
  tmp <- paste0(path,".partial-",Sys.getpid())
  if(file.exists(tmp))stop("Unfinished temporary file: ",tmp)
  writer(tmp)
  if(!file.rename(tmp,path))stop("Atomic rename failed: ",path)
  invisible(path)
}
cpk_mapping <- function(units,batches=10L) {
  stopifnot(length(units)>0L,!anyDuplicated(units),all(units>0),batches==10L)
  n <- length(units); B <- min(batches,n)
  data.frame(position=seq_len(n),unit=as.integer(units),batch=as.integer(ceiling(seq_len(n)*B/n)))
}
cpk_row <- function(out,file,kind,batch=0L,records=0L) data.frame(file=file,kind=kind,
  batch=as.integer(batch),records=as.integer(records),bytes=as.numeric(file.info(file.path(out,file))$size),
  sha256=cpk_hash(file.path(out,file)),validated=TRUE,stringsAsFactors=FALSE)
cpk_manifest <- function(out) {
  x <- read.delim(file.path(out,"manifest.tsv"),stringsAsFactors=FALSE,check.names=FALSE)
  stopifnot(identical(names(x),c("file","kind","batch","records","bytes","sha256","validated")),
    !anyDuplicated(x$file),all(x$validated),!any(grepl("(^/|(^|/)\\.\\.(/|$))",x$file)))
  for(i in seq_len(nrow(x))) {
    p <- file.path(out,x$file[i])
    if(!file.exists(p) || file.info(p)$size!=x$bytes[i] || cpk_hash(p)!=x$sha256[i])
      stop("Checkpoint checksum/size mismatch: ",x$file[i])
  }
  x
}
cpk_register <- function(out,row) {
  mf <- cpk_manifest(out)
  stopifnot(!row$file %in% mf$file)
  cpk_atomic(file.path(out,"manifest.tsv"),function(p)
    write.table(rbind(mf,row),p,sep="\t",row.names=FALSE,quote=TRUE),replace=TRUE)
}
cpk_create <- function(out,metadata) {
  stopifnot(!file.exists(out),is.list(metadata),is.data.frame(metadata$mapping),
    !anyDuplicated(metadata$mapping$unit),nzchar(metadata$scientific_signature),nzchar(metadata$execution_signature))
  stopifnot(dir.create(out,recursive=TRUE),dir.create(file.path(out,"checkpoints")))
  cpk_atomic(file.path(out,"metadata.rds"),function(p)saveRDS(metadata,p,version=3))
  write.table(cpk_row(out,"metadata.rds","metadata"),file.path(out,"manifest.tsv"),
    sep="\t",row.names=FALSE,quote=TRUE)
  invisible(metadata)
}
cpk_metadata <- function(out,scientific=NULL,execution=NULL) {
  mf <- cpk_manifest(out); stopifnot(sum(mf$kind=="metadata")==1L,mf$file[mf$kind=="metadata"]=="metadata.rds")
  m <- readRDS(file.path(out,"metadata.rds"))
  if(!is.null(scientific) && !identical(m$scientific_signature,scientific))stop("Scientific signature mismatch")
  if(!is.null(execution) && !identical(m$execution_signature,execution))stop("Execution signature mismatch")
  m
}
cpk_pack_file <- function(batch) file.path("checkpoints",sprintf("part_%02d.rds",batch))
cpk_receipt <- function(record,unit) data.frame(unit=as.integer(unit),record_sha256=cpk_object_hash(record),accepted=TRUE)
cpk_check_pack <- function(pack,metadata,batch,validate=function(...)invisible(TRUE),deep=TRUE) {
  ids <- metadata$mapping$unit[metadata$mapping$batch==batch]
  stopifnot(identical(pack$schema,"checkpoint_pack_v1"),identical(pack$batch,as.integer(batch)),
    identical(pack$scientific_signature,metadata$scientific_signature),
    identical(pack$execution_signature,metadata$execution_signature),identical(pack$units,ids),
    identical(names(pack$records),as.character(ids)),nrow(pack$receipts)==length(ids),
    identical(pack$receipts$unit,ids),all(pack$receipts$accepted),
    identical(pack$receipts$record_sha256,unname(vapply(pack$records,cpk_object_hash,""))))
  if(deep)for(i in seq_along(ids))validate(pack$records[[i]],ids[i])
  invisible(TRUE)
}
cpk_commit <- function(out,metadata,batch,records,receipts,validate=function(...)invisible(TRUE)) {
  ids <- metadata$mapping$unit[metadata$mapping$batch==batch]
  stopifnot(length(ids)>0L,setequal(names(records),as.character(ids)))
  records <- records[as.character(ids)]; receipts <- receipts[match(ids,receipts$unit),,drop=FALSE]
  rownames(receipts) <- NULL
  pack <- list(schema="checkpoint_pack_v1",batch=as.integer(batch),units=ids,
    scientific_signature=metadata$scientific_signature,execution_signature=metadata$execution_signature,
    records=records,receipts=receipts)
  cpk_check_pack(pack,metadata,batch,validate,deep=FALSE)
  file <- cpk_pack_file(batch); path <- file.path(out,file)
  cpk_atomic(path,function(p)saveRDS(pack,p,version=3,compress="gzip"))
  # Readback proves recoverability; scientific acceptance is bound by receipts.
  cpk_check_pack(readRDS(path),metadata,batch,validate,deep=FALSE)
  cpk_register(out,cpk_row(out,file,"pack",batch,length(ids)))
  invisible(cpk_hash(path))
}
cpk_scan <- function(out,scientific=NULL,execution=NULL,validate=function(...)invisible(TRUE),recover=FALSE,deep=TRUE) {
  if(any(grepl("[.]partial-",list.files(out,all.files=TRUE,no..=TRUE))))
    stop("Incomplete shared metadata write; preserve and inspect")
  m <- cpk_metadata(out,scientific,execution); mf <- cpk_manifest(out)
  names_on_disk <- list.files(file.path(out,"checkpoints"),all.files=TRUE,no..=TRUE)
  if(any(!grepl("^part_[0-9]{2}[.]rds$",names_on_disk)))stop("Incomplete/unexpected checkpoint file; preserve and inspect")
  known <- basename(mf$file[mf$kind=="pack"]); orphan <- setdiff(names_on_disk,known)
  if(length(orphan) && !recover)stop("Unregistered checkpoint pack; explicit recovery validation required")
  for(file in orphan) {
    batch <- as.integer(sub("^part_([0-9]+)[.]rds$","\\1",file))
    stopifnot(batch %in% m$mapping$batch)
    pack <- readRDS(file.path(out,"checkpoints",file))
    cpk_check_pack(pack,m,batch,validate,deep=TRUE)
    cpk_register(out,cpk_row(out,file.path("checkpoints",file),"pack",batch,length(pack$units)))
  }
  mf <- cpk_manifest(out); packs <- mf[mf$kind=="pack",]; packs <- packs[order(packs$batch),]
  if(nrow(packs))stopifnot(!anyDuplicated(packs$batch),identical(as.integer(packs$batch),seq_len(nrow(packs))))
  completed <- integer()
  for(i in seq_len(nrow(packs))) {
    p <- readRDS(file.path(out,packs$file[i])); cpk_check_pack(p,m,packs$batch[i],validate,deep)
    stopifnot(packs$records[i]==length(p$units)); completed <- c(completed,p$units)
  }
  list(metadata=m,manifest=mf,completed=completed,batches=as.integer(packs$batch))
}
cpk_dispatch <- function(units,workers,run,accept,pulse=function(...)NULL,can_launch=function()TRUE) {
  jobs <- list(); queue <- units; launched <- completed <- integer(); failed <- NULL
  on.exit({for(job in jobs)try(tools::pskill(job$pid),silent=TRUE)
    if(length(jobs))try(parallel::mccollect(jobs,wait=FALSE),silent=TRUE)},add=TRUE)
  while(length(queue)||length(jobs)) {
    while(is.null(failed) && length(queue) && length(jobs)<workers && can_launch()) {
      i <- queue[1]; queue <- queue[-1]; launched <- c(launched,i)
      job <- parallel::mcparallel(local({ii<-i;tryCatch(list(ok=TRUE,result=run(ii)),
        error=function(err)list(ok=FALSE,message=conditionMessage(err)))}),mc.set.seed=FALSE,silent=TRUE)
      attr(job,"unit") <- i; jobs[[as.character(job$pid)]] <- job
      pulse("start",i,length(jobs))
    }
    if(!length(jobs))break
    done <- parallel::mccollect(jobs,wait=FALSE)
    for(pid in names(done)) {
      i <- attr(jobs[[pid]],"unit"); jobs[[pid]] <- NULL; item <- done[[pid]]
      if(!is.list(item)||!isTRUE(item$ok)) {
        failed <- c(failed,list(list(unit=i,message=if(is.list(item))item$message else as.character(item),captured=item)))
        pulse("failed",i,length(jobs))
      } else {
        err <- tryCatch({accept(item$result,i);NULL},error=function(err)err)
        if(inherits(err,"error")) {
          failed <- c(failed,list(list(unit=i,message=conditionMessage(err),captured=item)))
          pulse("failed",i,length(jobs))
        }
        else {completed <- c(completed,i);pulse("accepted",i,length(jobs))}
      }
    }
    pulse("heartbeat",NA_integer_,length(jobs))
    if(length(jobs))Sys.sleep(.1)
  }
  list(completed=completed,launched=launched,failed=failed,pending=queue)
}
