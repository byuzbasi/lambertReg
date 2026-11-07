# Active-set cache used by the MAGIC maize study.
sch_cache_key <- function(active)paste0("active_",cpk_object_hash(as.integer(active)+0L))
sch_df <- function(Z,active,lambda,alpha,cache=NULL) {
  stopifnot(is.finite(lambda),lambda>0,alpha>=0,alpha<1)
  if(!length(active))return(0)
  active<-as.integer(active);key<-sch_cache_key(active)
  found<-if(!is.null(cache))get0(key,envir=cache,inherits=FALSE) else NULL
  ev<-NULL
  if(!is.null(found)) {stopifnot(identical(found$active,active));ev<-found$values}
  if(is.null(ev)) {
    A<-Z[,active,drop=FALSE]
    G<-if(ncol(A)<=nrow(A))crossprod(A)/nrow(A) else tcrossprod(A)/nrow(A)
    ev<-eigen(G,symmetric=TRUE,only.values=TRUE)$values
    stopifnot(min(ev)>= -1e-10*max(1,max(abs(ev))))
    ev<-pmax(ev,0)
    if(!is.null(cache))assign(key,list(active=active,values=ev),envir=cache)
  }
  sum(ev/(ev+lambda*(1-alpha)))
}
# Bind the replacement only inside this adapter's scorer; all formulas stay identical.
sch_candidates <- local({
  candidate<-lbc_candidates;private_df<-new.env(parent=environment(candidate));private_df$lbc_df<-sch_df
  environment(candidate)<-private_df
  scorer<-sch_candidates;private_score<-new.env(parent=environment(scorer));private_score$lbc_candidates<-candidate
  environment(scorer)<-private_score;scorer
})
