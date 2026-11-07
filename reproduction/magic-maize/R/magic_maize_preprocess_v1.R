# Training-only genotype preprocessing; deliberately accepts no response.
mmp_settings <- function()list(call_rate=.95,maf=.05,ld_r2=.4,window_bp=500000)

mmp_source <- function(root=getwd()) {
  raw<-"data/maize"
  prior<-"data/maize-marker-map"
  by_manifest(root,raw);by_manifest(root,prior)
  m<-read.csv(file.path(root,prior,"marker_inventory.csv"),stringsAsFactors=FALSE)
  con<-gzfile(file.path(root,raw,"MMlines.geno.gz"),"rt");on.exit(close(con))
  txt<-readLines(con);z<-strsplit(txt,"\t",fixed=TRUE)
  if(length(z[[1]])==1L)z<-strsplit(trimws(txt),"[[:space:]]+")
  stopifnot(identical(z[[1]][-1],m$marker),length(z)==530L,all(lengths(z)==54235L))
  ids<-vapply(z[-1],`[`,"",1L);C<-do.call(rbind,lapply(z[-1],`[`,-1L))
  X<-matrix(NA_real_,nrow(C),ncol(C),dimnames=list(ids,m$marker))
  for(j in which(m$homozygous_bases>0L & m$homozygous_bases<=2L)) {
    X[C[,j]==m$allele_0[j],j]<-0
    if(!is.na(m$allele_2[j]))X[C[,j]==m$allele_2[j],j]<-2
    X[C[,j]=="H",j]<-1
  }
  ph<-as.data.frame(readxl::read_excel(file.path(root,raw,"Table_S6.xlsx")))
  stopifnot(!anyDuplicated(ids),!anyNA(match(ids,ph$Line)))
  list(X=X,annotation=data.frame(marker=m$marker,chromosome=m$chromosome,
    position_bp=round(m$position_mb*1e6)),sample_id=ids,
    y=ph$PH[match(ids,ph$Line)],raw=raw,prior=prior)
}

mmp_pair_cor <- function(x,B) {
  if(!ncol(B))return(list(r=numeric(),n=integer()))
  stats::setNames(list(as.numeric(suppressWarnings(stats::cor(x,B,use="pairwise.complete.obs"))),
    colSums(is.finite(B)&is.finite(x))),c("r","n"))
}

mmp_pattern <- function(v) {
  ok<-is.finite(v);levels<-sort(unique(v[ok]));z<-rep.int(255L,length(v))
  z[ok]<-match(v[ok],levels)-1L
  first<-which(ok & z!=(length(levels)-1)/2)[1]
  if(z[first]>(length(levels)-1)/2)z[ok]<-length(levels)-1L-z[ok]
  digest::digest(as.raw(z),algo="sha256",serialize=FALSE)
}

mmp_fit <- function(X,annotation,settings=mmp_settings(),pulse=function(...)NULL,
                    purpose="training") {
  stopifnot(is.matrix(X),is.numeric(X),nrow(X)>=3L,ncol(X)>0L,
    !anyDuplicated(colnames(X)),!anyDuplicated(annotation$marker),
    all(is.na(X)|X%in%c(0,1,2)),all(is.finite(annotation$position_bp)),
    identical(settings,mmp_settings()))
  ann<-annotation[match(colnames(X),annotation$marker),,drop=FALSE]
  stopifnot(!anyNA(ann$marker),all(is.finite(ann$chromosome)))
  n<-nrow(X);p<-ncol(X);obs<-colSums(!is.na(X));Y<-X;Y[is.na(Y)]<-0
  mu<-colSums(Y)/obs;vv<-colSums(Y*Y)/obs-mu^2;maf<-pmin(mu/2,1-mu/2);rm(Y)
  reason<-rep("eligible",p);reason[obs==0L]<-"no_observed_calls"
  reason[obs>0L & vv==0]<-"constant_observed"
  reason[reason=="eligible" & obs/n<settings$call_rate]<-"low_call_rate"
  reason[reason=="eligible" & maf<settings$maf]<-"low_MAF"
  eligible<-which(reason=="eligible")
  stopifnot(length(eligible)>0L)
  priority<-eligible[order(-obs[eligible],-maf[eligible],ann$chromosome[eligible],
    ann$position_bp[eligible],ann$marker[eligible])]
  representative<-rep(NA_character_,p);blocked_r<-rep(NA_real_,p);blocked_n<-rep(NA_integer_,p)
  seen<-new.env(hash=TRUE,parent=emptyenv());unique_order<-integer()
  for(k in seq_along(priority)) {
    j<-priority[k];key<-mmp_pattern(X[,j])
    if(exists(key,seen,inherits=FALSE)) {
      b<-get(key,seen,inherits=FALSE);reason[j]<-"duplicate_pattern"
      representative[j]<-ann$marker[b]
    } else {assign(key,j,seen);unique_order<-c(unique_order,j)}
    if(k%%2000L==0L)pulse("duplicate_scan",k,length(priority),length(unique_order))
  }
  kept<-integer();kept_chr<-split(integer(),integer())
  for(k in seq_along(unique_order)) {
    j<-unique_order[k];chr<-as.character(ann$chromosome[j]);pool<-kept_chr[[chr]]
    near<-pool[abs(ann$position_bp[pool]-ann$position_bp[j])<=settings$window_bp]
    rr<-mmp_pair_cor(X[,j],X[,near,drop=FALSE]);bad<-which(rr$r^2>settings$ld_r2)
    if(length(bad)) {
      b<-bad[1];reason[j]<-"LD_pruned";representative[j]<-ann$marker[near[b]]
      blocked_r[j]<-rr$r[b];blocked_n[j]<-rr$n[b]
    } else {kept<-c(kept,j);kept_chr[[chr]]<-c(pool,j);reason[j]<-"retained"}
    if(k%%1000L==0L)pulse("LD_pruning",k,length(unique_order),length(kept))
  }
  # Genomic output order, independent of quality-priority visitation order.
  kept<-kept[order(ann$chromosome[kept],ann$position_bp[kept],ann$marker[kept])]
  M<-X[,kept,drop=FALSE];for(k in seq_along(kept))M[is.na(M[,k]),k]<-mu[kept[k]]
  center<-mu[kept];Z<-sweep(M,2,center,"-");scale<-sqrt(colMeans(Z^2))
  stopifnot(all(is.finite(scale)),all(scale>0),all(reason!="eligible"))
  ledger<-data.frame(ann,observed=obs,call_rate=obs/n,observed_mean=mu,MAF=maf,
    reason=reason,representative=representative,blocker_r=blocked_r,blocker_n=blocked_n)
  list(schema="magic_maize_preprocess_v1",purpose=purpose,settings=settings,
    training_n=n,training_ids=rownames(X),source_markers=colnames(X),
    retained=ann$marker[kept],center=unname(center),scale=unname(scale),
    ledger=ledger,eligible=ann$marker[eligible],quality_order=ann$marker[priority],
    ld_visit_order=ann$marker[unique_order])
}

mmp_transform <- function(context,X_new,standardize=TRUE) {
  stopifnot(is.matrix(X_new),is.numeric(X_new),!anyDuplicated(colnames(X_new)),
    all(is.na(X_new)|X_new%in%c(0,1,2)))
  idx<-match(context$retained,colnames(X_new));stopifnot(!anyNA(idx))
  A<-X_new[,idx,drop=FALSE]
  for(j in seq_len(ncol(A)))A[is.na(A[,j]),j]<-context$center[j]
  if(standardize)A<-sweep(sweep(A,2,context$center,"-"),2,context$scale,"/")
  stopifnot(all(is.finite(A)));A
}

mmp_pairs <- function(p,limit=100000L) {
  target<-min(limit,p*(p-1)/2);if(target==0)return(data.frame(j=integer(),k=integer()))
  if(p*(p-1)/2<=limit){z<-as.data.frame(t(utils::combn(p,2)));names(z)<-c("j","k");return(z)}
  set.seed(2026102701L,kind="Mersenne-Twister",normal.kind="Inversion",sample.kind="Rejection")
  z<-data.frame(j=integer(),k=integer())
  while(nrow(z)<target){a<-sample.int(p,120000L,TRUE);b<-sample.int(p,120000L,TRUE)
    z<-unique(rbind(z,data.frame(j=pmin(a,b),k=pmax(a,b))));z<-z[z$j!=z$k,]}
  head(z,target)
}

mmp_geometry <- function(A,pairs) {
  z<-pairs;z$r<-NA_real_;z$n_pair<-0L
  for(k in seq_len(nrow(z))) {r<-mmp_pair_cor(A[,z$j[k]],A[,z$k[k],drop=FALSE])
    z$r[k]<-r$r;z$n_pair[k]<-r$n}
  v<-abs(z$r[is.finite(z$r)]);q<-if(length(v))quantile(v,c(.5,.9,.95,.99,1),names=FALSE) else rep(NA_real_,5)
  list(pairs=z,summary=data.frame(pairs=nrow(z),available=length(v),median=q[1],
    q90=q[2],q95=q[3],q99=q[4],max=q[5],fraction_ge_05=mean(v>=.5),fraction_ge_09=mean(v>=.9)))
}
