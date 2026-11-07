# Data integrity check shared by the real-data engines.
by_manifest <- function(root,dir) {
  p<-file.path(root,dir);m<-read.csv(file.path(p,"manifest.csv"),stringsAsFactors=FALSE)
  stopifnot(!anyDuplicated(m$file),all(file.exists(file.path(p,m$file))))
  for(i in seq_len(nrow(m)))stopifnot(file.info(file.path(p,m$file[i]))$size==m$bytes[i],
    digest::digest(file=file.path(p,m$file[i]),algo="sha256",serialize=FALSE)==m$sha256[i])
  invisible(m)
}
