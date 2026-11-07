# Package interface used by the MAGIC maize study.
rfh_compile <- function(e) {
  stopifnot(as.character(utils::packageVersion("lambertReg"))=="0.1.11")
  ns<-asNamespace("lambertReg")
  solver<-get("lambert_cd_guarded_arma_v2",ns);context<-get("lambert_context_guarded_arma_v2",ns)
  e$lambert_values_cpp<-get("lambert_values_cpp",ns)
  e$lambert_cd_cpp<-function(X,y,init,lambda,shape,max_sweeps,kkt_tol,pulse)
    solver(X,y,init,lambda,shape,max_sweeps,kkt_tol,pulse)
  prior<-e$perf_path
  e$perf_path<-function(prep,setting,cfg,pulse=function(...)NULL,last=cfg$lambda_points) {
    if(setting$family!="lambert")return(prior(prep,setting,cfg,pulse,last))
    began<-proc.time()[[3]];cache<-context(prep$Z,prep$y);previous<-e$lambert_cd_cpp
    on.exit(e$lambert_cd_cpp<-previous,add=TRUE)
    e$lambert_cd_cpp<-function(X,y,init,lambda,shape,max_sweeps,kkt_tol,pulse)
      solver(X,y,init,lambda,shape,max_sweeps,kkt_tol,pulse,FALSE,cache)
    ans<-e$pacc_path(prep,setting,cfg,pulse,last,e)
    ans$elapsed<-proc.time()[[3]]-began;ans
  }
  invisible(e)
}
