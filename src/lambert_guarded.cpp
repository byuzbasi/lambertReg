// Opt-in ordered solver: per-path validation and per-coordinate zero-score cache.
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <chrono>
#include <cmath>
#include <limits>
#include <memory>
using namespace Rcpp;
namespace {
using Clock2=std::chrono::steady_clock;
double sec2(Clock2::time_point t){return std::chrono::duration<double>(Clock2::now()-t).count();}
double w2(double t,double lambda,double shape,double cutoff) {
  if(t==0)return 0;if(t>=cutoff)return shape;
  double x=shape*t/lambda,w=std::min(shape,std::log1p(x));
  for(int k=0;k<30;++k){double ew=std::exp(w),f=w*ew-x,d=f/(ew*(w+1)-(w+2)*f/(2*w+2));
    w-=d;if(std::abs(d)<=2e-15*(1+w))return std::max(0.0,std::min(shape,w));}
  stop("Lambert root did not converge");return NA_REAL;
}
double pen2(double w,double lambda,double shape) {
  return lambda*lambda*(std::expm1(2*w)/(4*shape)+std::exp(2*w)*w*(shape-w)/(2*shape*shape));
}
struct Context2 {
  NumericMatrix X;NumericVector y;int n,p;double score_max=0,seconds=0;
  Context2(NumericMatrix xx,NumericVector yy):X(xx),y(yy),n(X.nrow()),p(X.ncol()) {
    auto started=Clock2::now();
    if(n<1||p<1||y.size()!=n)stop("Invalid dimensions");
    const arma::mat view(X.begin(),n,p,false,true);
    if(!view.is_finite())stop("Nonfinite input");
    for(double value:y)if(!std::isfinite(value))stop("Nonfinite input");
    const double* yp=y.begin();
    for(int j=0;j<p;++j) {
      const double* col=view.colptr(j);double norm=0,score=0;
      for(int i=0;i<n;++i){norm+=col[i]*col[i]/n;score+=col[i]*yp[i]/n;}
      if(std::abs(norm-1)>1e-10)stop("Coordinate curvature must equal 1");
      score_max=std::max(score_max,std::abs(score));
    }
    seconds=sec2(started);
  }
};
}
// [[Rcpp::export]]
SEXP lambert_context_guarded_arma_v2(NumericMatrix X,NumericVector y) {
  XPtr<Context2> cache(new Context2(X,y),true);
  cache.attr("class")="lambert_guarded_context_v2";
  cache.attr("preparation_seconds")=cache->seconds;
  return cache;
}
// Same leading solver arguments; a cache is optional and scoped to one path.
// [[Rcpp::export]]
List lambert_cd_guarded_arma_v2(NumericMatrix X,NumericVector y,NumericVector init,
    double lambda,double shape,int max_sweeps,double kkt_tol,Function pulse,
    bool profile=false,SEXP cache=R_NilValue) {
  auto entered=Clock2::now();int n=X.nrow(),p=X.ncol();
  if(n<1||p<1||y.size()!=n||init.size()!=p||max_sweeps<1)stop("Invalid dimensions");
  if(!(lambda>0)||!std::isfinite(lambda)||!(shape>=1)||!std::isfinite(shape)||
    !(kkt_tol>0)||!std::isfinite(kkt_tol))stop("Invalid controls");
  for(double value:init)if(!std::isfinite(value))stop("Nonfinite input");
  std::unique_ptr<Context2> owned;const Context2* context;
  bool cached=cache!=R_NilValue;
  if(cached) {
    if(TYPEOF(cache)!=EXTPTRSXP||!Rf_inherits(cache,"lambert_guarded_context_v2")||
      R_ExternalPtrAddr(cache)==nullptr)stop("Invalid path context");
    context=static_cast<Context2*>(R_ExternalPtrAddr(cache));
    if(context->n!=n||context->p!=p||context->X.begin()!=X.begin()||context->y.begin()!=y.begin())
      stop("Path context input mismatch");
  } else {owned.reset(new Context2(X,y));context=owned.get();}
  const double* xp=X.begin();NumericVector b=clone(init),r=clone(y);
  double* bp=b.begin();double* rp=r.begin();
  arma::vec margin(p,arma::fill::zeros),anchor(p,arma::fill::zeros),gradient(p,arma::fill::zeros);
  arma::vec previous_b(p),roots(p,arma::fill::zeros);
  previous_b.fill(NA_REAL);arma::uvec ever(p,arma::fill::zeros);
  double score_max=context->score_max,cutoff=lambda*std::exp(shape);
  int active_count=0;
  for(int j=0;j<p;++j)if(bp[j]!=0) {
    const double* col=xp+n*j;for(int i=0;i<n;++i)rp[i]-=col[i]*bp[j];
    ever[j]=1;++active_count;
  }
  auto refresh_root=[&](int j) {
    if(bp[j]==previous_b[j])return false;
    roots[j]=w2(std::abs(bp[j]),lambda,shape,cutoff);
    previous_b[j]=bp[j];return true;
  };
  double q=0;for(int i=0;i<n;++i)q+=rp[i]*rp[i]/(2*n);
  for(int j=0;j<p;++j)if(bp[j]!=0){refresh_root(j);q+=pen2(roots[j],lambda,shape);}
  int initial_active=active_count,max_active=active_count,checks=0,sweep=0;
  double previous=q,max_increase=0,kkt=R_PosInf,movement=0,guard=0;
  bool cache_valid=false,converged=false,increase_failed=false;
  double updates=0,skipped=0,refreshed=0,root_refreshes=0,root_reuses=0;
  double update_seconds=0,residual_seconds=0,gradient_seconds=0,penalty_seconds=0;
  double initialization_seconds=sec2(entered);
  auto began=Clock2::now(),last_pulse=began,tick=began;
  for(sweep=1;sweep<=max_sweeps;++sweep) {
    if(profile)tick=Clock2::now();
    for(int j=0;j<p;++j) {
      if(sweep>3&&cache_valid&&bp[j]==0) {
        double distance=std::nextafter(std::max(0.0,movement-anchor[j])+
          8*std::numeric_limits<double>::epsilon()*(1+std::abs(movement)+std::abs(anchor[j])),R_PosInf);
        if(margin[j]>distance+guard){bp[j]=0;++skipped;continue;}
      }
      const double* col=xp+n*j;double old=bp[j],z=old;
      for(int i=0;i<n;++i)z+=col[i]*rp[i]/n;
      double a=std::abs(z),next=a<=lambda?0:a>=cutoff?z:z*std::log1p((a-lambda)/lambda)/shape;
      double delta=next-old;bp[j]=next;++updates;
      if(delta!=0) {
        for(int i=0;i<n;++i)rp[i]-=col[i]*delta;
        movement=std::nextafter(movement+std::abs(delta)*(1+1e-8),R_PosInf);
      } else if(cache_valid&&old==0&&next==0) {
        margin[j]=lambda-std::abs(z);anchor[j]=movement;++refreshed;
      }
      if(next!=0&&!ever[j]){ever[j]=1;++active_count;max_active=std::max(max_active,active_count);}
    }
    if(profile)update_seconds+=sec2(tick);
    if(sweep<=3||sweep%25==0||sweep==max_sweeps) {
      if(profile)tick=Clock2::now();++checks;
      std::copy(y.begin(),y.end(),rp);
      for(int j=0;j<p;++j)if(bp[j]!=0) {
        const double* col=xp+n*j;for(int i=0;i<n;++i)rp[i]-=col[i]*bp[j];
      }
      q=0;double rmax=0;for(int i=0;i<n;++i){q+=rp[i]*rp[i]/(2*n);rmax=std::max(rmax,std::abs(rp[i]));}
      if(profile){residual_seconds+=sec2(tick);tick=Clock2::now();}
      for(int j=0;j<p;++j) {
        const double* col=xp+n*j;double g=0;for(int i=0;i<n;++i)g-=col[i]*rp[i]/n;
        gradient[j]=g;margin[j]=lambda-std::abs(g);anchor[j]=0;
      }
      if(profile){gradient_seconds+=sec2(tick);tick=Clock2::now();}
      kkt=0;
      for(int j=0;j<p;++j) {
        double error;
        if(bp[j]==0)error=std::max(0.0,std::abs(gradient[j])-lambda);
        else {
          if(refresh_root(j))++root_refreshes;else ++root_reuses;
          double w=roots[j];
          // Keep the derivative/addition expression intact: caching its product
          // changes the final rounding when the compiler contracts operations.
          error=std::abs(gradient[j]+(bp[j]>0?1:-1)*lambda*std::exp(w)*(1-w/shape));
          q+=pen2(w,lambda,shape);
        }
        kkt=std::max(kkt,error);
      }
      cache_valid=true;movement=0;
      guard=256*std::numeric_limits<double>::epsilon()*(n+p)*(1+score_max+lambda+rmax);
      kkt/=1+score_max;
      max_increase=std::max(max_increase,(q-previous)/(1+std::abs(previous)));previous=q;
      if(profile)penalty_seconds+=sec2(tick);
      if(!std::isfinite(q)||max_increase>1e-9){increase_failed=true;break;}
      if(kkt<=kkt_tol){converged=true;break;}
      checkUserInterrupt();
      if(sec2(last_pulse)>=20){pulse("Lambert ordered cached updates");last_pulse=Clock2::now();}
    }
  }
  return List::create(_["beta"]=b,_["iterations"]=std::min(sweep,max_sweeps),
    _["converged"]=converged,_["objective"]=q,_["kkt"]=kkt,
    _["max_increase"]=max_increase,_["increase_failed"]=increase_failed,
    _["profile_enabled"]=profile,_["update_seconds"]=profile?update_seconds:NA_REAL,
    _["diagnostic_seconds"]=profile?residual_seconds+gradient_seconds+penalty_seconds:NA_REAL,
    _["residual_seconds"]=profile?residual_seconds:NA_REAL,_["gradient_seconds"]=profile?gradient_seconds:NA_REAL,
    _["penalty_seconds"]=profile?penalty_seconds:NA_REAL,_["initialization_seconds"]=initialization_seconds,
    _["diagnostic_checks"]=checks,_["full_scans"]=checks,_["coordinate_updates"]=updates,
    _["skipped_zero_updates"]=skipped,_["zero_cache_refreshes"]=refreshed,
    _["penalty_refreshes"]=root_refreshes,_["penalty_reuses"]=root_reuses,
    _["initial_active"]=initial_active,_["max_active"]=max_active,_["active_additions"]=active_count-initial_active,
    _["context_cached"]=cached,_["solver_seconds"]=sec2(began));
}
