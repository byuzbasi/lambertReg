// Adapted from the validated fast kernel; see inst/PROVENANCE.json.
#include <Rcpp.h>
#include <cmath>
#include <chrono>
using namespace Rcpp;
namespace {
// Same arithmetic and root tolerance as the frozen implementation.
double fast_w(double t,double lambda,double shape,double cutoff) {
  if(t==0) return 0;
  if(t>=cutoff) return shape;
  double x=shape*t/lambda,w=std::min(shape,std::log1p(x));
  for(int k=0;k<30;++k) {
    double ew=std::exp(w),f=w*ew-x,d=f/(ew*(w+1)-(w+2)*f/(2*w+2));
    w-=d;
    if(std::abs(d)<=2e-15*(1+w)) return std::max(0.0,std::min(shape,w));
  }
  stop("Lambert root did not converge"); return NA_REAL;
}
double fast_pen(double w,double lambda,double shape) {
  return lambda*lambda*(std::expm1(2*w)/(4*shape)+std::exp(2*w)*w*(shape-w)/(2*shape*shape));
}
using Clock=std::chrono::steady_clock;
double seconds(Clock::time_point start) {
  return std::chrono::duration<double>(Clock::now()-start).count();
}
}
// [[Rcpp::export]]
List lambert_values_cpp(NumericVector t,double lambda,double shape) {
  if(!(lambda>0)||!std::isfinite(lambda)||!(shape>=1)||!std::isfinite(shape))
    stop("Invalid Lambert parameters");
  int p=t.size(); NumericVector pv(p),ps(p),pc(p);
  double cutoff=lambda*std::exp(shape);
  for(int j=0;j<p;++j) {
    if(!std::isfinite(t[j])||t[j]<0) stop("Invalid magnitude");
    double w=fast_w(t[j],lambda,shape,cutoff);
    pv[j]=fast_pen(w,lambda,shape); ps[j]=lambda*std::exp(w)*(1-w/shape);
    pc[j]=t[j]==cutoff ? NA_REAL : t[j]>cutoff ? 0 : shape/(1+w)-1;
  }
  return List::create(_["penalty"]=pv,_["slope"]=ps,_["curvature"]=pc);
}
// [[Rcpp::export]]
List lambert_cd_cpp(NumericMatrix X,NumericVector y,NumericVector init,
    double lambda,double shape,int max_sweeps,double kkt_tol,Function pulse) {
  int n=X.nrow(),p=X.ncol();
  if(n<1||y.size()!=n||init.size()!=p||max_sweeps<1) stop("Invalid dimensions");
  if(!(lambda>0)||!std::isfinite(lambda)||!(shape>=1)||!std::isfinite(shape)||
     !std::isfinite(kkt_tol)||kkt_tol<=0) stop("Invalid controls");
  NumericVector b=clone(init),r=clone(y); double score_max=0,cutoff=lambda*std::exp(shape);
  for(int j=0;j<p;++j) {
    const double* col=X.begin()+n*j; double norm=0,score=0;
    for(int i=0;i<n;++i) {
      norm+=col[i]*col[i]/n; score+=col[i]*y[i]/n; r[i]-=col[i]*b[j];
    }
    if(std::abs(norm-1)>1e-10) stop("Coordinate curvature must equal 1");
    score_max=std::max(score_max,std::abs(score));
  }
  double q=0;
  for(int i=0;i<n;++i) q+=r[i]*r[i]/(2*n);
  for(int j=0;j<p;++j) q+=fast_pen(fast_w(std::abs(b[j]),lambda,shape,cutoff),lambda,shape);
  double previous=q,max_increase=0,kkt=R_PosInf;
  bool converged=false,increase_failed=false; int sweep=0;
  auto last_pulse=Clock::now();
  for(sweep=1;sweep<=max_sweeps;++sweep) {
    for(int j=0;j<p;++j) {
      const double* col=X.begin()+n*j; double z=b[j];
      for(int i=0;i<n;++i) z+=col[i]*r[i]/n;
      double a=std::abs(z),next=a<=lambda ? 0 : a>=cutoff ? z : z*std::log1p((a-lambda)/lambda)/shape;
      double delta=next-b[j]; b[j]=next;
      if(delta!=0) for(int i=0;i<n;++i) r[i]-=col[i]*delta;
    }
    if(sweep<=3||sweep%25==0||sweep==max_sweeps) {
      // Column-contiguous traversal; each residual keeps the same j summation order.
      std::copy(y.begin(),y.end(),r.begin());
      for(int j=0;j<p;++j) {
        const double* col=X.begin()+n*j;
        for(int i=0;i<n;++i) r[i]-=col[i]*b[j];
      }
      q=0; for(int i=0;i<n;++i) q+=r[i]*r[i]/(2*n);
      kkt=0;
      for(int j=0;j<p;++j) {
        const double* col=X.begin()+n*j; double g=0;
        for(int i=0;i<n;++i) g-=col[i]*r[i]/n;
        double w=fast_w(std::abs(b[j]),lambda,shape,cutoff);
        double e=b[j]==0 ? std::max(0.0,std::abs(g)-lambda) :
          std::abs(g+(b[j]>0 ? 1 : -1)*lambda*std::exp(w)*(1-w/shape));
        kkt=std::max(kkt,e); q+=fast_pen(w,lambda,shape);
      }
      kkt/=1+score_max;
      max_increase=std::max(max_increase,(q-previous)/(1+std::abs(previous))); previous=q;
      if(!std::isfinite(q)||max_increase>1e-9) { increase_failed=true; break; }
      if(kkt<=kkt_tol) { converged=true; break; }
      checkUserInterrupt();
      if(seconds(last_pulse)>=20) { pulse("Lambert coordinate updates"); last_pulse=Clock::now(); }
    }
  }
  return List::create(_["beta"]=b,_["iterations"]=std::min(sweep,max_sweeps),
    _["converged"]=converged,_["objective"]=q,_["kkt"]=kkt,
    _["max_increase"]=max_increase,_["increase_failed"]=increase_failed);
}
