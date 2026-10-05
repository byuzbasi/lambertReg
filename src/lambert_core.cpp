#include <Rcpp.h>
#include <cmath>
#include <chrono>
using namespace Rcpp;
// Numerical implementation of the frozen rule in R/lambert_family.R.
double capped_w(double t, double lambda, double shape) {
  if (t == 0) return 0;
  if (t >= lambda * std::exp(shape)) return shape;
  double x=shape*t/lambda, w=std::min(shape,std::log1p(x));
  for(int k=0;k<30;++k) {
    double ew=std::exp(w), f=w*ew-x;
    double d=f/(ew*(w+1)-(w+2)*f/(2*w+2));
    w-=d;
    if(std::abs(d)<=2e-15*(1+w)) return std::max(0.0,std::min(shape,w));
  }
  stop("Lambert root did not converge"); return NA_REAL;
}
double pen(double t,double lambda,double shape) {
  double w=capped_w(t,lambda,shape);
  return lambda*lambda*(std::expm1(2*w)/(4*shape)+std::exp(2*w)*w*(shape-w)/(2*shape*shape));
}
double slope(double t,double lambda,double shape) {
  double w=capped_w(t,lambda,shape);
  return lambda*std::exp(w)*(1-w/shape);
}
double threshold(double z,double lambda,double shape) {
  double a=std::abs(z);
  if(a<=lambda) return 0;
  if(a>=lambda*std::exp(shape)) return z;
  return z*std::log1p((a-lambda)/lambda)/shape;
}
// [[Rcpp::export]]
List lambert_values_cpp(NumericVector t,double lambda,double shape) {
  if(!(lambda>0)||!std::isfinite(lambda)||!(shape>=1)||!std::isfinite(shape))
    stop("Invalid Lambert parameters");
  int p=t.size(); NumericVector pv(p),ps(p),pc(p);
  double cutoff=lambda*std::exp(shape);
  for(int j=0;j<p;++j) {
    if(!std::isfinite(t[j])||t[j]<0) stop("Invalid magnitude");
    double w=capped_w(t[j],lambda,shape);
    pv[j]=pen(t[j],lambda,shape); ps[j]=slope(t[j],lambda,shape);
    pc[j]=t[j]==cutoff ? NA_REAL : t[j]>cutoff ? 0 : shape/(1+w)-1;
  }
  return List::create(_["penalty"]=pv,_["slope"]=ps,_["curvature"]=pc);
}
// [[Rcpp::export]]
List lambert_cd_cpp(NumericMatrix X,NumericVector y,NumericVector init,
                    double lambda,double shape,int max_sweeps,double kkt_tol,Function pulse) {
  int n=X.nrow(),p=X.ncol();
  if(n<1||y.size()!=n||init.size()!=p||max_sweeps<1) stop("Invalid dimensions");
  NumericVector b=clone(init),r=clone(y); double score_max=0;
  for(int j=0;j<p;++j) {
    double norm=0,score=0;
    for(int i=0;i<n;++i) {
      norm+=X(i,j)*X(i,j)/n; score+=X(i,j)*y[i]/n; r[i]-=X(i,j)*b[j];
    }
    if(std::abs(norm-1)>1e-10) stop("Coordinate curvature must equal 1");
    score_max=std::max(score_max,std::abs(score));
  }
  auto objective=[&]() {
    double q=0;
    for(int i=0;i<n;++i) q+=r[i]*r[i]/(2*n);
    for(int j=0;j<p;++j) q+=pen(std::abs(b[j]),lambda,shape);
    return q;
  };
  double previous=objective(),max_increase=0,kkt=R_PosInf,q=previous;
  bool converged=false,increase_failed=false; int sweep=0;
  auto last_pulse=std::chrono::steady_clock::now();
  for(sweep=1;sweep<=max_sweeps;++sweep) {
    for(int j=0;j<p;++j) {
      double z=b[j];
      for(int i=0;i<n;++i) z+=X(i,j)*r[i]/n;
      double next=threshold(z,lambda,shape),delta=next-b[j]; b[j]=next;
      if(delta!=0) for(int i=0;i<n;++i) r[i]-=X(i,j)*delta;
    }
    if(sweep<=3||sweep%25==0||sweep==max_sweeps) {
      for(int i=0;i<n;++i) {
        r[i]=y[i]; for(int j=0;j<p;++j) r[i]-=X(i,j)*b[j];
      }
      kkt=0;
      for(int j=0;j<p;++j) {
        double g=0;
        for(int i=0;i<n;++i) g-=X(i,j)*r[i]/n;
        double e=b[j]==0 ? std::max(0.0,std::abs(g)-lambda) :
          std::abs(g+(b[j]>0 ? 1 : -1)*slope(std::abs(b[j]),lambda,shape));
        kkt=std::max(kkt,e);
      }
      kkt/=1+score_max; q=objective();
      max_increase=std::max(max_increase,(q-previous)/(1+std::abs(previous))); previous=q;
      if(!std::isfinite(q)||max_increase>1e-9) { increase_failed=true; break; }
      if(kkt<=kkt_tol) { converged=true; break; }
      checkUserInterrupt();
      auto now=std::chrono::steady_clock::now();
      if(std::chrono::duration<double>(now-last_pulse).count()>=20) {
        pulse("Lambert coordinate updates"); last_pulse=now;
      }
    }
  }
  return List::create(_["beta"]=b,_["iterations"]=std::min(sweep,max_sweeps),
    _["converged"]=converged,_["objective"]=q,_["kkt"]=kkt,
    _["max_increase"]=max_increase,_["increase_failed"]=increase_failed);
}
