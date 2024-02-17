#include <cmath>  // std::pow
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <roptim.h>
// [[Rcpp::depends(roptim)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
// [[Rcpp::depends(RcppProgress)]]

/* Library functions of the paper:
 Classical and Bayesian inference of C_{pc} 
 for  Wilson-Hilferty Distribution Under Progressively
 First Failure  Type-II Censoring Samples
 by: Sanku Dey and  Riyadh Al-Mosawi
*/

//#include<iostream>
#include<algorithm>
#include"post.h"

using namespace Rcpp;
using namespace std;
using namespace Numer;
using namespace arma;
using namespace roptim;

// These codes for the progress bar inside the loop
#include <progress.hpp>
#include <progress_bar.hpp>

#include<stdio.h>
#include<ctype.h>

// Compute gamma function
double Gamma(double x){
  Function gam_rcpp("gamma");
  NumericVector res=gam_rcpp(x);
  return res[0];
};

// Compute upper incomplete gamma function
double Incgamma(double x, double y){
  Function incgam_rcpp("incgam");
  NumericVector res=incgam_rcpp(std::max(0.0,x), y);
  return res[0];
};

// Compute digamma (the firts derivative of log(Gamma)) function
double Digamma(double x){
  Function digam_rcpp("digamma");
  NumericVector res=digam_rcpp(x);
  return res[0];
};

// Compute trigamma (the second derivative of log(Gamma)) function
double Trigamma(double x){
  Function trigam_rcpp("trigamma");
  NumericVector res=trigam_rcpp(x);
  return res[0];
};

// Compute trigamma (the second derivative of log(Gamma)) function
double Psigamma(double x){
  Function psigam_rcpp("psigamma");
  NumericVector res=psigam_rcpp(x);
  return res[0];
};


// The pdf of WH Dist
// [[Rcpp::export]]
double pdf(double x, Rcpp::NumericVector para){
  return 3.0*(pow(x,3.0*para[0]-1)/Gamma(para[0]))*pow(para[0]/para[1],para[0])*std::exp(-para[0]*pow(x,3.0)/para[1]);
}

// The cdf of WH Dist
// [[Rcpp::export]]
double cdf(double x, Rcpp::NumericVector para){
  return 1.0-Incgamma(para[0]*pow(x,3)/para[1],para[0])/Gamma(para[0]);
}
// The survival of WH Dist
// [[Rcpp::export]]
double sur(double x, Rcpp::NumericVector para){
  return 1.0-cdf(x, para);
}

// This function is used to compute quantile of WH dist
double fun1(double x, double p,Rcpp::NumericVector para){
  return cdf(x,para)-p;
}


// The quantile of WH Dist
// [[Rcpp::export]]
double quan(Rcpp::NumericVector init, double p, Rcpp::NumericVector para){
  Rcpp::Environment base("package:stats"); 
  Rcpp::Function nr = base["uniroot"];
  Rcpp::List res = nr(Rcpp::_["f"]     = Rcpp::InternalFunction(&fun1),
                      Rcpp::_["interval"]= init,
                      Rcpp::_["p"]  = p,
                      Rcpp::_["para"]=para);
  return (double) res[0];
}

NumericVector my_as_mcmc(NumericVector x){
  Rcpp::Function as_mcmc_rcpp("as.mcmc");
  Rcpp::NumericVector res=as_mcmc_rcpp(x);
  return res;
}


// This function is used to compute highest posterior density interval
NumericVector my_HPD(NumericVector x){
  Rcpp::Function HPD_rcpp("HPDinterval");
  Rcpp::NumericVector res=HPD_rcpp(x);
  return res;
}


// These classes are used to compute derivative of incomplete gamma
class P1: public Func {
private:
  double alp;
  double lam;
public:
  P1(double a_, double b_) : alp(a_), lam(b_) {}
  
  double operator()(const double& x) const {
    return std::log(x)*std::pow(x,alp-1)*exp(-x);
  };
};

// Definition of Cpc function
// [[Rcpp::export]]
double cpc_fun(Rcpp::NumericVector para, double L, double U, double P0){
  return((1.0-P0)/(1.0-cdf(U,para)+cdf(L,para)));
}

// Gradient of Cpc function
// [[Rcpp::export]]
Rcpp::NumericVector cpc_grad(Rcpp::NumericVector para, double L, double U, double P0){
  
  double alp=para[0];
  double lam=para[1];
  
  double err_est1;
  int err_code1;
  double err_est2;
  int err_code2;
  
  P1 f1(alp, lam);
  
  double I_l = integrate(f1, alp*L*L*L/lam, 3000, err_est1, err_code1);
  double I_u = integrate(f1, alp*U*U*U/lam, 3000, err_est2, err_code2);
  
  double psi_l=std::log(Incgamma(alp*pow(L,3)/lam,alp));
  double psi_u=std::log(Incgamma(alp*pow(U,3)/lam,alp));
  double psi_alp_l=(-pow(alp,alp-1)*exp(-alp*L*L*L/lam)*pow(L,3*alp)*pow(lam,-alp)+I_l)/Incgamma(alp*L*L*L/lam,alp);
  double psi_alp_u=(-pow(alp,alp-1)*exp(-alp*U*U*U/lam)*pow(U,3*alp)*pow(lam,-alp)+I_u)/Incgamma(alp*U*U*U/lam,alp);
  double psi_lam_l=pow(alp,alp)*exp(-alp*L*L*L/lam)*pow(L,3*alp)*pow(lam,-alp-1)/Incgamma(alp*L*L*L/lam,alp);
  double psi_lam_u=pow(alp,alp)*exp(-alp*U*U*U/lam)*pow(U,3*alp)*pow(lam,-alp-1)/Incgamma(alp*U*U*U/lam,alp);
  double P=(Incgamma(alp*L*L*L/lam,alp)-Incgamma(alp*U*U*U/lam,alp))/Gamma(alp);
  double g1=(1.0-P0)/((1.0-P)*(1.0-P));
  double g2=Incgamma(alp*L*L*L/lam,alp)*(psi_alp_l-Digamma(alp))-Incgamma(alp*U*U*U/lam,alp)*(psi_alp_u-Digamma(alp));
  double g3=Incgamma(alp*L*L*L/lam,alp)*psi_lam_l-Incgamma(alp*U*U*U/lam,alp)*psi_lam_u;
  double G1=g1*g2/Gamma(alp);
  double G2=g1*g3/Gamma(alp);
  return {G1,G2};
} 

// Generating progressive first failure time data from WH distribution 
// [[Rcpp::export]]
Rcpp::NumericVector GenDataRcpp(Rcpp::NumericVector init, 
                                Rcpp::NumericVector para,
                            Rcpp::NumericVector R, int k){  
  int m=R.size();
  Rcpp::NumericVector w(m);  
  Rcpp::NumericVector u(m);  
  Rcpp::NumericVector v(m);  
  Rcpp::NumericVector t(m);  
  Rcpp::NumericVector vv(m);  
  Rcpp::NumericVector RR(m);  
  w=Rcpp::runif(m,0.0,1.0);
  for(int i=0;i<m;i++){
    RR(i)=0;
    for(int j=m-i-1;j<m;j++){
      RR(i)+=R(j);
    }
    v(i)=std::pow(w(i),1/(k*(i+RR(i))));
  }
  for(int i=0;i<m;i++){
    vv(i)=1;
    for(int j=m-i-1;j<m;j++){
      vv(i)*=v(j);
    }
    u(i)=1.0-vv(i);
    t(i)=quan(init,u(i),para);
  };
  return t;
};


// likelihood function
double like(Rcpp::NumericVector para, 
            Rcpp::NumericVector X, 
            Rcpp::NumericVector R,
            int k){
  
  // Compute objective value
  double lk=1;
  int m=X.size();
  for(int i=0;i<m;i++){
    lk *= pdf(X(i),para)*pow(sur(X(i),para),k*(R(i)+1)-1);
  }
  return lk;
}  

// log-likelihood function
// [[Rcpp::export]]
double loglike(Rcpp::NumericVector para, 
               Rcpp::NumericVector X, 
               Rcpp::NumericVector R,
               int k){
  
  return -std::log(like(para,X,R,k));
}

class LogLike : public Functor {
  private: arma::vec X;
    arma::vec R;
    int k;
public:
  LogLike(arma::vec xx_, arma::vec rr_, int kk_) : X(xx_), R(rr_), k(kk_) {}
  
  double operator()(const arma::vec &y) override {
    return loglike(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)) ,Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)),k);
  }
};

double MPS(Rcpp::NumericVector para, 
           Rcpp::NumericVector X, 
           Rcpp::NumericVector R,
           int k){
  
  // Compute objective value
  int m=X.size();
  double mp=sur(X(m-1),para)*cdf(X(0),para)*pow(sur(X(0),para),k*(R(0)+1)-1);
  for(int i=1;i<m;i++){
    mp *= (cdf(X(i),para)-cdf(X(i-1),para))*pow(sur(X(i),para),k*(R(i)+1)-1);
  }
  return mp;
}

// [[Rcpp::export]]
double logMPS(Rcpp::NumericVector para, 
              Rcpp::NumericVector X, 
              Rcpp::NumericVector R,
              int k){
  
  return -std::log(MPS(para,X,R,k));
}

class LogMPS : public Functor {
  private: arma::vec X;
    arma::vec R;
    int k;
public:
  LogMPS(arma::vec xx_, arma::vec rr_, int kk_) : X(xx_), R(rr_), k(kk_) {}
  
  double operator()(const arma::vec &y) override {
    return logMPS(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                  Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                  Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)),k);
  };
};

class LogObjective : public Functor {
  private: arma::vec X;
    arma::vec R;
    int k;
    std::string Type;
public:
  LogObjective(arma::vec xx_, arma::vec rr_, int kk_, std::string type_) : X(xx_), R(rr_), k(kk_), Type(type_){}
  
  double operator()(const arma::vec &y) override {
    if(Type=="LF")
      return loglike(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                     Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                     Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)),k);
    if(Type=="SPF")
      return logMPS(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                    Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                    Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)),k);
  }
};

// Non informative Jeffrey's prior
// [[Rcpp::export]]
double prior(Rcpp::NumericVector para){
  return (1.0/para[1])*std::sqrt(Trigamma(para[0])-1.0/para[0]);
}
/*
// Gradient of prior function
// [[Rcpp::export]]
Rcpp::NumericVector prior_grad(double a1, double b1,double a2,double b2,Rcpp::NumericVector para){
  
  double alp=para[0];
  double lam=para[1];
  
  double G1=psig1*g2/Gamma(alp);
  double G2=g1*g3/Gamma(alp);
  return {G1,G2};
} 
*/
// Posterior distribution 
struct Post{
  double postlike(NumericVector para, 
                  NumericVector X,
                  NumericVector R,
                  int k){
    return like(para,X,R,k)*prior(para);
  }  
  double postmps(NumericVector para, 
                 NumericVector X, 
                 NumericVector R,
                 int k){
    return MPS(para,X,R,k)*prior(para);
  }
};


// This function is used to compute the estimate the MLE and MPS
// [[Rcpp::export]]
Rcpp::List Estim(arma::vec para, 
                 arma::vec X, 
                 arma::vec R, 
                 double L, double U, double P0,int k,
                 std::string type){
  
  Rcpp::NumericVector para_estim(2);
  Rcpp::NumericVector para_se(2);
  Rcpp::NumericVector para_ci(4);
  double obj_val=0.0;
  double cpc_estim;
  double cpc_se;
  Rcpp::NumericVector cpc_gr(2);
  Rcpp::NumericVector cpc_ci(2);
  
  LogObjective rb(X,R,k,type);
  Roptim<LogObjective> opt("Nelder-Mead");
  opt.control.trace = 0;
  opt.set_hessian(true);
  opt.minimize(rb, para);
  Rcpp::NumericVector estim=Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(opt.par()));
  para_estim[0]=estim[0];
  para_estim[1]=estim[1];
  obj_val=opt.value();
  arma::mat inf1=arma::inv(as<arma::mat>(wrap(opt.hessian())));
  Rcpp::NumericMatrix inf=Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(inf1));
  para_se[0] =std::sqrt(inf(0,0));
  para_se[1] =std::sqrt(inf(1,1));
  
  para_ci[0]=std::max(0.0,para_estim[0]-1.96*para_se[0]);
  para_ci[1]=para_estim[0]+1.96*para_se[0];
  para_ci[2]=std::max(0.0,para_estim[1]-1.96*para_se[1]);
  para_ci[3]=para_estim[1]+1.96*para_se[1];
  cpc_estim=cpc_fun(para_estim,L,U,P0);
  cpc_gr=cpc_grad(para_estim,L,U,P0);
  cpc_se=std::sqrt(cpc_gr(0)*cpc_gr(0)*inf(0,0)+cpc_gr(0)*cpc_gr(1)*(inf(0,1)+inf(1,0))+cpc_gr(1)*cpc_gr(1)*inf(1,1));
  cpc_ci[0]=std::max(0.0,cpc_estim-1.96*cpc_se);
  cpc_ci[1]=cpc_estim+1.96*cpc_se;
  return Rcpp::List::create(
    Rcpp::Named("value") = obj_val,
    Rcpp::Named("Estim_para") = para_estim,
    Rcpp::Named("SE_para") = para_se,
    Rcpp::Named("CI_para") = para_ci,
    Rcpp::Named("Estim_Cpc") = cpc_estim,
    Rcpp::Named("SE_Cpc") = cpc_se,
    Rcpp::Named("CI_Cpc") = cpc_ci,
    Rcpp::Named("Hess") = opt.hessian());
};  

// Tierney-Kadane approximation method
class TK_Obj : public Functor {
  private: arma::vec X;
    arma::vec R;
    int k;
    int N;
    std::string Type;
public:
  TK_Obj(arma::vec xx_, arma::vec rr_, int kk_, int NN_,std::string type_) : X(xx_), R(rr_), k(kk_), N(NN_), Type(type_){}
  
  double operator()(const arma::vec &y) override {
    if(Type=="LF")
      return (loglike(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                      Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                      Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)),k)-std::log(prior(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)))))/N;
    if(Type=="SPF")
      return (logMPS(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                     Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                     Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)),k)-std::log(prior(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)))))/N;
  };
};

// Tierney-Kadane approximation method of computing Bayesian of alp for SEl and GEL
class TK_Obj_alp_gel : public Functor {
  private: arma::vec X;
    arma::vec R;
    int k;
    int N;
    double q;
    std::string Type;
public:
  TK_Obj_alp_gel(arma::vec xx_, arma::vec rr_, int kk_, int NN_, double qq_, std::string type_) : X(xx_), R(rr_), k(kk_), N(NN_), q(qq_), Type(type_){}
  
  double operator()(const arma::vec &y) override {
    if(Type=="LF")
      return (loglike(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                      Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                      Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)),k)-std::log(prior(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y))))+q*std::log(y[0]))/N;
    if(Type=="SPF")
      return (logMPS(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                     Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                     Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)),k)-std::log(prior(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y))))+q*std::log(y[0]))/N;
  };
};

// Tierney-Kadane approximation method of computing Bayesian of lam for SEl and GEL
class TK_Obj_lam_gel : public Functor {
  private: arma::vec X;
    arma::vec R;
    int k;
    int N;
    double q;
    std::string Type;
public:
  TK_Obj_lam_gel(arma::vec xx_, arma::vec rr_, int kk_, int NN_, double qq_,std::string type_) : X(xx_), R(rr_), k(kk_), N(NN_), q(qq_), Type(type_){}
  
  double operator()(const arma::vec &y) override {
    if(Type=="LF")
      return (loglike(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                      Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                      Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)),k)-std::log(prior(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y))))+q*std::log(y[1]))/N;
    if(Type=="SPF")
      return (logMPS(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                     Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                     Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)),k)-std::log(prior(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y))))+q*std::log(y[1]))/N;
  };
};

// Tierney-Kadane approximation method of computing Bayesian of cpc for SEL and GEL
class TK_Obj_cpc_gel : public Functor {
  private: arma::vec X;
    arma::vec R;
    int k;
    int N;
    double q;
    double L;
    double U;
    double P0;
    std::string Type;
public:
  TK_Obj_cpc_gel(arma::vec xx_, arma::vec rr_, int kk_, int NN_, double qq_, double LL_, double UU_, double PO_,std::string type_) : X(xx_), R(rr_), k(kk_), N(NN_), q(qq_), L(LL_), U(UU_), P0(PO_),Type(type_){}
  
  double operator()(const arma::vec &y) override {
    if(Type=="LF")
      return (loglike(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                      Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                      Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)),k)-std::log(prior(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y))))+q*std::log(cpc_fun(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),L,U,P0)))/N;
    if(Type=="SPF")
      return (logMPS(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                     Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                     Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)),k)-std::log(prior(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y))))+q*std::log(cpc_fun(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),L,U,P0)))/N;
  };
};


// Tierney-Kadane approximation method of computing Bayesian of alp for LINEX
class TK_Obj_alp_Linex : public Functor {
  private: arma::vec X;
    arma::vec R;
    int k;
    int N;
    double c;
    std::string Type;
public:
  TK_Obj_alp_Linex(arma::vec xx_, arma::vec rr_, int kk_, int NN_, double cc_, std::string type_) : X(xx_), R(rr_), k(kk_), N(NN_), c(cc_), Type(type_){}
  
  double operator()(const arma::vec &y) override {
    if(Type=="LF")
      return (loglike(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                      Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                      Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)),k)-std::log(prior(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y))))+c*y[0])/N;
    if(Type=="SPF")
      return (logMPS(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                     Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                     Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)),k)-std::log(prior(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y))))+c*y[0])/N;
  };
};

// Tierney-Kadane approximation method of computing Bayesian of lam for LINEX
class TK_Obj_lam_Linex : public Functor {
  private: arma::vec X;
    arma::vec R;
    int k;
    int N;
    double c;
    std::string Type;
public:
  TK_Obj_lam_Linex(arma::vec xx_, arma::vec rr_, int kk_, int NN_, double cc_,std::string type_) : X(xx_), R(rr_), k(kk_), N(NN_), c(cc_), Type(type_){}
  
  double operator()(const arma::vec &y) override {
    if(Type=="LF")
      return (loglike(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                      Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                      Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)),k)-std::log(prior(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y))))+c*y[1])/N;
    if(Type=="SPF")
      return (logMPS(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                     Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                     Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)),k)-std::log(prior(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y))))+c*y[1])/N;
  };
};

// Tierney-Kadane approximation method of computing Bayesian of cpc for LINEX
class TK_Obj_cpc_Linex : public Functor {
  private: arma::vec X;
    arma::vec R;
    int k;
    int N;
    double c;
    double L;
    double U;
    double P0;
    std::string Type;
public:
  TK_Obj_cpc_Linex(arma::vec xx_, arma::vec rr_, int kk_, int NN_, double cc_, double LL_, double UU_, double PO_,std::string type_) : X(xx_), R(rr_), k(kk_), N(NN_), c(cc_), L(LL_), U(UU_), P0(PO_),Type(type_){}
  
  double operator()(const arma::vec &y) override {
    if(Type=="LF")
      return (loglike(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                      Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                      Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)),k)-std::log(prior(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y))))+c*cpc_fun(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),L,U,P0))/N;
    if(Type=="SPF")
      return (logMPS(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                     Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                     Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)),k)-std::log(prior(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y))))+c*cpc_fun(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),L,U,P0))/N;
  };
};

// This is the main function of TK method
// [[Rcpp::export]]
Rcpp::List TK(arma::vec para, 
                   arma::vec X, 
                   arma::vec R, 
                   double q, double c,
                   double L, double U, double P0,int k,
                   std::string type){
    
  Rcpp::NumericVector para_est_sel(3);
  Rcpp::NumericVector para_est_gel(3);
  Rcpp::NumericVector para_est_linex(3);
  
  int m=R.size();
  int N=k*(sum(R)+m);
  
  double obj_val,obj_val1,obj_val2,obj_val3;
  arma::mat I, I1, I2, I3;
  double det,det1,det2,det3;
  Rcpp::NumericMatrix para_inf,para_inf1,para_inf2,para_inf3;
  
  // Bayesian estimation under SEL
  TK_Obj rb(X,R,k,N,type);
  Roptim<TK_Obj> opt("Nelder-Mead");
  opt.control.trace = 0;
  opt.set_hessian(true);
  opt.minimize(rb, para);
  obj_val=-opt.value();
  I=arma::inv(as<arma::mat>(wrap(opt.hessian())));
  para_inf=Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(I));
  det =arma::det(I);
  
  TK_Obj_alp_gel rb_alp_sel(X,R,k,N,-1,type);
  Roptim<TK_Obj_alp_gel> opt_alp_sel("Nelder-Mead");  
  opt_alp_sel.control.trace = 0;
  opt_alp_sel.set_hessian(true);
  opt_alp_sel.minimize(rb_alp_sel, para);
  obj_val1=-opt_alp_sel.value();
  I1=arma::inv(as<arma::mat>(wrap(opt_alp_sel.hessian())));
  para_inf1=Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(I1));
  det1=arma::det(I1);
  
  TK_Obj_lam_gel rb_lam_sel(X,R,k,N,-1,type);
  Roptim<TK_Obj_lam_gel> opt_lam_sel("Nelder-Mead");
  opt_lam_sel.control.trace = 0;
  opt_lam_sel.set_hessian(true);
  opt_lam_sel.minimize(rb_lam_sel, para);
  obj_val2=-opt_lam_sel.value();
  I2=arma::inv(as<arma::mat>(wrap(opt_lam_sel.hessian())));
  para_inf2=Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(I2));
  det2=arma::det(I2);
  

  TK_Obj_cpc_gel rb_cpc_sel(X,R,k,N,-1,L,U,P0,type);
  Roptim<TK_Obj_cpc_gel> opt_cpc_sel("Nelder-Mead");
  opt_cpc_sel.control.trace = 0;
  opt_cpc_sel.set_hessian(true);
  opt_cpc_sel.minimize(rb_cpc_sel, para);
  obj_val3=-opt_cpc_sel.value();
  I3=arma::inv(as<arma::mat>(wrap(opt_cpc_sel.hessian())));
  para_inf3=Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(I3));
  det3=arma::det(I3);
  
  para_est_sel={std::sqrt(det1/det)*std::exp(N*(obj_val1-obj_val)),std::sqrt(det2/det)*std::exp(N*(obj_val2-obj_val)),std::sqrt(det3/det)*std::exp(N*(obj_val3-obj_val))};
  // Bayesian estimation under GEL
  
  TK_Obj_alp_gel rb_alp_gel(X,R,k,N,q,type);
  Roptim<TK_Obj_alp_gel> opt_alp_gel("Nelder-Mead");
  opt_alp_gel.control.trace = 0;
  opt_alp_gel.set_hessian(true);
  opt_alp_gel.minimize(rb_alp_gel,para);
  obj_val1=-opt_alp_gel.value();
  I1=arma::inv(as<arma::mat>(wrap(opt_alp_gel.hessian())));
  para_inf1=Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(I1));
  det1=arma::det(I1);
  
  TK_Obj_lam_gel rb_lam_gel(X,R,k,N,q,type);
  Roptim<TK_Obj_lam_gel> opt_lam_gel("Nelder-Mead");
  opt_lam_gel.control.trace = 0;
  opt_lam_gel.set_hessian(true);
  opt_lam_gel.minimize(rb_lam_gel, para);
  obj_val2=-opt_lam_gel.value();
  I2=arma::inv(as<arma::mat>(wrap(opt_lam_gel.hessian())));
  para_inf2=Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(I2));
  det2=arma::det(I2);
  
  TK_Obj_cpc_gel rb_cpc_gel(X,R,k,N,q,L,U,P0,type);
  Roptim<TK_Obj_cpc_gel> opt_cpc_gel("Nelder-Mead");
  opt_cpc_gel.control.trace = 0;
  opt_cpc_gel.set_hessian(true);
  opt_cpc_gel.minimize(rb_cpc_gel, para);
  obj_val3=-opt_cpc_gel.value();
  I3=arma::inv(as<arma::mat>(wrap(opt_cpc_gel.hessian())));
  para_inf3=Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(I3));
  det3=arma::det(I3);
  para_est_gel={pow(std::sqrt(det1/det)*std::exp(N*(obj_val1-obj_val)),-1.0/q),pow(std::sqrt(det2/det)*std::exp(N*(obj_val2-obj_val)),-1.0/q),pow(std::sqrt(det3/det)*std::exp(N*(obj_val3-obj_val)),-1.0/q)};
  
  // Bayesian estimation under Linex
  TK_Obj_alp_Linex rb_alp_linex(X,R,k,N,c,type);
  Roptim<TK_Obj_alp_Linex> opt_alp_linex("Nelder-Mead");
  opt_alp_linex.control.trace = 0;
  opt_alp_linex.set_hessian(true);
  opt_alp_linex.minimize(rb_alp_linex, para);
  obj_val1=-opt_alp_linex.value();
  I1=arma::inv(as<arma::mat>(wrap(opt_alp_linex.hessian())));
  para_inf1=Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(I1));
  det1=arma::det(I1);
  
  TK_Obj_lam_Linex rb_lam_linex(X,R,k,N,c,type);
  Roptim<TK_Obj_lam_Linex> opt_lam_linex("Nelder-Mead");
  opt_lam_linex.control.trace = 0;
  opt_lam_linex.set_hessian(true);
  opt_lam_linex.minimize(rb_lam_linex, para);
  obj_val2=-opt_lam_linex.value();
  I2=arma::inv(as<arma::mat>(wrap(opt_lam_linex.hessian())));
  para_inf2=Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(I2));
  det2=arma::det(I2);
  
  TK_Obj_cpc_Linex rb_cpc_linex(X,R,k,N,c,L,U,P0,type);
  Roptim<TK_Obj_cpc_Linex> opt_cpc_linex("Nelder-Mead");
  opt_cpc_linex.control.trace = 0;
  opt_cpc_linex.set_hessian(true);
  opt_cpc_linex.minimize(rb_cpc_linex, para);
  obj_val3=-opt_cpc_linex.value();
  I3=arma::inv(as<arma::mat>(wrap(opt_cpc_linex.hessian())));
  para_inf3=Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(I3));
  det3=arma::det(I3);

  para_est_linex={(-1.0/c)*std::log(std::sqrt(det1/det)*std::exp(N*(obj_val1-obj_val))),(-1.0/c)*log(std::sqrt(det2/det)*std::exp(N*(obj_val2-obj_val))),(-1.0/c)*std::log(std::sqrt(det3/det)*std::exp(N*(obj_val3-obj_val)))};

   return Rcpp::List::create(
     Rcpp::Named("SEL") = para_est_sel,
     Rcpp::Named("GEL") = para_est_gel,
     Rcpp::Named("LINEX") = para_est_linex);
};  


// Function to compute Bayesian estimation and HPD using MCMC technique v=based on MH samples
// and LF (or PSF) functions
// [[Rcpp::export]]
Rcpp::List MH_sample(std::string type,
                       NumericVector para,NumericVector se,
                       NumericVector R,NumericVector X,
                       int k,double L, double U, double P0, 
                       int MC_size,int MC_burn,
                       double q, double c,int verbose=0, bool display_progress=true){
    
    NumericVector MH_alp_sel(MC_size),MH_lam_sel(MC_size),MH_cpc_sel(MC_size);
    NumericVector MH_alp_gel(MC_size),MH_lam_gel(MC_size),MH_cpc_gel(MC_size);
    NumericVector MH_alp_linex(MC_size),MH_lam_linex(MC_size),MH_cpc_linex(MC_size);
    
    MH_alp_sel(0)=para[0];
    MH_lam_sel(0)=para[1];
    MH_cpc_sel(0)=cpc_fun(para,L,U,P0);
    MH_alp_gel(0)=pow(MH_alp_sel(0),-q);
    MH_lam_gel(0)=pow(MH_lam_sel(0),-q);
    MH_cpc_gel(0)=pow(MH_cpc_sel(0),-q);
    MH_alp_linex(0)=std::exp(-c*MH_alp_sel(0));
    MH_lam_linex(0)=std::exp(-c*MH_lam_sel(0));
    MH_cpc_linex(0)=std::exp(-c*MH_cpc_sel(0));
    
    double alp1,lam1,del,dad;
    int rho=0;
    
    Post post;
    Functor_MH f(&post,&Post::postlike,&Post::postmps);
    
    Progress p(MC_size, display_progress);
    
    int i=1; 
    while(i<=MC_size){
      alp1=std::exp(rnorm(1,std::log(MH_alp_sel[i-1]),se[0])[0]);
      if(alp1==R_NaN) continue;
      if(alp1<=0.5) continue;
      if(alp1>3.0*para[0]) continue;
      lam1=std::exp(rnorm(1,std::log(MH_lam_sel[i-1]),se[1])[0]);
      if(lam1==R_NaN) continue;
      if(lam1>3.0*para[1]) continue;
      double per=f(type,{alp1,lam1},X,R,k)*alp1*lam1/(f(type,{MH_alp_sel[i-1],MH_lam_sel[i-1]},X,R,k)*MH_alp_sel[i-1]*MH_lam_sel[i-1]);
      if(per==R_NaN) continue;
      del=std::min(1.0,per);
      dad=runif(1)[0];
      rho=0;
      if(dad<del) rho=1;
      MH_alp_sel[i] =alp1*rho+(1-rho)*MH_alp_sel[i-1]; 
      MH_lam_sel[i] =lam1*rho+(1-rho)*MH_lam_sel[i-1]; 
      MH_cpc_sel[i]=cpc_fun({MH_alp_sel[i],MH_lam_sel[i]},L,U,P0);
      
      MH_alp_gel[i] =pow(MH_alp_sel[i],-q);
      MH_lam_gel[i] =pow(MH_lam_sel[i],-q);
      MH_cpc_gel[i]=pow(MH_cpc_sel[i],-q);
      MH_alp_linex[i] =std::exp(-c*MH_alp_sel[i]);
      MH_lam_linex[i] =std::exp(-c*MH_lam_sel[i]);
      MH_cpc_linex[i]=std::exp(-c*MH_cpc_sel[i]);
      
      p.increment(); 
      if(verbose>0)
        std::cout<<"Estimate SEL= "<<MH_alp_sel[i]<<" "<<MH_lam_sel[i]<<" "<<MH_cpc_sel[i]<<std::endl; 
      i+=1;
    };
    
    MH_alp_sel.erase(0,MC_burn-1);
    MH_lam_sel.erase(0,MC_burn-1);
    MH_cpc_sel.erase(0,MC_burn-1);
    
    MH_alp_gel.erase(0,MC_burn-1);
    MH_lam_gel.erase(0,MC_burn-1);
    MH_cpc_gel.erase(0,MC_burn-1);
    
    MH_alp_linex.erase(0,MC_burn-1);
    MH_lam_linex.erase(0,MC_burn-1);
    MH_cpc_linex.erase(0,MC_burn-1);
    
    NumericVector est_sel={mean(MH_alp_sel),mean(MH_lam_sel),mean(MH_cpc_sel)};
    NumericVector est_gel={pow(mean(MH_alp_gel),-1/q),pow(mean(MH_lam_gel),-1/q),pow(mean(MH_cpc_gel),-1/q)};
    NumericVector est_linex={(-1/c)*log(mean(MH_alp_linex)),(-1/c)*log(mean(MH_lam_linex)),(-1/c)*log(mean(MH_cpc_linex))};
    NumericVector hpd_alp={my_HPD(my_as_mcmc(MH_alp_sel))[0],my_HPD(my_as_mcmc(MH_alp_sel))[1],my_HPD(my_as_mcmc(MH_alp_sel))[1]-my_HPD(my_as_mcmc(MH_alp_sel))[0]};
    NumericVector hpd_lam={my_HPD(my_as_mcmc(MH_lam_sel))[0],my_HPD(my_as_mcmc(MH_lam_sel))[1],my_HPD(my_as_mcmc(MH_lam_sel))[1]-my_HPD(my_as_mcmc(MH_lam_sel))[0]};
    NumericVector hpd_cpc={my_HPD(my_as_mcmc(MH_cpc_sel))[0],my_HPD(my_as_mcmc(MH_cpc_sel))[1],my_HPD(my_as_mcmc(MH_cpc_sel))[1]-my_HPD(my_as_mcmc(MH_cpc_sel))[0]};
    
    return Rcpp::List::create(
      Rcpp::Named("est_SEL")  =est_sel,
      Rcpp::Named("est_GEL")  =est_gel,
      Rcpp::Named("est_LINEX")=est_linex,
      Rcpp::Named("HPD_alp")=hpd_alp,
      Rcpp::Named("HPD_lam")=hpd_lam,
      Rcpp::Named("HPD_cpc")=hpd_cpc,
      Rcpp::Named("alp_sample")=MH_alp_sel,
      Rcpp::Named("lam_sample")=MH_lam_sel,
      Rcpp::Named("cpc_sample")=MH_cpc_sel
    );
  };
  
  