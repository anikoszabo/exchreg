#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

// Computes log(sum(exp(x))) with better precision
double logSumExp(arma::vec x){
  NumericVector xa = wrap(x);
  int i = which_max(xa);
  double m = xa[i];
  xa.erase(i);
  double lse = log1p(sum(exp(xa-m))) + m;
  return lse;
}

// g function (logit transformation from appendix)
arma::vec g(arma::vec mu, double m, double M){
  arma::vec res = log(mu-m) - log(M-mu);
  return res;
}  
  

// [[Rcpp::export]]
List getThetaC(arma::vec spt, arma::vec f0, arma::vec mu, arma::vec weights, 
               arma::vec thetaStart, List thetaControl){
 // bool logit  = thetaControl["logit"];
//  double eps = thetaControl["eps"];
//  int maxiter = thetaControl["maxiter"];
//  double maxhalf = thetaControl["maxhalf"];
//  double maxtheta = thetaControl["maxtheta"];
 // bool logsumexp = thetaControl["logsumexp"];
  
  int sptN = spt.size();
//  double m = min(spt);
//  double M = max(spt);
  int n = mu.size();

  //initialize values
  arma::vec theta = thetaStart;
  arma::vec thetaOld(n);
  arma::vec bPrimeErrOld(n);
  LogicalVector conv(n);
  conv.fill(FALSE);
  LogicalVector maxedOut(n);
  maxedOut.fill(FALSE);
  
  arma::rowvec oo = ones(sptN);
    
  mat fUnstd = f0 * exp(spt * theta.t()); // |spt| x n matrix of tilted f0 values
  rowvec b = oo * fUnstd;  //column sums;
  mat fTilt = fUnstd.each_row() / b;  // normalized
  
  colvec bPrime = fTilt.t() * spt;  // mean as a function of theta
  colvec bPrime2(sptN); // variance as a function of theta
  for (int j=0; j<sptN; j++){
    bPrime2(j) = 0;
    for (int i=0; i<n; i++){
      bPrime2(j) += pow(spt(j) - bPrime(i), 2.0) * fTilt(i,j);
    }
  }
  colvec bPrimeErr = bPrime - mu;  // used to assess convergence
    
  int iter = 0;    
  List res;
  res["theta"] = theta;
  res["fTilt"] = fTilt;
  res["bPrime"] = bPrime;
  res["bPrime2"] = bPrime2;
  res["conv"] = conv;
  res["iter"] = iter;
  return res;
}
