
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

// Computes log(sum(exp(x))) with better precision
double logSumExp(arma::vec x){
  uword i = x.index_max();
  double m = x(i);
  x.shed_row(i);
  double lse = log1p(sum(exp(x-m))) + m;
  return lse;
}

// g function (logit transformation from appendix)
arma::vec g(arma::vec mu, double m, double M){
  arma::vec res = log(mu-m) - log(M-mu);
  return res;
}  
  

// [[Rcpp::export]]
List getThetaC(arma::vec spt, arma::vec f0, arma::vec mu,  
               arma::vec thetaStart, List thetaControl){
 // bool logit  = thetaControl["logit"];
  double eps = thetaControl["eps"];
  int maxiter = thetaControl["maxiter"];
  double maxhalf = thetaControl["maxhalf"];
  double maxtheta = thetaControl["maxtheta"];
 // bool logsumexp = thetaControl["logsumexp"];
  
  int sptN = spt.size();
//  double m = min(spt);
//  double M = max(spt);
  int n = mu.size();

  //initialize values
  arma::vec theta = thetaStart;
  arma::vec thetaOld(n);
  arma::vec bPrimeErrOld(n);
  uvec conv(n, fill::zeros);
  uvec maxedOut(n, fill::zeros);
  
  arma::rowvec oo(sptN, fill::ones);

  mat fUnstd = exp(spt * theta.t());
  fUnstd.each_col() %= f0;  // |spt| x n matrix of tilted f0 values
  rowvec b = oo * fUnstd;  //column sums;
  mat fTilt = fUnstd.each_row() / b;  // normalized

  colvec bPrime = fTilt.t() * spt;  // mean as a function of theta
  colvec bPrime2(n); // variance as a function of theta
  for (int j=0; j<n; j++){       //iterate over the columns of fTilt
    bPrime2(j) = 0;
    for (int i=0; i<sptN; i++){  //iterate over the rows of fTilt
      bPrime2(j) += pow(spt(i) - bPrime(j), 2.0) * fTilt(i,j);
    }
  }
  colvec bPrimeErr = bPrime - mu;  // used to assess convergence
    
  conv = (abs(bPrimeErr) < eps) || (theta==maxtheta && bPrimeErr<0) ||
      (theta==-maxtheta && bPrimeErr>0);
  uvec s = find(conv == 0);
  int iter = 0;
    

 while(s.size() > 0 && iter < maxiter) {
    iter++;
    bPrimeErrOld(s) = bPrimeErr(s);  // used to assess convergence
    
// 1) Update theta
    thetaOld(s) = theta(s);
    colvec thetaS = theta(s) - bPrimeErr(s) / bPrime2(s);
    thetaS(find(thetaS > maxtheta)).fill(maxtheta);
    thetaS(find(thetaS < -maxtheta)).fill(-maxtheta);
    theta(s)= thetaS;
    
// 2) Update fTilt, bPrime, and bPrime2 and take half steps if bPrimeErr not improved
    uvec ss = s;
    int nhalf = 0;
    while(ss.size() > 0 && nhalf < maxhalf) {
// 2a) Update fTilt, bPrime, and bPrime2
      fUnstd.cols(ss) = exp(spt * theta(ss).t());
      fUnstd.each_col(ss) %= f0;  // |spt| x n matrix of tilted f0 values
      b = oo * fUnstd.cols(ss);  //column sums;
      mat tmp = fUnstd.cols(ss);
      fTilt.cols(ss) = tmp.each_row() / b;  // normalized
      
      bPrime(ss) = fTilt.cols(ss).t() * spt;  // mean as a function of theta
      for (uword j=0; j<ss.size(); j++){       //iterate over the columns of fTilt[,ss]
        bPrime2(ss(j)) = 0;
        for (int i=0; i<sptN; i++){  //iterate over the rows of fTilt[,ss]
          bPrime2(ss(j)) += pow(spt(i) - bPrime(ss(j)), 2.0) * fTilt(i,ss(j));
        }
      }
      bPrimeErr(ss) = bPrime(ss) - mu(ss); 
      
// 2b) Take half steps if necessary
      ss = ss(find(abs(bPrimeErr(ss)) > abs(bPrimeErrOld(ss))));
      if (ss.size() > 0){
        nhalf++;
      } 
      theta(ss) = (theta(ss) + thetaOld(ss)) / 2;
    }
// If maximum half steps are exceeded, set theta to previous value
    maxedOut(ss).fill(1);
    theta(ss) = thetaOld(ss);
    
// 3) Check convergence
    conv(s) = (abs(bPrimeErr(s)) < eps) ||  (theta(s)==maxtheta && bPrimeErr(s) < 0) ||
      (theta(s)==-maxtheta && bPrimeErr(s) > 0);
    s = s(find(conv(s) == 0 && maxedOut(s) == 0));
 }
      
  List res;
  res["theta"] = theta;
  res["fTilt"] = fTilt;
  res["bPrime"] = bPrime;
  res["bPrime2"] = bPrime2;
  res["conv"] = conv;
  res["iter"] = iter;
  return res;
}

