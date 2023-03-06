#include <Rcpp.h>
#include "focus_new_implementation.h"
using namespace Rcpp;

long get_changepoint (const Cost& Q, const CUSUM& cs, const double& theta0) {
  
  auto max = -INFINITY;
  long cp = 0;
  
  // std::cout << m0val << std::endl;
  
  for (auto i = 0; i <= Q.k; i++) {
    
    auto q_max = get_max(Q.ps[i], cs, theta0);
    
    if (q_max > max) {
      max = q_max;
      cp = Q.ps[i]->tau;
    }
    
  }
  // std::cout << std::endl;
  
  return cp;
}


Rcpp::List unpack_cost(const Cost& cost) {
  Rcpp::List pieces;
  for (int i = 0; i < cost.k; i++) {
    const auto& p = cost.ps[i];
    pieces.push_back(Rcpp::NumericVector::create(Rcpp::Named("St") = p->St, Rcpp::Named("tau") = p->tau, Rcpp::Named("m0") = p->m0, Rcpp::Named("Mdiff") = p->Mdiff));
  }
  return pieces;
}

// [[Rcpp::export(.focus_offline_new_imp)]]
List focus_offline_new_imp (NumericVector Z, double threshold, double theta0) {
  
  auto Y = clone(Z);
  
  // here we define the function to initialize a new piece
  std::function<std::unique_ptr<Piece>(double, int, double)> newP;
  bool adp_max_check = 0;
  auto cp = -1;
  
  newP = [](double St, int tau, double m0){
    std::unique_ptr<Piece> p = std::make_unique<PieceGau>();
    p->St = St;
    p->tau = tau;
    p->m0 = m0;
    
    return p;
  };
  
  // new init with constructor
  Info info(newP);
  
  // gaussian cost, in case pre-change mean is known, center the observations on zero
  if (!std::isnan(theta0)) {
    Y = Y - theta0;
    theta0 = 0;
  }
  
  std::vector<double> stat(Y.size());

  
  for (auto& y:Y) {
    info.update(y, newP, threshold, theta0, adp_max_check);
    //stat.push_back(std::max(info.Ql.opt, info.Qr.opt));
    stat[info.cs.n - 1] = std::max(info.Ql.opt, info.Qr.opt);

    
    if (stat[info.cs.n - 1] >= threshold) {
      if (info.Ql.opt > info.Qr.opt) {
        cp = get_changepoint(info.Ql, info.cs, theta0);
      } else {
        cp = get_changepoint(info.Qr, info.cs, theta0);
      }
      break;
    }
      
  }
  
  auto cs = List::create(Rcpp::Named("Sn") = info.cs.Sn,
                         Rcpp::Named("n") = info.cs.n);
  
  return List::create(Rcpp::Named("maxs") = stat,
                      Rcpp::Named("t") = info.cs.n,
                      Rcpp::Named("changepoint") = cp,
                      Rcpp::Named("Ql") = unpack_cost(info.Ql),
                      Rcpp::Named("Qr") = unpack_cost(info.Qr),
                      Rcpp::Named("cs") = cs);
  
}