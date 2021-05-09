#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;
#include <list>
#include <algorithm>


/* ------------------------------------------------------------
 
 OFFLINE VERSIONS OF YU CUSUM
 
 -------------------------------------------------------------- */

// [[Rcpp::export]]
List YuCUSUM_offline(NumericVector Y, const double thres) {
  
  long t = 0;
  long cp = -1;
  
  std::list<double> cusums;
  std::list<double> max_at_time_t;
  

  for (auto& y:Y) {
    if (cusums.size() == 0)
      cusums.push_back(y);
    else
      cusums.push_back(cusums.back() + y);
    
    auto maximum = -INFINITY;
    t++;
    auto s = 0;
    
    
    for (auto par_sum:cusums) {
      s++;
      
      double test;
      
      if (t == s) 
        test = 0;
      else 
        test = sqrt((t)/(s * (t - s))) * std::abs(((double) s / (double) t) * cusums.back() - par_sum);
      

      if (test > maximum)
        maximum = test;
    }
    
    max_at_time_t.push_back(maximum);

    if (max_at_time_t.back() >= thres) {
      cp = t;
      break;
    }
  }
  
  
  // std::cout << info.Q1.size() << std::endl;
  // std::cout << info.global_max << std::endl;
  
  return List::create(Rcpp::Named("t") = cp,
                      Rcpp::Named("maxs") = max_at_time_t, 
                      Rcpp::Named("cusums") = cusums);
}