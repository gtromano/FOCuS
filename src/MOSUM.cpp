#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;
using namespace std;
#include <list>
#include <algorithm>
#include <vector>

/* ------------------------------------------------------------
 
 OFFLINE VERSIONS OF MOSUM
 
 -------------------------------------------------------------- */

/*
// [[Rcpp::export]]
List MOSUM_offline_kirch(NumericVector Y, const double thres, std::vector<int> W) {
  
  long cp = -1;


  std::list<double> stats;
  std::list<double> maxs;


  for (auto t = 0; t <= Y.size(); t++) {

    auto max = 0.0;
    
    //cout << "value of t: " << t << endl;
    
    for (auto& w:W) {
      if (w <= t) {
        auto stat = std::accumulate(Y.begin() + (t - w), Y.begin() + t, 0.0);

        stat = std::abs(stat / sqrt(w));
        if (stat > max)
          max = stat;
      }
    }
    
    maxs.push_back(max);
    
    if (maxs.back() >= thres) {
      cp = t;
      break;
    }
  }
  
  maxs.pop_front();
  
  return List::create(Rcpp::Named("t") = cp,
                      Rcpp::Named("maxs") = maxs);
  
}
*/


// [[Rcpp::export]]
List MOSUM_offline_kirch(NumericVector Y, const double thres, std::vector<int> W) {

  long cp = -1;


  std::list<double> maxs;
  std::vector<double> cusum(W.size()); // container for cusums
  std::vector<double> stat(W.size()); // container for stats

  for (auto t = 0; t < Y.size(); t++) {

    auto max = 0.0;

    //cout << "value of t: " << t << endl;
    auto j = 0; // index for all the W different statistics
    for (auto& w:W) {

      if (w <= t) {
        cusum[j] = cusum[j] + Y[t] - Y[t-w];

        stat[j] = 2 * std::abs(cusum[j]) / sqrt(w);
        //stat[j] = stat[j] * stat[j];

        if (stat[j] > max)
          max = stat[j];
      } else {
        cusum[j] += Y[t];
      }
      j++; // incrementing the j index (for the various window sizes)
    }

    maxs.push_back(max);

    if (maxs.back() >= thres) {
      cp = t;
      break;
    }
  }

  //maxs.pop_front();

  return List::create(Rcpp::Named("t") = cp,
                      Rcpp::Named("maxs") = maxs);

}
