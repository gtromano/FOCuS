#include <Rcpp.h>
using namespace Rcpp;
#include "FOCuS.h"

// This is the script for exporting the C++ functions to R.


// [[Rcpp::export]]
int FOCuS_offline(NumericVector Y, double thres) {
  long t = 0;
  long cp = -1;
  
  Quadratic Q0, q1;
  Info info = {Q0, {q1}, 0};
  
  for (auto& y:Y) {
    t += 1;
    info = FOCuS_step(std::move(info), y);
    
    if (info.global_max >= thres) {
      cp = t;
      break;
    }
  }
  std::cout << info.Q1.size() << std::endl;
  std::cout << info.global_max << std::endl;
  
  // last Q1
  std::cout << "last Q1"<< std::endl;
  for (auto& q:info.Q1)
    print(q);
  std::cout << std::endl;
  
  
  return cp;
}



/*** R
set.seed(42)
Y <- rnorm(10)
FOCuS_offline(Y, 10)
*/
