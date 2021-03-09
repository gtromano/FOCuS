#include <Rcpp.h>
#include <cmath>
#include <iostream>
#include <vector>

using namespace Rcpp;


void oneStepUpdate(int n, double sumAllX,
		std::vector<int> &tau,
		std::vector<double> &bound,
		std::vector<double> &sum);


double getBestCost(	
		int n, double sumAllX,
		std::vector<int> &tau, 
		std::vector<double> &sum);

// [[Rcpp::export]]
List simpleMelkman(NumericVector x, bool onlyPrune, bool exportInR){
  List out;
  /*.................................................................................*/
  /* ((INITIALIZE */
  /*.................................................................................*/
  int initSize = 10*log(x.size()) + 10;
  // TAUS
  std::vector<int> tauUp(0);
  std::vector<int> tauDw(0);
  tauUp.reserve(initSize); 
  tauDw.reserve(initSize);

  // BOUND INTERVALS
  std::vector<double> boundUp(0); 
  std::vector<double> boundDw(0); 
  boundUp.reserve(initSize);
  boundDw.reserve(initSize);

  // SUMS
  std::vector<double> sumUp(0);
  std::vector<double> sumDw(0);
  sumUp.reserve(initSize);
  sumDw.reserve(initSize);
  /*.................................................................................*/
  /* INITIALIZE)) */
  /*.................................................................................*/

  /*.................................................................................*/
  /* ((FIRST DATA POINT */
  /*.................................................................................*/
  tauUp.push_back(0);
  tauDw.push_back(0);

  sumUp.push_back(x[0]);
  sumDw.push_back(-x[0]);
  
  boundUp.push_back(x[0]);
  boundDw.push_back(-x[0]);
  
  /*.................................................................................*/
  /* FIRST DATA POINT)) */
  /*.................................................................................*/

  /*.................................................................................*/
  /* ((LOOP OVER OTHER POINTS */
  /*.................................................................................*/
  double sumAllX = x[0];
  double minCurrent = 10;
  if(onlyPrune){ // does not compute best at each step
    for(int i=1; i < x.size(); i++){
      sumAllX = sumAllX + x[i];
      oneStepUpdate(i, sumAllX, tauUp, boundUp, sumUp);
      oneStepUpdate(i, -sumAllX, tauDw, boundDw, sumDw);
    }
  } else {  // recover best at each step
    
    for(int i=1; i < x.size(); i++){
      sumAllX = sumAllX + x[i];
      oneStepUpdate(i+1, sumAllX, tauUp, boundUp, sumUp);
      oneStepUpdate(i+1, -sumAllX, tauDw, boundDw, sumDw);

      minCurrent = std::min(getBestCost(i+1, sumAllX, tauUp, sumUp), getBestCost(i+1, -sumAllX, tauDw, sumDw));
    }

  }

  /*.................................................................................*/
  /* LOOP OVER OTHER POINTS)) */
  /*.................................................................................*/
  
  // OUTPUT IN R 
  if(exportInR){
    out["tauUp"] = tauUp;
    out["tauDw"] = tauDw;
    out["boundUp"] = boundUp;
    out["boundDw"] = boundDw;
    out["sumUp"] = sumUp;
    out["sumDw"] = sumDw;
    out["minCurrent"] = minCurrent;
  }
  return(out);
}


// One-step interval update of the Melkman 2d algo for Up-changes
void oneStepUpdate(	
		int n, double sumAllX,
		std::vector<int> &tau,
		std::vector<double> &bound,
		std::vector<double> &sum){

   int cur = tau.size();
   double muNew = (sumAllX - sum[cur-1]) / (n - tau[cur-1] - 1);
  
   // loop
    while(cur > 0 & (muNew <= bound[cur-1]) ){
      cur--;
      muNew = (sumAllX - sum[cur-1]) / (n - tau[cur-1] - 1);
      
   } 

   // update
   if(cur == 0) { // All have been removed
     muNew = sumAllX / n;
   }
 
   // Resize if necessary
   tau.resize(cur);
   bound.resize(cur);
   sum.resize(cur);
   
   // insert
   tau.push_back(n-1);
   bound.push_back(muNew);
   sum.push_back(sumAllX);
}

// One-step cost calculation
double getBestCost(	
		int n, double sumAllX,
		std::vector<int> &tau, 
		std::vector<double> &sum){


   int cur = tau.size();
   double sumAf=0;
   double minCost = -sum[cur-1]*sum[cur-1]/n;
   //std::cout << "Init: " << sum[cur-1] << "^2/" << n <<", " << sumAf << "^2/" << (n - tau[cur-1]-1) << " : " << minCost << std::endl;
   double tmp = 0;
   // loop
    while(cur > 0 ){
      cur--;
      sumAf = sumAllX - sum[cur-1];
      tmp = -sum[cur-1]*sum[cur-1]/(tau[cur-1]+1) - sumAf*sumAf / (n - tau[cur-1] -1);
     // std::cout << "Next: " << sum[cur-1] << "^2/" << tau[cur-1]+1 <<", " << sumAf << "^2/" << (n - tau[cur-1] - 1) << " : " << tmp << std::endl;
      if(tmp < minCost) minCost = tmp;
   } 
   return(minCost);

}


