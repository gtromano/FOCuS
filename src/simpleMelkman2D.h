#include <Rcpp.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <list>

using namespace Rcpp;

//' Function to update list of intervals (up or down) 
//'
//' @param n number of datapoints
//' @param sumAllX sum of all x_i up to n (For down change should give -sumAllX)   
//' @param tau vector of candidate changes
//' @param bound vector of interval bound m = (\mu1+\mu2)/2
//' @param sum vector of sum of x_i up to tau	 
//' @return void
void oneStepUpdate(int n, double sumAllX,
		std::vector<int> &tau,
		std::vector<double> &bound,
		std::vector<double> &sum);


//' Function to get the best cost given a set of candidate and there sum
//' 
//' @param n number of datapoints
//' @param sumAllX sum of all x_i up to n (For down change should give -sumAllX)   
//' @param tau vector of candidate changes
//' @param sum vector of sum of x_i up to tau	  
//' @return double
double getBestCost(	
		int n, double sumAllX,
		std::vector<int> &tau, 
		std::vector<double> &sum);

//' Function to run a Melkman-like algorithm for unknown first and second segment mean
//'
//' @param x datapoints
//' @param onlyPrune if TRUE only update intervals and does not compute the best cost at each step
//' @param exportInR if TRUE results (tau, bound and sum) are exported in R
//' @return Product of v1 and v2
// [[Rcpp::export(.simpleMelkman)]]
List simpleMelkman(NumericVector x, bool onlyPrune, bool exportInR);
