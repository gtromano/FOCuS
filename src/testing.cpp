#include <Rcpp.h>
#include "quadratic.h"
using namespace Rcpp;

// [[Rcpp::export]]
void test() {
  
  Quadratic ex1, ex2, ex3;
  
  ex1.a = 13;
  ex1.c = 2;
  ex2.c = 5;
  
  double Q0_max = 0;
  std::cout << (Q0_max) <<  std::endl;
  
  for(const auto& i:ex1.ints) {
    auto m = std::get<0>(get_minimum(ex1, i));
    if (m > Q0_max)
      Q0_max = m;
  }
  std::cout << "Max is " << Q0_max <<  std::endl;
  
  
  
  auto res = get_intersections(ex1, ex2);

  std::cout << !std::isnan(std::get<0>(res))  << " " << !std::isnan(std::get<1>(res)) << std::endl;
  
  get_min_of_two_quadratics(ex1, ex2);
  
  print(ex1);
  std::cout << std::endl;
  print(ex2);
  
  ex3.a = 1;
  ex3.b = 20;
  ex3.c = 3;
  
  get_min_of_two_quadratics(ex3, ex2);
  
  print(ex3);
  std::cout << std::endl;
  print(ex2);
  
  
  
  

  /*
  // std::list<int> mylist1 = {2, 3, 4, 7, 5, 6};
  std::list<int> mylist1 = {7, 3};
  std::list<int> mylist2 = {2, 3, 4, 6};


  for (std::list<int>::iterator it=mylist1.begin() ; it != mylist1.end(); ++it)
    std::cout << ' ' << *it;

  std::cout << std::endl;

  for (std::list<int>::iterator it1=mylist1.begin() ; it1 != mylist1.end(); ++it1) {
    for (std::list<int>::iterator it2=mylist2.begin() ; it2 != mylist2.end(); ++it2) {
      if ((*it1) == (*it2))
        it1 = mylist1.erase(it1);
    }
  }



  for (std::list<int>::iterator it=mylist1.begin() ; it != mylist1.end(); ++it)
    std::cout << ' ' << *it;
  */
  
  std::list<int> mylist1 = {2, 3, 4, 7, 5, 6};
  
  auto testf = [](auto& n) {n++;};
  
  for (auto i = mylist1.begin(); i != mylist1.end();) {
    
    testf(*i);
    
    if (*i > 5)
      i = mylist1.erase(i);
    else
      ++i;
  }
  
  
  for (std::list<int>::iterator it=mylist1.begin() ; it != mylist1.end(); ++it)
    std::cout << ' ' << *it;
    
  return;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
test()
*/
