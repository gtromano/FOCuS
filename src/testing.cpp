#include <Rcpp.h>
#include "quadratic.h"
using namespace Rcpp;

// [[Rcpp::export]]
void test() {
  

  // std::list<int> mylist1 = {2, 3, 4, 7, 5, 6};
  std::list<int> mylist1 = {1, 2, 3, 4, 5, 6};
  std::list<int> mylist2 = {2, 3, 4, 6};

  
  std::unique_ptr<int> foo;
  
  std::cout << static_cast<bool>(foo) << std::endl;
  
  

  // auto testf = [](auto& n) {n++;};
  // 
  // for (auto i = mylist1.begin(); i != mylist1.end();) {
  //   
  //   testf(*i);
  //   
  //   if (*i > 5)
  //     i = mylist1.erase(i);
  //   else
  //     ++i;
  // }
  // 
  // 
  // for (std::list<int>::iterator it=mylist1.begin() ; it != mylist1.end(); ++it)
  //   std::cout << ' ' << *it;
    
  return;
}
