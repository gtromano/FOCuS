#include "quadratic.h"

Interval I(const double& l, const double& u) {
  Interval I = {l, u};
  return I;
}

void print(const Quadratic& q) {
  std::cout << q.a << "x^2 + " << q.b << "x + " << q.c << std::endl;
  std::cout << "quadratic domain:"<< std::endl;
  for(auto i:q.ints)
    std::cout << "(" << i.l << ", " << i.u << ")   ";
   std::cout << std::endl;
}

auto evaluate_quadratic(const Quadratic& q, const double& x) {
  return q.a * (x * x) + q.b * x + q.c;
}


// to use with std::move() in order to speed up times
void invert_quadratic(Quadratic& q) {
  q.a = - q.a;
  q.b = - q.b;
  q.c = - q.c;
}


// finds the minimum of a quadratic on an interval
std::tuple<double, double> get_minimum(const Quadratic& q, const Interval& i){

  if (q.a == 0) {
    return std::make_tuple(q.c, i.l);
  }
  
  auto at = - q.b / (2.0 * q.a);
  if (at <= i.l) {
    at = i.l;
  } else if (at >= i.u) {
    at = i.u;
  }

  
  auto minim = evaluate_quadratic(q, at);

  return std::make_tuple(minim, at);
}


// check if a value is in range of an interval
bool inRange(const double& x, const Interval& i) { return i.l <= x && i.u >= x; }


// finds the intersection of two quadratics
// note: should return both the intersections, or none in case the two quadratic
// has no intersections at all
std::tuple<double, double> get_intersections (const Quadratic& q1, const Quadratic& q2) {
  auto a = q1.a - q2.a;
  auto b = q1.b - q2.b;
  auto c = q1.c - q2.c;
  
  if (c == 0) { // to avoid the sqrt
    auto inter = - b / a;
    return std::make_tuple(std::min(0.0, inter), std::max(0.0, inter));
  }

  auto sqrt_z = sqrt((b * b) - (4 * a * c));
  return std::make_tuple((- b - sqrt_z) / (2 * a), (- b + sqrt_z) / (2 * a));
}


// careful, this function takes by reference and actually modifies the intervals
// by finding on which domain the quadratics are smallest
void get_min_of_two_quadratics (Quadratic& q1, Quadratic& q2) {
  
  auto inters = get_intersections(q1, q2);
  auto n_inters = !std::isnan(std::get<0>(inters)) + !std::isnan(std::get<0>(inters)); // sum of non null intersections
  if (n_inters == 0 || (std::get<0>(inters) == std::get<1>(inters))) {
    
    if (std::get<0>(get_minimum(q1, q1.ints.front())) >= std::get<0>(get_minimum(q2, q1.ints.front())))
      q1.ints = {}; // deleting the quadratic since the line is better
    else
      q2.ints = {}; // deleting the line since the quad is better
    return;
  }
  
  // check for data races. Shouldn't be an issue since iterator skips to the end if list empty
  for (auto& i1:q1.ints) {
    for (auto& i2:q2.ints) {
      // run only if the domain of the first is contained in the second
      if ((i2.l <= i1.l) && (i2.u >= i1.u)) {
        // check whether the left or right conditions are in range
        auto lCond = !std::isnan(std::get<0>(inters)) &&
          std::get<0>(inters) != i1.u &&
          inRange(std::get<0>(inters), i1) &&
          inRange(std::get<0>(inters), i2);
        auto rCond = !std::isnan(std::get<1>(inters)) &&
          std::get<1>(inters) != i1.l &&
          inRange(std::get<1>(inters), i1) &&
          inRange(std::get<1>(inters), i2);
        
        if ((lCond + rCond) == 2) {
          // both in range, we cut the coefficient of the line in two parts
          q2.ints.push_back(I(std::get<1>(inters), i2.u));
          i2.u = std::get<0>(inters);
          i1 = I(std::get<0>(inters), std::get<1>(inters));
        } else if (lCond) {
          // left in range, we cut first the line and then the quad
          if (i1.u < i2.u)
            q2.ints.push_back(I(i1.u, i2.u));
          i2.u = std::get<0>(inters);
          i1.l = std::get<0>(inters);
        } else if (rCond) {
          // right in range, we cut first the quad and then the line
          if (i2.l < i1.l)
            q2.ints.push_back(I(i2.l, i1.l));
          i2.l = std::get<1>(inters);
          i1.u = std::get<1>(inters);
        } else {
          // here we don't have intersections and we have to figure out 
          // whether the line is highest, or the quadratic
          auto interval = i1;
          //if (std::get<0>(get_minimum(q1, interval)) >= std::get<0>(get_minimum(q2, interval))) {
          if (std::get<0>(get_minimum(q1, interval)) >= 0.0) { // here this is a particular case that only works in this scenario ***********************************
            // std::cout<<"******* erasing interval **********"<<std::endl;
            //i1 = q1.ints.erase(i1); // if line is highest we prune the quadratic
            i1 = I(std::nanf(""), std::nanf(""));
          } else {
            // otherwise we have to trim the line
            if (i2.l == interval.l)
              i2.l = std::move(interval.u);
            else if (i2.u == interval.u)
              i2.u = std::move(interval.l);
            else {
              q2.ints.push_back(std::move(I(interval.u, i2.u)));
              i2.u = std::move(interval.l);
            }
          }
        }
      } // end domain condition
    } // end q2 for
  } // end q1 for
  
  q1.ints.remove_if([](auto& i){return std::isnan(i.l);});
  q2.ints.remove_if([](auto& i){return i.l == i.u;});
  
} // end function



// this function takes in a list of quadratics (cost) and a new line (newq)
// and try to return the updated cost with a new line added in it
// to be used with std::move() after testing
void get_min_of_cost(std::list<Quadratic>& cost, Quadratic& newq) {
  for (auto& q:cost)
    get_min_of_two_quadratics(q, newq);
  
  cost.remove_if([](auto& q){
    return q.ints.size() == 0;
  });
    
  cost.push_back(std::move(newq));
}


// add std::move()
// this is simply get min_of_cost but with the inverted coefficients
std::list<Quadratic> get_max_of_cost(std::list<Quadratic> cost, Quadratic newq) {

  
  for (auto& q:cost)
    invert_quadratic(q);
  
  get_min_of_cost(cost, newq); // add std::move
  
  for (auto& q:cost)
    invert_quadratic(q);

  
  return cost;
}


// MELKMANS LIKE ALGORITHM ****   experimental stuff   *****
/*
void get_max_of_cost_melk_right(std::list<Quadratic>& cost, Quadratic newq) {
  
  bool found = false;
  
  for(auto it = cost.begin(); it != cost.end(); ++it) {
    if (found) {
      (*it).ints = {}; // if we have found an intersection already we prune the remaining quadratics
    } else {
      auto nx = std::next(it, 1);
      
      //std::cout << "Break 4 " << (nx == cost.end()) << std::endl;
      if (nx == cost.end()) {
        auto inter = - (*it).b / (*it).a; // in this case this only happens to the last one
        if (inter > 0)
          (*it).ints.front().u = newq.ints.front().l = inter;
        else
          (*it).ints = {};
        found = true;
      } else {
        auto a = (*it).a - (*nx).a;
        auto b = (*it).b - (*nx).b;
        auto inter = - b / a;
        
        if (evaluate_quadratic((*it), inter) > 0.0 && inter > 0.0) {
          (*it).ints.front().u = (*nx).ints.front().l = inter; // in this case we intercept the quadratic
        } else {
          inter = - (*it).b / (*it).a; // in this case we have intercepted the line the line
          if (inter > 0)
            (*it).ints.front().u = newq.ints.front().l = inter;
          else
            (*it).ints = {};
          found = true;
        } //end else 
      } // end else 
    } // end else
  } // end for
  
  cost.remove_if([](auto& q){
    return q.ints.size() == 0;
  });
  cost.push_back(newq);
  
}


void get_max_of_cost_melk_left(std::list<Quadratic>& cost, Quadratic newq) {
  
  bool found = false;
  
  for(auto it = cost.begin(); it != cost.end(); ++it) {
    
    if (found) {
      (*it).ints = {}; // if we have found an intersection already we prune the remaining quadratics
    } else {
      auto nx = std::next(it, 1);
      
      if (nx == cost.end()) {
        //print(newq);
        //print((*it));
        auto inter = - (*it).b / (*it).a; // in this case this only happens to the last one
        if (inter < 0)
          (*it).ints.front().l = newq.ints.front().u = inter;
        else
          (*it).ints = {};
        found = true;
      } else {
        auto a = (*it).a - (*nx).a;
        auto b = (*it).b - (*nx).b;
        auto inter = - b / a;
        
        if (evaluate_quadratic((*it), inter) > 0.0 && inter < 0.0) {
          (*it).ints.front().l = (*nx).ints.front().u = inter; // in this case we intercept the quadratic
        } else {
          inter = - (*it).b / (*it).a; // in this case we have intercepted the line
          if (inter < 0)
            (*it).ints.front().l = newq.ints.front().u = inter;
          else
            (*it).ints = {};
          found = true;
        } //end else 
      } // end else 
    } // end else
    
  } // end for
  
  cost.remove_if([](auto& q){
    return q.ints.size() == 0;
  });
  cost.push_back(newq);
}
*/


// so this is the Melkmans algo using only the linked list features. 
// here we actually break from the for cycle as soon as an intersection is detected
// the pruning is done with the pop_back function, that only touches the back of the 
// list and it is constant in time. 
 
void get_max_of_cost_melk_right(std::list<Quadratic>& cost, Quadratic newq) {
  
  auto n_to_retain = 0; // number of quadratics to retain (needed for the pruning later)
  
  for(auto it = cost.begin(); it != cost.end(); ++it) {
    
    auto nx = std::next(it, 1); // this is the pointer to the next quadratic
    
    if (nx == cost.end()) {
      // in this case we are at the end of the list
      auto inter = - (*it).b / (*it).a; 
      if (inter > 0) {
        (*it).ints.front().u = newq.ints.front().l = inter;
        n_to_retain++;
      }
      break;
      
    } else {
      auto a = (*it).a - (*nx).a;
      auto b = (*it).b - (*nx).b;
      auto inter = - b / a;
      
      if (evaluate_quadratic((*it), inter) > 0.0 && inter > 0.0) {
        (*it).ints.front().u = (*nx).ints.front().l = inter; // in this case we intercept the quadratic
        n_to_retain++;
      } else {
        // in this case we do not have any valid intersection so we need to cut with the line (newq)
        inter = - (*it).b / (*it).a;
        if (inter > 0) {
          (*it).ints.front().u = newq.ints.front().l = inter;
          n_to_retain++;
        }
        break; // and we break
      } //end else 
    } // end else 
  } // end for
  
  //std::cout << "\nMy list is long: " <<  cost.size() << std::endl;
  auto to_prune =  (cost.size() - n_to_retain);
  for (int a = 0; a < to_prune; a++) {
    cost.pop_back(); // please dear god have mercy on our souls
  }
  
  cost.push_back(newq);
}


void get_max_of_cost_melk_left(std::list<Quadratic>& cost, Quadratic newq) {
  
  auto n_to_retain = 0;
  
  for(auto it = cost.begin(); it != cost.end(); ++it) {
    
    auto nx = std::next(it, 1);
    if (nx == cost.end()) {
      
      auto inter = - (*it).b / (*it).a;
      if (inter < 0) {
        (*it).ints.front().l = newq.ints.front().u = inter;
        n_to_retain++;
      }
      break;
      
    } else {
      auto a = (*it).a - (*nx).a;
      auto b = (*it).b - (*nx).b;
      auto inter = - b / a;
      
      if (evaluate_quadratic((*it), inter) > 0.0 && inter < 0.0) {
        (*it).ints.front().l = (*nx).ints.front().u = inter;
        n_to_retain++;
      } else {
        inter = - (*it).b / (*it).a; // in this case we have intercepted the line
        if (inter < 0) {
          (*it).ints.front().l = newq.ints.front().u = inter;
          n_to_retain++;
        }
        break;
      } //end else 
    } // end else 
  } // end for
  
  //std::cout << "\nMy list is long: " <<  cost.size() << std::endl;
  auto to_prune =  (cost.size() - n_to_retain);
  for (int a = 0; a < to_prune; a++) {
    cost.pop_back(); // please dear god have mercy on our souls
  }
  
  cost.push_back(newq);
}



// ********* end of the experimental stuff *************

// function for the approximation
void approximation_grid (std::list<Quadratic>& Q, const std::list<double>& grid) {
  
  if (int(Q.size()) - int(grid.size()) < 1) {return;}
  
  //std::cout << "running" << std::endl;
  
  for (auto& q:Q) { // iterate for each quadratic
    for (auto& i:q.ints) { // and for each interval
      
      // check if we have to retain the quadratic (i.e. skip the process)
      auto point_in_range = false;
      for (auto p:grid) {
        if (inRange(p, i)) {
          point_in_range = true;
          break;
        } 
      }
      if (point_in_range)
        continue; // go to next iteration since we do not have to do anything
      
      auto found = false;
      for (auto& left:Q) { // left here search for the left interval

        if (found)
          break;
        
        if (left.ints.front().u == i.l) { // if the upper of left matches the lower of i, we found it
          
          // and now we look for the right
          for (auto& right:Q) {
            for (auto& ri:right.ints) {
              
              if (found)
                break;
              
              if (i.u == ri.l) { // if the upper of i matches the lower of right
                
                // here we do the magic :)
                auto inters = get_intersections(left, right); // we find the intersections of left and right
                
                if (std::get<1>(inters) > left.ints.front().l && std::get<1>(inters) < ri.u)
                  left.ints.front().u = ri.l = std::get<1>(inters); // if the right intersection is in between you get that
                else
                  left.ints.front().u = ri.l = std::get<0>(inters); // otherwise you put the left
                
                // all is left is to mark as to delete the interval in the middle, so our i
                i = I(std::nanf(""), std::nanf(""));
                
              } // end if for the right match
            } // end right quad interval iteration
          } // end right quad iteration
        } // end if for the left match
      } // end left quad iteration
    } // end interval iteration
  } // end quad iteration
  
  Q.remove_if([](auto& q){
    for(auto& i:q.ints) {
      if(std::isnan(i.l))
        return true;
    }
    return false;
  });

} // end function
