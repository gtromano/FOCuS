#include "biweight_loss.h"
//#include "quadratic.h"

void update_cost_biweight(std::list<Quadratic>& Q, const double& x, const double& K, const double& m0 = 0.0) {
  auto k = sqrt(2 * K);
  //auto scaling = std::max(-.5 * (x * x), - k * k);
  auto scaling = std::max(-.5 * (m0 * m0 - 2 * m0 * x + x * x), - k * k);
  
  Quadratic new_q;
  new_q.a = - .5; new_q.b = x; new_q.c = -.5 * (x * x) - scaling;
  
  // this is some experimental shenanigans (handling of the entirety of the update in update_cost_biweight)
  if (std::isinf(K)) {
    // if K is infinity then we simply update with the quadratic
    for(auto& q:Q) {
      q.a += new_q.a;
      q.b += new_q.b;
      q.c += new_q.c;
    }
    return; // this is the end here you have the update and that's all folks no intersections are ever expected
  }
  // end of some experimental crap
  
  Quadratic new_l;
  new_l.c = - (k * k) - scaling;
  
  auto inters = get_intersections(new_l, new_q); // those are the intersection for the updated cost
  
  auto left = std::get<0>(inters);
  auto right = std::get<1>(inters);
  
  // here we actually figure out which intervals need to be updated with the quadratic, which need to be updated 
  // with the line. For those we create a new quadratic.
  
  std::list<Quadratic> new_quads;
  
  // for every quadratic
  for (auto& q : Q) {
    
    bool update_w_q = false;
    std::list<Interval> new_intervals;
    
    // for every interval of that quadratic
    for (auto& i:q.ints) {
      if (i.l >= left && i.u <= right) { // we update with the quadratic
        update_w_q = true;
        
      } else if (i.l > right || i.u < left) { // we update with the line (and mark the relative interval for deletion on the quad)
        new_intervals.push_back(i);
        i.l = std::nanf(""); i.u = std::nanf("");
        
      } else if ( inRange(left, i) &&  inRange(right, i)){ // we split in the middle
        new_intervals.push_back(I(i.l, left));
        new_intervals.push_back(I(right, i.u));
        i.l = left;
        i.u = right;
        update_w_q = true;
        
      } else if (i.l < right && i.u > right) {  // we split on the right side (quad, line)
        new_intervals.push_back(I(right, i.u));
        i.u = right;
        update_w_q = true;
        
      } else if (i.l < left && i.u > left) {  // we split on the left (line quad)
        new_intervals.push_back(I(i.l, left));
        i.l = left;
        update_w_q = true;
      }
    }
    
    // here we're done with setting up the new intervals, and we actually update
    
    // removing the pruned intervals
    q.ints.remove_if([](auto& interval){return std::isnan(interval.l);});
    
    
    if (!new_intervals.empty()) {
      // here if we have new intervals we create a new quad
      Quadratic extra_q;
      extra_q.a = q.a + new_l.a; extra_q.b = q.b + new_l.b; extra_q.c = q.c + new_l.c;
      extra_q.ints = new_intervals;
      
      new_quads.push_back(extra_q); // which we add to our new quads list
    }
    
    
    // here you do the standard update with the quad
    if (update_w_q) {
      q.a += new_q.a;
      q.b += new_q.b;
      q.c += new_q.c;
    }
      
  }
  
  // after all the quadratics have gone trough remove the empty quads from the list
  Q.remove_if([](auto& q){
    return q.ints.size() == 0;
  });
  
  // and add the new quadratics to finish
  
  
  for (auto& q:new_quads)
    Q.push_back(std::move(q));
  
  
}