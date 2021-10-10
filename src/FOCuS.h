#ifndef ___FOCuS_H___
#define ___FOCuS_H___

//#include <vector>
//#include <deque>
#include <algorithm>
#include "quadratic.h"
#include "biweight_loss.h"

typedef struct {
  Quadratic Q0;
  std::list<Quadratic> Q1;
  double global_max;
} Info;


typedef struct {
  std::list<Quadratic> Qright;
  std::list<Quadratic> Qleft;
  double global_max;
} mInfo;


Info FOCuS_step(Info, const double&, const std::list<double>&, const double&);
Info FOCuS_training_step(Info, const double&, const std::list<double>&, const double&);
Info FOCuS_step_sim(Info, const double&, const std::list<double>&, const double&);
mInfo FOCuS_step_melk(mInfo, const double&, const std::list<double>&,  const std::list<double>&, const double&);
#endif
