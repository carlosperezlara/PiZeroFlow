#ifndef MODEL3_H
#define MODEL3_H
#include "DoubleFitter.h"

using std::string;

class Model3: public DoubleFitter {
  virtual double Sgn(double x, double *p);
  virtual double Bgr(double x, double *p);
  virtual double FBgr(double x, double *p);

 public:
  Model3();
  virtual void Init();
};

#endif
