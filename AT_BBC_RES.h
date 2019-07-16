#ifndef __AT_BBC_RES_HH__
#define __AT_BBC_RES_HH__

#include "AT_ReadTree.h"

class TProfile2D;

class AT_BBC_RES : public AT_ReadTree {
 public:
  AT_BBC_RES();
  virtual ~AT_BBC_RES();
  virtual void MyInit();
  virtual void MyExec();
  virtual void MyFinish();

 private:
  TProfile2D *hRes[4]; // ord
};

#endif
