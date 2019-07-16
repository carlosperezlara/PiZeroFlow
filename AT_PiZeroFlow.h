#ifndef __AT_PIZEROFLOW_HH__
#define __AT_PIZEROFLOW_HH__

#include "AT_PiZero.h"

class AT_PiZeroFlow : public AT_PiZero {
 public:
  AT_PiZeroFlow();
  virtual ~AT_PiZeroFlow();
  virtual void MyInit();
  virtual void MyExec();
  virtual void MyFinish();

};

#endif
