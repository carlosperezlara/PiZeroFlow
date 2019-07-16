#ifndef __AT_CHARGED_HH__
#define __AT_CHARGED_HH__

#include <vector>
#include <TLorentzVector.h>
#include "AT_ReadTree.h"

class TH1F;

class AT_Charged : public AT_ReadTree {
 public:
  AT_Charged();
  virtual ~AT_Charged();
  virtual void MyInit();
  virtual void MyExec();
  virtual void MyFinish();

 private:
  TH1F *hPt;
  TH1F *hNTrk;

};

#endif
