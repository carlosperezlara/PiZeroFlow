#ifndef __AT_PIDFLOW_HH__
#define __AT_PIDFLOW_HH__

#include "AT_ReadTree.h"

class TH1F;
class TH2F;

class AT_PIDFlow : public AT_ReadTree {
 public:
  AT_PIDFlow();
  virtual ~AT_PIDFlow();
  virtual void MyInit();
  virtual void MyExec();
  virtual void MyFinish();

 private:
  TH1F *hPt;
  TH1F *hNTrk;
  TH2F *hPtDPhi[4];
  TH2F *hPtDPhiME[4];
  TH1F *hEP_BBC[4];
  float fPsi1_BBC_PE;
  float fPsi2_BBC_PE;
  float fPsi3_BBC_PE;
  float fPsi4_BBC_PE;

};

#endif
