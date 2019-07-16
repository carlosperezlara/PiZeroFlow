#ifndef __AT_EVENTCHECKER_HH__
#define __AT_EVENTCHECKER_HH__

#include "AT_ReadTree.h"

class TH1F;
class TH2F;
class TProfile;

class AT_EventChecker : public AT_ReadTree {
 public:
  AT_EventChecker();
  virtual ~AT_EventChecker();
  virtual void MyInit();
  virtual void MyExec();
  virtual void MyFinish();

 private:
  TH1F *hPsi2[10];
  TH1F *hTrackMultiplicity[2];
  TH1F *hBBCMultiplicity[2];
  TH1F *hBBCTmean[2];
  TH1F *hBBCTrms[2];
  TH1F *hBBCTrmsS[2];
  TH1F *hBBCTnrmsN[2];
  TH1F *hBBCTnrmsS[2];
  TH1F *hBBCTrmsN[2];
  TH2F *hBBCmeanT[2];
  TH2F *hBBCrmsT[2];
  TH2F *hBBCsgn[2];
  TH1F *hBBCsgnD[2];
  TH1F *hBBCvtx[2];  
  TH1F *hFrac[2];
};

#endif
