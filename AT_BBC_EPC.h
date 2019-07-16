#ifndef __AT_BBC_EPC_HH__
#define __AT_BBC_EPC_HH__

#include "AT_ReadTree.h"

class TList;
class TH1D;
class TH2D;
class TProfile2D;
class TProfile;

class AT_BBC_EPC : public AT_ReadTree {
 public:
  AT_BBC_EPC();
  virtual ~AT_BBC_EPC();
  virtual void MyInit();
  virtual void MyExec();
  virtual void MyFinish();

 private:
  TList *fListBBCEPC;

  TList *fListCalib;
  TH2D *hQxVtx[kBBC_NQ][3][kBinsCen][3]; // ord se cbin step
  TH2D *hQyVtx[kBBC_NQ][3][kBinsCen][3]; // ord se cbin step
  TProfile2D *hPsiC[3][kBBC_NQR][kBinsCen]; // se ord cbin
  TProfile2D *hPsiS[3][kBBC_NQR][kBinsCen]; // se ord cbin

  TList *fListDelta;
  TH2D *hDeltaPsi[3][kBBC_NQR][kBinsCen]; // se ord cbin

  TList *fListReso;
  TProfile *hRes[kBBC_NQR][kBinsCen]; // ord cbin
};

#endif
