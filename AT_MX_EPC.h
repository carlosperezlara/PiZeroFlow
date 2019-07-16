#ifndef __AT_MX_EPC_HH__
#define __AT_MX_EPC_HH__

#include "AT_ReadTree.h"

class TH1F;
class TH2F;
class TProfile2D;

class AT_MX_EPC : public AT_ReadTree {
 public:
  AT_MX_EPC();
  virtual ~AT_MX_EPC();
  virtual void MyInit();
  virtual void MyExec();
  virtual void MyFinish();

 private:
  TH2F *hQxVtx[6][8][60]; // ord se cbin
  TH2F *hQyVtx[6][8][60]; // ord se cbin

  TH1F *hQxC[3][4][8][60]; // step ord se cbin
  TH1F *hQyC[3][4][8][60]; // step ord se cbin

  TProfile2D *hPsiC[4][60]; // ord cbin
  TProfile2D *hPsiS[4][60]; // ord cbin
  TH2F *hDeltaPsi[4][60]; // ord cbin

  TProfile *hRes[4][60]; // ord cbin
};

#endif
