#ifndef __AT_EP_HH__
#define __AT_EP_HH__

#include "QC_EventPlane.h"
#include "QC_ScalarProduct.h"
#include "QC_Cumulants.h"
#include "AnalysisTask.h"

class TH1D;
class TList;
class TProfile;

class AT_EP : public AnalysisTask {
 public:
  AT_EP();
  virtual ~AT_EP();
  virtual void Init();
  virtual Bool_t Exec();
  virtual void Finish();

 private:
  TList *fListEP;
  QC_EventPlane    *fEP[4];
  QC_ScalarProduct *fSP[4];
  QC_Cumulants     *fQC[4];

  int BinPt(float);
  int BinMass(float);

  Int_t fNpt;
  Double_t fPtBins[100];
  Int_t fNma;
  Double_t fMassBins[300];

  TH1D *hStats;
  TH1D *hEta;
  TH1D *hEta2;
  TH1D *hMass[100];
  TH1D *hMass2[100];
  TProfile *hMass3[100];
};

#endif
