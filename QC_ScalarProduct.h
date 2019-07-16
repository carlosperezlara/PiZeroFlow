#ifndef __QC_ScalarProduct_HH__
#define __QC_ScalarProduct_HH__

#include <map>
#include <TString.h>
#include "QCorrelation.h"
#include "qcQ.h"

class TProfile;
class TProfile2D;
class TList;
class TH1D;

class QC_ScalarProduct : public QCorrelation {
 public:
  QC_ScalarProduct();
  QC_ScalarProduct(TString,qcQ*,qcQ*);
  virtual ~QC_ScalarProduct();
  virtual void SetList(TList*);
  virtual void Reset();
  virtual void FillCandidate(Double_t pt, Double_t ma, Double_t phi);
  virtual void FillEvent();
  void Results();
  void Init(TString,qcQ*,qcQ*);
  
 private:
  qcQ *fQa; // not owned
  qcQ *fQb; // not owned
  qcQ *fNULL;
  std::map<Int_t,qcQ> fP2;
  //
  TProfile *fQPC;
  //
  TProfile2D *fPQa_2PC;
  TProfile2D *fPQb_2PC;
  TProfile2D *fPQaxQa_2PC;
  TProfile2D *fPQbxQb_2PC;
  TProfile2D *fPQaxQb_2PC;
  TProfile2D *fPQbxQa_2PC;
  TProfile2D *fPQaxPQb_2PC;
  //
  TH1D *fQA_Qa;
  TH1D *fQA_Qb;
};

#endif
