#ifndef __QC_Cumulants_HH__
#define __QC_Cumulants_HH__

#include <map>
#include <TString.h>
#include "QCorrelation.h"
#include "qcQ.h"

class TProfile;
class TProfile2D;
class TList;
class TH1D;

class QC_Cumulants : public QCorrelation {
 public:
  QC_Cumulants();
  QC_Cumulants(TString,qcQ*,qcQ*);
  virtual ~QC_Cumulants();
  virtual void SetList(TList*);
  virtual void Reset();
  virtual void FillCandidate(Double_t pt, Double_t ma, Double_t phi);
  virtual void FillEvent();
  void Results();
  void Init(TString,qcQ*,qcQ*);
  
 private:
  qcQ *fQ2; // not owned
  qcQ *fQ4; // not owned
  qcQ *fNULL;
  std::map<Int_t,qcQ> fP2;
  //
  TProfile *fQPC;
  //
  TProfile2D *fP2PC;
  TProfile2D *fP4PC;
  TProfile2D *fP4P2PC;
  TProfile2D *fP4Q2PC;
  TProfile2D *fP4Q4PC;
  TProfile2D *fP2Q2PC;
  TProfile2D *fP2Q4PC;
  //
  TH1D *fQA_Q2;
  TH1D *fQA_Q4;
};

#endif
