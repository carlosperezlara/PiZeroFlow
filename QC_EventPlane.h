#ifndef __QC_EventPlane_HH__
#define __QC_EventPlane_HH__

#include <map>
#include <TString.h>
#include "QCorrelation.h"
#include "qcQ.h"

class TProfile;
class TProfile2D;
class TList;
class TH1D;

class QC_EventPlane : public QCorrelation {
 public:
  QC_EventPlane();
  QC_EventPlane(TString,qcQ*,qcQ*,qcQ*);
  virtual ~QC_EventPlane();
  virtual void SetList(TList*);
  virtual void Reset();
  virtual void FillCandidate(Double_t pt, Double_t ma, Double_t phi);
  virtual void FillEvent();
  void Results();
  void Init(TString,qcQ*,qcQ*,qcQ*);
  
 private:
  qcQ *fQa; // not owned
  qcQ *fQb; // not owned
  qcQ *fQ;  // not owned
  //
  TProfile *fQC;
  //
  TProfile2D *fPQC;
  //
  TH1D *fQA_Qa;
  TH1D *fQA_Qb;
  TH1D *fQA_Q;
};

#endif
