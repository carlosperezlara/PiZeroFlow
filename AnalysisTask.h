#ifndef __ANALYSIS_TASK_HH__
#define __ANALYSIS_TASK_HH__

#include <iostream>
#include <vector>
#include <TObject.h>
#include <TString.h>
#include <TLorentzVector.h>
#include "qcQ.h"

const Int_t kBBC_NQ = 9;
const Int_t kBBC_NQR = 6;

class AnalysisTask : public TObject { // needed to add to TLists
 public:
  AnalysisTask() {
    fCandidates = NULL;
    fCandidates2 = NULL;
    for(int i=0; i!=kBBC_NQR; ++i) {
      fQ[i]=NULL;
      fQab[i][0]=NULL;
      fQab[i][1]=NULL;
    }
  }
  virtual ~AnalysisTask() {}
  virtual void Init() {std::cout << "AT::INIT" << std::endl;}
  virtual Bool_t Exec() {std::cout << "AT::EXEC" << std::endl; return kTRUE;}
  virtual void Finish() {std::cout << "AT::FINISH" << std::endl;}
  void LinkCandidates( std::vector<TLorentzVector> *ref ) { fCandidates = ref; }
  void LinkCandidates2( std::vector<TLorentzVector> *ref ) { fCandidates2 = ref; }
  void LinkQ( Int_t idx, qcQ *qt, qcQ *qa, qcQ *qb )
  { fQ[idx]=qt; fQab[idx][0]=qa; fQab[idx][1]=qb; }

 protected:
  std::vector<TLorentzVector> *fCandidates;
  std::vector<TLorentzVector> *fCandidates2;
  qcQ *fQ[kBBC_NQR];
  qcQ *fQab[kBBC_NQR][2];
};

#endif
