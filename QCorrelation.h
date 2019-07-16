#ifndef __QCORRELATION_HH__
#define __QCORRELATION_HH__
#include <vector>
#include <TString.h>
#include <TList.h>
#include <TH2D.h>
#include "qcQ.h"

class QCorrelation {
 public:
  QCorrelation();
  virtual ~QCorrelation();
  virtual void FillCandidate(Double_t pt, Double_t ma, Double_t phi) {}
  virtual void FillEvent() {}
  void SaveResults(TString name);
  void Save(TString name);
  virtual void SetList(TList*);
  TList* GetList() {return fList;}
  TList* GetListResults() {return fListResults;}
  //
  Double_t GetMaMin() {return fMaList.size()<1?-1:fMaList[0];}
  Double_t GetMaMax() {return fMaList.size()<1?-1:fMaList[fMaList.size()-1];}
  std::vector<Double_t> GetPt() {return fPtList;}
  std::vector<Double_t> GetMa() {return fMaList;}
  std::vector< std::vector<Double_t> > GetNewMa() {return fNewMaList;}
  void SetPt( std::vector<Double_t> val ) {fPtList=val;}
  void SetMa( std::vector<Double_t> val ) {fMaList=val;}
  void SetNewMa( std::vector< std::vector<Double_t> > val ) {fNewMaList=val;}

 protected:
  void ClearResults();
  void Init(TString);
  Double_t Make2PC(qcQ*, qcQ*, qcQ*, Double_t&); // a b ab err
  Double_t Make4PC(qcQ*, qcQ*, qcQ*, qcQ*, qcQ*, Double_t&); // a b ab b2 ab2 err
  Double_t MakeQC2(Double_t, Double_t, Double_t&);
  Double_t MakeQC4(Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t&);
  Double_t MakeV2(Double_t,Double_t,Double_t,Double_t,Double_t,Double_t&);
  Double_t MakeV4(Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t,Double_t&);

  TList *fList;
  TList *fListResults;
  TH2D *fBinMap;
  // used to create binmap during init
  std::vector<Double_t> fPtList;
  std::vector<Double_t> fMaList;
  // used to create results during closure
  std::vector< std::vector<Double_t> > fNewMaList;

};
#endif
