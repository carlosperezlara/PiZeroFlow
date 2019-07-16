#ifndef __Analysis_HH__
#define __Analysis_HH__

#include <vector>
#include <TString.h>
#include <TList.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include "qcQ.h"
#include "AnalysisTask.h"

class TFile;
class TTree;

class Analysis {
 public:
  static Analysis* Instance() {
    if(!fAnalysis) fAnalysis = new Analysis();
    return fAnalysis;
  }
  virtual ~Analysis();
  void Run();
  void Init();
  void Exec();
  void Finish();
  void AddTask(AnalysisTask *tsk) {fListOfTasks->Add(tsk);}
  void InputFileName(TString name) {fInputFileName = name;}
  void OutputFileName(TString name) {fOutputFileName = name;}
  void DataSetTag(TString name) {fDSTag =name;}
  void NumberOfEventsToSkipAtBeginning(Long64_t skp) {fNoSkipEventsAtBeginning = skp;}
  void NumberOfEventsToAnalyze(Long64_t nev) {fNoEventsAnalyzed = nev;}
  TTree* GetTree() {return fTree;}
  TString GetInputFileName() {return fInputFileName;} // used for calibration purposes
  std::vector<TLorentzVector>* GetCandidates() {return fCandidates;}
  std::vector<TLorentzVector>* GetCandidates2() {return fCandidates2;}
  qcQ* GetQ(int n) {return fQ[n];}
  qcQ* GetQab(int n,int se) {return fQab[n][se];}
  int RunNumber();
  int SegmentNumber();

 protected:
  Analysis();
  
 private:
  static Analysis *fAnalysis;
  TString fInputFileName;
  TString fOutputFileName;
  TString fDSTag;
  Long64_t fNoSkipEventsAtBeginning;
  Long64_t fNoEventsAnalyzed;
  TList *fListOfTasks;
  TFile *fInputFile;
  TTree *fTree;
  std::vector<TLorentzVector> *fCandidates;
  std::vector<TLorentzVector> *fCandidates2;
  qcQ *fQ[kBBC_NQR];
  qcQ *fQab[kBBC_NQR][2];
};

#endif
