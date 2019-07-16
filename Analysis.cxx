#include <iostream>
#include <vector>

#include <TString.h>
#include <TList.h>
#include <TFile.h>
#include <TH2F.h>
#include <TTree.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TLorentzVector.h>

#include "Analysis.h"
#include "AnalysisTask.h"

Analysis *Analysis::fAnalysis = NULL;


Analysis::Analysis() {
  fInputFileName = "input.root";
  fOutputFileName = "output.root";
  fDSTag="";
  fNoSkipEventsAtBeginning=0;
  fNoEventsAnalyzed=-1;
  fListOfTasks = new TList();
  fListOfTasks->SetOwner();
  fInputFile = NULL;
  fTree = NULL;
  fCandidates = new std::vector<TLorentzVector>;
  fCandidates2 = new std::vector<TLorentzVector>;
  for(int i=0; i!=kBBC_NQR; ++i) {
    fQ[i] = new qcQ(i+1);
    fQab[i][0] = new qcQ(i+1);
    fQab[i][1] = new qcQ(i+1);
  }
}
//=====
Analysis::~Analysis() {
  delete fListOfTasks;
  if(fInputFile) delete fInputFile;
  delete fCandidates;
  delete fCandidates2;
  for(int i=0; i!=kBBC_NQR; ++i) {
    delete fQ[i];
    delete fQab[i][0];
    delete fQab[i][1];
  }
}
//=====
void Analysis::Init() {
  std::cout << "** Analysis::Init() **" << std::endl;
  std::cout << " Reading from file " << fInputFileName.Data() << std::endl;
  fInputFile = new TFile(fInputFileName.Data(),"READ");
  if(!fInputFile) return;
  fTree = (TTree*) fInputFile->Get("TOP");
  if(!fTree) {
    std::cout << " No Tree found!!" << std::endl;
    return;
  }
  //---
  int ntsk = fListOfTasks->GetEntries();
  for(int i=0; i!=ntsk; ++i) {
    AnalysisTask *tsk = (AnalysisTask*) fListOfTasks->At(i);
    //=== linking
    tsk->LinkCandidates( fCandidates );
    tsk->LinkCandidates2( fCandidates2 );
    for(int i=0; i!=kBBC_NQR; ++i)
      tsk->LinkQ( i, fQ[i], fQab[i][0], fQab[i][1] );
    //=======
    tsk->Init();
  }
}
//=====
void Analysis::Finish() {
  std::cout << "** Analysis::Finish() **" << std::endl;
  TFile *fOutputFile = new TFile(fOutputFileName.Data(),"RECREATE");
  fOutputFile->cd();
  //---
  int ntsk = fListOfTasks->GetEntries();
  for(int i=0; i!=ntsk; ++i) {
    AnalysisTask *tsk = (AnalysisTask*) fListOfTasks->At(i);
    tsk->Finish();
  }
  //  hEvents->Write();
  fOutputFile->Close();
  delete fOutputFile;
  std::cout << "Results saved into " << fOutputFileName.Data() << std::endl;
  fInputFile->Close();
}
//=====
void Analysis::Exec() {
  std::cout << "** Analysis::Exec() **" << std::endl;
  if(!fTree) return;
  Long64_t EndOfLoop = fTree->GetEntries();
  if(fNoEventsAnalyzed>0) {
    Long64_t sum = fNoSkipEventsAtBeginning + fNoEventsAnalyzed;
    if(sum<EndOfLoop) EndOfLoop = sum;
  }
  if(fNoSkipEventsAtBeginning<0) fNoSkipEventsAtBeginning=0;
  for(Long64_t i1=fNoSkipEventsAtBeginning;
      i1<EndOfLoop; ++i1) {
    if(i1%50000 == 0) {
      std::cout << " Executing event number :  " << i1 << "/" << EndOfLoop;
      std::cout << Form(" (%.1f)",i1*100.0/EndOfLoop) << std::endl;
    }
    //std::cout << " LOADTREE " << fTree->LoadTree(i1) << std::endl;
    //std::cout << " SIZE " << fTree->GetEntry(i1) << std::endl;
    fTree->GetEntry(i1);
    //---
    int ntsk = fListOfTasks->GetEntries();
    for(int i=0; i!=ntsk; ++i) {
      AnalysisTask *tsk = (AnalysisTask*) fListOfTasks->At(i);
      if( !tsk->Exec() ) break; // if task return false event is discarded
    }
    //---
  }
}
//=====
int Analysis::RunNumber() {
  TObjArray *arr = fDSTag.Tokenize("_");
  TObjString *obj = (TObjString*) arr->At(0);
  TString str = obj->GetString();
  return str.Atoi();
}
//=====
int Analysis::SegmentNumber() {
  TObjArray *arr = fDSTag.Tokenize("_");
  TObjString *obj = (TObjString*) arr->At(1);
  TString str = obj->GetString();
  return str.Atoi();
}
//=====
void Analysis::Run() {
  Init();
  Exec();
  Finish();
}
