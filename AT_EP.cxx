#include <iostream>
#include <TTree.h>
#include <TMath.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TList.h>
#include "Analysis.h"
#include "AT_EP.h"
#include "QC_EventPlane.h"
#include "QC_ScalarProduct.h"
#include "QC_Cumulants.h"

AT_EP::AT_EP() {
}

AT_EP::~AT_EP() {
  if(fListEP) delete fListEP;
}

void AT_EP::Init() {
  fListEP = new TList();
  fListEP->SetOwner();
  fListEP->SetName("AT_EP");

  hStats = new TH1D("hStats","Stats",4,-0.5,3.5);
  hStats->GetXaxis()->SetBinLabel(1,"Reached");
  hStats->GetXaxis()->SetBinLabel(2,"QM>=1");
  hStats->GetXaxis()->SetBinLabel(3,"Candidates>0");
  fListEP->Add( hStats );
  //==
  fEP[0] = new QC_EventPlane( Form("EP_Ord%d_Psi%d",0,0), fQ[0], fQab[0][0], fQab[0][1]); // v1
  fEP[1] = new QC_EventPlane( Form("EP_Ord%d_Psi%d",1,1), fQ[1], fQab[1][0], fQab[1][1]); // v2
  fEP[2] = new QC_EventPlane( Form("EP_Ord%d_Psi%d",2,2), fQ[2], fQab[2][0], fQab[2][1]); // v3
  fEP[3] = new QC_EventPlane( Form("EP_Ord%d_Psi%d",3,3), fQ[3], fQab[3][0], fQab[3][1]); // v4
  //==
  fSP[0] = new QC_ScalarProduct( Form("SP_Ord%d_Psi%d",0,0), fQab[0][0], fQab[0][1]); // v1
  fSP[1] = new QC_ScalarProduct( Form("SP_Ord%d_Psi%d",1,1), fQab[1][0], fQab[1][1]); // v2
  fSP[2] = new QC_ScalarProduct( Form("SP_Ord%d_Psi%d",2,2), fQab[2][0], fQab[2][1]); // v3
  fSP[3] = new QC_ScalarProduct( Form("SP_Ord%d_Psi%d",3,3), fQab[3][0], fQab[3][1]); // v4
  //==
  fQC[0] = new QC_Cumulants( Form("QC_Ord%d_Psi%d",0,0), fQ[0], fQ[1]); // v1,2
  fQC[1] = new QC_Cumulants( Form("QC_Ord%d_Psi%d",1,1), fQ[1], fQ[3]); // v2,4
  fQC[2] = new QC_Cumulants( Form("QC_Ord%d_Psi%d",2,2), fQ[2], fQ[5]); // v3,6
  fQC[3] = new QC_Cumulants( Form("QC_Ord%d_Psi%d",1,3), fQ[3], fQ[2]); // v4,3

  std::vector<Double_t> pts = fEP[0]->GetPt();
  Int_t nma = 110;
  fNma = 1;
  fMassBins[0] = fEP[0]->GetMaMin();
  fMassBins[1] = fEP[1]->GetMaMax();

  for(int i=0; i!=pts.size(); ++i) {
    fPtBins[i] = pts.at(i);
  }
  fNpt = pts.size()-1;
  hEta = new TH1D("hEta","hEta",100,-1,+1);
  fListEP->Add( hEta );
  hEta2 = new TH1D("hEta2","hEta2",100,-1,+1);
  fListEP->Add( hEta2 );
  for(int p=0; p!=fNpt; ++p) {
    hMass[p] = new TH1D( Form("hMass_PB%d",p),
			 Form("hMass_PB%d;Mass",p),
			 nma, fMassBins[0], fMassBins[1] );
    fListEP->Add( hMass[p] );
    hMass2[p] = new TH1D( Form("hMass2_PB%d",p),
			  Form("hMass2_PB%d;Mass",p),
			  nma, fMassBins[0], fMassBins[1] );
    fListEP->Add( hMass2[p] );
    hMass3[p] = new TProfile( Form("hMass3_PB%d",p),
			      Form("hMass3_PB%d;Mass",p),
			      nma, fMassBins[0], fMassBins[1] );
    fListEP->Add( hMass3[p] );
  }
}

void AT_EP::Finish() {
  fListEP->Write(fListEP->GetName(),kSingleKey);

  for(int ord=0; ord!=4; ++ord) {
    TList *tmp = fEP[ord]->GetList();
    tmp->Write(tmp->GetName(),kSingleKey);
    tmp = fSP[ord]->GetList();
    tmp->Write(tmp->GetName(),kSingleKey);
    tmp = fQC[ord]->GetList();
    tmp->Write(tmp->GetName(),kSingleKey);
  }
}

Bool_t AT_EP::Exec() {
  hStats->Fill(0);
  if(fQ[0]->M()<1) return kTRUE; // exit here but do not break loop
  //for(int i=0; i!=kBBC_NQR; ++har) {
    //std::cout << " AT::EP => " << i << " PSIT : " << fQ[i]->Psi() << " ";
    //std::cout << " PSIA : " << fQab[i][0]->Psi() << " ";
    //std::cout << " PSIB : " << fQab[i][1]->Psi() << std::endl;
  //}

  for(int ord=0; ord!=4; ++ord) {
    fEP[ord]->Reset();
    fSP[ord]->Reset();
    fQC[ord]->Reset();
  }
  hStats->Fill(1);

  // CANDIDATES 1
  uint npa = fCandidates->size();
  if(npa>0) hStats->Fill(2);
  for(int i=0; i!=npa; ++i) {
    TLorentzVector a = fCandidates->at(i);
    double ma = a.M();
    double pt = a.Pt();
    int mb = BinMass( ma );
    int pb = BinPt( pt );
    if(mb<0||pb<0) continue;
    /// recording
    hEta->Fill( a.Eta() );
    hMass[pb]->Fill(ma);
    hMass3[pb]->Fill(ma,pt);
    for(int ord=0; ord!=4; ++ord) {
      fEP[ord]->FillCandidate(pt,ma,a.Phi());
      fSP[ord]->FillCandidate(pt,ma,a.Phi());
      fQC[ord]->FillCandidate(pt,ma,a.Phi());
    }
  }

  for(int ord=0; ord!=4; ++ord) {
    fEP[ord]->FillEvent();
    fSP[ord]->FillEvent();
    fQC[ord]->FillEvent();
  }
  hStats->Fill(1);

  // CANDIDATES 2
  npa = fCandidates2->size();
  for(int i=0; i!=npa; ++i) {
    TLorentzVector a = fCandidates2->at(i);
    double ma = a.M();
    double pt = a.Pt();
    int mb = BinMass( ma );
    int pb = BinPt( pt );
    if(mb<0||pb<0) continue;
    /// recording
    hEta2->Fill( a.Eta() );
    hMass2[pb]->Fill(ma);
  }
  return kTRUE;
}

int AT_EP::BinPt(float pt) {
  int ret = -1;
  if( pt<fPtBins[0] || pt>fPtBins[fNpt] ) return ret;
  for(int p=0; p!=fNpt+1; ++p) 
    if(pt>fPtBins[p]) ret = p;
  return ret;
}

int AT_EP::BinMass(float ma) {
  int ret = -1;
  if( ma<fMassBins[0] || ma>fMassBins[fNma] ) return ret;
  for(int m=0; m!=fNma+1; ++m) 
    if(ma>fMassBins[m]) ret = m;
  return ret;
}
