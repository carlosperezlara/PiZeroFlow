#include <iostream>
#include <TTree.h>
#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>
#include "Analysis.h"
#include "AT_PIDFlow.h"

AT_PIDFlow::AT_PIDFlow() : AT_ReadTree() {
}

AT_PIDFlow::~AT_PIDFlow() {
  if(hPt) delete hPt;
  if(hNTrk) delete hNTrk;
  for(int i=0; i!=4; ++i) {
    if(hPtDPhi[i]) delete hPtDPhi[i];
    if(hPtDPhiME[i]) delete hPtDPhiME[i];
    if(hEP_BBC[i]) delete hEP_BBC[i];
  }
}

void AT_PIDFlow::MyInit() {
  hPt =   new TH1F("hPt",  "hPt", 100,0.0,5.0);
  hNTrk = new TH1F("hNTrk","hNTrk",100,-0.5,99.5);
  double twopi = TMath::TwoPi();
  int nbin[4] = {50,30,20,20};
  for(int ord=0; ord!=4; ++ord) {
    int nh = ord+1;
    hPtDPhi[ord] = new TH2F(Form("hPtDPhi%d",nh), Form(";pt;phi-Psi%d",nh),    100,0.0,5.0, nbin[ord],-twopi,+twopi);
    hPtDPhiME[ord] = new TH2F(Form("hPtDPhi%dME",nh),Form(";pt;phi-Psi%d",nh), 100,0.0,5.0, nbin[ord],-twopi,+twopi);
    hEP_BBC[ord] = new TH1F(Form("hPsi%d_BBC",nh),Form(";Psi%d",nh),100,0,7);
  }
}

void AT_PIDFlow::MyFinish() {
  hPt->Write();
  hNTrk->Write();
  for(int i=0; i!=4; ++i) {
    hPtDPhi[i]->Write();
    hPtDPhiME[i]->Write();
    hEP_BBC[i]->Write();
  }
}

void AT_PIDFlow::MyExec() {
  float vtxz = fGLB.vtxZ;
  float cent = fGLB.cent;
  if(TMath::Abs(vtxz)>20) return;
  if(cent<0||cent>5) return;
  if(ReferenceTracks()<2) return;
  if(!Psi_BBC) return;
  hEP_BBC[0]->Fill(Psi1_BBC);
  hEP_BBC[1]->Fill(Psi2_BBC);
  hEP_BBC[2]->Fill(Psi3_BBC);
  hEP_BBC[3]->Fill(Psi4_BBC);
  int ntrk=0;
  uint ntrks = pTRKpt->size();
  for(uint itrk=0; itrk!=ntrks; ++itrk) {
    float pt   = TMath::Abs( pTRKpt->at(itrk) );
    float pz   = pTRKpz->at(itrk);
    float phi  = pTRKphi->at(itrk);
    float ep   = pTRKecore->at(itrk) / TMath::Sqrt(pt*pt+pz*pz);
    float tof  = pTRKetof->at(itrk);
    float chi2 = pTRKchisq->at(itrk);
    float dphi = pTRKpc3sdphi->at(itrk);
    float dz   = pTRKpc3sdz->at(itrk);
    float zed  = pTRKzed->at(itrk);
    int qua = pTRKqua->at(itrk);
    if(qua!=63) continue;
    if(TMath::Abs(zed)<3||TMath::Abs(zed)>70) continue;
    if(TMath::Abs(dphi)>3) continue;
    if(TMath::Abs(dz)>3) continue;
    ntrk++;
    hPt->Fill(pt);
    float dphi1 = phi - Psi1_BBC;
    float dphi2 = phi - Psi2_BBC;
    float dphi3 = phi - Psi3_BBC;
    float dphi4 = phi - Psi4_BBC;
    hPtDPhi[0]->Fill(pt,dphi1);
    hPtDPhi[1]->Fill(pt,dphi2);
    hPtDPhi[2]->Fill(pt,dphi3);
    hPtDPhi[3]->Fill(pt,dphi4);
    float dphi1me = phi - fPsi1_BBC_PE;
    float dphi2me = phi - fPsi2_BBC_PE;
    float dphi3me = phi - fPsi3_BBC_PE;
    float dphi4me = phi - fPsi4_BBC_PE;
    hPtDPhiME[0]->Fill(pt,dphi1me);
    hPtDPhiME[1]->Fill(pt,dphi2me);
    hPtDPhiME[2]->Fill(pt,dphi3me);
    hPtDPhiME[3]->Fill(pt,dphi4me);
  }
  hNTrk->Fill(ntrk);
  fPsi1_BBC_PE = Psi1_BBC;
  fPsi2_BBC_PE = Psi2_BBC;
  fPsi3_BBC_PE = Psi3_BBC;
  fPsi4_BBC_PE = Psi4_BBC;
}
