#include <iostream>
#include <TTree.h>
#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>
#include "Analysis.h"
#include "AT_Charged.h"

AT_Charged::AT_Charged() : AT_ReadTree() {
}

AT_Charged::~AT_Charged() {
  if(hPt) delete hPt;
  if(hNTrk) delete hNTrk;
}

void AT_Charged::MyInit() {
  hPt =   new TH1F("hPt",  "hPt", 100,0.0,5.0);
  hNTrk = new TH1F("hNTrk","hNTrk",100,-0.5,99.5);
}

void AT_Charged::MyFinish() {
  hPt->Write();
  hNTrk->Write();
}

void AT_Charged::MyExec() {
  fCandidates->clear();
  float vtxz = fGLB.vtxZ;
  float cent = fGLB.cent;
  if(TMath::Abs(vtxz)>20) return;
  if(cent<0||cent>5) return;
  if(ReferenceTracks()<2) return;
  int ntrk=0;
  uint ntrks = pTRKpt->size();
  for(uint itrk=0; itrk!=ntrks; ++itrk) {
    float pt   = TMath::Abs( pTRKpt->at(itrk) );
    float pz   = pTRKpz->at(itrk);
    float pp = TMath::Sqrt( pt*pt + pz*pz );
    float eta = TMath::ATanH( pz / pp );
    float phi  = pTRKphi->at(itrk);
    //float ep   = pTRKecore->at(itrk) / TMath::Sqrt(pt*pt+pz*pz);
    //float tof  = pTRKetof->at(itrk);
    //float chi2 = pTRKchisq->at(itrk);
    float dphi = pTRKpc3sdphi->at(itrk);
    float dz   = pTRKpc3sdz->at(itrk);
    float zed  = pTRKzed->at(itrk);
    int qua = pTRKqua->at(itrk);
    if(qua!=63) continue;
    if(TMath::Abs(zed)<3||TMath::Abs(zed)>70) continue;
    if(TMath::Abs(dphi)>3) continue;
    if(TMath::Abs(dz)>3) continue;
    TLorentzVector lv;
    lv.SetPtEtaPhiM(pt, eta, phi, 0);
    fCandidates->push_back(lv);
    ntrk++;
    hPt->Fill(pt);
  }
  hNTrk->Fill(ntrk);
}
