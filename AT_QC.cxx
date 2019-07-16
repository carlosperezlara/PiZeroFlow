#include <iostream>
#include <TTree.h>
#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>
#include "Analysis.h"
#include "AT_QC.h"

AT_QC::AT_QC() {
}

AT_QC::~AT_QC() {
}

void AT_QC::Init() {
  Analysis *ana = Analysis::Instance();
  fCandidates = ana->GetCandidates();
  for(int i=0; i!=4; ++i) fQ[i] = ana->GetQ(i);

  float ptbins[11] = {0.4, 0.6, 0.8, 1.0, 1.2,
                      1.4, 1.6, 1.8, 2.0, 2.5,
                      3.0};

  TProfile *fPQ_Cos[4];
  for(int n=0; n!=4; ++n) {
    fPxQx[n] = new TProfile( Form("P%dxQ%dx",i,i),
			     Form("P_{%d,x}Q_{%d,x}",i,i),
			     10,ptbins,"s");
    fPyQy[n] = new TProfile( Form("P%dyQ%dy",i,i),
			     Form("P_{%d,y}Q_{%d,y}",i,i),
			     10,ptbins,"s");
    fQxQx[n] = new TProfile( Form("Q%dxQ%dx",i,i),
			     Form("Q_{%d,x}Q_{%d,x}",i,i),
			     10,ptbins,"s");
    fQyQy[n] = new TProfile( Form("Q%dyQ%dy",i,i),
			     Form("Q_{%d,y}Q_{%d,y}",i,i),
			     10,ptbins,"s");
  }
}

void AT_QC::Finish() {
}

void AT_QC::Exec() {
  uint npa = fCandidates->size();
  //std::cout << "CANDIDATES " << npa <<std::endl;
  //std::cout << "Q2.M " << fQ[1]->M() << std::endl;
  if(npa==0) return;
  for(int i=0; i!=4; ++i)
    hQxQx[i]->Fill( fQ[i].X() * fQ[i].X() );

  qcQ u0(0);
  qcQ u1(1);
  qcQ u2(2);
  qcQ u3(3);
  for(int i=0; i!=npa; ++i) {
    TLorentzVector a = fCandidates->at(i);
    u0->Fill( a.Phi(), 1 );
    u1->Fill( a.Phi(), 1 );
    u2->Fill( a.Phi(), 1 );
    u3->Fill( a.Phi(), 1 );
  }


}
