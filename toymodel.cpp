#include <vector>
#include <TRandom3.h>
#include <TFile.h>
#include <TList.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TMath.h>
#include <TProfile.h>
#include "qcQ.h"
#include "QCorrelation.h"
#include "QC_EventPlane.h"
#include "QC_ScalarProduct.h"
#include "QC_Cumulants.h"

int main() {
  qcQ *fQ = new qcQ(2);
  qcQ *fQ2 = new qcQ(4);
  qcQ *fQa = new qcQ(2);
  qcQ *fQb = new qcQ(2);
  //QC_EventPlane *epc = new QC_EventPlane();
  //QC_ScalarProduct *spc = new QC_ScalarProduct();
  QC_EventPlane *epc = new QC_EventPlane();
  QC_ScalarProduct *spc = new QC_ScalarProduct();
  QC_Cumulants *qcc = new QC_Cumulants();
  QCorrelation *all[3] = {epc,spc,qcc};
  for(int i=0; i!=3; ++i) {
    std::vector<Double_t> ptlist = all[i]->GetPt();
    ptlist.clear();
    ptlist.push_back(0);
    ptlist.push_back(1);
    all[i]->SetPt(ptlist);
  }
  epc->Init("epc",fQ,fQa,fQb);
  spc->Init("spc",fQa,fQb);
  qcc->Init("qcc",fQ,fQ2);

  TH1D *histQa = new TH1D("histQa","histQa",100,-TMath::TwoPi(),+TMath::TwoPi());
  TH1D *histQb = new TH1D("histQb","histQb",100,-TMath::TwoPi(),+TMath::TwoPi());
  TH1D *histQT = new TH1D("histQT","histQT",100,-TMath::TwoPi(),+TMath::TwoPi());
  Double_t mamin = epc->GetMaMin();
  Double_t mamax = epc->GetMaMax();
  TString sSgn = "TMath::Gaus(x,0.135,0.01,1)";
  TString sBgr = "1";
  TF1 *sgn = new TF1("sgn",sSgn.Data());
  TF1 *bgr = new TF1("bgr",sBgr.Data());
  std::cout << "v2ref 0.1 | v2sgn 0.15 | v2bgr 0.09" << std::endl;
  TF1 *v2ref = new TF1("v2ref","1-2*0.10*TMath::Cos(2*x)",0,TMath::TwoPi());
  TF1 *v2sgn = new TF1("v2sgn","1-2*0.15*TMath::Cos(2*x)",0,TMath::TwoPi());
  TF1 *v2bgr = new TF1("v2bgr","1-2*0.09*TMath::Cos(2*x)",0,TMath::TwoPi());
  //======
  for(int ev=0; ev!=10000; ++ev) {
    for(int rt=0; rt!=50; ++rt) {
      Double_t phi = v2ref->GetRandom(0,TMath::TwoPi());
      histQT->Fill( phi );
      fQ->Fill(phi,1.2);
      fQ2->Fill(phi,1.2);
      if(rt%2==0) {
	fQa->Fill(phi,1.2);
	histQa->Fill( phi );
      } else {
	fQb->Fill(phi,1.2);
	histQb->Fill( phi );
      }
    }
    epc->Reset();
    spc->Reset();
    qcc->Reset();
    for(int tt=0; tt!=100; ++tt) {
      Double_t phi = 0;
      Double_t mass = 0;
      if( gRandom->Rndm() < 0.1 ) { // sgn
	phi = v2sgn->GetRandom(0,TMath::TwoPi());
	mass = sgn->GetRandom(mamin,mamax);
      } else { // bgr
	phi = v2bgr->GetRandom(0,TMath::TwoPi());
	mass = bgr->GetRandom(mamin,mamax);
      }
      epc->FillCandidate(0.5,mass,phi);
      spc->FillCandidate(0.5,mass,phi);
      qcc->FillCandidate(0.5,mass,phi);
    }
    epc->FillEvent();
    spc->FillEvent();
    qcc->FillEvent();
  }
  //======
  std::cout << "data generated" << std::endl;
  epc->Results();
  spc->Results();
  qcc->Results();
  std::cout << "results processed" << std::endl;
  epc->Save("tm_ep.root");
  spc->Save("tm_sp.root");
  qcc->Save("tm_qc.root");
  std::cout << "all saved" << std::endl;

  TCanvas *main = new TCanvas("main","main");
  main->SaveAs("qaplots.pdf[","PDF");
  histQa->SetLineColor(kRed);
  histQb->SetLineColor(kBlue);
  histQT->SetLineColor(kBlack);
  histQa->DrawNormalized();
  histQb->DrawNormalized("same");
  histQT->DrawNormalized("same");
  main->SaveAs("qaplots.pdf","PDF");
  
  TFile *file;
  TList *listA;
  file = new TFile("tm_qc.root");
  listA = (TList*) file->Get("qcc");
  TH1D *yie = ( (TH2D*) listA->FindObject("binmap_qcc") )->ProjectionY("yiel",1,1);
  listA = (TList*) listA->FindObject("results");
  TH1D *qc2 = (TH1D*) listA->FindObject("PV2_0_qcc");
  TH1D *qc4 = (TH1D*) listA->FindObject("PV4_0_qcc");

  main->Divide(1,2);
  qc2->SetLineColor(kRed-3);
  qc4->SetLineColor(kBlue-3);
  main->cd(1);
  yie->Draw();
  main->cd(2);
  qc4->Draw("");
  qc2->Draw("same");
  main->SaveAs("qaplots.pdf","PDF");

  main->SaveAs("qaplots.pdf]","PDF");

  return 0;
}
