#include <iostream>
#include <map>
#include <TMath.h>
#include <TH2D.h>
#include <TF1.h>
#include <TAxis.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TArray.h>
#include "qcQ.h"
#include "QC_EventPlane.h"

QC_EventPlane::QC_EventPlane() :
  QCorrelation() {
  fQ = NULL;
  fQa = NULL;
  fQb = NULL;
}

QC_EventPlane::QC_EventPlane(TString name,qcQ *q1,qcQ *q2,qcQ *q3) :
  QCorrelation() {
  Init(name,q1,q2,q3);
}

void QC_EventPlane::Init(TString name,qcQ *q1,qcQ *q2,qcQ *q3) {
  QCorrelation::Init(name);
  fQ = q1;
  fQa = q2;
  fQb = q3;
  fQC = new TProfile(Form("QC_%s",name.Data()),"QC",3,-0.5,2.5,-1,+1,"s"); fList->Add(fQC);
  fQC->Sumw2();
  fQC->GetXaxis()->SetBinLabel(1,"<Cos2(ab_2pi)>");
  fQC->GetXaxis()->SetBinLabel(2,"<Cos2(ab)>");
  fQC->GetXaxis()->SetBinLabel(3,"<Cos2abs(ab)>");
  //
  const TArrayD *arrpt = fBinMap->GetXaxis()->GetXbins();
  Int_t npt = arrpt->GetSize();
  Double_t npts[1000];
  for(int i=0; i!=npt; ++i) {
    npts[i] = arrpt->GetAt(i);
  }
  const TArrayD *arrma = fBinMap->GetYaxis()->GetXbins();
  Int_t nma = arrma->GetSize();
  Double_t nmas[1000];
  for(int i=0; i!=nma; ++i) {
    nmas[i] = arrma->GetAt(i);
  }
  //
  fPQC = new TProfile2D(Form("PQC_%s",name.Data()),"PQC",npt-1,npts,nma-1,nmas,"s"); fList->Add(fPQC);
  fPQC->Sumw2();
  //
  fQA_Qa = new TH1D(Form("QA_Qa_%s",name.Data()),"QA_Qa",200,-0.1,+TMath::TwoPi()); fList->Add(fQA_Qa);
  fQA_Qb = new TH1D(Form("QA_Qb_%s",name.Data()),"QA_Qb",200,-0.1,+TMath::TwoPi()); fList->Add(fQA_Qb);
  fQA_Q  = new TH1D(Form("QA_Q_%s",name.Data()), "QA_Q", 200,-0.1,+TMath::TwoPi()); fList->Add(fQA_Q);
}

void QC_EventPlane::SetList(TList *list) {
  QCorrelation::SetList(list);
  fQC  = (TProfile*)   fList->FindObject( Form("QC_%s",fList->GetName()) );
  fPQC = (TProfile2D*) fList->FindObject( Form("PQC_%s",fList->GetName()) );
}

QC_EventPlane::~QC_EventPlane() {
  if(fList) delete fList;
}

void QC_EventPlane::Results() {
  if(!fList) return;
  //std::cout << "fList found!" << std::endl;
  ClearResults();
  //std::cout << "fListResults created!" << std::endl;

  TH1D *tRAB = new TH1D(Form("RAB_%s",fList->GetName()),"RAB;rab\{2}",1,-0.5,0.5);
  fListResults->Add(tRAB);
  TH1D *tRQ = new TH1D(Form("RQ_%s",fList->GetName()),"RQ;r\{2}",1,-0.5,0.5);
  fListResults->Add(tRQ);
  //std::cout << "histograms created!" << std::endl;
  //====
  Double_t qab = fQC->GetBinContent(1);
  Double_t qab_sw = fQC->GetBinEntries(1);
  Double_t qab_sww = (fQC->GetBinSumw2())->At(1);
  Double_t qab_err = TMath::Sqrt(qab_sww)/qab_sw*fQC->GetBinError(1);
  Double_t rab2 = TMath::Sqrt( qab );
  Double_t rab2_err = TMath::Sqrt( qab_err );
  tRAB->SetBinContent(1,rab2);
  tRAB->SetBinError(1,rab2_err);
  std::cout << "QC_EventPlane::Results("<< fList->GetName() <<") => Rab2  ";
  std::cout << rab2 << "(" << rab2_err << ")" << std::endl;
  //====
  Double_t coef = TMath::Sqrt(TMath::Pi())/2;
  TF1 res("res", Form("%f*x*TMath::Exp(-0.5*x*x)*(TMath::BesselI0(0.5*x*x)+TMath::BesselI1(0.5*x*x))",coef), 0, 10 );
  Double_t xx = res.GetX( rab2 );
  Double_t r2 = res.Eval( xx*TMath::Sqrt(2.0) );
  Double_t r2_err = rab2_err/TMath::Sqrt(2.0);
  //reading from published v2 paper
  r2 = 1.07e-1;
  r2_err = 3.03e-4;
  //
  tRQ->SetBinContent(1,r2);
  tRQ->SetBinError(1,r2_err);
  std::cout << "QC_EventPlane::Results("<< fList->GetName() <<") => R2  ";
  std::cout << r2 << "(" << r2_err << ")" << std::endl;
  //====
  TH1D *tD2[1000];
  TH1D *tPV2[1000];
  TH1D *tFD2[1000];
  TH1D *tFPV2[1000];
  const TArrayD *arrma = fBinMap->GetYaxis()->GetXbins();
  Int_t nma = arrma->GetSize();
  Double_t nmas[1000];
  Int_t nma0=0;
  Double_t nmas0[1000];
  for(int i=0; i!=nma; ++i) {
    nmas[i] = arrma->GetAt(i);
    if(i<36) {
      if(i%6==0) {
	nmas0[nma0] = arrma->GetAt(i);
	nma0++;
      }
    } else if(i>=78) {
      if(i%6==0) {
	nmas0[nma0] = arrma->GetAt(i);
	nma0++;
      }
    } else {
      if(i%3==0) {
	nmas0[nma0] = arrma->GetAt(i);
	nma0++;
      }
    }
  }
  for(int pt=0; pt!=fPQC->GetXaxis()->GetNbins(); ++pt) { // as good as any
    tD2[pt]  = new TH1D( Form("D2_%d_%s_UNBINNED",pt,fList->GetName()), "UNBINNED_D2;d\{2}", nma-1,nmas);
    tPV2[pt] = new TH1D( Form("PV2_%d_%s_UNBINNED",pt,fList->GetName()),"UNBINNED_PV2;v\{2}",nma-1,nmas);
    tFD2[pt]  = new TH1D( Form("D2_%d_%s_BINNED",pt,fList->GetName()), "BINNED_D2;d\{2}", nma0-1,nmas0);
    tFPV2[pt] = new TH1D( Form("PV2_%d_%s_BINNED",pt,fList->GetName()),"BINNED_PV2;v\{2}",nma0-1,nmas0);
    fListResults->Add(tD2[pt]);
    fListResults->Add(tFD2[pt]);
    fListResults->Add(tPV2[pt]);
    fListResults->Add(tFPV2[pt]);
    //UNBINNED
    for(int ma=0; ma!=fPQC->GetYaxis()->GetNbins(); ++ma) {
      Int_t binx = pt+1;
      Int_t biny = ma+1;
      Int_t bin = fPQC->GetBin(binx,biny,0);
      //std::cout << "BIN " << bin << "  BINX " << binx << "  BINY " << biny << std::endl;
      Double_t pq2pc = fPQC->GetBinContent( bin );
      Double_t pq2pc_sw = fPQC->GetBinEntries( bin );
      Double_t pq2pc_sww = (fPQC->GetBinSumw2())->At( bin );
      Double_t pq2pc_err = TMath::Sqrt(pq2pc_sww)/pq2pc_sw*fPQC->GetBinError( bin );
      tD2[pt]->SetBinContent(ma+1,pq2pc);
      tD2[pt]->SetBinError(ma+1,pq2pc_err);
      Double_t pv2 = pq2pc/r2;
      Double_t pv2_err = pq2pc_err/r2;
      tPV2[pt]->SetBinContent(ma+1,pv2);
      tPV2[pt]->SetBinError(ma+1,pv2_err);
    }
    //BINNED
    TProfile *tempBin0 = fPQC->ProfileY(Form("tempBin_PT%d",pt),pt+1,pt+1,"s");
    TProfile *tempBin = (TProfile*) tempBin0->Rebin(nma0-1,Form("tempBin_PT%d_BINNED",pt),nmas0);
    for(int ma=0; ma!=tempBin->GetXaxis()->GetNbins(); ++ma) {
      Int_t bin = ma+1;
      //std::cout << "BIN " << bin << "  BINX " << binx << "  BINY " << biny << std::endl;
      Double_t pq2pc = tempBin->GetBinContent( bin );
      Double_t pq2pc_sw = tempBin->GetBinEntries( bin );
      Double_t pq2pc_sww = (tempBin->GetBinSumw2())->At( bin );
      Double_t pq2pc_err = TMath::Sqrt(pq2pc_sww)/pq2pc_sw*tempBin->GetBinError( bin );
      tFD2[pt]->SetBinContent(ma+1,pq2pc);
      tFD2[pt]->SetBinError(ma+1,pq2pc_err);
      Double_t pv2 = pq2pc/r2;
      Double_t pv2_err = pq2pc_err/r2;
      tFPV2[pt]->SetBinContent(ma+1,pv2);
      tFPV2[pt]->SetBinError(ma+1,pv2_err);
    }
  }
}

void QC_EventPlane::Reset() {
}
//==========
void QC_EventPlane::FillCandidate(Double_t pt, Double_t ma, Double_t phi) {
  fBinMap->Fill(pt,ma);
  Double_t psi  = fQ->Psi2Pi();
  Int_t n = fQ->Order();
  Double_t cos = TMath::Cos( n*(phi-psi) );
  fPQC->Fill(pt,ma,cos);
}

void QC_EventPlane::FillEvent() {
  //===
  Double_t psiA2pi = fQa->Psi2Pi();
  Double_t psiB2pi = fQb->Psi2Pi();
  Double_t psiA = fQa->Psi();
  Double_t psiB = fQb->Psi();
  Int_t n = fQ->Order();
  Double_t cosAB2pi = TMath::Cos( n*(psiA2pi-psiB2pi) );
  Double_t cosAB    = TMath::Cos( n*(psiA-psiB) );
  Double_t cosabsAB = TMath::Cos( n*TMath::Abs(psiA-psiB) );
  fQC->Fill(0.0, cosAB2pi);
  fQC->Fill(1.0, cosAB);
  fQC->Fill(2.0, cosabsAB);

  Double_t psiR = fQ->Psi2Pi();
  fQA_Qa->Fill(fQa->Psi());
  fQA_Qb->Fill(fQb->Psi());
  fQA_Q->Fill(fQ->Psi());

  //std::cout << " QC_EventPlane => " << n << " PSIT : " << fQ->Psi() << " ";
  //std::cout << " PSIA : " << fQa->Psi() << " ";
  //std::cout << " PSIB : " << fQb->Psi() << std::endl;

}
