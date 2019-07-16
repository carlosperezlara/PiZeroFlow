#include <iostream>
#include <map>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TAxis.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TArray.h>
#include "qcQ.h"
#include "QC_Cumulants.h"

QC_Cumulants::QC_Cumulants() :
  QCorrelation() {
  fQ2 = NULL;
  fQ4 = NULL;
  fNULL = new qcQ();
  fQPC = NULL;
  fP2PC = NULL;
  fP4PC = NULL;
  fP4P2PC = NULL;
  fP4Q2PC = NULL;
  fP4Q4PC= NULL;
  fP2Q2PC = NULL;
  fP2Q4PC = NULL;
  fQA_Q2 = NULL;
  fQA_Q4 = NULL;
}

QC_Cumulants::QC_Cumulants(TString name,qcQ *q1,qcQ *q2) :
  QCorrelation() {
  fNULL = new qcQ();
  Init(name,q1,q2);
}

void QC_Cumulants::Init(TString name,qcQ *q1,qcQ *q2) {
  QCorrelation::Init(name);
  fQ2 = q1;
  fQ4 = q2;
  fQPC = new TProfile(Form("QPC_%s",name.Data()),"QPC",3,-0.5,2.5,"s"); fList->Add(fQPC);
  fQPC->Sumw2();
  fQPC->GetXaxis()->SetBinLabel(1,"<<2>>");
  fQPC->GetXaxis()->SetBinLabel(2,"<<4>>");
  fQPC->GetXaxis()->SetBinLabel(3,"<<2><4>>");
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
  fP2PC = new TProfile2D(Form("P2PC_%s",name.Data()),"P2PC",npt-1,npts,nma-1,nmas,"s"); fList->Add(fP2PC);
  fP2PC->Sumw2();
  fP4PC = new TProfile2D(Form("P4PC_%s",name.Data()),"P4PC",npt-1,npts,nma-1,nmas,"s"); fList->Add(fP4PC);
  fP4PC->Sumw2();

  fP4P2PC = new TProfile2D(Form("P4P2PC_%s",name.Data()),"P4P2PC",npt-1,npts,nma-1,nmas,"s"); fList->Add(fP4P2PC);
  fP4P2PC->Sumw2();
  fP4Q2PC = new TProfile2D(Form("P4Q2PC_%s",name.Data()),"P4Q2PC",npt-1,npts,nma-1,nmas,"s"); fList->Add(fP4Q2PC);
  fP4Q2PC->Sumw2();
  fP4Q4PC = new TProfile2D(Form("P4Q4PC_%s",name.Data()),"P4Q4PC",npt-1,npts,nma-1,nmas,"s"); fList->Add(fP4Q4PC);
  fP4Q4PC->Sumw2();
  fP2Q2PC = new TProfile2D(Form("P2Q2PC_%s",name.Data()),"P2Q2PC",npt-1,npts,nma-1,nmas,"s"); fList->Add(fP2Q2PC);
  fP2Q2PC->Sumw2();
  fP2Q4PC = new TProfile2D(Form("P2Q4PC_%s",name.Data()),"P2Q4PC",npt-1,npts,nma-1,nmas,"s"); fList->Add(fP2Q4PC);
  fP2Q4PC->Sumw2();
  //==
  fQA_Q2 = new TH1D(Form("QA_Q2_%s",name.Data()),"QA_Q2",200,-0.1,+TMath::TwoPi()); fList->Add(fQA_Q2);
  fQA_Q4 = new TH1D(Form("QA_Q4_%s",name.Data()),"QA_Q4",200,-0.1,+TMath::TwoPi()); fList->Add(fQA_Q4);
}

void QC_Cumulants::SetList(TList *list) {
  QCorrelation::SetList(list);
  fQPC    = (TProfile*)   fList->At(1);
  fP2PC   = (TProfile2D*) fList->At(2);
  fP4PC   = (TProfile2D*) fList->At(3);
  fP4P2PC = (TProfile2D*) fList->At(4);
  fP4Q2PC = (TProfile2D*) fList->At(5);
  fP4Q4PC = (TProfile2D*) fList->At(6);
  fP2Q2PC = (TProfile2D*) fList->At(7);
  fP2Q4PC = (TProfile2D*) fList->At(8);
}

QC_Cumulants::~QC_Cumulants() {
  delete fNULL;
  if(fList) delete fList;
}

void QC_Cumulants::Results() {
  if(!fList) return;
  ClearResults();

  TH1D *tR2 = new TH1D(Form("R2_%s",fList->GetName()),"R2;r\{2}",1,-0.5,0.5); fListResults->Add(tR2);
  TH1D *tR4 = new TH1D(Form("R4_%s",fList->GetName()),"R4;r\{4}",1,-0.5,0.5); fListResults->Add(tR4);
  TH1D *tQV2 = new TH1D(Form("QV2_%s",fList->GetName()),"QV2;V\{2}",1,-0.5,0.5); fListResults->Add(tQV2);
  TH1D *tQV4 = new TH1D(Form("QV4_%s",fList->GetName()),"QV4;V\{4}",1,-0.5,0.5); fListResults->Add(tQV4);
  //====
  Double_t q2pc = fQPC->GetBinContent(1);
  Double_t q2pc_sw = fQPC->GetBinEntries(1);
  Double_t q2pc_sww = (fQPC->GetBinSumw2())->At(1);
  Double_t q2pc_err = TMath::Sqrt(q2pc_sww)/q2pc_sw*fQPC->GetBinError(1);
  Double_t r2_err;
  Double_t r2 = MakeQC2( q2pc, q2pc_err, r2_err );
  tR2->SetBinContent( 1, r2 );
  tR2->SetBinError( 1, r2_err );
  std::cout << "QC_Cumulants::Results("<< fList->GetName() <<") => R2  " << r2 << "(" << r2_err << ")" << std::endl;
  if(r2>0) {
    Double_t v2_err;
    Double_t v2 = MakeV2(q2pc,q2pc,0,q2pc_err,q2pc_err,v2_err);
    tQV2->SetBinContent( 1, v2 );
    tQV2->SetBinError( 1, v2_err );
    std::cout << "QC_Cumulants::Results("<< fList->GetName() <<") => QV2 " << v2 << "(" << v2_err << ")" << std::endl;
  }
  //====
  Double_t q4pc = fQPC->GetBinContent(2);
  Double_t q4pc_sww = (fQPC->GetBinSumw2())->At(2);
  Double_t q4pc_sw = fQPC->GetBinEntries(2);
  Double_t q4pc_err = TMath::Sqrt(q4pc_sww)/q4pc_sw*fQPC->GetBinError(2);
  Double_t r4_err;
  Double_t r4 = MakeQC4( q4pc, q2pc, q2pc, q4pc_err, q2pc_err, q2pc_err, 0, 0, 0, r4_err );
  tR4->SetBinContent( 1, r4 );
  tR4->SetBinError( 1, r4_err );
  std::cout << "QC_Cumulants::Results("<< fList->GetName() <<") => R4  " << r4 << "(" << r4_err << ")" << std::endl;
  if(r4<0) {
    Double_t v4_err;
    Double_t v4 = MakeV4( r4, r4, q4pc, q4pc, q2pc, q2pc,
			  q4pc_err, q4pc_err, q2pc_err, q2pc_err,
			  0, 0, 0, 0, 0, 0, v4_err );
    tQV4->SetBinContent( 1, v4 );
    tQV4->SetBinError( 1, v4_err );
    std::cout << "QC_Cumulants::Results("<< fList->GetName() <<") => QV4 " << v4 << "(" << v4_err << ")" << std::endl;
  }
  //====
  Double_t q4pcq2pc = fQPC->GetBinContent( 3 );
  Double_t q4pcq2pc_sw = fQPC->GetBinEntries( 3 );
  Double_t cov_q4pcq2pc = q4pcq2pc_sw / q2pc_sw/ q4pc_sw * (q4pcq2pc - q2pc*q4pc) / (1-q4pcq2pc_sw/q2pc_sw/q4pc_sw);
  //====


  TH1D *tD2[1000];
  TH1D *tD4[1000];
  TH1D *tPV2[1000];
  TH1D *tPV4[1000];
  const TArrayD *arrma = fBinMap->GetYaxis()->GetXbins();
  Int_t nma = arrma->GetSize();
  Double_t nmas[1000];
  for(int i=0; i!=nma; ++i) {
    nmas[i] = arrma->GetAt(i);
  }
  for(int pt=0; pt!=fP2PC->GetXaxis()->GetNbins(); ++pt) {
    tD2[pt]  = new TH1D( Form("D2_%d_%s",pt,fList->GetName()), "D2;d\{2}", nma-1,nmas); fListResults->Add(tD2[pt]);
    tD4[pt]  = new TH1D( Form("D4_%d_%s",pt,fList->GetName()), "D4;d\{4}", nma-1,nmas); fListResults->Add(tD4[pt]);
    tPV2[pt] = new TH1D( Form("PV2_%d_%s",pt,fList->GetName()),"PV2;v\{2}",nma-1,nmas); fListResults->Add(tPV2[pt]);
    tPV4[pt] = new TH1D( Form("PV4_%d_%s",pt,fList->GetName()),"PV4;v\{4}",nma-1,nmas); fListResults->Add(tPV4[pt]);
  }
  //====
  //std::map<Int_t,qcQ>::iterator it;
  //for(it=fP2.begin(); it!=fP2.end(); ++it) {
  for(Int_t binx=1; binx<fBinMap->GetXaxis()->GetNbins()+1; ++binx) for(Int_t biny=1; biny<fBinMap->GetYaxis()->GetNbins()+1; ++biny) {
      //std::cout << " AAAA " << std::endl;
      Double_t ptval = fBinMap->GetXaxis()->GetBinCenter(binx);
      Double_t maval = fBinMap->GetYaxis()->GetBinCenter(biny);
      Int_t bin = fBinMap->FindBin(ptval,maval);
      Double_t p2pc = fP2PC->GetBinContent( bin );
      Double_t p2pc_sw = fP2PC->GetBinEntries( bin );
      Double_t p2pc_sww = (fP2PC->GetBinSumw2())->At( bin );
      Double_t p2pc_err = TMath::Sqrt(p2pc_sww)/p2pc_sw*fP2PC->GetBinError( bin );
      Double_t d2_err;
      Double_t d2 = MakeQC2( p2pc, p2pc_err, d2_err );
      Double_t p2pcq2pc = fP2Q2PC->GetBinContent( bin );
      Double_t p2pcq2pc_sw = fP2Q2PC->GetBinEntries( bin );
      Double_t cov_p2pcq2pc = p2pcq2pc_sw / q2pc_sw/ p2pc_sw * (p2pcq2pc - q2pc*p2pc) / (1-p2pcq2pc_sw/q2pc_sw/p2pc_sw);
      tD2[binx-1]->SetBinContent( biny, d2 );
      tD2[binx-1]->SetBinError( biny, d2_err );
      //std::cout << "  D2_" << ptval << " " << d2 << "(" << d2_err << ")" << std::endl;
      if(r2>0) {
	Double_t v2_err;
	Double_t v2 = MakeV2(p2pc,q2pc,cov_p2pcq2pc,p2pc_err,q2pc_err,v2_err);
	tPV2[binx-1]->SetBinContent( biny, v2 );
	tPV2[binx-1]->SetBinError( biny, v2_err );
	//std::cout << "   => V2_" << ptval << " " << v2 << "(" << v2_err << ")" << std::endl;
      }
      //====
      //std::cout << " BBBB " << std::endl;
      Double_t p4pc = fP4PC->GetBinContent( bin );
      Double_t p4pc_sww = (fP4PC->GetBinSumw2())->At(bin);
      Double_t p4pc_sw = fP4PC->GetBinEntries(bin);
      if(p4pc_sww<0||
	 TMath::Abs(p4pc_sw)<1e-10||
	 TMath::Abs(q4pc_sw)<1e-10) {
	//std::cout << " cannot process QC4 " << std::endl;
	continue;
      }
      //std::cout << " BBBB0 " << std::endl;
      Double_t p4pc_err = TMath::Sqrt(p4pc_sww)/p4pc_sw*fP4PC->GetBinError( bin );
      //std::cout << " BBBB0 " << std::endl;
      Double_t p4pcq2pc = fP4Q2PC->GetBinContent( bin );
      //std::cout << " BBBB0 " << std::endl;
      Double_t p4pcq2pc_sw = fP4Q2PC->GetBinEntries( bin );
      //std::cout << " BBBB0 " << std::endl;
      Double_t cov_p4pcq2pc = p4pcq2pc_sw / q2pc_sw/ p4pc_sw * (p4pcq2pc - q2pc*p4pc) / (1-p4pcq2pc_sw/q2pc_sw/p4pc_sw);
      //std::cout << " BBBB1 " << std::endl;
      Double_t p4pcp2pc = fP4P2PC->GetBinContent( bin );
      Double_t p4pcp2pc_sw = fP4P2PC->GetBinEntries( binx );
      Double_t cov_p4pcp2pc = p4pcp2pc_sw / p2pc_sw/ p4pc_sw * (p4pcp2pc - p2pc*p4pc) / (1-p4pcp2pc_sw/p2pc_sw/p4pc_sw);
      Double_t p4pcq4pc = fP4Q4PC->GetBinContent( bin );
      Double_t p4pcq4pc_sw = fP4Q4PC->GetBinEntries( bin );
      Double_t cov_p4pcq4pc = p4pcq4pc_sw / q4pc_sw/ p4pc_sw * (p4pcq4pc - q4pc*p4pc) / (1-p4pcq4pc_sw/q4pc_sw/p4pc_sw);
      //std::cout << " BBBB2 " << std::endl;
      Double_t p2pcq4pc = fP2Q4PC->GetBinContent( bin );
      Double_t p2pcq4pc_sw = fP2Q4PC->GetBinEntries( bin );
      Double_t cov_p2pcq4pc = p2pcq4pc_sw / q4pc_sw/ p2pc_sw * (p2pcq4pc - q4pc*p2pc) / (1-p2pcq4pc_sw/q4pc_sw/p2pc_sw);
      Double_t d4_err;
      //std::cout << " BBBB3 " << std::endl;
      Double_t d4 = MakeQC4( p4pc, p2pc, q2pc, p4pc_err, p2pc_err, q2pc_err,
			     cov_p4pcp2pc, cov_p4pcq2pc, cov_p2pcq2pc, d4_err );
      tD4[binx-1]->SetBinContent( biny, d4 );
      tD4[binx-1]->SetBinError( biny, d4_err );
      //std::cout << "  D4_" << ptval << " " << d4 << "(" << d4_err << ")" << std::endl;
      //std::cout << " CCCC " << std::endl;
      if(r4<0) {
	Double_t v4_err;
	Double_t v4 = MakeV4( d4, r4, p4pc, q4pc, p2pc, q2pc,
			      p4pc_err, q4pc_err, p2pc_err, q2pc_err,
			      cov_p4pcq4pc, cov_p4pcq2pc, cov_p2pcq4pc, cov_p4pcp2pc, cov_p2pcq2pc, cov_q4pcq2pc, v4_err );
	tPV4[binx-1]->SetBinContent( biny, v4 );
	tPV4[binx-1]->SetBinError( biny, v4_err );
	//std::cout << "   => V4_" << ptval << " " << v4 << "(" << v4_err << ")" << std::endl;
      }
      //std::cout << " DDDD " << std::endl;
    }
}

void QC_Cumulants::Reset() {
  fP2.clear();
}

void QC_Cumulants::FillCandidate(Double_t pt, Double_t ma, Double_t phi) {
  fBinMap->Fill(pt,ma);
  Int_t bin = fBinMap->FindBin(pt,ma);
  std::map<Int_t,qcQ>::iterator it = fP2.find(bin);
  if(it == fP2.end()) {
    fP2[bin] = qcQ( fQ2->Order() );
  }
  fP2[bin].Fill( phi );
}

void QC_Cumulants::FillEvent() {
  Double_t q2pcw;
  Double_t q4pcw;
  Double_t q2pc = Make2PC( fQ2, fQ2, fQ2, q2pcw );
  Double_t q4pc = Make4PC( fQ2, fQ2, fQ2, fQ4, fQ4, q4pcw );
  if( !TMath::AreEqualAbs( q2pc, 37767, 1 ) ) {
    fQPC->Fill( 0.0, q2pc, q2pcw );
  } else {
    std::cout << "QC_Cumulants::MakePC ERROR IN Q2PC! " << q2pc << std::endl;
  }
  if( !TMath::AreEqualAbs( q4pc, 37767, 1 ) ) {
    fQPC->Fill( 1.0, q4pc, q4pcw );
    fQPC->Fill( 2.0, q2pc*q4pc, q2pcw*q4pcw );
  } else {
    std::cout << "QC_Cumulants::MakePC ERROR IN Q4PC! " << q4pc << std::endl;
  }

  std::map<Int_t,qcQ>::iterator it;
  for(it=fP2.begin(); it!=fP2.end(); ++it) {
    Int_t bin = it->first;
    qcQ *pn = &it->second;
    Int_t binx, biny, binz;
    fP2PC->GetBinXYZ(bin,binx,biny,binz);
    Double_t ptval = fP2PC->GetXaxis()->GetBinCenter(binx);
    Double_t maval = fP2PC->GetYaxis()->GetBinCenter(biny);
    Double_t p2pcw, p4pcw;
    Double_t p2pc = Make2PC( pn, fQ2, fNULL, p2pcw );
    Double_t p4pc = Make4PC( pn, fQ2, fNULL, fQ4, fNULL, p4pcw );
    if( !TMath::AreEqualAbs( p2pc, 37767, 1 ) ) {
      fP2PC->Fill( ptval, maval, p2pc, p2pcw );
      fP2Q2PC->Fill( ptval, maval, p2pc*q2pc, p2pcw*q2pcw );
    } else
      std::cout << "qcAnalysisQC::MakePC ERROR IN P2PC! " << ptval << ":" << p2pc << std::endl;
    if( !TMath::AreEqualAbs( p4pc, 37767, 1 ) ) {
      fP4PC->Fill( ptval, maval, p4pc, p4pcw );
      fP4P2PC->Fill( ptval, maval, p2pc*p4pc, p2pcw*p4pcw );
      fP4Q2PC->Fill( ptval, maval, q2pc*p4pc, q2pcw*p4pcw );
      fP4Q4PC->Fill( ptval, maval, q4pc*p4pc, q4pcw*p4pcw );
      fP2Q4PC->Fill( ptval, maval, q4pc*p2pc, q4pcw*p2pcw );
    } else 
      std::cout << "qcAnalysisQC::MakePC ERROR IN P4PC! " << ptval << ":" << p4pc << std::endl;
  }  
  //==
  Double_t psi2 = fQ2->Psi();
  Double_t psi4 = fQ4->Psi();
  fQA_Q2->Fill(psi2);
  fQA_Q4->Fill(psi4);
}
