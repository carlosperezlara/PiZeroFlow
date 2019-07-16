#include <iostream>
#include <map>
#include <TMath.h>
#include <TH2D.h>
#include <TAxis.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TArray.h>
#include "qcQ.h"
#include "QC_ScalarProduct.h"

QC_ScalarProduct::QC_ScalarProduct() :
  QCorrelation() {
  fQa = NULL;
  fQb = NULL;
  fNULL = new qcQ();
}

QC_ScalarProduct::QC_ScalarProduct(TString name,qcQ *q1,qcQ *q2) :
  QCorrelation() {
  fNULL = new qcQ();
  Init(name,q1,q2);
}

void QC_ScalarProduct::Init(TString name,qcQ *q1,qcQ *q2) {
  QCorrelation::Init(name);
  fQa = q1;
  fQb = q2;
  fQPC = new TProfile(Form("QPC_%s",name.Data()),"QPC",3,-0.5,2.5,"s"); fList->Add(fQPC);
  fQPC->Sumw2();
  fQPC->GetXaxis()->SetBinLabel(1,"<<2a>>");
  fQPC->GetXaxis()->SetBinLabel(2,"<<2b>>");
  fQPC->GetXaxis()->SetBinLabel(3,"<<2a><2b>>");
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
  fPQa_2PC = new TProfile2D(Form("PQa_2PC_%s",name.Data()),"PQa_2PC",npt-1,npts,nma-1,nmas,"s"); fList->Add(fPQa_2PC);
  fPQa_2PC->Sumw2();
  fPQb_2PC = new TProfile2D(Form("PQb_2PC_%s",name.Data()),"PQb_2PC",npt-1,npts,nma-1,nmas,"s"); fList->Add(fPQb_2PC);
  fPQb_2PC->Sumw2();
  fPQaxQa_2PC  = new TProfile2D(Form("PQaxQa_2PC_%s",name.Data()), "PQaxQa_2PC",npt-1,npts,nma-1,nmas,"s");  fList->Add(fPQaxQa_2PC);
  fPQaxQa_2PC->Sumw2();
  fPQbxQb_2PC  = new TProfile2D(Form("PQbxQb_2PC_%s",name.Data()), "PQbxQb_2PC",npt-1,npts,nma-1,nmas,"s");  fList->Add(fPQbxQb_2PC);
  fPQbxQb_2PC->Sumw2();
  //==
  fPQaxQb_2PC  = new TProfile2D(Form("PQaxQb_2PC_%s",name.Data()), "PQaxQb_2PC",npt-1,npts,nma-1,nmas,"s");  fList->Add(fPQaxQb_2PC);
  fPQaxQb_2PC->Sumw2();
  fPQbxQa_2PC  = new TProfile2D(Form("PQbxQa_2PC_%s",name.Data()), "PQbxQa_2PC",npt-1,npts,nma-1,nmas,"s");  fList->Add(fPQbxQa_2PC);
  fPQbxQa_2PC->Sumw2();
  fPQaxPQb_2PC = new TProfile2D(Form("PQaxPQb_2PC_%s",name.Data()),"PQaxPQb_2PC",npt-1,npts,nma-1,nmas,"s"); fList->Add(fPQaxPQb_2PC);
  fPQaxPQb_2PC->Sumw2();
  //==
  fQA_Qa = new TH1D(Form("QA_Qa_%s",name.Data()),"QA_Qa",200,-0.1,+TMath::TwoPi()); fList->Add(fQA_Qa);
  fQA_Qb = new TH1D(Form("QA_Qb_%s",name.Data()),"QA_Qb",200,-0.1,+TMath::TwoPi()); fList->Add(fQA_Qb);
}

void QC_ScalarProduct::SetList(TList *list) {
  QCorrelation::SetList(list);
  fQPC         = (TProfile*)   fList->FindObject( Form("QPC_%s",fList->GetName()) );
  fPQa_2PC     = (TProfile2D*) fList->FindObject( Form("PQa_2PC_%s",fList->GetName()) );
  fPQb_2PC     = (TProfile2D*) fList->FindObject( Form("PQb_2PC_%s",fList->GetName()) );
  fPQaxQa_2PC  = (TProfile2D*) fList->FindObject( Form("PQaxQa_2PC_%s",fList->GetName()) );
  fPQbxQb_2PC  = (TProfile2D*) fList->FindObject( Form("PQbxQb_2PC_%s",fList->GetName()) );
  fPQaxQb_2PC  = (TProfile2D*) fList->FindObject( Form("PQaxQb_2PC_%s",fList->GetName()) );
  fPQbxQa_2PC  = (TProfile2D*) fList->FindObject( Form("PQbxQa_2PC_%s",fList->GetName()) );
  fPQaxPQb_2PC = (TProfile2D*) fList->FindObject( Form("PQaxPQb_2PC_%s",fList->GetName()) );
}

QC_ScalarProduct::~QC_ScalarProduct() {
  delete fNULL;
  if(fList) delete fList;
}

void QC_ScalarProduct::Results() {
  if(!fList) return;
  ClearResults();

  TH1D *tRa2 = new TH1D(Form("Ra2_%s",fList->GetName()),"Ra2;r\{2}",1,-0.5,0.5); fListResults->Add(tRa2);
  TH1D *tRb2 = new TH1D(Form("Rb2_%s",fList->GetName()),"Rb2;r\{4}",1,-0.5,0.5); fListResults->Add(tRb2);
  TH1D *tQVa2 = new TH1D(Form("QVa2_%s",fList->GetName()),"QVa2;V\{2}",1,-0.5,0.5); fListResults->Add(tQVa2);
  TH1D *tQVb2 = new TH1D(Form("QVb2_%s",fList->GetName()),"QVb2;V\{4}",1,-0.5,0.5); fListResults->Add(tQVb2);
  //====
  Double_t qa2pc = fQPC->GetBinContent(1);
  Double_t qa2pc_sw = fQPC->GetBinEntries(1);
  Double_t qa2pc_sww = (fQPC->GetBinSumw2())->At(1);
  Double_t qa2pc_err = TMath::Sqrt(qa2pc_sww)/qa2pc_sw*fQPC->GetBinError(1);
  Double_t ra2_err;
  Double_t ra2 = MakeQC2( qa2pc, qa2pc_err, ra2_err );
  tRa2->SetBinContent( 1, ra2 );
  tRa2->SetBinError( 1, ra2_err );
  std::cout << "QC_ScalarProduct::Results("<< fList->GetName() <<") => Ra2  " << ra2 << "(" << ra2_err << ")" << std::endl;
  if(ra2>0) {
    Double_t va2_err;
    Double_t va2 = MakeV2(qa2pc,qa2pc,0,qa2pc_err,qa2pc_err,va2_err);
    tQVa2->SetBinContent( 1, va2 );
    tQVa2->SetBinError( 1, va2_err );
    std::cout << "QC_ScalarProduct::Results("<< fList->GetName() <<") => QVa2 " << va2 << "(" << va2_err << ")" << std::endl;
  }
  //====
  Double_t qb2pc = fQPC->GetBinContent(2);
  Double_t qb2pc_sw = fQPC->GetBinEntries(2);
  Double_t qb2pc_sww = (fQPC->GetBinSumw2())->At(2);
  Double_t qb2pc_err = TMath::Sqrt(qb2pc_sww)/qb2pc_sw*fQPC->GetBinError(2);
  Double_t rb2_err;
  Double_t rb2 = MakeQC2( qb2pc, qb2pc_err, rb2_err );
  tRb2->SetBinContent( 1, rb2 );
  tRb2->SetBinError( 1, rb2_err );
  std::cout << "QC_ScalarProduct::Results("<< fList->GetName() <<") => Rb2  " << rb2 << "(" << rb2_err << ")" << std::endl;
  if(rb2>0) {
    Double_t vb2_err;
    Double_t vb2 = MakeV2(qb2pc,qb2pc,0,qb2pc_err,qb2pc_err,vb2_err);
    tQVb2->SetBinContent( 1, vb2 );
    tQVb2->SetBinError( 1, vb2_err );
    std::cout << "QC_ScalarProduct::Results("<< fList->GetName() <<") => QVb2 " << vb2 << "(" << vb2_err << ")" << std::endl;
  }
  //====
  Double_t qaqb2pc = fQPC->GetBinContent(3);
  Double_t qaqb2pc_sw = fQPC->GetBinEntries(3);
  Double_t cov_qaqb2pc = qaqb2pc_sw / qa2pc_sw/ qb2pc_sw * (qaqb2pc - qa2pc*qb2pc) / (1-qaqb2pc_sw/qa2pc_sw/qb2pc_sw);
  //====
  TH1D *tDa2[1000];
  TH1D *tDb2[1000];
  TH1D *tPVa2[1000];
  TH1D *tPVb2[1000];
  TH1D *tPV2[1000];
  TH1D *tFDa2[1000];
  TH1D *tFDb2[1000];
  TH1D *tFPVa2[1000];
  TH1D *tFPVb2[1000];
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
  for(int pt=0; pt!=fPQa_2PC->GetXaxis()->GetNbins(); ++pt) { // as good as any
    tDa2[pt]  = new TH1D( Form("Da2_%d_%s_UNBINNED",pt,fList->GetName()), "UNBINNED_Da2;d\{2}", nma-1,nmas);
    tDb2[pt]  = new TH1D( Form("Db2_%d_%s_UNBINNED",pt,fList->GetName()), "UNBINNED_Db2;d\{4}", nma-1,nmas);
    tPVa2[pt] = new TH1D( Form("PVa2_%d_%s_UNBINNED",pt,fList->GetName()),"UNBINNED_PVa2;v\{2}",nma-1,nmas);
    tPVb2[pt] = new TH1D( Form("PVb2_%d_%s_UNBINNED",pt,fList->GetName()),"UNBINNED_PVb2;v\{4}",nma-1,nmas);
    tPV2[pt] = new TH1D( Form("PV2_%d_%s_UNBINNED",pt,fList->GetName()),"UNBINNED_PV2;v\{4}",nma-1,nmas);
    tFDa2[pt]  = new TH1D( Form("Da2_%d_%s_BINNED",pt,fList->GetName()), "BINNED_Da2;d\{2}", nma0-1,nmas0);
    tFDb2[pt]  = new TH1D( Form("Db2_%d_%s_BINNED",pt,fList->GetName()), "BINNED_Db2;d\{4}", nma0-1,nmas0);
    tFPVa2[pt] = new TH1D( Form("PVa2_%d_%s_BINNED",pt,fList->GetName()),"BINNED_PVa2;v\{2}",nma0-1,nmas0);
    tFPVb2[pt] = new TH1D( Form("PVb2_%d_%s_BINNED",pt,fList->GetName()),"BINNED_PVb2;v\{4}",nma0-1,nmas0);
    tFPV2[pt] = new TH1D( Form("PV2_%d_%s_BINNED",pt,fList->GetName()),"BINNED_PV2;v\{4}",nma0-1,nmas0);
    fListResults->Add(tDa2[pt]);
    fListResults->Add(tFDa2[pt]);
    fListResults->Add(tDb2[pt]);
    fListResults->Add(tFDb2[pt]);
    fListResults->Add(tPVa2[pt]);
    fListResults->Add(tFPVa2[pt]);
    fListResults->Add(tPVb2[pt]);
    fListResults->Add(tFPVb2[pt]);
    fListResults->Add(tPV2[pt]);
    fListResults->Add(tFPV2[pt]);
    //UNBINNED
    for(int ma=0; ma!=fPQa_2PC->GetYaxis()->GetNbins(); ++ma) {
      Int_t binx = pt+1;
      Int_t biny = ma+1;
      Int_t bin = fPQa_2PC->GetBin(binx,biny,0);
      //std::cout << "binx " << binx << " biny " << biny << " => bin " << bin << std::endl;
      //====
      Double_t pqa2pc = fPQa_2PC->GetBinContent( bin );
      Double_t pqa2pc_sw = fPQa_2PC->GetBinEntries( bin );
      Double_t pqa2pc_sww = (fPQa_2PC->GetBinSumw2())->At( bin );
      Double_t pqa2pc_err = TMath::Sqrt(pqa2pc_sww)/pqa2pc_sw*fPQa_2PC->GetBinError( bin );
      Double_t da2_err;
      Double_t da2 = MakeQC2( pqa2pc, pqa2pc_err, da2_err );
      Double_t pqaqa2pc = fPQaxQa_2PC->GetBinContent( bin );
      Double_t pqaqa2pc_sw = fPQaxQa_2PC->GetBinEntries( bin );
      Double_t cov_pqaqa2pc = pqaqa2pc_sw / qa2pc_sw/ pqa2pc_sw * (pqaqa2pc - qa2pc*pqa2pc) / (1-pqaqa2pc_sw/qa2pc_sw/pqa2pc_sw);
      tDa2[pt]->SetBinContent( ma+1, da2 );
      tDa2[pt]->SetBinError( ma+1, da2_err );
      Double_t va2_err, vb2_err;
      Double_t va2=-999;
      Double_t vb2=-999;
      if(ra2>0) {
	va2 = MakeV2(pqa2pc,qa2pc,cov_pqaqa2pc,pqa2pc_err,qa2pc_err,va2_err);
	tPVa2[pt]->SetBinContent( ma+1, va2 );
	tPVa2[pt]->SetBinError( ma+1, va2_err );
      }
      //====
      Double_t pqb2pc = fPQb_2PC->GetBinContent( bin );
      Double_t pqb2pc_sw = fPQb_2PC->GetBinEntries( bin );
      Double_t pqb2pc_sww = (fPQb_2PC->GetBinSumw2())->At( bin );
      Double_t pqb2pc_err = TMath::Sqrt(pqb2pc_sww)/pqb2pc_sw*fPQb_2PC->GetBinError( bin );
      Double_t db2_err;
      Double_t db2 = MakeQC2( pqb2pc, pqb2pc_err, db2_err );
      Double_t pqbqb2pc = fPQbxQb_2PC->GetBinContent( bin );
      Double_t pqbqb2pc_sw = fPQbxQb_2PC->GetBinEntries( bin );
      Double_t cov_pqbqb2pc = pqbqb2pc_sw / qb2pc_sw/ pqb2pc_sw * (pqbqb2pc - qb2pc*pqb2pc) / (1-pqbqb2pc_sw/qb2pc_sw/pqb2pc_sw);
      tDb2[pt]->SetBinContent( ma+1, db2 );
      tDb2[pt]->SetBinError( ma+1, db2_err );
      if(rb2>0) {
	vb2 = MakeV2(pqb2pc,qb2pc,cov_pqbqb2pc,pqb2pc_err,qb2pc_err,vb2_err);
	tPVb2[pt]->SetBinContent( ma+1, vb2 );
	tPVb2[pt]->SetBinError( ma+1, vb2_err );
      }
      //====
      if(va2>-1 && vb2>-1) {
	Double_t v2 = TMath::Sqrt( TMath::Abs( va2*vb2 ) );
	if(va2<0&&vb2<0) {
	  v2 = -v2;
	}
	Double_t a = pqa2pc;
	Double_t b = pqb2pc;
	Double_t c = TMath::Sqrt(qa2pc);
	Double_t d = TMath::Sqrt(qb2pc);
	Double_t f = TMath::Sqrt(a*b/c/d);
	Double_t ea = pqa2pc_err;
	Double_t eb = pqb2pc_err;
	Double_t ec = qa2pc_err;
	Double_t ed = qb2pc_err;
	Double_t eab = 0;
	Double_t eac = pqaqa2pc;
	Double_t ead = 0;
	Double_t ebc = 0;
	Double_t ebd = pqbqb2pc;
	Double_t ecd = cov_qaqb2pc;
	//Double_t v2_err = TMath::Sqrt( f*f/4*(ea*ea/a/a + eb*eb/b/b + ec*ec/c/c ) + eta/f/c/c*( c*eab - a*ebc - b*eac ) ); // three terms error
	//Double_t v2_err = v2*TMath::Sqrt( va2_err*va2_err/va2/va2 + vb2_err*vb2_err/vb2/vb2 + 2*0/va2/vb2  ); // two terms error (no cov)
	Double_t v2_err = TMath::Sqrt( f*f/4*( ea*ea/a/a + eb*eb/b/b + ec*ec/c/c + ed*ed/d/d )
				       /*+ TMath::Abs(f/2*(eab/a/b - eac/a/c - ead/a/d - ebc/b/c - ebd/b/d +ecd/c/d))*/ );
	tPV2[pt]->SetBinContent( ma+1, v2 );
	tPV2[pt]->SetBinError( ma+1, v2_err );
      }
    }
    //BINNED
    TProfile *t0PQa = fPQa_2PC->ProfileY(Form("t0PQa_PT%d",pt),pt+1,pt+1,"s");
    TProfile *t0PQb = fPQb_2PC->ProfileY(Form("t0PQb_PT%d",pt),pt+1,pt+1,"s");
    TProfile *t0PQaQa = fPQaxQa_2PC->ProfileY(Form("t0PQaQa_PT%d",pt),pt+1,pt+1,"s");
    TProfile *t0PQbQb = fPQbxQb_2PC->ProfileY(Form("t0PQbQb_PT%d",pt),pt+1,pt+1,"s");
    TProfile *tPQa = (TProfile*) t0PQa->Rebin(nma0-1,Form("tPQa_PT%d_BINNED",pt),nmas0);
    TProfile *tPQb = (TProfile*) t0PQb->Rebin(nma0-1,Form("tPQb_PT%d_BINNED",pt),nmas0);
    TProfile *tPQaQa = (TProfile*) t0PQaQa->Rebin(nma0-1,Form("tPQaQa_PT%d_BINNED",pt),nmas0);
    TProfile *tPQbQb = (TProfile*) t0PQbQb->Rebin(nma0-1,Form("tPQbQb_PT%d_BINNED",pt),nmas0);
    for(int ma=0; ma!=tPQa->GetXaxis()->GetNbins(); ++ma) {
      Int_t bin = ma+1;
      //====
      Double_t pqa2pc = tPQa->GetBinContent( bin );
      Double_t pqa2pc_sw = tPQa->GetBinEntries( bin );
      Double_t pqa2pc_sww = (tPQa->GetBinSumw2())->At( bin );
      Double_t pqa2pc_err = TMath::Sqrt(pqa2pc_sww)/pqa2pc_sw*tPQa->GetBinError( bin );
      Double_t da2_err;
      Double_t da2 = MakeQC2( pqa2pc, pqa2pc_err, da2_err );
      Double_t pqaqa2pc = tPQaQa->GetBinContent( bin );
      Double_t pqaqa2pc_sw = tPQaQa->GetBinEntries( bin );
      Double_t cov_pqaqa2pc = pqaqa2pc_sw / qa2pc_sw/ pqa2pc_sw * (pqaqa2pc - qa2pc*pqa2pc) / (1-pqaqa2pc_sw/qa2pc_sw/pqa2pc_sw);
      tFDa2[pt]->SetBinContent( ma+1, da2 );
      tFDa2[pt]->SetBinError( ma+1, da2_err );
      Double_t va2_err, vb2_err;
      Double_t va2=-999;
      Double_t vb2=-999;
      if(ra2>0) {
        va2 = MakeV2(pqa2pc,qa2pc,cov_pqaqa2pc,pqa2pc_err,qa2pc_err,va2_err);
        tFPVa2[pt]->SetBinContent( ma+1, va2 );
        tFPVa2[pt]->SetBinError( ma+1, va2_err );
      }
      //====
      Double_t pqb2pc = tPQb->GetBinContent( bin );
      Double_t pqb2pc_sw = tPQb->GetBinEntries( bin );
      Double_t pqb2pc_sww = (tPQb->GetBinSumw2())->At( bin );
      Double_t pqb2pc_err = TMath::Sqrt(pqb2pc_sww)/pqb2pc_sw*tPQb->GetBinError( bin );
      Double_t db2_err;
      Double_t db2 = MakeQC2( pqb2pc, pqb2pc_err, db2_err );
      Double_t pqbqb2pc = tPQbQb->GetBinContent( bin );
      Double_t pqbqb2pc_sw = tPQbQb->GetBinEntries( bin );
      Double_t cov_pqbqb2pc = pqbqb2pc_sw / qb2pc_sw/ pqb2pc_sw * (pqbqb2pc - qb2pc*pqb2pc) / (1-pqbqb2pc_sw/qb2pc_sw/pqb2pc_sw);
      tFDb2[pt]->SetBinContent( ma+1, db2 );
      tFDb2[pt]->SetBinError( ma+1, db2_err );
      if(rb2>0) {
	vb2 = MakeV2(pqb2pc,qb2pc,cov_pqbqb2pc,pqb2pc_err,qb2pc_err,vb2_err);
	tFPVb2[pt]->SetBinContent( ma+1, vb2 );
	tFPVb2[pt]->SetBinError( ma+1, vb2_err );
      }
      //====
      if(va2>-1 && vb2>-1) {
	Double_t v2 = TMath::Sqrt( TMath::Abs( va2*vb2 ) );
	if(va2<0&&vb2<0) {
	  v2 = -v2;
	}
	Double_t a = pqa2pc;
	Double_t b = pqb2pc;
	Double_t c = TMath::Sqrt(qa2pc);
	Double_t d = TMath::Sqrt(qb2pc);
	Double_t f = TMath::Sqrt(a*b/c/d);
	Double_t ea = pqa2pc_err;
	Double_t eb = pqb2pc_err;
	Double_t ec = qa2pc_err;
	Double_t ed = qb2pc_err;
	Double_t eab = 0;
	Double_t eac = pqaqa2pc;
	Double_t ead = 0;
	Double_t ebc = 0;
	Double_t ebd = pqbqb2pc;
	Double_t ecd = cov_qaqb2pc;
	//Double_t v2_err = TMath::Sqrt( f*f/4*(ea*ea/a/a + eb*eb/b/b + ec*ec/c/c ) + eta/f/c/c*( c*eab - a*ebc - b*eac ) ); // three terms error
	//Double_t v2_err = v2*TMath::Sqrt( va2_err*va2_err/va2/va2 + vb2_err*vb2_err/vb2/vb2 + 2*0/va2/vb2  ); // two terms error (no cov)
	Double_t v2_err = TMath::Sqrt( f*f/4*( ea*ea/a/a + eb*eb/b/b + ec*ec/c/c + ed*ed/d/d )
				       /*+ TMath::Abs(f/2*(eab/a/b - eac/a/c - ead/a/d - ebc/b/c - ebd/b/d +ecd/c/d))*/ );
	tFPV2[pt]->SetBinContent( ma+1, v2 );
	tFPV2[pt]->SetBinError( ma+1, v2_err );
      }
    }
  }
  
}
//==========
void QC_ScalarProduct::Reset() {
  fP2.clear();
}
//==========
void QC_ScalarProduct::FillCandidate(Double_t pt, Double_t ma, Double_t phi) {
  fBinMap->Fill(pt,ma);
  Int_t bin = fBinMap->FindBin(pt,ma);
  std::map<Int_t,qcQ>::iterator it = fP2.find(bin);
  if(it == fP2.end()) {
    fP2[bin] = qcQ( fQa->Order() );
  }
  fP2[bin].Fill( phi );
}
//==========
void QC_ScalarProduct::FillEvent() {
  if(fQa->M()<2) return; // rejecting questionable
  if(fQb->M()<2) return; // events

  Double_t qa2pcw;
  Double_t qb2pcw;
  Double_t qa2pc = Make2PC( fQa, fQa, fQa, qa2pcw );
  Double_t qb2pc = Make2PC( fQb, fQb, fQb, qb2pcw );
  if( !TMath::AreEqualAbs( qa2pc, 37767, 1 ) ) {
    fQPC->Fill( 0.0, qa2pc, qa2pcw );
  } else {
    std::cout << "QC_ScalarProduct::MakePC ERROR IN Qa2PC! " << qa2pc << std::endl;
    std::cout << fQPC->GetName() << std::endl;
    return;
  }
  if( !TMath::AreEqualAbs( qb2pc, 37767, 1 ) ) {
    fQPC->Fill( 1.0, qb2pc, qb2pcw );
  } else {
    std::cout << "QC_ScalarProduct::MakePC ERROR IN Qb2PC! " << qb2pc << std::endl;
    std::cout << fQPC->GetName() << std::endl;
    return;
  }
  if( !TMath::AreEqualAbs( qb2pc, 37767, 1 ) && !TMath::AreEqualAbs( qb2pc, 37767, 1 ) ) {
    fQPC->Fill( 2.0, qa2pc*qb2pc, qa2pcw*qb2pcw );
  }
  std::map<Int_t,qcQ>::iterator it;
  for(it=fP2.begin(); it!=fP2.end(); ++it) {
    Int_t bin = it->first;
    qcQ *pn = &it->second;
    Int_t binx, biny, binz;
    fPQa_2PC->GetBinXYZ(bin,binx,biny,binz);
    Double_t ptval = fPQa_2PC->GetXaxis()->GetBinCenter(binx);
    Double_t maval = fPQa_2PC->GetYaxis()->GetBinCenter(biny);
    Double_t pqa2pcw, pqb2pcw;
    Double_t pqa2pc = Make2PC( pn, fQa, fNULL, pqa2pcw );
    Double_t pqb2pc = Make2PC( pn, fQb, fNULL, pqb2pcw );
    if( !TMath::AreEqualAbs( pqa2pc, 37767, 1 ) ) {
      fPQa_2PC->Fill( ptval, maval, pqa2pc, pqa2pcw );
      fPQaxQa_2PC->Fill( ptval, maval, pqa2pc*qa2pc, pqa2pcw*qa2pcw );
    } else {
      std::cout << "qcAnalysisQC::MakePC ERROR IN P2PC! " << ptval << ":" << pqa2pc << std::endl;
      std::cout << fPQa_2PC->GetName() << std::endl;
    }
    if( !TMath::AreEqualAbs( pqb2pc, 37767, 1 ) ) {
      fPQb_2PC->Fill( ptval, maval, pqb2pc, pqb2pcw );
      fPQbxQb_2PC->Fill( ptval, maval, pqb2pc*qb2pc, pqb2pcw*qb2pcw );
    } else {
      std::cout << "qcAnalysisQC::MakePC ERROR IN P2PC! " << ptval << ":" << pqb2pc << std::endl;
      std::cout << fPQb_2PC->GetName() << std::endl;
    }
  }
  //==
  Double_t psiA = fQa->Psi();
  Double_t psiB = fQb->Psi();
  fQA_Qa->Fill(psiA);
  fQA_Qb->Fill(psiB);
}
