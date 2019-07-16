#include <TString.h>
#include <TFile.h>
#include <TList.h>
#include <TMath.h>
#include <TH2D.h>
#include "qcQ.h"
#include "QCorrelation.h"

QCorrelation::QCorrelation() {
  //default binning
  /*
  Double_t ptbins[27] = {0.5, 0.6, 0.7, 0.8, 0.9,
			 1.0, 1.1, 1.2, 1.3, 1.4,
			 1.5, 1.6, 1.7, 1.8, 2.0,
			 2.2, 2.4, 2.6, 2.8, 3.0,
			 3.5, 4.0, 5.0, 6.0, 8.0,
			 10., 20.};
  }
  */
  const int nptbins = 20;
  Double_t ptbins[21] = {0.40, 0.60, 0.80, 1.00, 1.20,
			 1.40, 1.60, 1.80, 2.00, 2.20,
			 2.40, 2.75, 3.00, 3.25, 3.50,
			 4.00, 5.00, 6.00, 8.00, 10.0,
			 20.0};
  Double_t mabins[111];
  std::vector<Double_t> cfg;
  for(int i=0; i!=111; ++i) {
    fMaList.push_back( 0.040 + i*0.002 );
    cfg.push_back( 1 );
  }
  for(int i=0; i!=nptbins+1; ++i) {
    fPtList.push_back( ptbins[i] );
    if(i<nptbins)
      fNewMaList.push_back( cfg );
  }
  //
  fList = NULL;
  fListResults = NULL;
  fBinMap = NULL;
}
//=====
void QCorrelation::SetList(TList *list) {
  if(fList) delete fList;
  fList = list;
  fBinMap = (TH2D*) fList->At(0);
}
//=====
void QCorrelation::Init(TString name) {
  fList = new TList();
  fList->SetOwner();
  fList->SetName( name.Data() );
  Int_t npt = fPtList.size();
  Double_t npts[1000];
  for(int i=0; i!=npt; ++i) {
    npts[i] = fPtList.at(i);
  }
  Int_t nma = fMaList.size();
  Double_t nmas[1000];
  for(int i=0; i!=nma; ++i) {
    nmas[i] = fMaList.at(i);
  }
  std::cout << npt << " " << nma << std::endl;
  fBinMap = new TH2D(Form("binmap_%s",name.Data()),"binmap;pt,mass",npt-1,npts,nma-1,nmas);
  fList->Add(fBinMap);
}
//=====
void QCorrelation::ClearResults() {
  if(!fList) return;
  if(fListResults) delete fListResults;
  fListResults = new TList();
  fListResults->SetName("results");
  fListResults->SetOwner();
  fList->Add(fListResults);
}
//=====
QCorrelation::~QCorrelation() {
  if(fList) delete fList;
  if(fListResults) delete fListResults;
}
//=====
void QCorrelation::SaveResults(TString name) {
  TFile *file = new TFile(name.Data(),"RECREATE");
  if(fListResults&&fList) {
    fListResults->Write(fList->GetName(),TObject::kSingleKey);
  }
  file->Close();
  delete file;
  return;
}
//=====
void QCorrelation::Save(TString name) {
  TFile *file = new TFile(name.Data(),"RECREATE");
  if(fList) {
    fList->Write(fList->GetName(),TObject::kSingleKey);
  }
  file->Close();
  delete file;
  return;
}
//========
Double_t QCorrelation::Make2PC( qcQ *a, qcQ *b, qcQ *ab, Double_t &w ) {
  // maker
  if(!a || !b || !ab) {
    std::cout << "QCorrelation::Make2PC ERROR null pointers " << std::endl;
    return 37767;
  }
  if(a->M()<1) {
    std::cout << "QCorrelation::Make2PC ERROR MultiplictyA " << std::endl;
    return 37767;
  }
  if(b->M()<1) {
    std::cout << "QCorrelation::Make2PC ERROR MultiplictyB " << std::endl;
    return 37767;
  }
  w = a->M()*b->M() - ab->M();
  if( w < 1 ) {
    std::cout << "QCorrelation::Make2PC ERROR Weight " << std::endl;
    return 37767;
  }
  Double_t num = (a->X()*b->X() + a->Y()*b->Y() - ab->M());
  return num/w;
}
//========
Double_t QCorrelation::Make4PC( qcQ *a, qcQ *b, qcQ *ab, qcQ *b2, qcQ *ab2, Double_t &w ) {
  // maker
  if(!a || !b || !ab || !ab2) {
    std::cout << "QCorrelation::Make4PC ERROR null pointers " << std::endl;
    return 37767;
  }
  if(a->M()<1) {
    std::cout << "QCorrelation::Make4PC ERROR MultiplicityA " << std::endl;
    return 37767;
  }
  if(b->M()<4) {
    std::cout << "QCorrelation::Make4PC ERROR MultiplicityB " << std::endl;
    return 37767;
  }
  w = (a->M()*b->M() - 3*ab->M())*(b->M()-2)*(b->M()-3);
  if( w < 1 ) {
    std::cout << "QCorrelation::Make4PC ERROR Weight " << std::endl;
    return 37767;
  }
  Double_t num = 0.0;
  num += (a->X()*b->X() + a->Y()*b->Y())*b->ModulusSquared(); // 1
  num -= (ab2->X()*b->X() + ab->Y()*b->Y())*b->X() - (ab2->X()*b->Y() + ab2->Y()*b->X())*b->Y(); // 2
  num -= (b->X()*b2->X() + b->Y()*b2->Y())*a->X() - (a->X()*b2->Y() + a->Y()*b2->X())*a->Y(); // 3
  num += ab2->X()*b2->X() + ab2->Y()*b2->Y(); // 4
  num -= 2*b->M()*(a->X()*b->X() + a->Y()*b->Y()); // 5
  num -= 2*ab->M()*b->ModulusSquared(); // 6
  num += 7*(ab->X()*b->X() + ab->Y()*b->Y()); // 7
  num -= b->X()*ab->X() + b->Y()*ab->Y(); // 8
  num += 2*(a->X()*b->X() + a->Y()*b->Y()); // 9
  num += 2*ab->M()*b->M(); // 10
  num -= 6*ab->M(); // 11
  return num/w;
}
//========
Double_t QCorrelation::MakeQC2( Double_t p2, Double_t er, Double_t &err ) {
  err = er;
  return p2;
}
//========
Double_t QCorrelation::MakeQC4( Double_t p4, Double_t p2, Double_t q2,
				 Double_t p4_er, Double_t p2_er, Double_t q2_er,
				 Double_t c_p4p2, Double_t c_p4q2, Double_t c_p2q2,
				 Double_t &err) {
  err= ( +4.0*p2*p2*q2_er*q2_er
	 +4.0*q2*q2*p2_er*p2_er
	 +p4_er*p4_er
	 +8.0*q2*p2*c_p2q2
	 -4.0*p2*c_p4q2
	 -4.0*q2*c_p4p2 );
  err = TMath::Sqrt(err);
  return p4 - 2*p2*q2;
}
//========
Double_t QCorrelation::MakeV2( Double_t p2, Double_t q2, Double_t c_p2q2,
				Double_t p2_err, Double_t q2_err,
				Double_t &err ) {
  Double_t v = p2 / TMath::Sqrt(q2);
  err = 0.25/q2/q2/q2*(p2*p2*q2_err*q2_err + 4*q2*q2*p2_err*p2_err - 4*q2*p2*c_p2q2);
  err = TMath::Sqrt(err);
  return v;
}
//========
Double_t QCorrelation::MakeV4( Double_t d4, Double_t r4,
				Double_t p4, Double_t q4, Double_t p2, Double_t q2,
				Double_t p4_er, Double_t q4_er, Double_t p2_er, Double_t q2_er,
				Double_t c_p4q4, Double_t c_p4q2, Double_t c_p2q4, Double_t c_p2p4, Double_t c_p2q2, Double_t c_q2q4,
				Double_t &err ) {
  Double_t v = -d4 / TMath::Power(-r4,0.75);
  Double_t dterm1 = 2*q2*q2*p2 - 3*q2*p4 + 2*q4*p2;
  Double_t dterm2 = 9.0/16.0*d4*d4;
  Double_t dterm3 = 4.0*q2*q2*r4*r4;
  Double_t dterm4 = r4*r4;
  Double_t dterm5 = -3.0/2.0*d4*dterm1;
  Double_t dterm6 = -4.0*q2*r4*dterm1;
  Double_t dterm7 = -2.0*r4*dterm1;
  Double_t dterm8 = 3.0*q2*r4*d4;
  Double_t dterm9 = 3.0/2.0*r4*d4;
  Double_t dterm10= 4*q2*r4*r4;
  err = 1.0/TMath::Power(-r4,3.5)*(+dterm1*dterm1*q2_er*q2_er
				   +dterm2*q4_er*q4_er
				   +dterm3*p2_er*p2_er
				   +dterm4*p4_er*p4_er
				   -dterm5*c_q2q4
				   -dterm6*c_p2q2
				   +dterm7*c_p4q2
				   +dterm8*c_p2q4
				   -dterm9*c_p4q4
				   -dterm10*c_p2p4);
  err = TMath::Sqrt(err);
  return v;
}
