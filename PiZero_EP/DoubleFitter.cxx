#include <string.h>
#include <TH1D.h>
#include <TF1.h>
#include <TVirtualFitter.h>
#include <TMinuit.h>
#include <fstream>

#include "DoubleFitter.h"

using std::string;

//=======================================================
DoubleFitter::DoubleFitter() {
}

//=======================================================
DoubleFitter::~DoubleFitter() {
}

//============================================================
// Default Parametrisation
// FSGN[0] + FBGR[1,2] + SGN[3,4,5] + BGR[6,7,8]
void DoubleFitter::Init() {
  SetNPars( 9 );
  SetParameter(0, "v2sgn"  , +0.10, -0.01, +0.40);
  SetParameter(1, "vsbgr_a", +0.10, -5.00, +5.00);
  SetParameter(2, "v2bgr_b", +0.10, -1.00, +1.00);
  SetParameter(3, "sgn_raw", +1e+2, +1e+1, +5e+3);
  SetParameter(4, "sgn_mu" , +1.86, +1.85, +1.87);
  SetParameter(5, "sgn_sig", +0.01, +0.005, +0.10);
  SetParameter(6, "bgr_c"  , -1.00, -1000, +1000);
  SetParameter(7, "bgr_b"  , -1.00, -1000, +1000);
  SetParameter(8, "bgr_a"  , +1e+4, +1e+2, +1e+5);
  SetFitRangeFromMassHistogram();
  for(int i=0; i!=fNPar; ++i)
    fParFix[i] = false;
}
double DoubleFitter::FBgr(double x, double *p) {
  return double( p[1]*(x-p[4])+p[2] );
}
double DoubleFitter::Sgn(double x, double *p) {
  return double( p[3]/TMath::Sqrt(2*TMath::Pi())/p[5]*TMath::Gaus(x,p[4],p[5]) );
}
double DoubleFitter::Bgr(double x, double *p) {
  return double( p[8]*( 1 + p[6]*(x-p[4])*(x-p[4]) + p[7]*(x-p[4]) ) );
}
//============================================================

void DoubleFitter::SetFitRangeFromMassHistogram() {
  fMinFit = fMass->GetXaxis()->GetBinLowEdge( 1 );
  fMaxFit = fMass->GetXaxis()->GetBinUpEdge( fMass->GetNbinsX() );
}
//=======================================================
double DoubleFitter::Mass(double x, double *p, double xmin, double xmax) {
  double yield = Sgn(x,p) + Bgr(x,p);
  if((xmin+xmax)>0){
    TF1 sgn("sgn",this,&DoubleFitter::SgnPointer,fMinFit,fMaxFit,fNPar,"DoubleFitter","SgnPointer");
    TF1 bgr("bgr",this,&DoubleFitter::BgrPointer,fMinFit,fMaxFit,fNPar,"DoubleFitter","BgrPointer");
    sgn.SetParameters(p);
    bgr.SetParameters(p);
    yield = ( sgn.Integral(xmin,xmax) + bgr.Integral(xmin,xmax) )/(xmax-xmin);
  }
  return yield;
}
//=======================================================
double DoubleFitter::Flow(double x, double *p, double xmin, double xmax) {
  double fracSgn = Sgn(x,p);
  double fracBgr = Bgr(x,p);
  double total = fracSgn+fracBgr;
  fracSgn/=total;
  fracBgr/=total;
  if((xmin+xmax)>0){
    TF1 sgnf("sgnfrac",this,&DoubleFitter::SgnFracPointer,fMinFit,fMaxFit,fNPar,"DoubleFitter","SgnFracPointer");
    TF1 bgrf("bgrfrac",this,&DoubleFitter::BgrFracPointer,fMinFit,fMaxFit,fNPar,"DoubleFitter","BgrFracPointer");
    sgnf.SetParameters(p);
    bgrf.SetParameters(p);
    fracSgn = sgnf.Integral(xmin,xmax)/(xmax-xmin);
    fracBgr = bgrf.Integral(xmin,xmax)/(xmax-xmin);
  }
  return double( fracSgn*p[0]+fracBgr*FBgr(x,p) ); // assume v2bgr is constant within mass bin
}
//=======================================================
void DoubleFitter::FCN(Int_t& , Double_t* , Double_t &val, Double_t *par, Int_t ) {
  cout << "**************** CALLED FCN!!!" << endl;
  fNDF=0;
  fLastChi2=0;
  // MASS FIT

  int minBinFit = fMass->FindBin( fMinFit + 1e-6 );
  int maxBinFit = fMass->FindBin( fMaxFit - 1e-6 );
  double tmp;
  for(int mb=minBinFit; mb!=maxBinFit+1; ++mb)
    if(fMass->GetBinContent(mb)>0) {
      tmp = fMass->GetBinContent( mb );
      tmp-= Mass( fMass->GetBinCenter(mb), par, fMass->GetBinLowEdge(mb), fMass->GetBinLowEdge(mb+1) );
      tmp/= fMass->GetBinError( mb );
      fLastChi2 += tmp*tmp;
      ++fNDF;
    }
  // FLOW FIT
  minBinFit = fFlow->FindBin( fMinFit + 1e-6 );
  maxBinFit = fFlow->FindBin( fMaxFit - 1e-6 );
  for(int mb=minBinFit; mb!=maxBinFit+1; ++mb) {
    tmp = fFlow->GetBinContent( mb );
    tmp-= Flow( fFlow->GetBinCenter(mb), par, fFlow->GetBinLowEdge(mb), fFlow->GetBinLowEdge(mb+1) );
    tmp/= fFlow->GetBinError( mb );
    fLastChi2 += tmp*tmp;
    ++fNDF;
  }
  val = fLastChi2;
  fNDF -= fNPar;
}
//=======================================================
void DoubleFitter::Fit( TVirtualFitter *minuit ) {
  for(int i=0; i!=fNPar; ++i)
    minuit->SetParameter( i, fParNam[i].data(), fParIni[i], 1e-6, fParMin[i], fParMax[i] );
  for(int i=0; i!=fNPar; ++i)
    if( fParFix[i] && !minuit->IsFixed(i) ) minuit->FixParameter( i );
    else if ( !fParFix[i] && minuit->IsFixed(i) ) minuit->ReleaseParameter( i );
  double argList[100];
  argList[0] = 5000; // FUNCTION CALLS
  argList[1] = 1e-6; // TOLERANCE
  //fStatus = minuit->ExecuteCommand("DEBUG");
  fStatus = minuit->ExecuteCommand("MIGRAD",argList,2);
  for(int i=0; i!=fNPar; ++i) {
    fParVal[i] = minuit->GetParameter(i);
    fParErr[i] = minuit->GetParError(i);
  }
}
//=======================================================
void DoubleFitter::SaveAs( string out ) {
  std::ofstream fout;
  fout.open( out.data() );
  fout << Form( "%s\n", out.data() );
  fout << Form( "minFit  %f\n", fMinFit );
  fout << Form( "maxFit  %f\n", fMaxFit );
  fout << Form( "nPars   %d\n", fNPar );
  for(int i=0; i!=fNPar; ++i) {
    fout << Form( "%s  %f  %f\n", fParNam[i].data(), fParVal[i], fParErr[i] );
  }
  fout << Form( "Chi2NDF  %f  %d\n", fLastChi2, fNDF );
  fout.close();
}
//=======================================================
void DoubleFitter::LoadFrom( string in ) {
  std::ifstream fin;
  fin.open( in.data() );
  string name;
  double val, err;
  fin >> name; // filename
  fin >> name; // fMinFit
  fin >> fMinFit;
  fin >> name; // fMaxFit
  fin >> fMaxFit;
  fin >> name; // fNPar
  fin >> val;
  fNPar = (int) val;
  for(int i=0; i!=fNPar; ++i) {
    fin >> name;
    fin >> val;
    fin >> err;
    fParNam[i] = name;
    fParVal[i] = val;
    fParErr[i] = err;
  }
  fin >> name; // fLastCh2 fNDF
  fin >> fLastChi2;
  fin >> val;
  fNDF = (int) val;
  fin.close();
}
//=======================================================
void DoubleFitter::LoadIni( string in ) {
  std::ifstream fin;
  fin.open( in.data() );
  string name;
  double val1, val2, val3;
  int val4;
  fin >> name; // fMinFit
  if(!fin.good()) return;
  fin >> fMinFit;
  fin >> name; // fMaxFit
  fin >> fMaxFit;
  fin >> name; // fNPar
  fin >> val4;
  fNPar = (int) val4;
  for(int i=0; i!=fNPar; ++i) {
    fin >> name;
    fin >> val1;
    fin >> val2;
    fin >> val3;
    fin >> val4;
    fParNam[i] = name;
    fParIni[i] = val1;
    fParMin[i] = val2;
    fParMax[i] = val3;
    fParFix[i] = val4;
  }
  fin.close();
}
//=======================================================
void DoubleFitter::Print( string out ) {
  std::ofstream fout;
  fout.open( out.data() );
  fout << Form( "fMinFit  %f\n", fMinFit );
  fout << Form( "fMaxFit  %f\n", fMaxFit );
  fout << Form( "fNPar    %d\n", fNPar );
  for(int i=0; i!=fNPar; ++i) {
    fout << Form( "%s  %f  %f  %f  %d\n", fParNam[i].data(), fParIni[i], fParMin[i], fParMax[i], fParFix[i]?1:0 );
  }
  fout.close();
}
//=======================================================
void DoubleFitter::ExportTF1( TF1 **mass, TF1 **flow, TF1 **msgn, TF1 **mbgr, int color ) {
  TF1 *tmass = new TF1("tf1Mass",this,&DoubleFitter::MassPointer,fMinFit,fMaxFit,fNPar,"DoubleFitter","MassPointer");
  TF1 *tflow = new TF1("tf1Flow",this,&DoubleFitter::FlowPointer,fMinFit,fMaxFit,fNPar,"DoubleFitter","FlowPointer");
  TF1 *tmsgn = new TF1("tf1MSgn",this,&DoubleFitter::SgnPointer, fMinFit,fMaxFit,fNPar,"DoubleFitter","SgnPointer");
  TF1 *tmbgr = new TF1("tf1MBgr",this,&DoubleFitter::BgrPointer, fMinFit,fMaxFit,fNPar,"DoubleFitter","BgrPointer");
  tmass->SetLineColor( color );
  tflow->SetLineColor( color );
  tmass->SetParameters(fParVal);
  tflow->SetParameters(fParVal);
  tmsgn->SetParameters(fParVal);
  tmbgr->SetParameters(fParVal);

  *mass = tmass;
  *flow = tflow;
  *msgn = tmsgn;
  *mbgr = tmbgr;
}
//=======================================================
void DoubleFitter::ExportTF1( TF1 **mass, TF1 **flow, TF1 **mbgr, int color ) {
  TF1 *tmass = new TF1("tf1Mass",this,&DoubleFitter::MassPointer,fMinFit,fMaxFit,fNPar,"DoubleFitter","MassPointer");
  TF1 *tflow = new TF1("tf1Flow",this,&DoubleFitter::FlowPointer,fMinFit,fMaxFit,fNPar,"DoubleFitter","FlowPointer");
  TF1 *tmbgr = new TF1("tf1MBgr",this,&DoubleFitter::BgrPointer, fMinFit,fMaxFit,fNPar,"DoubleFitter","BgrPointer");
  tmass->SetLineColor( color );
  tflow->SetLineColor( color );
  tmass->SetParameters(fParVal);
  tflow->SetParameters(fParVal);
  tmbgr->SetParameters(fParVal);

  *mass = tmass;
  *flow = tflow;
  *mbgr = tmbgr;
}
//=======================================================
void DoubleFitter::ExportTF1( TF1 **mass, TF1 **flow, int color ) {
  TF1 *tmass = new TF1("tf1Mass",this,&DoubleFitter::MassPointer,fMinFit,fMaxFit,fNPar,"DoubleFitter","MassPointer");
  TF1 *tflow = new TF1("tf1Flow",this,&DoubleFitter::FlowPointer,fMinFit,fMaxFit,fNPar,"DoubleFitter","FlowPointer");
  tmass->SetLineColor( color );
  tflow->SetLineColor( color );
  tmass->SetParameters(fParVal);
  tflow->SetParameters(fParVal);

  *mass = tmass;
  *flow = tflow;
}
