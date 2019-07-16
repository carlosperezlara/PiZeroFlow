#ifndef DOUBLEFITTER_H
#define DOUBLEFITTER_H

#include <string.h>
#include <TH1D.h>
#include <TMath.h>
#include <TF1.h>
#include <TVirtualFitter.h>

using std::string;

class DoubleFitter {
 protected:
  TH1D *fMass;
  TH1D *fFlow;

  int fNPar;
  string fParNam[30];
  double fParIni[30];
  double fParMin[30];
  double fParMax[30];
  double fParVal[30];
  double fParErr[30];
  bool   fParFix[30];

  double fMinFit;
  double fMaxFit;

  double fLastChi2;
  int fNDF;
  int fStatus;

  // function parametrisation
  virtual double Sgn(double x, double *p);
  virtual double Bgr(double x, double *p);
  virtual double FBgr(double x, double *p);

  // main functions
  double Mass(double x, double *p, double xmin=-1, double xmax=-1); // SGN(x) + BGR(x)
  double Flow(double x, double *p, double xmin=-1, double xmax=-1); // ( SGN(x)*FSGN(x) + BGR(x)*FBGR(x) ) / TOT(x)

  // helpers
  double MassPointer(double *x, double *p); // MASS(x)
  double FlowPointer(double *x, double *p); // FLOW(x)
  double SgnFracPointer(double *x, double *p); // SGN(x)/TOT(x)
  double BgrFracPointer(double *x, double *p); // BGR(x)/TOT(x)
  double SgnPointer(double *x, double *p); // SGN(x)
  double BgrPointer(double *x, double *p); // BGR(x)

  void SetFitRangeFromMassHistogram();
 public:
  DoubleFitter();
  virtual ~DoubleFitter();
  virtual void Init();
  
  void Fit( TVirtualFitter *minuit ); // main
  void SaveAs( string output ); // save results
  void LoadFrom( string input ); // load results
  void LoadIni( string input ); // load init parameters
  void Print( string output ); // print data members
  void ExportTF1( TF1 **mass, TF1 **flow, TF1 **msgn, TF1 **mbgr, int color=0 ); // tf1
  void ExportTF1( TF1 **mass, TF1 **flow, TF1 **mbgr, int color=0 ); // tf1
  void ExportTF1( TF1 **mass, TF1 **flow, int color=0 ); // tf1

  // core
  virtual void FCN(Int_t& , Double_t* , Double_t &val, Double_t *par, Int_t ); // error function

  // setters
  void SetMassHistogram( TH1D *p1 ) { fMass = p1; }
  void SetFlowHistogram( TH1D *p1 ) { fFlow = p1; }
  void SetNPars( int value ) { fNPar = value; }
  void SetParameter( int i1, string s1, double d1, double d2, double d3 ) {
    fParIni[i1] = d1;
    fParNam[i1] = s1;
    fParMin[i1] = d2;
    fParMax[i1] = d3;
  }
  void SetParameterValue( int i1, double v1, bool rearrangelimits=false, double per=0.0001 ) {
    fParIni[i1]=fParVal[i1]=v1;
    if(rearrangelimits) {
      fParMin[i1] = TMath::Min( v1*(1-per), v1*(1+per) );
      fParMax[i1] = TMath::Max( v1*(1-per), v1*(1+per) );
    }
  }
  void SetMinFit( double val ) { fMinFit = val; }
  void SetMaxFit( double val ) { fMaxFit = val; }

  // getters
  int GetNumberOfParameters() { return fNPar; }
  string GetParameterName(  int i1 ) { return fParNam[i1]; }
  double GetParameterValue( int i1 ) { return fParVal[i1]; }
  double GetParameterError( int i1 ) { return fParErr[i1]; }
  double GetChiSquared() { return fLastChi2; }
  int GetNDF() { return fNDF; }
  int GetSatus() { return fStatus; }
  double GetMinFit() { return fMinFit; }
  double GetMaxFit() { return fMaxFit; }
  double GetParMin(int v) { return fParMin[v]; }
  double GetParMax(int v) { return fParMax[v]; }
};


//============================================================
// inline definitions
inline double DoubleFitter::SgnPointer(double *x, double *p) {
  return double( Sgn(x[0],p) );
}
inline double DoubleFitter::BgrPointer(double *x, double *p) {
  return double( Bgr(x[0],p) );
}
inline double DoubleFitter::SgnFracPointer(double *x, double *p) {
  return double( Sgn(x[0],p)/(Sgn(x[0],p)+Bgr(x[0],p)) );
}
inline double DoubleFitter::BgrFracPointer(double *x, double *p) {
  return double( Bgr(x[0],p)/(Sgn(x[0],p)+Bgr(x[0],p)) );
}
inline double DoubleFitter::MassPointer(double *x, double *p) {
  return double( Mass(x[0],p) );
}
inline double DoubleFitter::FlowPointer(double *x, double *p) {
  return double( Flow(x[0],p) );
}

#endif
