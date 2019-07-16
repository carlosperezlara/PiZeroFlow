//=========================
// written by Carlos Perez 
// 2015-2016               
//=========================
#include "TMath.h"
#include "TRandom3.h"
#include "qcQ.h"

//========
qcQ::qcQ(int n) :
  fN(n),
  fNP(0),
  fM(0.),
  fX(0.),
  fY(0.) {
  // ctor
}
//========
qcQ::~qcQ() {
  // dtor
}
//========
void qcQ::Reset() {
  fX = 0.0;
  fY = 0.0;
  fM = 0.0;
  fNP = 0;
}
//========
Double_t qcQ::Psi() {
  if(fNP==0) return -99999;
  return (TMath::Pi() + TMath::ATan2(-fY,-fX) )/fN;
}
//========
Double_t qcQ::Psi2Pi() {
  if(fNP==0) return -99999;
  Double_t psi = Psi();
  psi += gRandom->Integer( fN )*TMath::TwoPi()/fN;
  return psi;
}
//========
void qcQ::Fill(Double_t phi, Double_t w) {
  // filler
  fX += w*TMath::Cos(fN*phi);
  fY += w*TMath::Sin(fN*phi);
  fM += w;
  fNP++;
}
//========
Double_t qcQ::Xhat() {
  Double_t x = fX;
  Double_t y = fY;
  Double_t mod = x*x + y*y;
  if(mod<0) return 0;
  return x/TMath::Sqrt(mod);
}
//========
Double_t qcQ::Yhat() {
  Double_t x = fX;
  Double_t y = fY;
  Double_t mod = x*x+y*y;
  if(mod<0) return 0;
  return y/TMath::Sqrt(mod);
}
//========
void qcQ::SetXY(Double_t x, Double_t y, Int_t np, Double_t m) {
  // filler
  fX = x;
  fY = y;
  fM = m;
  fNP = np;
}
//========
Double_t qcQ::Reduced() {
  if(fM>0) return TMath::Sqrt(ModulusSquared()/M());
  return 0;
}
