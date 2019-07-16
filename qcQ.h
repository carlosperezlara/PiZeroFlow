//=========================
// written by Carlos Perez
// 2015-2016
//=========================
#ifndef __qcQ_HH__
#define __qcQ_HH__

#include <iostream>

class qcQ {
 public:
  qcQ(Int_t n=0);
  virtual ~qcQ();
  Int_t Order() {return fN;}
  Int_t NP() {return fNP;}
  Double_t M() {return fM;}
  Double_t X() {return fX;}
  Double_t Y() {return fY;}
  Double_t Xhat();
  Double_t Yhat();
  Double_t Psi();
  Double_t Psi2Pi();
  void Fill(Double_t p, Double_t w=1.0);
  void SetXY(Double_t x, Double_t y, Int_t n=1, Double_t m=1);
  Double_t ModulusSquared() {return (double(fX)*fX+double(fY)*fY);}
  Double_t Reduced();
  void Reset();
  void SetOrder(Int_t val) {fN=val;}
  qcQ operator+(qcQ const &obj) {
    qcQ res;
    if(fN!=obj.fN)
      std::cout << "Error: qcQ adding two vectors of different order!" << std::endl;
    res.fN = fN;
    res.fNP = fNP + obj.fNP;
    res.fM = fM + obj.fM;
    res.fX = fX + obj.fX;
    res.fY = fY + obj.fY;
    return res;
  }
  void CopyFrom(qcQ q) {
    fN = q.fN;
    fNP = q.fNP;
    fM = q.fM;
    fX = q.fX;
    fY = q.fY;
  }

 protected:
  Int_t fN;
  Int_t fNP;
  Double_t fM;
  Double_t fX;
  Double_t fY;
};
#endif /* __qcQ_H__ */
