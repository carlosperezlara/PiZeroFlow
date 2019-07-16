//=========================
// written by Carlos Perez
// 2015-2016
//=========================
#ifndef __qcCorrelations_HH__
#define __qcCorrelations_HH__

class qcQ;

class qcCorrelations {
 public:
  qcCorrelations();
  virtual ~qcCorrelations();
  float Make2PC(qcQ*, qcQ*, qcQ*, float&); // a b ab err
  float Make4PC(qcQ*, qcQ*, qcQ*, qcQ*, qcQ*, float&); // a b ab b2 ab2 err
  float MakeQC2(float, float, float&);
  float MakeQC4(float,float,float,float,float,float,float,float,float,float&);
  float MakeV2(float,float,float,float,float,float&);
  float MakeV4(float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float&);
  void Debug(int v) {fDebug=v;}

 protected:
  int fDebug;
};
#endif /* __qcCorrelations_H__ */
