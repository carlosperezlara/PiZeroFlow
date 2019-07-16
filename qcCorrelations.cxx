//=========================
// written by Carlos Perez 
// 2015-2016               
//=========================
#include <iostream>
#include "TMath.h"
#include "qcQ.h"
#include "qcCorrelations.h"

//========
qcCorrelations::qcCorrelations() : fDebug(0) {
  // ctor
}
//========
qcCorrelations::~qcCorrelations() {
  // dtor
}
//========
float qcCorrelations::Make2PC( qcQ *a, qcQ *b, qcQ *ab, float &w ) {
  // maker
  if(!a || !b || !ab) {
    if(fDebug>0) std::cout << "qcCorrelations::Make2PC ERROR null pointers " << std::endl;
    return 37767;
  }
  if(a->M()<1) {
    if(fDebug>0) std::cout << "qcCorrelations::Make2PC ERROR MultiplictyA " << std::endl;
    return 37767;
  }
  if(b->M()<1) {
    if(fDebug>0) std::cout << "qcCorrelations::Make2PC ERROR MultiplictyB " << std::endl;
    return 37767;
  }
  w = a->M()*b->M() - ab->M();
  if( w < 1 ) {
    if(fDebug>0) std::cout << "qcCorrelations::Make2PC ERROR Weight " << std::endl;
    return 37767;
  }
  float num = (a->X()*b->X() + a->Y()*b->Y() - ab->M());
  return num/w;
}
//========
float qcCorrelations::Make4PC( qcQ *a, qcQ *b, qcQ *ab, qcQ *b2, qcQ *ab2, float &w ) {
  // maker
  if(!a || !b || !ab || !ab2) {
    if(fDebug>0) std::cout << "qcCorrelations::Make4PC ERROR null pointers " << std::endl;
    return 37767;
  }
  if(a->M()<1) {
    if(fDebug>0) std::cout << "qcCorrelations::Make4PC ERROR MultiplicityA " << std::endl;
    return 37767;
  }
  if(b->M()<4) {
    if(fDebug>0) std::cout << "qcCorrelations::Make4PC ERROR MultiplicityB " << std::endl;
    return 37767;
  }
  w = (a->M()*b->M() - 3*ab->M())*(b->M()-2)*(b->M()-3);
  if( w < 1 ) {
    if(fDebug>0) std::cout << "qcCorrelations::Make4PC ERROR Weight " << std::endl;
    return 37767;
  }
  float num = 0.0;
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
float qcCorrelations::MakeQC2( float p2, float er, float &err ) {
  err = er;
  return p2;
}
//========
float qcCorrelations::MakeQC4( float p4, float p2, float q2,
			       float p4_er, float p2_er, float q2_er,
			       float c_p4p2, float c_p4q2, float c_p2q2,
			       float &err) {
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
float qcCorrelations::MakeV2( float p2, float q2, float c_p2q2,
			      float p2_err, float q2_err,
			      float &err ) {
  float v = p2 / TMath::Sqrt(q2);
  err = 0.25/q2/q2/q2*(p2*p2*q2_err*q2_err + 4*q2*q2*p2_err*p2_err - 4*q2*p2*c_p2q2);
  err = TMath::Sqrt(err);
  return v;
}
//========
float qcCorrelations::MakeV4( float d4, float r4,
			      float p4, float q4, float p2, float q2,
			      float p4_er, float q4_er, float p2_er, float q2_er,
			      float c_p4q4, float c_p4q2, float c_p2q4, float c_p2p4, float c_p2q2, float c_q2q4,
			      float &err ) {
  float v = -d4 / TMath::Power(-r4,0.75);
  double dterm1 = 2*q2*q2*p2 - 3*q2*p4 + 2*q4*p2;
  double dterm2 = 9.0/16.0*d4*d4;
  double dterm3 = 4.0*q2*q2*r4*r4;
  double dterm4 = r4*r4;
  double dterm5 = -3.0/2.0*d4*dterm1;
  double dterm6 = -4.0*q2*r4*dterm1;
  double dterm7 = -2.0*r4*dterm1;
  double dterm8 = 3.0*q2*r4*d4;
  double dterm9 = 3.0/2.0*r4*d4;
  double dterm10= 4*q2*r4*r4;
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


