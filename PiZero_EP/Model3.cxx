#include <string.h>
#include "Model3.h"

using std::string;

//=======================================================
Model3::Model3() {
}

//============================================================
// Default Parametrisation
// FSGN[0] + FBGR[1,2] + SGN[3,4,5,6,7] + BGR[8,9,10]
void Model3::Init() {
  SetNPars( 18 );
  SetParameter( 0, "v2sgn"  , +0.10, -0.20, +0.40);
  SetParameter( 1, "vsbgr_0", +0.10, -5.00, +5.00);
  SetParameter( 2, "v2bgr_1", +0.10, -1.00, +1.00);
  SetParameter( 3, "v2bgr_2", +0.10, -10.0, +10.0);
  SetParameter( 4, "v2bgr_3", +0.10, -10.0, +10.0);
  ///
  SetParameter( 5, "sgn_raw", +1e+3, +1e+2, +1e+10);
  SetParameter( 6, "sgn_mu" , 0.138, 0.135, 0.139);
  SetParameter( 7, "sgn_sg" , 0.010, 0.001, 0.100);
  ///
  for(int i=0; i!=8; ++i)
    fParFix[i] = false;
  SetParameter( 8, "bgr_0", 0.00, -1.0, +1.0);
  SetParameter( 9, "bgr_1", 0.00, -1.0, +1.0);
  SetParameter(10, "bgr_2", 0.00, -1.0, +1.0);
  SetParameter(11, "bgr_3", 0.00, -1.0, +1.0);
  SetParameter(12, "bgr_4", 0.00, -1.0, +1.0);
  SetParameter(13, "bgr_5", 0.00, -1.0, +1.0);
  SetParameter(14, "bgr_6", 0.00, -1.0, +1.0);
  SetParameter(15, "bgr_7", 0.00, -1.0, +1.0);
  SetParameter(16, "bgr_8", 0.00, -1.0, +1.0);
  SetParameter(17, "bgr_9", 0.00, -1.0, +1.0);
  for(int i=8; i!=fNPar; ++i)
    fParFix[i] = true;
  SetFitRangeFromMassHistogram();
}
double Model3::FBgr(double x, double *p) {
  return double( p[1] + p[2]*(x-p[6]) + p[3]*(x-p[6])*(x-p[6]) );
}
double Model3::Sgn(double x, double *p) {
  return double( p[5]*TMath::Gaus(x,p[6],p[7],1) );
}
double Model3::Bgr(double x, double *p) {
  return double( p[8]
		 + p[9]*x
		 + p[10]*(2*x*x-1)
		 + p[11]*(4*x*x*x-3*x)
		 + p[12]*(8*pow(x,4)-8*x*x+1)
		 + p[13]*(16*pow(x,5)-20*pow(x,3)+5*x)
		 + p[14]*(32*pow(x,6)-48*pow(x,4)+18*x*x-1)
		 + p[15]*(64*pow(x,7)-112*pow(x,5)+56*x*x*x-7*x)
		 + p[16]*(128*pow(x,8)-256*pow(x,6)+160*pow(x,4)-32*x*x+1)
		 + p[17]*(256*pow(x,9)-576*pow(x,7)+432*pow(x,5)-120*x*x*x+9*x )
		 );
}
