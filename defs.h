//====C++
#include <iostream>
#include <fstream>
#include <vector>
//===ROOT
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TList.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "TRandom3.h"
//===OTHERS
#include "PbScIndexer.h"
#include "PbScIndexer.C"
#include "PbGlIndexer.h"
#include "PbGlIndexer.C"
#include "EmcIndexer.h"
#include "EmcIndexer.C"
#include "qcQ.h"
#include "qcQ.cxx"
#include "mxGeometry.h"
#include "mxGeometry.cxx"

using namespace std;

const double kSpeedOfLightinNS = 29.979245829979; //[cm/ns]
const int kNTWRS = 24768;
mxGeometry *fGeo = mxGeometry::Instance();

unsigned int kBBCnc = 0x00000008;
unsigned int kBBCn  = 0x00000010;

TTree *fTree;

typedef struct MyTreeRegister {
  Float_t vtxZ;
  Float_t cent;
  Float_t bbcs;
  Float_t frac;
  UInt_t  trig;
} MyTreeRegister_t;
MyTreeRegister_t fGLB;

std::vector<qcQ> *pQ1ex;
std::vector<qcQ> *pQ2ex;
std::vector<qcQ> *pQ3ex;
std::vector<qcQ> *pQ4ex;
std::vector<qcQ> *pQ6ex;
std::vector<qcQ> *pQ1fv;
std::vector<qcQ> *pQ2fv;
std::vector<qcQ> *pQ1bb;
std::vector<qcQ> *pQ2bb;

std::vector<Int_t>   *pEMCid;
std::vector<Int_t>   *pEMCtwrid;
std::vector<Float_t> *pEMCx;
std::vector<Float_t> *pEMCy;
std::vector<Float_t> *pEMCz;
std::vector<Float_t> *pEMCecore;
std::vector<Float_t> *pEMCecent;
std::vector<Float_t> *pEMCchisq;
std::vector<Float_t> *pEMCtimef;

//  0 (1)   X1 used
//  1 (2)   X2 used
//  2 (4)   UV found
//  3 (8)   UV unique
//  4 (16)  PC1 found
//  5 (48)  PC1 unique
std::vector<Int_t>   *pTRKqua;
std::vector<Float_t> *pTRKpt;
std::vector<Float_t> *pTRKphi;
std::vector<Float_t> *pTRKpz;
std::vector<Float_t> *pTRKecore;
std::vector<Float_t> *pTRKetof;
std::vector<Int_t>   *pTRKtwrid;
std::vector<Float_t> *pTRKplemc;
std::vector<Float_t> *pTRKchisq;// chi2/npe0 -> seme-normalized chi2/dof.
std::vector<Float_t> *pTRKdphi; // diff in rads btw track model projection and hit in EMC
std::vector<Float_t> *pTRKdz;   // diff in cm btw track model projection and hit in EMC
std::vector<Float_t> *pTRKpc3sdphi;
std::vector<Float_t> *pTRKpc3sdz;
std::vector<Float_t> *pTRKzed;  // z coord at which track crosses PC1
std::vector<Float_t> *pTRKdisp; // displacement of ring center wrt Track projection in the RICH PMT array
std::vector<Float_t> *pTRKprob; // probability that particle shower is electromagnetic
std::vector<Int_t>   *pTRKcid;

std::vector<Float_t> *pMXSempccent;
std::vector<Float_t> *pMXSempc3x3;
std::vector<Float_t> *pMXSpt;
std::vector<Float_t> *pMXSpz;
std::vector<Float_t> *pMXSeta;
std::vector<Float_t> *pMXSphi;
std::vector<Int_t>   *pMXSflyr;
std::vector<Float_t> *pMXSsingleD;
std::vector<Int_t>   *pMXSsingleP;
