#ifndef __AT_PIZERO_HH__
#define __AT_PIZERO_HH__

#include <vector>
#include "AT_ReadTree.h"
#include "TLorentzVector.h"

class TH1D;
class TH2D;
class TList;
class TProfile;

class AT_PiZero : public AT_ReadTree {
 public:
  AT_PiZero();
  virtual ~AT_PiZero();
  virtual void MyInit();
  virtual void MyExec();
  virtual void MyFinish();
  void SetPt(float m, float M) {fCuts.minPt=m;fCuts.maxPt=M;}
  void SetDist(float val) {fCuts.dist=val;}
  void SetAlpha(float val) {fCuts.alpha=val;}
  void SetTime(float val) {fCuts.time=val;}

 private:
  bool IsBad(int sc, int y, int z);
  int  P0_VertexBin(float vtx);
  int  EMCMAP[8][48][96];

  TList *fListPZ;
  TH1D *hVertex;
  TH1D *hCentrality;
  TProfile *hNClu0;
  TProfile *hNClu1;
  TH2D *hPizeroMass[4][8]; // Step Section
  TH2D *hPizeroMixMass[8]; // Section

  struct FASTCLU {
    float ecore;
    int idx;
    float x;
    float y;
    float z;
    float t;
  };
  struct PI0CUTS {
    float minPt;
    float maxPt;
    float dist;
    float alpha;
    float time;
  };

  PI0CUTS fCuts;
  std::vector<FASTCLU> fPrevious[20]; //!
  std::vector<FASTCLU> fBuffer; //!
};

#endif
