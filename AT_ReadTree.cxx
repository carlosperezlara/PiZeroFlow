#include <iostream>
#include <fstream>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TGraph.h>
#include "Analysis.h"
#include "AT_ReadTree.h"

AT_ReadTree::AT_ReadTree() : AnalysisTask() {
  // -20.0 ==> +20.0 (40+1)
  // 0.5 ==> 60.5 (60+1)
  fListRT = NULL;
  hEvents = NULL;
  hTriggers0 = NULL;
  hCentrality0 = NULL;
  hVertex0 = NULL;
  hPileUpRejectionS = NULL;
  hPileUpRejectionN = NULL;
  hCentralitySelection = NULL;
  for(int i=0; i!=kBBC_NQ; ++i) {
    for(int j=0; j!=3; ++j) {
      for(int k=0; k!=5; ++k) {
	hPsi[i][j][k] = NULL;
      }
    }
  }
  for(int i=0; i!=kBBC_NQR; ++i) {
    hPsi_Q[i] = NULL;
    hPsi_Qab[i][0] = NULL;
    hPsi_Qab[i][1] = NULL;
  }
  fBBCQCal = true;
  for(int bce=0; bce!=kBinsCen; ++bce) {
    for(int bvt=0; bvt!=kBinsVtx; ++bvt) {
      for(int ord=0; ord!=kBBC_NQ; ++ord) {
        for(int se=0; se!=2; ++se) {
	  for(int st=0; st!=3; ++st) {
	    fBBCm[se][ord][0][bce][bvt][st] = 0.0;
	    fBBCm[se][ord][1][bce][bvt][st] = 0.0;
	  }
        }
      }
      for(int ord=0; ord!=kBBC_NQR; ++ord) {
        for(int ior=0; ior!=kFlatten; ++ior) {
          for(int se=0; se!=3; ++se) {
            fBBCc[se][ior][ord][bce][bvt] = 0.0;
            fBBCs[se][ior][ord][bce][bvt] = 0.0;
          }
        }
      }
    }
  }
  fTime[0] = fTime[1] = fTime[2] = 0;

  pMXSempc3x3 = new std::vector<Float_t>;
  pQ1ex = new std::vector<qcQ>;
  pQ2ex = new std::vector<qcQ>;
  pQ3ex = new std::vector<qcQ>;
  pQ4ex = new std::vector<qcQ>;
  pQ6ex = new std::vector<qcQ>;
  pQ8ex = new std::vector<qcQ>;
  pQ1fv = new std::vector<qcQ>;
  pQ2fv = new std::vector<qcQ>;
  pQ3fv = new std::vector<qcQ>;
  
  pQ1bb = new std::vector<qcQ>;
  pQ2bb = new std::vector<qcQ>;
  pQ3bb = new std::vector<qcQ>;
  pQ4bb = new std::vector<qcQ>;
  pQ5bb = new std::vector<qcQ>;
  pQ6bb = new std::vector<qcQ>;
  pQ8bb = new std::vector<qcQ>;
  pQ10bb = new std::vector<qcQ>;
  pQ12bb = new std::vector<qcQ>;

  pEMCid = new std::vector<Int_t>;
  pEMCtwrid = new std::vector<Int_t>;
  pEMCx = new std::vector<Float_t>;
  pEMCy = new std::vector<Float_t>;
  pEMCz = new std::vector<Float_t>;
  pEMCecore = new std::vector<Float_t>;
  pEMCecent = new std::vector<Float_t>;
  pEMCchisq = new std::vector<Float_t>;
  pEMCtimef = new std::vector<Float_t>;

  pTRKqua = new std::vector<Int_t>;
  pTRKpt = new std::vector<Float_t>;
  pTRKphi = new std::vector<Float_t>;
  pTRKpz = new std::vector<Float_t>;
  pTRKecore = new std::vector<Float_t>;
  pTRKetof = new std::vector<Float_t>;
  pTRKtwrid = new std::vector<Int_t>;
  pTRKplemc = new std::vector<Float_t>;
  pTRKchisq = new std::vector<Float_t>;
  pTRKdphi = new std::vector<Float_t>;
  pTRKdz = new std::vector<Float_t>;
  pTRKpc3sdphi = new std::vector<Float_t>;
  pTRKpc3sdz = new std::vector<Float_t>;
  pTRKzed = new std::vector<Float_t>;
  pTRKdisp = new std::vector<Float_t>;
  pTRKprob = new std::vector<Float_t>;
  pTRKcid = new std::vector<Int_t>;

  pMXSempccent = new std::vector<Float_t>;
  pMXSempc3x3 = new std::vector<Float_t>;
  pMXSpt = new std::vector<Float_t>;
  pMXSpz = new std::vector<Float_t>;
  pMXSeta = new std::vector<Float_t>;
  pMXSphi = new std::vector<Float_t>;
  pMXSflyr = new std::vector<Int_t>;
  pMXSsingleD = new std::vector<Float_t>;
  pMXSsingleP = new std::vector<Int_t>;

  fSkipDetails = false;
  fSkipClusters = false;
  fSkipTracks = false;
  fSkipShowers = false;
  fSkipPileUpCuts = false;
}

void AT_ReadTree::Init() {
  std::cout << "AT_ReadTree::Init()" << std::endl;
  Analysis *ana = Analysis::Instance();

  fListRT = new TList();
  fListRT->SetName("AT_ReadTree");
  fListRT->SetOwner();
  hEvents = new TH1D("hEvents","hEvents",15,-0.5,14.5);
  hEvents->GetXaxis()->SetBinLabel(1,"AllEvents");
  hEvents->GetXaxis()->SetBinLabel(2,"AT_Trigger");
  hEvents->GetXaxis()->SetBinLabel(3,"AT_Vertex");
  hEvents->GetXaxis()->SetBinLabel(4,"AT_PileUp");
  hEvents->GetXaxis()->SetBinLabel(5,"AT_Centrality");
  hEvents->GetXaxis()->SetBinLabel(6,"AT_VTXCEN_BIN");
  hEvents->GetXaxis()->SetBinLabel(10,"AT_X");
  fListRT->Add( hEvents );
  hTriggers0 = new TH1D("hTriggers","hTriggers",10,-0.5,9.5);
  hTriggers0->GetXaxis()->SetBinLabel(1,"0x00000008");
  hTriggers0->GetXaxis()->SetBinLabel(2,"0x00000010");
  hTriggers0->GetXaxis()->SetBinLabel(3,"0x00000040");
  hTriggers0->GetXaxis()->SetBinLabel(4,"0x00010000");
  hTriggers0->GetXaxis()->SetBinLabel(5,"0x00020000");
  hTriggers0->GetXaxis()->SetBinLabel(6,"0x00100000");
  hTriggers0->GetXaxis()->SetBinLabel(7,"0x00400000");
  fListRT->Add( hTriggers0 );
  hCentrality0 = new TH1D("hCentrality","hCentrality",100,0.,100.);
  fListRT->Add( hCentrality0 );
  hVertex0 = new TH1D("hVertex","hVertex",100,-40,+40);
  fListRT->Add( hVertex0 );
  hPileUpRejectionS = new TH2D("hPileUpRejectionS","PileUpRejection;;BBC south",
                               2,-0.5,1.5,300,0.,300.);
  fListRT->Add( hPileUpRejectionS );
  hPileUpRejectionN = new TH2D("hPileUpRejectionN","PileUpRejection;;BBC north",
                               2,-0.5,1.5,300,0.,300.);
  fListRT->Add( hPileUpRejectionN );
  hCentralitySelection = new TH2D("hCentralitySelection",";Cent;BBCs",
                                  100,0.,100.,1000,0.,200.);
  fListRT->Add( hCentralitySelection );
  TList *listHar = new TList();
  listHar->SetName("PsiSteps");
  listHar->SetOwner();
  for(int i=0; i!=5; ++i) {
    for(int se=0; se!=3; ++se) {
      for(int ord=0; ord!=kBBC_NQ; ++ord) {
	int har = ord+1;
        hPsi[ord][se][i] = new TH1D(Form("hPsi_HAR%d_SE%d_STEP%d",har,se,i),
                                    Form("Psi_HAR%d_SE%d_STEP%d",har,se,i),
                                    200, -0.1, TMath::TwoPi()+0.1 );
        listHar->Add( hPsi[ord][se][i] );
      }
    }
  }
  fListRT->Add( listHar );

  //==
  TList *hPsiList = new TList();
  hPsiList->SetName("PsiQualityAssurance");
  hPsiList->SetOwner();
  for(int i=0; i!=kBBC_NQR; ++i) {
    Int_t har = fQ[i]->Order();
    hPsi_Q[i] =      new TH1D( Form("PsiQ%d",har),    Form("PsiQ%d",har),    200, -0.1, TMath::TwoPi() );
    hPsiList->Add( hPsi_Q[i] );
    hPsi_Qab[i][0] = new TH1D( Form("PsiQ%dSE0",har), Form("PsiQ%dSE0",har), 200, -0.1, TMath::TwoPi() );
    hPsiList->Add( hPsi_Qab[i][0] );
    hPsi_Qab[i][1] = new TH1D( Form("PsiQ%dSE1",har), Form("PsiQ%dSE1",har), 200, -0.1, TMath::TwoPi() );
    hPsiList->Add( hPsi_Qab[i][1] );
  }
  fListRT->Add( hPsiList );

  TTree *tree = ana->GetTree();
  if(!tree) {
    std::cout << "AT_ReadTree:Init says: Tree not found." << std::endl;
    return;
  }
  //Opening assigning branches
  tree->SetBranchAddress("Event",&fGLB);
  tree->SetBranchAddress("EventC",&fGLB2);
  if(!fSkipDetails) {
    //=
    tree->SetBranchAddress("Q1ex",&pQ1ex);
    tree->SetBranchAddress("Q2ex",&pQ2ex);
    tree->SetBranchAddress("Q3ex",&pQ3ex);
    tree->SetBranchAddress("Q4ex",&pQ4ex);
    tree->SetBranchAddress("Q6ex",&pQ6ex);
    tree->SetBranchAddress("Q8ex",&pQ8ex);
    tree->SetBranchAddress("Q1fv",&pQ1fv);
    tree->SetBranchAddress("Q2fv",&pQ2fv);
    tree->SetBranchAddress("Q3fv",&pQ3fv);
    tree->SetBranchAddress("Q1bb",&pQ1bb);
    tree->SetBranchAddress("Q2bb",&pQ2bb);
    tree->SetBranchAddress("Q3bb",&pQ3bb);
    tree->SetBranchAddress("Q4bb",&pQ4bb);
    tree->SetBranchAddress("Q5bb",&pQ5bb);
    tree->SetBranchAddress("Q6bb",&pQ6bb);
    tree->SetBranchAddress("Q8bb",&pQ8bb);
    tree->SetBranchAddress("Q10bb",&pQ10bb);
    tree->SetBranchAddress("Q12bb",&pQ12bb);
    //=
    if(!fSkipClusters) {
      tree->SetBranchAddress("EMCid",   &pEMCid);
      tree->SetBranchAddress("EMCtwrid",&pEMCtwrid);
      tree->SetBranchAddress("EMCx",    &pEMCx);
      tree->SetBranchAddress("EMCy",    &pEMCy);
      tree->SetBranchAddress("EMCz",    &pEMCz);
      tree->SetBranchAddress("EMCecore",&pEMCecore);
      tree->SetBranchAddress("EMCecent",&pEMCecent);
      tree->SetBranchAddress("EMCchisq",&pEMCchisq);
      tree->SetBranchAddress("EMCtimef",&pEMCtimef);
    }
    //=
    if(!fSkipTracks) {
      tree->SetBranchAddress("TRKqua",  &pTRKqua);
      tree->SetBranchAddress("TRKpt",   &pTRKpt);
      tree->SetBranchAddress("TRKphi",  &pTRKphi);
      tree->SetBranchAddress("TRKpz",   &pTRKpz);
      tree->SetBranchAddress("TRKecore",&pTRKecore);
      tree->SetBranchAddress("TRKetof", &pTRKetof);
      tree->SetBranchAddress("TRKplemc",&pTRKplemc);
      tree->SetBranchAddress("TRKtwrid",&pTRKtwrid);
      tree->SetBranchAddress("TRKchisq",&pTRKchisq);
      tree->SetBranchAddress("TRKdphi", &pTRKdphi);
      tree->SetBranchAddress("TRKdz",   &pTRKdz);
      tree->SetBranchAddress("TRKpc3sdphi",&pTRKpc3sdphi);
      tree->SetBranchAddress("TRKpc3sdz",  &pTRKpc3sdz);
      tree->SetBranchAddress("TRKzed",  &pTRKzed);
      tree->SetBranchAddress("TRKdisp", &pTRKdisp);
      tree->SetBranchAddress("TRKprob", &pTRKprob);
      tree->SetBranchAddress("TRKcid",  &pTRKcid);
    }
    //=
    if(!fSkipShowers) {
      tree->SetBranchAddress("MXSpt",  &pMXSpt);
      tree->SetBranchAddress("MXSpz",  &pMXSpz);
      tree->SetBranchAddress("MXSphi", &pMXSphi);
      tree->SetBranchAddress("MXSflyr",&pMXSflyr);
      tree->SetBranchAddress("MXSsingleD", &pMXSsingleD);
      tree->SetBranchAddress("MXSsingleP", &pMXSsingleP);
      tree->SetBranchAddress("MXSempccent",&pMXSempccent);
      tree->SetBranchAddress("MXSempc3x3", &pMXSempc3x3);
    }
  }
  LoadTableEP();
  LoadTableTime();

  MyInit();
}

void AT_ReadTree::Finish() {
  fListRT->Write( fListRT->GetName(), kSingleKey );
  MyFinish();
}

AT_ReadTree::~AT_ReadTree() {
  if(fListRT) delete fListRT;
}

Bool_t AT_ReadTree::Exec() {
  hEvents->Fill(0); // =====> 0
  float vtx = fGLB.vtxZ;
  float cen = fGLB.cent;
  unsigned int trigger = fGLB.trig;
  bool trig = false;
  if(0x00000008 | trigger) hTriggers0->Fill(0);
  if(0x00000010 | trigger) hTriggers0->Fill(1);
  if(0x00000040 | trigger) hTriggers0->Fill(2);
  if(0x00010000 | trigger) hTriggers0->Fill(3);
  if(0x00020000 | trigger) hTriggers0->Fill(4);
  if(0x00100000 | trigger) hTriggers0->Fill(5);
  if(0x00400000 | trigger) hTriggers0->Fill(6);
  if(trigger & kTriggerMask) trig = true;
  if(!trig) return kFALSE;
  hEvents->Fill(1); // =====> 1
  int bvtx = BinVertex( vtx );
  if(bvtx<0) return kFALSE;
  hEvents->Fill(2); // =====> 2
  hVertex0->Fill(vtx);
  float sgnS = fGLB.bbcs;
  float sgnN = fGLB2.bbcn;
  hPileUpRejectionS->Fill(0.,sgnS);
  hPileUpRejectionN->Fill(0.,sgnN);
  if(!fSkipPileUpCuts) {
    float meanS = fGLB2.bbcsTmean;
    float meanN = fGLB2.bbcnTmean;
    if( TMath::Abs(meanS-meanN+fTime[0]) > 5*0.60 ) return kFALSE;
    double rmsS = fGLB2.bbcsTrms;
    double rmsN = fGLB2.bbcnTrms;
    if( rmsS/fTime[1] + rmsN/fTime[2] > 1 ) return kFALSE;
  }
  hEvents->Fill(3); // =====> 3
  int bcen = BinCentrality( cen );
  hPileUpRejectionS->Fill(1.,sgnS);
  hPileUpRejectionN->Fill(1.,sgnN);
  hCentralitySelection->Fill(cen,sgnS);
  if(bcen<0) return kFALSE;
  hEvents->Fill(4); // =====> 4
  hCentrality0->Fill(cen);
  if(fBBCQCal&&!fSkipDetails) MakeBBCEventPlanes(bcen,bvtx);
  for(int i=0; i!=kBBC_NQR; ++i) {
    hPsi_Qab[i][0]->Fill( fQab[i][0]->Psi() );
    hPsi_Qab[i][1]->Fill( fQab[i][1]->Psi() );
    hPsi_Q[i]->Fill( fQ[i]->Psi() );
    //std::cout << " AT_ReadTree::MakeBBC => " << i << " PSIT : " << fQ[i]->Psi() << " ";
    //std::cout << " PSIA : " << fQab[i][0]->Psi() << " ";
    //std::cout << " PSIB : " << fQab[i][1]->Psi() << std::endl;
  }
  hEvents->Fill(5); // =====> 5
  
  MyExec();
  return kTRUE;
}

void AT_ReadTree::MakeBBCEventPlanes(int bcen, int bvtx) {
  qcQ qvec[kBBC_NQ][3];
  for(int se=0; se!=2; ++se) {
    qvec[0][se] = pQ1bb->at(se); // Q1
    qvec[1][se] = pQ2bb->at(se); // Q2
    qvec[2][se] = pQ3bb->at(se); // Q3
    qvec[3][se] = pQ4bb->at(se); // Q4
    qvec[4][se] = pQ5bb->at(se); // Q5
    qvec[5][se] = pQ6bb->at(se); // Q6
    qvec[6][se] = pQ8bb->at(se); // Q8
    qvec[7][se] = pQ10bb->at(se); // Q10
    qvec[8][se] = pQ12bb->at(se); // Q12
    if(qvec[0][se].M()<1) {
      return;
    }
  }
  for(int k=0; k!=kBBC_NQ; ++k) { // order
    qvec[k][2] = qvec[k][0] + qvec[k][1];
  }

  int snapshot;
  /// SNAPSHOT 0
  snapshot=0;
  for(int ord=0; ord!=kBBC_NQ; ++ord) {
    for(int se=0; se!=3; ++se) {
      hPsi[ord][se][snapshot]->Fill(qvec[ord][se].Psi2Pi());
    }
  }

  // ======= STAGE 2: Recentering SubEvents (STEP1)  =======
  // Qx => Qx - <Qx>
  // Qy => Qy - <Qy>
 for(int k=0; k!=kBBC_NQ; ++k) { // order
    for(int se=0; se!=3; ++se) { // subevent
      double x = qvec[k][se].X();
      double y = qvec[k][se].Y();
      double cn = fBBCm[se][k][0][bcen][bvtx][0];
      double sn = fBBCm[se][k][1][bcen][bvtx][0];
      qvec[k][se].SetXY( x - cn, y - sn, qvec[k][se].NP(), qvec[k][se].M() );
    }
  }

  /// SNAPSHOT 1
  snapshot=1;
  for(int ord=0; ord!=kBBC_NQ; ++ord) {
    for(int se=0; se!=3; ++se) {
      hPsi[ord][se][snapshot]->Fill(qvec[ord][se].Psi2Pi());
    }
  }

  int twon[6] = {1,3,5,6,7,8}; // [1,2,3,4,5,6] => [2,4,6,8,10,12]
  // ======= STAGE 4: Twisting SubEvents (STEP2)  =======
  // LdaSm = s2n/(1+c2n)
  // LdaSp = s2n/(1-c2n)
  // Qx => (Qx-LdaSm*Qy)/(1-LdaSm*ldaSp)
  // Qy => (Qy-LdaSm*Qx)/(1-LdaSm*ldaSp)
  for(int k=0; k!=kBBC_NQR; ++k) { // order
    for(int se=0; se!=3; ++se) { // subevent
      double x = qvec[k][se].X();
      double y = qvec[k][se].Y();
      double c2n = fBBCm[se][twon[k]][0][bcen][bvtx][1] / qvec[k][se].M();
      double s2n = fBBCm[se][twon[k]][1][bcen][bvtx][1] / qvec[k][se].M();
      double ldaSm = s2n/(1.0+c2n);
      double ldaSp = s2n/(1.0-c2n);
      double den = 1.0 - ldaSm*ldaSp;
      qvec[k][se].SetXY( (x-ldaSm*y) / den,
                        (y-ldaSp*x) / den,
                        qvec[k][se].NP(),
                        qvec[k][se].M() );
      if( TMath::IsNaN( qvec[k][se].X() ) || TMath::IsNaN( qvec[k][se].Y() ) ) {
	std::cout << "Error building coefficient ";
	std::cout << " | qvec.M: " << qvec[k][se].M();
	std::cout << " | c2n: " << c2n;
	std::cout << " | s2n: " << s2n;
	std::cout << " | ldaSm: " << ldaSm;
	std::cout << " | ldaSp: " << ldaSp;
	std::cout << " | den: " << den << std::endl;
      }
    }
  }

  /// SNAPSHOT 2
  snapshot=2;
  for(int ord=0; ord!=kBBC_NQ; ++ord) {
    for(int se=0; se!=3; ++se) {
      hPsi[ord][se][snapshot]->Fill(qvec[ord][se].Psi2Pi());
    }
  }

  // ======= STAGE 6: Rescaling SubEvents (STEP3)  =======
  // Qx => Qx/(1+c2n)
  // Qy => Qy/(1-c2n)
  for(int k=0; k!=kBBC_NQR; ++k) { // order
    for(int se=0; se!=3; ++se) { // subevent
      double x = qvec[k][se].X();
      double y = qvec[k][se].Y();
      double c2n = fBBCm[se][twon[k]][0][bcen][bvtx][2] / qvec[k][se].M();
      double a2np = 1.0+c2n/2; //1.0+c2n;
      double a2nm = 1.0-c2n/2; //1.0-c2n;
      qvec[k][se].SetXY(x / a2np,
                        y / a2nm,
                        qvec[k][se].NP(),
                        qvec[k][se].M() );
      if( TMath::IsNaN( qvec[k][se].X() ) || TMath::IsNaN( qvec[k][se].Y() ) ) {
	std::cout << "Error building coefficient [2] ";
	std::cout << " | a2np: " << a2np;
	std::cout << " | a2nm: " << a2nm << std::endl;
      }
    }
  }

  /// SNAPSHOT 3
  snapshot=3;
  for(int ord=0; ord!=kBBC_NQ; ++ord) {
    for(int se=0; se!=3; ++se) {
      hPsi[ord][se][snapshot]->Fill(qvec[ord][se].Psi2Pi());
    }
  }

  // ======= STAGE 8: Bulding Full Q and Storing Flattening Coeficients  =======
  for(int k=0; k!=kBBC_NQR; ++k) { // order
    for(int se=0; se!=3; ++se) { // sub events
      double psi = qvec[k][se].Psi2Pi();
      double delta = 0;      
      for(int ik=0; ik!=32; ++ik) { // correction order
	int nn = ik+1;
	delta -= TMath::Cos(nn*psi)*fBBCs[se][ik][k][bcen][bvtx]*2.0/nn;
	delta += TMath::Sin(nn*psi)*fBBCc[se][ik][k][bcen][bvtx]*2.0/nn;
      }
      double cn = TMath::Cos( (k+1)*delta );
      double sn = TMath::Sin( (k+1)*delta );
      double xprime = qvec[k][se].X()*cn - qvec[k][se].Y()*sn;
      double yprime = qvec[k][se].X()*sn + qvec[k][se].Y()*cn;
      qvec[k][se].SetXY( xprime, yprime, qvec[k][se].NP(),qvec[k][se].M() );
    }
    //moving to analysis manager
    //saving

    fQab[k][0]->CopyFrom( qvec[k][0] );
    fQab[k][1]->CopyFrom( qvec[k][1] );
    fQ[k]->CopyFrom( qvec[k][2] );
    if( TMath::Abs( fQab[k][0]->Psi() - qvec[k][0].Psi() ) > 1e-6 ) {
      std::cout << " ERROR!!! " << k << " A " << fQab[k][0]->Psi() << " " << qvec[k][0].Psi() << std::endl;
    }
    if( TMath::Abs( fQab[k][1]->Psi() - qvec[k][1].Psi() ) > 1e-6 ) {
      std::cout << " ERROR!!! " << k << " B " << fQab[k][1]->Psi() << " " << qvec[k][1].Psi() << std::endl;
    }
    if( TMath::Abs( fQ[k]->Psi() - qvec[k][2].Psi() ) > 1e-6 ) {
      std::cout << " ERROR!!! " << k << " T " << fQ[k]->Psi() << " " << qvec[k][2].Psi() << std::endl;
    }
  }

  /// SNAPSHOT 4
  snapshot=4;
  for(int ord=0; ord!=kBBC_NQ; ++ord) {
    for(int se=0; se!=3; ++se) {
      hPsi[ord][se][snapshot]->Fill(qvec[ord][se].Psi2Pi());
    }
  }
}

int AT_ReadTree::ReferenceTracks() {
  int ntrk=0;
  uint ntrks = pTRKpt->size();
  for(uint itrk=0; itrk!=ntrks; ++itrk) {
    float zed  = pTRKzed->at(itrk);
    int qua = pTRKqua->at(itrk);
    float dphi = pTRKpc3sdphi->at(itrk);
    float dz   = pTRKpc3sdz->at(itrk);
    if(qua!=63) continue;
    if(TMath::Abs(zed)<3||TMath::Abs(zed)>70) continue;
    if(TMath::Abs(dphi)>3) continue;
    if(TMath::Abs(dz)>3) continue;
    ntrk++;
  }
  return ntrk;
}

int AT_ReadTree::BinVertex(float vtx) {
  int ret=-1;
  for(int i=0; i!=kBinsVtx+1; ++i) {
    if(vtx<kMinBinVtx+i) {
      ret = i-1;
      break;
    }
  }
  return ret;
}

int AT_ReadTree::BinCentrality(float cen) {
  int ret=-1;
  for(int i=0; i!=kBinsCen+1; ++i) {
    if(cen<kMinBinCen+i) {
      ret = i-1;
      break;
    }
  }
  return ret;
}

void AT_ReadTree::LoadTableTime( int run ) {
  if(run<0) {
    Analysis *ana = Analysis::Instance();
    run = ana->RunNumber();
    std::cout << " SEGMENT " << ana->SegmentNumber() << std::endl;
  }
  std::cout << " RUN " << run << std::endl;

  std::ifstream fin;
  int irun;
  float tmp, mean, rmsS, rmsN;
  fin.open( "EventChecker/dat/TimeConstants.dat" );
  int nn=0;
  fTime[0] = -9999;
  fTime[1] = 1;
  fTime[2] = 1;
  for(;;++nn) {
    fin >> irun >> mean >> tmp >> tmp >> tmp >> rmsS >> tmp >> tmp >> rmsN;
    if( irun == run ) {
      fTime[0] = mean;
      fTime[1] = rmsS;
      fTime[2] = rmsN;
      break;
    }
    if(!fin.good()) break;
  }
  fin.close();
  std::cout << " Time Constants " << fTime[0] << " " << fTime[1];
  std::cout << " " << fTime[2] << std::endl;
  return;
}

void AT_ReadTree::LoadTableEP( int run ) {
  if(run<0) {
    Analysis *ana = Analysis::Instance();
    run = ana->RunNumber();
    std::cout << " SEGMENT " << ana->SegmentNumber() << std::endl;
  }
  std::cout << " RUN " << run << std::endl;

  int cenXvtx = kBinsVtx*kBinsCen;
  int flaXvtx = kBinsVtx*kFlatten;
  int cenXflaXvtx = kBinsCen*kBinsVtx*kFlatten;

  std::ifstream fin;
  int se, ord, xy, bce, bvt;
  double tmp;

  fin.open( Form("BBC_EPC/tables/BBC_%d.dat",run) );
  int nn=0;
  for(;;++nn) {
    fin >> tmp;
    if(!fin.good()) break;
    int ord = (nn/(cenXvtx*6))%kBBC_NQ;
    int xy = (nn/(cenXvtx*3))%2;
    int se = (nn/cenXvtx)%3;
    int bce = (nn/kBinsVtx)%kBinsCen;
    int bvt = nn%kBinsVtx;
    fBBCm[se][ord][xy][bce][bvt][0] = tmp;
    fBBCm[se][ord][xy][bce][bvt][1] = tmp;
    fBBCm[se][ord][xy][bce][bvt][2] = tmp;
  }
  fin.close();
  std::cout << "   BBC ReCenter coefficients loaded: " << nn << std::endl;

  /*
  fin.open( Form("BBC_EPC/tables/BBC_P_%d.dat",run) );
  nn=0;
  for(;;++nn) {
    fin >> tmp;
    if(!fin.good()) break;
    int ord = (nn/(cenXvtx*6))%kBBC_NQ;
    int xy = (nn/(cenXvtx*3))%2;
    int se = (nn/cenXvtx)%3;
    int bce = (nn/kBinsVtx)%kBinsCen;
    int bvt = nn%kBinsVtx;
    //fBBCm[se][ord][xy][bce][bvt][1] = tmp;
  }
  fin.close();
  std::cout << "   BBC primed coefficients loaded: " << nn << std::endl;

  fin.open( Form("BBC_EPC/tables/BBC_Q_%d.dat",run) );
  nn=0;
  for(;;++nn) {
    fin >> tmp;
    if(!fin.good()) break;
    int ord = (nn/(cenXvtx*6))%kBBC_NQ;
    int xy = (nn/(cenXvtx*3))%2;
    int se = (nn/cenXvtx)%3;
    int bce = (nn/kBinsVtx)%kBinsCen;
    int bvt = nn%kBinsVtx;
    //fBBCm[se][ord][xy][bce][bvt][2] = tmp;
  }
  fin.close();
  std::cout << "   BBC pri-primed coefficients loaded: " << nn << std::endl;
  */

  char cse[3] = {'A','B','C'};
  for(int se=0; se!=3; ++se) {
    fin.open( Form("BBC_EPC/tables/BBC_%c_%d.dat",cse[se],run) );
    std::cout << Form("BBC_EPC/tables/BBC_%c_%d.dat",cse[se],run) << std::endl;
    nn=0;
    for(;;++nn) {
      fin >> tmp;
      if(!fin.good()) break;
      int ord = (nn/(cenXflaXvtx*2))%kBBC_NQR;
      int bce = (nn/(flaXvtx*2))%kBinsCen;
      int bcs = (nn/flaXvtx)%2;
      int bor = (nn/kBinsVtx)%kFlatten;
      int bvt = nn%kBinsVtx;
      if(bcs==0) fBBCc[se][bor][ord][bce][bvt] = tmp;
      else fBBCs[se][bor][ord][bce][bvt] = tmp;
    }
    std::cout << "   BBC Flattening coefficients loaded: " << cse[se] << " " << nn << std::endl;
    fin.close();
  }
}
