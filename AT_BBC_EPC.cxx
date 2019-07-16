#include <iostream>
#include <TString.h>
#include <TTree.h>
#include <TList.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TMath.h>

#include "Analysis.h"
#include "AT_BBC_EPC.h"

AT_BBC_EPC::AT_BBC_EPC() : AT_ReadTree() {
  fListBBCEPC = NULL;
  for(int nq=0; nq!=kBBC_NQ; ++nq) {
    for(int ce=0; ce!=kBinsCen; ++ce) {
      for(int se=0; se!=3; ++se) {
	for(int st=0; st!=3; ++st) {
	  hQxVtx[nq][se][ce][st] = NULL;
	  hQyVtx[nq][se][ce][st] = NULL;
	}
      }
    }
  }
  for(int nq=0; nq!=kBBC_NQR; ++nq) {
    for(int ce=0; ce!=kBinsCen; ++ce) {
      for(int se=0; se!=3; ++se) {
        hPsiC[se][nq][ce] = NULL;
        hPsiS[se][nq][ce] = NULL;
        hDeltaPsi[se][nq][ce] = NULL;
      }
      hRes[nq][ce] = NULL;
    }
  }
}

AT_BBC_EPC::~AT_BBC_EPC() {
  if(fListBBCEPC) delete fListBBCEPC;
}

void AT_BBC_EPC::MyInit() {
  std::cout << "AT_BBC_EPC::MyInit()" << std::endl;
  hEvents->GetXaxis()->SetBinLabel(10,"AT_BBC_EPC init");
  hEvents->GetXaxis()->SetBinLabel(11,"AT_BBC_EPC_good");
  //std::cout << "BinsCen " << kBinsCen << " MinBinCen " << kMinBinCen << std::endl;
  //std::cout << "BinsVtx " << kBinsVtx << " MinBinVtx " << kMinBinVtx << std::endl;
  fListBBCEPC = new TList();
  fListBBCEPC->SetOwner();
  fListBBCEPC->SetName("BBCEPC");

  fListCalib = new TList();
  fListCalib->SetOwner();
  fListCalib->SetName("Calib");
  fListBBCEPC->Add( fListCalib );

  fListDelta = new TList();
  fListDelta->SetOwner();
  fListDelta->SetName("Delta");
  fListBBCEPC->Add( fListDelta );

  fListReso = new TList();
  fListReso->SetOwner();
  fListReso->SetName("Resolution");
  fListBBCEPC->Add( fListReso );

  Analysis *ana = Analysis::Instance();
  int run = ana->RunNumber();
  for(int ice=0; ice!=kBinsCen; ++ice) { //centrality
    for(int ord=0; ord!=kBBC_NQ; ++ord) { // order
      for(int se=0; se!=3; ++se) { // subevent
	for(int st=0; st!=3; ++st) { //step
	  hQxVtx[ord][se][ice][st] = new TH2D(Form("QxVtx_ORD%d_SE%d_CB%02d_ST%d",ord,se,ice,st),
					      "BBCQx;VtxBin",
					      kBinsVtx, -0.5, kBinsVtx-0.5, 100, -80., +80. );
	  fListCalib->Add(hQxVtx[ord][se][ice][st]);
	  hQyVtx[ord][se][ice][st] = new TH2D(Form("QyVtx_ORD%d_SE%d_CB%02d_ST%d",ord,se,ice,st),
					      "BBCQy;VtxBin",
					      kBinsVtx, -0.5, kBinsVtx-0.5, 100, -80., +80. );
	  fListCalib->Add(hQyVtx[ord][se][ice][st]);
	}
      }
    }

    for(int ord=0; ord!=kBBC_NQR; ++ord) { // order
      for(int se=0; se!=3; ++se) { // subevent
        hPsiC[se][ord][ice] = new TProfile2D(Form("PsiC_SE%d_ORD%d_CB%02d",se,ord,ice),
					     "PsiC;VtxBin;Order",
                                             kBinsVtx, -0.5, kBinsVtx-0.5,
					     kFlatten, 0.5, kFlatten+0.5 );
        fListCalib->Add(hPsiC[se][ord][ice]);
        hPsiS[se][ord][ice] = new TProfile2D(Form("PsiS_SE%d_ORD%d_CB%02d",se,ord,ice),
					     "PsiS;VtxBin;Order",
                                             kBinsVtx, -0.5, kBinsVtx-0.5, 
					     kFlatten, 0.5, kFlatten+0.5 );
        fListCalib->Add(hPsiS[se][ord][ice]);
        hDeltaPsi[se][ord][ice] = new TH2D(Form("DeltaPsi_SE%d_ORD%d_CB%02d",se,ord,ice),
					   "DeltaPsi;Order;Delta",
                                           kFlatten, -0.5, kFlatten-0.5, 100, -0.2, +0.2 );
        fListDelta->Add(hDeltaPsi[se][ord][ice]);
      }

      hRes[ord][ice] = new TProfile(Form("Res_ORD%d_CB%02d",ord,ice),
				    "Res;VtxBin;Res",
                                    kBinsVtx, -0.5, kBinsVtx-0.5, -1, +1 );
      fListReso->Add(hRes[ord][ice]);
    }
  }
}

void AT_BBC_EPC::MyFinish() {
  fListBBCEPC->Write( fListBBCEPC->GetName(), kSingleKey );
}

void AT_BBC_EPC::MyExec() {
  hEvents->Fill(9);
  float vtx = fGLB.vtxZ;
  float cen = fGLB.cent;
  int bvtx = BinVertex( vtx );
  int bcen = BinCentrality( cen );

  qcQ qvec[kBBC_NQ][3];
  for(int se=0; se!=2; ++se) {
    qvec[0][se] = pQ1bb->at(se);
    qvec[1][se] = pQ2bb->at(se);
    qvec[2][se] = pQ3bb->at(se);
    qvec[3][se] = pQ4bb->at(se);
    qvec[4][se] = pQ5bb->at(se);
    qvec[5][se] = pQ6bb->at(se);
    qvec[6][se] = pQ8bb->at(se);
    qvec[7][se] = pQ10bb->at(se);
    qvec[8][se] = pQ12bb->at(se);
    if(qvec[0][se].M()<1) return;
  }
  hEvents->Fill(10);
  for(int k=0; k!=kBBC_NQ; ++k) { // order
    qvec[k][2] = qvec[k][0] + qvec[k][1];
  }

  // ======= STAGE 1: Storing Raw Centroids =======
  for(int se=0; se!=3; ++se) { // subevent
    for(int k=0; k!=kBBC_NQ; ++k) { //order
      hQxVtx[k][se][bcen][0]->Fill( bvtx, qvec[k][se].X() );
      hQyVtx[k][se][bcen][0]->Fill( bvtx, qvec[k][se].Y() );
    }
  }

  // ======= STAGE 2: Recentering SubEvents (STEP1)  =======
  for(int k=0; k!=kBBC_NQ; ++k) { // order
    for(int se=0; se!=3; ++se) { // subevent
      double x = qvec[k][se].X();
      double y = qvec[k][se].Y();
      double cn = fBBCm[se][k][0][bcen][bvtx][0];
      double sn = fBBCm[se][k][1][bcen][bvtx][0];
      qvec[k][se].SetXY( x - cn, y - sn, qvec[k][se].NP(), qvec[k][se].M() );
    }
  }

  // ======= STAGE 3: Storing Prime Centroids =======
  for(int se=0; se!=3; ++se) { // subevent
    for(int k=0; k!=kBBC_NQ; ++k) { //order
      hQxVtx[k][se][bcen][1]->Fill( bvtx, qvec[k][se].X() );
      hQyVtx[k][se][bcen][1]->Fill( bvtx, qvec[k][se].Y() );
    }
  }

  int twon[6] = {1,3,5,6,7,8}; // [1,2,3,4,5,6] => [2,4,6,8,10,12]
  // ======= STAGE 4: Twisting SubEvents (STEP2)  =======
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

  // ======= STAGE 5: Storing Doubly-Prime Centroids =======
  for(int se=0; se!=3; ++se) { // subevent
    for(int k=0; k!=kBBC_NQ; ++k) { //order
      hQxVtx[k][se][bcen][2]->Fill( bvtx, qvec[k][se].X() );
      hQyVtx[k][se][bcen][2]->Fill( bvtx, qvec[k][se].Y() );
    }
  }

  // ======= STAGE 6: Rescaling SubEvents (STEP3)  =======
  for(int k=0; k!=kBBC_NQR; ++k) { // order
    for(int se=0; se!=3; ++se) { // subevent
      double x = qvec[k][se].X();
      double y = qvec[k][se].Y();
      double c2n = fBBCm[se][twon[k]][0][bcen][bvtx][2] / qvec[k][se].M();
      double a2np = 1.0+c2n/2; //1.0+c2n;
      double a2nm = 1.0-c2n/2; //1.0-c2n;
      qvec[k][se].SetXY( x / a2np,
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

  // ======= STAGE 7: Storing Triple-Prime Centroids =======


  // ======= STAGE 8: Bulding Full Q and Storing Flattening Coeficients  =======
  for(int k=0; k!=kBBC_NQR; ++k) { // order
    for(int se=0; se!=3; ++se) {
      for(int ik=0; ik!=32; ++ik) { // correction order
        int nn = 1+ik;
        hPsiC[se][k][bcen]->Fill(bvtx, nn, TMath::Cos(nn*qvec[k][se].Psi2Pi()) );
        hPsiS[se][k][bcen]->Fill(bvtx, nn, TMath::Sin(nn*qvec[k][se].Psi2Pi()) );
      }
      double psi = qvec[k][se].Psi2Pi();
      for(int ik=0; ik!=32; ++ik) { // correction order
        int nn = ik+1;
        double prime = 0.0;
        prime += TMath::Cos(nn*psi)*fBBCc[se][ik][k][bcen][bvtx]*2.0/nn;
        prime += TMath::Sin(nn*psi)*fBBCs[se][ik][k][bcen][bvtx]*2.0/nn;
        hDeltaPsi[se][k][bcen]->Fill(ik, prime);
      }
    }
    hRes[k][bcen]->Fill(bvtx, TMath::Cos( (k+1)*(qvec[k][0].Psi2Pi()-qvec[k][1].Psi2Pi()) ) );
  }

  //std::cout << "STEP5" << std::endl;
}
