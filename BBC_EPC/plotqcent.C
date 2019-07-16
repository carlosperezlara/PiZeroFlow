#include "../RunningConf.h"
const int kBBC_NQ = 9;
int cenXvtx = kBinsVtx*kBinsCen;
int flaXvtx = kBinsVtx*kFlatten;
int cenXflaXvtx = kBinsCen*kBinsVtx*kFlatten;
TH2D *data[3][2][kBBC_NQ];

void reader(int run) {
  ifstream fin( Form("tables/BBC_%d.dat",run) );
  Double_t tmp;
  for(int se=0; se!=3; ++se) {
    for(int xy=0; xy!=2; ++xy) {
      for(int nq=0; nq!=kBBC_NQ; ++nq) {
        if(data[se][xy][nq]) delete data[se][xy][nq];
        data[se][xy][nq] = new TH2D( Form("DATA_%d_%d_nq%d",se,xy,nq), Form("SE%d XY%d NQ%d",se,xy,nq), kBinsVtx, -0.5, kBinsVtx-0.5, kBinsCen, -0.5, kBinsCen-0.5 );
      }
    }
  }
  int n=0;
  int nn=0;
  for(;;++nn) {
    fin >> tmp;
    if(!fin.good()) break;
    //=== EXTRACTION CODE
    int ord = (nn/(cenXvtx*6))%kBBC_NQ;
    int xy = (nn/(cenXvtx*3))%2;
    int se = (nn/cenXvtx)%3;
    int bce = (nn/kBinsVtx)%kBinsCen;
    int bvt = nn%kBinsVtx;
    //===
    data[se][xy][ord]->SetBinContent(bvt+1,bce+1,tmp);
  }
  //cout << "Points " << nn << " from a total of " << nn << " values " << endl;
}

int plotqcent(int run=454777) {
  gStyle->SetOptStat(0);
  reader(run);
  TCanvas *main[kBBC_NQ];
  char xoy[2] = {'x','y'};
  int order[9] = {1,2,3,4,5,6,8,10,12};
  for(int nq=0; nq!=kBBC_NQ; ++nq) {
    main[nq]= new TCanvas();
    main[nq]->Divide(2,2);
    for(int se=0; se!=2; ++se) {
      for(int xy=0; xy!=2; ++xy) {
        main[nq]->cd(se*2+xy+1);
        data[se][xy][nq]->Draw("COLZ");
	data[se][xy][nq]->GetXaxis()->SetTitle( "Vertex  bin  index" );
	data[se][xy][nq]->GetYaxis()->SetTitle( "Centrality  bin  index" );
	data[se][xy][nq]->SetTitle( Form("Run %d  <Q%c>  Order %d  SubEvent %d",run,xoy[xy],order[nq],se) );
      }
    }
    main[nq]->SaveAs( Form("QC_%d_ORD%d.pdf",run,nq) );
  }
  return 0;
}

int plotallqcent() {
  ifstream fin("../runs.dat");
  Double_t sumxdata[200];
  Double_t sumydata[2][2][kBBC_NQ][200];
  int nn=0;
  for(;;++nn) {
    int irun;
    if(!fin.good()) break;
    fin >> irun;
    reader(irun);
    sumxdata[nn] = irun;
    for(int nq=0; nq!=kBBC_NQ; ++nq) {
      for(int se=0; se!=2; ++se) {
        for(int xy=0; xy!=2; ++xy) {
          sumydata[se][xy][nq][nn] = data[se][xy][nq]->GetBinContent(20,4);
        }
      }
    }
  }

  TCanvas *main[kBBC_NQ];
  TGraph *gr[2][2][kBBC_NQ];
  for(int nq=0; nq!=kBBC_NQ; ++nq) {
    main[nq]= new TCanvas();
    main[nq]->Divide(2,2);
    for(int se=0; se!=2; ++se) {
      for(int xy=0; xy!=2; ++xy) {
        main[nq]->cd(se*2+xy+1);
        gr[se][xy][nq] = new TGraph(nn,sumxdata,sumydata[se][xy][nq]);
        gr[se][xy][nq]->Draw("A*");
      }
    }
  }

  return 0;
}
