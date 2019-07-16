int nn;
float rr[150];
float zz[150];
float aa[6][60][150];
float mm[6][60][150];
float ea[6][60][150];
float em[6][60][150];

void loadres(int run=454777) {
  ifstream fin( Form("tables/BBC_R_%d.dat",run) );
  bool breakme = false;
  for(int ord=0; ord!=6; ++ord) {
    for(int i=0; i!=60; ++i) {
      fin >> aa[ord][i][ nn ];
      if(!fin.good()) {
	breakme = true;
	break;
      }
      fin >> ea[ord][i][ nn ];
      fin >> mm[ord][i][ nn ] >> em[ord][i][ nn ];
    }
    if(breakme) break;
  }
  fin.close();
  rr[nn] = run;
  zz[nn] = 0;
  nn++;
}

void plotres1() {
  nn=0;
  ifstream runs("../runs.dat");
  for(;;) {
    int run;
    runs >> run;
    if(!runs.good()) break;
    loadres(run);
  }
  runs.close();

  int ord=1;
  int cen=1;
  new TCanvas();
  TGraphErrors *grA = new TGraphErrors(nn,rr,aa[ord][cen],zz,ea[ord][cen]);
  grA->Draw("A*");

  new TCanvas();
  TGraphErrors *grM = new TGraphErrors(nn,rr,mm[ord][cen],zz,em[ord][cen]);
  grM->Draw("A*");
}

void plotres0() {
  nn=0;
  loadres(0);
  int color[4] = {kBlue-3, kCyan-3, kOrange-3, kRed-3};
  float xx[60], ex[60];
  float ZZaa[4][60];
  float ZZmm[4][60];
  float ZZea[4][60];
  float ZZem[4][60];
  TGraphErrors *grA[4];
  TGraphErrors *grM[4];
  for(int ord=0; ord!=4; ++ord) {
    for(int i=0; i!=60; ++i) {
      ZZaa[ord][i] = aa[ord][i][0];
      ZZea[ord][i] = ea[ord][i][0];
      ZZmm[ord][i] = mm[ord][i][0];
      ZZem[ord][i] = em[ord][i][0];
      xx[i] = i;
      ex[i] = 0;
    }
    grA[ord] = new TGraphErrors(60,xx,ZZaa[ord],ex,ZZea[ord]);
    grM[ord] = new TGraphErrors(60,xx,ZZmm[ord],ex,ZZem[ord]);
    grA[ord]->SetLineColor( color[ord] );
    grM[ord]->SetLineColor( color[ord] );
    grA[ord]->SetMarkerColor( color[ord] );
    grM[ord]->SetMarkerColor( color[ord] );
    grA[ord]->SetMarkerStyle( 20 );
    grM[ord]->SetMarkerStyle( 20 );
  }  
  new TCanvas();
  grA[1]->Draw("AP");
  grA[0]->Draw("PSAME");
  grA[2]->Draw("PSAME");
  grA[3]->Draw("PSAME");
  new TCanvas();
  grM[1]->Draw("AP");
  grM[0]->Draw("PSAME");
  grM[2]->Draw("PSAME");
  grM[3]->Draw("PSAME");

}

int plotres() {
  plotres1();
  return 0;
}
