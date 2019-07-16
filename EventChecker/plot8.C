int plot8() {
  Int_t run[10000];
  Double_t x[10000];
  Double_t mean[10000];
  Double_t sigma[10000];
  Double_t mpvS[10000];
  Double_t mpvN[10000];
  Double_t sgS[10000];
  Double_t sgN[10000];
  Double_t cdfS[10000];
  Double_t cdfN[10000];

  ifstream fin("TimeConstants.dat");
  int nruns = 0;
  for(;;) {
    fin >> run[nruns];
    if(!fin.good()) break;
    fin >> mean[nruns] >> sigma[nruns];
    fin >> mpvS[nruns] >> sgS[nruns] >> cdfS[nruns];
    fin >> mpvN[nruns] >> sgN[nruns] >> cdfN[nruns];
    x[nruns] = run[nruns];//nruns;
    nruns++;
  }
  
  TGraph *gm = new TGraph( nruns, x, mean );
  TGraph *gs = new TGraph( nruns, x, sigma );
  TGraph *gmS = new TGraph( nruns, x, mpvS );
  TGraph *gsS = new TGraph( nruns, x, sgS );
  TGraph *gcS = new TGraph( nruns, x, cdfS );
  TGraph *gmN = new TGraph( nruns, x, mpvN );
  TGraph *gsN = new TGraph( nruns, x, sgN );
  TGraph *gcN = new TGraph( nruns, x, cdfN );

  TCanvas *main = new TCanvas();
  main->Divide(3,3);
  main->cd(1); gm->Draw("A*L");
  main->cd(2); gs->Draw("A*L");
  main->cd(4); gmS->Draw("A*L");
  main->cd(5); gsS->Draw("A*L");
  main->cd(6); gcS->Draw("A*L");
  main->cd(7); gmN->Draw("A*L");
  main->cd(8); gsN->Draw("A*L");
  main->cd(9); gcN->Draw("A*L");

}
