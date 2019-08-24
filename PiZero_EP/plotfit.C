int plotfit(TString cut="NOM") {
  gStyle->SetOptStat(0);
  int pts=0;
  double pt[100]; double ept[100];
  double rw[100]; double erw[100];
  double mu[100]; double emu[100];
  double sg[100]; double esg[100];
  double vnEP[100]; double evnEP[100];
  double vnSP[100]; double evnSP[100];
  double vnQC2[100]; double evnQC2[100];
  double vnQC4[100]; double evnQC4[100];
  double vnEPs[100]; double evnEPs[100];
  double zero[100];
  for(int i=0; i!=100; ++i) {
    zero[i] = 0;
  }

  double ptchpi[8] = { 0.49179, 0.69111, 0.89019, 1.14896, 1.57488, 1.95116, 2.37066, 2.77607 };
  double v2chpi[8] = { 0.04727, 0.07161, 0.08611, 0.10316, 0.10815, 0.13724, 0.15455, 0.11416 };
  double stchpi[8] = { 0.00163, 0.00217, 0.00291, 0.00310, 0.00523, 0.00886, 0.01546, 0.02720 };
  double sychpi[8] = { 0.00284, 0.00430, 0.00517, 0.00619, 0.00649, 0.00823, 0.00927, 0.00685 };
  TGraphErrors *gstchpi = new TGraphErrors( 8, ptchpi, v2chpi, zero, stchpi );
  TGraphErrors *gsychpi = new TGraphErrors( 8, ptchpi, v2chpi, zero, sychpi );
  gstchpi->SetLineColor( kRed-3 );
  gstchpi->SetFillColor( kRed-3 );
  gstchpi->SetMarkerColor( kRed-3 );
  gstchpi->SetMarkerStyle( 21 );
  gsychpi->SetLineColor( kRed-3 );
  gsychpi->SetFillColor( kRed-3 );
  gsychpi->SetMarkerColor( kRed-3 );
  gsychpi->SetMarkerStyle( 21 );

  ifstream fin;
  for(int ptbin=2; ptbin!=19; ++ptbin) {
    double tmp;
    fin.open( Form("fit/res/%s_EP_PB%02d_DF.res",cut.Data(),ptbin) );
    fin >> tmp >> pt[pts] >> pt[pts+1];
    if(!fin.good()) continue;
    fin >> rw[pts] >> erw[pts];
    fin >> mu[pts] >> emu[pts];
    fin >> sg[pts] >> esg[pts];
    fin >> vnEP[pts] >> evnEP[pts];
    cout << Form("%.1f  %.1f   |   %f  %f",pt[pts], pt[pts+1], vnEP[pts], evnEP[pts]) << endl;
    rw[pts] /= (pt[pts+1]-pt[pts]);
    erw[pts] /= (pt[pts+1]-pt[pts]);
    ept[pts] = (pt[pts+1]-pt[pts])/2;
    pt[pts] = (pt[pts+1]+pt[pts])/2;
    vnEPs[pts] = vnEP[pts]*4;
    evnEPs[pts] = evnEP[pts]*4;
    fin.close();
    pts++;
  }
  /*
  int ptss=0;
  for(int ptbin=2; ptbin!=20; ++ptbin) {
    double tmp;
    fin.open( Form("fit/res/SP_PB%02d_DF.res",ptbin) );
    fin >> tmp >> tmp >> tmp;
    fin >> tmp >> tmp;
    fin >> tmp >> tmp;
    fin >> tmp >> tmp;
    fin >> vnSP[ptss] >> evnSP[ptss];
    fin.close();
    //cout << Form("%.1f  %.1f",pt[pts], pt[pts+1]) << endl;
    ptss++;
  }
  */
  TGraphErrors *grw = new TGraphErrors(pts,pt,rw,ept,erw);
  TGraphErrors *gmu = new TGraphErrors(pts,pt,mu,ept,emu);
  TGraphErrors *gsg = new TGraphErrors(pts,pt,sg,ept,esg);
  TGraphErrors *gvnEP = new TGraphErrors(pts,pt,vnEP,ept,evnEP);
  //TGraphErrors *gvnSP = new TGraphErrors(pts,pt,vnSP,ept,evnSP);
  //TGraphErrors *gvnEPs = new TGraphErrors(pts,pt,vnEPs,ept,evnEPs);

  TCanvas *mainSGN = new TCanvas();
  mainSGN->Divide(3,1);
  mainSGN->cd(1)->SetLogy(1);
  grw->Draw("A*");
  mainSGN->cd(2);
  gmu->Draw("A*");
  mainSGN->cd(3);
  gsg->Draw("A*");

  TCanvas *mainVN = new TCanvas();
  TH2F *axisvn = new TH2F("axisvn",";p_{T}  GeV/c;v_{2}",100,0,10,100,-0.15,0.30);
  axisvn->Draw();
  gsychpi->Draw("PSAME");
  gstchpi->Draw("PSAME");
  gvnEP->Draw("PSAME");
  //gvnSP->Draw("PSAME");
  //gvnEPs->Draw("PSAME");
  gvnEP->SetMarkerStyle(20);
  //gvnSP->SetMarkerStyle(20);
  //gvnEPs->SetMarkerStyle(24);
  gsychpi->SetFillColor(kWhite);
  gstchpi->SetFillColor(kWhite);
  gvnEP->SetFillColor(kWhite);
  //gvnSP->SetFillColor(kWhite);
  //gvnEPs->SetFillColor(kWhite);
  //gvnSP->SetLineColor( kGreen-3 );
  //gvnSP->SetMarkerColor( kGreen-3 );


  TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->AddEntry(gstchpi,"charged pions");
  leg->AddEntry(gvnEP,"neutral pions EP");
  //leg->AddEntry(gvnSP,"pi0 SP");
  //leg->AddEntry(gvnEPs,"pi0 EP (x4)");
  leg->Draw();

  return 0;
}
