TH2F *hist2[200];
TH2F *hist2AC[200];
TH1F *hist[200];
TH1F *histAC[200];

void LoadTH2Ftime(int idx, int run) {
  TFile *file = new TFile( Form("out/run%d.root",run) );
  hist2[idx] = (TH2F*) file->Get("hBBCTimeMean2D0");
  hist2AC[idx] = (TH2F*) file->Get("hBBCTimeMean2D1");
  if(!hist2[idx]) return;
  hist2[idx]->SetTitle( Form("%d",run) );
  hist2[idx]->GetXaxis()->SetRangeUser(-2,20);
  hist2[idx]->GetYaxis()->SetRangeUser(-2,20);
  hist2AC[idx]->SetTitle( Form("%d",run) );
  hist2AC[idx]->GetXaxis()->SetRangeUser(-2,20);
  hist2AC[idx]->GetYaxis()->SetRangeUser(-2,20);
}

void LoadTH1Ftime(int idx, int run) {
  TFile *file = new TFile( Form("out/run%d.root",run) );
  if(!file) {
    hist[idx] = NULL;
    histAC[idx] = NULL;
    return;
  }
  hist[idx] = (TH2F*) file->Get("hBBCTimeMean0");
  histAC[idx] = (TH2F*) file->Get("hBBCTimeMean1");
  if(!hist[idx]) return;
  hist[idx]->SetTitle( Form("%d",run) );
  hist[idx]->SetLineColor( kBlue-3 );
  histAC[idx]->SetTitle( Form("%d",run) );
  histAC[idx]->SetLineColor( kRed-3 );
}

void LoadTH2Ftimerms(int idx, int run) {
  TFile *file = new TFile( Form("out/run%d.root",run) );
  hist2[idx] = (TH2F*) file->Get("hBBCTimeRMS2D0");
  hist2AC[idx] = (TH2F*) file->Get("hBBCTimeRMS2D1");
  if(!hist2[idx]) return;
  hist2[idx]->SetTitle( Form("%d",run) );
  hist2[idx]->GetXaxis()->SetRangeUser(0,0.5);
  hist2[idx]->GetYaxis()->SetRangeUser(0,0.5);
  hist2AC[idx]->SetTitle( Form("%d",run) );
  hist2AC[idx]->GetXaxis()->SetRangeUser(0,0.5);
  hist2AC[idx]->GetYaxis()->SetRangeUser(0,0.5);
}

void LoadTH1Ftimerms(int idx, int run) {
  TFile *file = new TFile( Form("out/run%d.root",run) );
  hist[idx] = (TH2F*) file->Get("hBBCTimeRMS0");
  histAC[idx] = (TH2F*) file->Get("hBBCTimeRMS1");
  if(!hist[idx]) return;
  hist[idx]->SetTitle( Form("%d",run) );
  hist[idx]->SetLineColor( kBlue-3 );
  histAC[idx]->SetTitle( Form("%d",run) );
  histAC[idx]->SetLineColor( kRed-3 );
}

int plotQAtimevstime() {
  gStyle->SetOptStat(0);
  ifstream fin("../runs.anal.dat");
  int run, nruns=0;
  for(;;) {
    fin >> run;
    if(!fin.good()) break;
    LoadTH2Ftime(nruns,run);
    nruns++;
  }
  cout << nruns << endl;

  TCanvas *main = new TCanvas();
  main->Divide(11,10,0,0);
  for(int i=0; i!=nruns; ++i) {
    if(i>=11*10) break;
    main->cd(i+1);
    if(hist2[i]) hist2[i]->Draw("colz");
  }
  main->SaveAs("timevstime.pdf");

  TCanvas *main1 = new TCanvas();
  main1->Divide(11,10,0,0);
  for(int i=0; i!=nruns; ++i) {
    if(i>=11*10) break;
    main1->cd(i+1);
    if(hist2AC[i]) hist2AC[i]->Draw("colz");
  }
  main1->SaveAs("timevstimeAC.pdf");

  return 0;
}

int plotQAtime() {
  gStyle->SetOptStat(0);
  ifstream fin("../runs.anal.dat");
  int run, nruns=0;
  for(;;) {
    fin >> run;
    if(!fin.good()) break;
    LoadTH1Ftime(nruns,run);
    nruns++;
  }
  cout << nruns << endl;

  TCanvas *main = new TCanvas();
  main->Divide(11,10,0,0);
  for(int i=0; i!=nruns; ++i) {
    if(i>=11*10) break;
    main->cd(i+1)->SetLogy(1);
    if(hist[i]) hist[i]->Draw("");
    if(histAC[i]) histAC[i]->Draw("same");
  }
  main->SaveAs("time.pdf");

  return 0;
}

int plotQArmsvsrms() {
  gStyle->SetOptStat(0);
  ifstream fin("../runs.anal.dat");
  int run, nruns=0;
  for(;;) {
    fin >> run;
    if(!fin.good()) break;
    LoadTH2Ftimerms(nruns,run);
    nruns++;
  }
  cout << nruns << endl;

  TCanvas *main = new TCanvas();
  main->Divide(11,10,0,0);
  for(int i=0; i!=nruns; ++i) {
    if(i>=11*10) break;
    main->cd(i+1);
    if(hist2[i]) hist2[i]->Draw("colz");
  }
  main->SaveAs("rmsvsrms.pdf");

  TCanvas *main1 = new TCanvas();
  main1->Divide(11,10,0,0);
  for(int i=0; i!=nruns; ++i) {
    if(i>=11*10) break;
    main1->cd(i+1);
    if(hist2AC[i]) hist2AC[i]->Draw("colz");
  }
  main1->SaveAs("rmsvsrmsAC.pdf");

  return 0;
}

int plotQArms() {
  gStyle->SetOptStat(0);
  ifstream fin("../runs.anal.dat");
  int run, nruns=0;
  for(;;) {
    fin >> run;
    if(!fin.good()) break;
    LoadTH1Ftimerms(nruns,run);
    nruns++;
  }
  cout << nruns << endl;

  TCanvas *main = new TCanvas();
  main->Divide(11,10,0,0);
  for(int i=0; i!=nruns; ++i) {
    if(i>=11*10) break;
    main->cd(i+1)->SetLogy(1);
    if(hist[i]) hist[i]->Draw("");
    if(histAC[i]) histAC[i]->Draw("same");
  }
  main->SaveAs("rms.pdf");

  return 0;
}

