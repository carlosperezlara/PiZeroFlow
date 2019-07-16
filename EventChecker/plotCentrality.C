int plotCentrality() {
  gStyle->SetOptStat(0);
  TFile *file = new TFile("out/run455364.root");
  TH2F *mean2D_0 = (TH2F*) file->Get("hBBCTimeMean2D0");
  TH2F *mean2D_1 = (TH2F*) file->Get("hBBCTimeMean2D1");
  TH2F *rms2D_0 = (TH2F*) file->Get("hBBCTimeRMS2D0");
  TH2F *rms2D_1 = (TH2F*) file->Get("hBBCTimeRMS2D1");

  mean2D_0->GetXaxis()->SetRangeUser(-2,16);
  mean2D_0->GetYaxis()->SetRangeUser(-2,16);
  mean2D_1->GetXaxis()->SetRangeUser(-2,16);
  mean2D_1->GetYaxis()->SetRangeUser(-2,16);

  rms2D_1->GetXaxis()->SetRangeUser(0,0.35);
  rms2D_1->GetYaxis()->SetRangeUser(0,0.35);

  if(1) {
  TCanvas *canM0 = new TCanvas();
  canM0->SetLogz(1);
  mean2D_0->Draw("colz");
  canM0->SaveAs("dAu_mean2Draw.pdf","pdf");

  TCanvas *canM1 = new TCanvas();
  canM1->SetLogz(1);
  mean2D_1->Draw("colz");
  canM1->SaveAs("dAu_mean2Dac.pdf","pdf");

  TCanvas *canR0 = new TCanvas();
  canR0->SetLogz(1);
  rms2D_0->Draw("colz");
  canR0->SaveAs("dAu_rms2Draw.pdf","pdf");

  TCanvas *canR1 = new TCanvas();
  canR1->SetLogz(1);
  rms2D_1->Draw("colz");
  canR1->SaveAs("dAu_rms2Dac.pdf","pdf");
  }
  TH1F *mean1D_0 = (TH1F*) file->Get("hBBCTimeMean0");
  TH1F *mean1D_1 = (TH1F*) file->Get("hBBCTimeMean1");
  mean1D_0->SetLineColor(kBlue-3);
  mean1D_1->SetLineColor(kRed-3);
  mean1D_0->SetLineWidth(2);
  mean1D_1->SetLineWidth(2);

  TH1F *rms1D_0 = (TH1F*) file->Get("hBBCTimeRMS0");
  TH1F *rms1D_1 = (TH1F*) file->Get("hBBCTimeRMS1");
  rms1D_0->SetLineColor(kBlue-3);
  rms1D_1->SetLineColor(kRed-3);
  rms1D_0->SetLineWidth(2);
  rms1D_1->SetLineWidth(2);

  TH1F *frac1D_0 = (TH1F*) file->Get("hFrac0");
  TH1F *frac1D_1 = (TH1F*) file->Get("hFrac1");
  frac1D_0->SetLineColor(kBlue-3);
  frac1D_1->SetLineColor(kRed-3);
  frac1D_0->SetLineWidth(2);
  frac1D_1->SetLineWidth(2);

  TLegend *leg = new TLegend(0.12,0.72,0.27,0.88);
  leg->AddEntry(mean1D_0,"raw");
  leg->AddEntry(mean1D_1,"after cuts");


  TH2F *axisM = new TH2F("axisM",";mean  time  S - N  (a.u.);counts  (a.u.)",
			100,-5,+5,100,0,6e4);
  TCanvas *canM2 = new TCanvas();
  axisM->Draw();
  mean1D_0->Draw("SAME");
  mean1D_1->Draw("SAME");
  leg->Draw();
  canM2->RedrawAxis();
  canM2->SaveAs("dAu_mean1D.pdf","pdf");


  TH2F *axisR = new TH2F("axisR",";rms  time  S + N  (a.u.);counts  (a.u.)",
			100,0,+3,100,0,3.5e4);
  TCanvas *canR2 = new TCanvas();
  axisR->Draw();
  rms1D_0->Draw("SAME");
  rms1D_1->Draw("SAME");
  leg->Draw();
  canR2->RedrawAxis();
  canR2->SaveAs("dAu_rms1D.pdf","pdf");


  TH2F *axisF = new TH2F("axisF",";frac;counts  (a.u.)",
			100,0.5,+1,100,1,3e8);
  TCanvas *canF2 = new TCanvas();
  canF2->SetLogy(1);
  axisF->Draw();
  frac1D_0->Draw("SAME");
  frac1D_1->Draw("SAME");
  leg->Draw();
  canF2->RedrawAxis();
  canF2->SaveAs("dAu_frac1D.pdf","pdf");



  return 0;
}
