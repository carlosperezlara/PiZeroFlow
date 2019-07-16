TF1 *CreateFit(TString name) {
  TF1 *ret = new TF1(name.Data(),"[0]*TMath::Landau(x,[1],[2])",0,0.5);
  ret->SetParameter(1,0.06); ret->SetParLimits(1,0.01,0.1);
  ret->SetParameter(2,0.01); ret->SetParLimits(2,0.001,0.1);
  return ret;
}

int ComputeTimeConstants(int run=455547) {
  gStyle->SetOptStat(0);
  TFile *file = new TFile( Form("out/run%d.root",run) );
  file->ls();
  TH1F *hM0 = (TH1F*) file->Get("hBBCTimeMean0");
  TH1F *hS0 = (TH1F*) file->Get("hBBCTimeRMS_S0");
  TH1F *hN0 = (TH1F*) file->Get("hBBCTimeRMS_N0");
  TF1 *fS0 = CreateFit("fS0");
  TF1 *fN0 = CreateFit("fN0");
  TF1 *fS0cdf = new TF1("fS0cdf","[0]*ROOT::Math::landau_cdf(x,[2],[1])",0,1.0);
  TF1 *fN0cdf = new TF1("fN0cdf","[0]*ROOT::Math::landau_cdf(x,[2],[1])",0,1.0);
  TF1 *fM0 = new TF1("fM0",
		     "[0]*TMath::Gaus(x,[1],[2])",
		     -3,+3);
  fM0->SetParameter(0,1.0); fM0->SetParLimits(0,0,100000);
  fM0->SetParameter(1,0.0); fM0->SetParLimits(1,-2.,+2.);
  fM0->SetParameter(2,0.2); fM0->SetParLimits(2,0.1,1.0);
  new TCanvas();
  hM0->Draw();
  hM0->Fit("fM0","RN","",-3,+3);

  TLatex *tex = new TLatex();
  tex->SetTextSize(0.05);
  tex->SetTextColor(kGray);
  TLine *lin = new TLine();
  lin->SetLineColor(kGray);

  TCanvas *saveS = new TCanvas();
  hS0->Draw();
  hS0->Fit("fS0","RN","",0.035.,0.3);
  //fS0->Draw("same");
  fS0cdf->SetParameter(0,70000);//fS0->GetParameter(0));
  fS0cdf->SetParameter(1,fS0->GetParameter(1));
  fS0cdf->SetParameter(2,fS0->GetParameter(2));
  fS0cdf->Draw("same");
  lin->DrawLine( fS0cdf->GetX(70000*0.95),0,fS0cdf->GetX(70000*0.95),70000 );
  lin->DrawLine( 0,70000,0.5,70000 );
  tex->DrawLatex(0.51,70000, "100\%" );
  tex->DrawLatexNDC(0.6,0.5, Form("MPV  %.3f", fS0->GetParameter(1)) );
  tex->DrawLatexNDC(0.6,0.4, Form("Sigma  %.3f", fS0->GetParameter(2) ) );
  tex->DrawLatexNDC(0.6,0.3, Form("95per  %.3f", fS0cdf->GetX(70000*0.95) ) );
  saveS->SaveAs( Form("rmsSouth_%d.pdf",run), "pdf" );

  TCanvas *saveN = new TCanvas();
  hN0->Draw();
  hN0->Fit("fN0","RN","",0.035.,0.3);
  //fN0->Draw("same");
  fN0cdf->SetParameter(0,25000);//fN0->GetParameter(0));
  fN0cdf->SetParameter(1,fN0->GetParameter(1));
  fN0cdf->SetParameter(2,fN0->GetParameter(2));
  fN0cdf->Draw("same");
  lin->DrawLine( fN0cdf->GetX(25000*0.95),0,fN0cdf->GetX(25000*0.95),25000 );
  lin->DrawLine( 0,25000,0.5,25000 );
  tex->DrawLatex(0.51,25000, "100\%" );
  tex->DrawLatexNDC(0.6,0.5, Form("MPV  %.3f", fN0->GetParameter(1)) );
  tex->DrawLatexNDC(0.6,0.4, Form("Sigma  %.3f", fN0->GetParameter(2) ) );
  tex->DrawLatexNDC(0.6,0.3, Form("95per  %.3f", fN0cdf->GetX(25000*0.95) ) );
  saveN->SaveAs( Form("rmsNorth_%d.pdf",run),"pdf" );

  ofstream fout( Form("dat/TimeConstants%d.dat",run) );
  fout << run << " ";
  fout << Form("%0.2f",fM0->GetParameter(1)) << " ";
  fout << Form("%0.3f",fM0->GetParameter(2)) << " ";
  fout << Form("%0.4f",fS0->GetParameter(1)) << " ";
  fout << Form("%0.4f",fS0->GetParameter(2)) << " ";
  fout << Form("%0.3f",fS0cdf->GetX(70000*0.95)) << " ";
  fout << Form("%0.4f",fN0->GetParameter(1)) << " ";
  fout << Form("%0.4f",fN0->GetParameter(2)) << " ";
  fout << Form("%0.3f",fN0cdf->GetX(25000*0.95)) << " " << endl;
  fout.close();


  return 0;
}
