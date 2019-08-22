int plotMethod(TString method="EP", TString cut="NOM") {
  TFile *file = new TFile( Form("allfiles/all%s_%s_Ord1.root",cut.Data(),method.Data()) );
  TList *list = (TList*) file->Get( Form("%s_Ord1_Psi1",method.Data()) );
  TList *listR = (TList*) list->FindObject("results");
  TList *listY = (TList*) list->FindObject("yields");
  listY->ls();
  TH1D *binPT[100];
  TH1D *sgnPT[100];
  TH1D *mixPT[100];
  TH1D *mixPTL[100];
  TH1D *mixPTR[100];
  TH1D *D2BINNED[100];
  TH1D *D2UNBINNED[100];
  TH1D *PV2BINNED[100];
  TH1D *PV2UNBINNED[100];
  TCanvas *main = new TCanvas();
  main->Divide(5,5);
  TCanvas *main2 = new TCanvas();
  main2->Divide(5,5);
  TCanvas *main3 = new TCanvas();
  main3->Divide(5,5);
  TCanvas *mainPUB = new TCanvas("pub","pub");
  mainPUB->Divide(2,1);
  for(int pt=2; pt!=19; ++pt) {
    //Double_t minPt = map->GetXaxis()->GetBinLowEdge( pt+1 );
    //Double_t maxPt = map->GetXaxis()->GetBinLowEdge( pt+2 );
    binPT[pt] = (TH1D*) listY->FindObject( Form("Yield_PB%d",pt) );
    mixPT[pt] = (TH1D*) listY->FindObject( Form("hMass2_PB%d",pt) );
    mixPTL[pt] = (TH1D*) mixPT[pt]->Clone( Form("hMass2_PB%dL",pt) );
    mixPTR[pt] = (TH1D*) mixPT[pt]->Clone( Form("hMass2_PB%dR",pt) );
    sgnPT[pt] = (TH1D*) binPT[pt]->Clone( Form("SGN_PB%dR",pt) );
    //sgnPT[pt]->Reset();
    int bin10 = mixPT[pt]->GetXaxis()->FindBin( 0.050 );
    int bin20 = mixPT[pt]->GetXaxis()->FindBin( 0.100 );
    int bin11 = mixPT[pt]->GetXaxis()->FindBin( 0.176 );
    int bin21 = mixPT[pt]->GetXaxis()->FindBin( 0.250 );
    Double_t counts1L = mixPT[pt]->Integral(bin10,bin20);
    Double_t counts1R = mixPT[pt]->Integral(bin11,bin21);
    Double_t counts2L = binPT[pt]->Integral(bin10,bin20);
    Double_t counts2R = binPT[pt]->Integral(bin11,bin21);
    mixPTL[pt]->Scale(counts2L/counts1L);
    mixPTR[pt]->Scale(counts2R/counts1R);
    mixPT[pt]->Reset();
    int thL = binPT[pt]->GetXaxis()->FindBin( 0.110 );
    int thR = binPT[pt]->GetXaxis()->FindBin( 0.166 );
    for(int iii=0; iii!=mixPT[pt]->GetXaxis()->GetNbins()+2; ++iii) { //under and over flows
      float fracL, fracR;
      if(iii<thL) {
	fracL = 1.0;
	fracR = 0.0;
      } else if(iii>thR) {
	fracL = 0.0;
	fracR = 1.0;
      } else {
	fracL = (thR-iii)*1.0/(thR-thL);
	fracR = 1-fracL;
      }
      mixPT[pt]->SetBinContent( iii,
				fracL*mixPTL[pt]->GetBinContent(iii) +
				fracR*mixPTR[pt]->GetBinContent(iii) );
    }
    sgnPT[pt]->Add(mixPT[pt],-1.0);
    PV2UNBINNED[pt] = (TH1D*) listR->FindObject( Form("PV2_%d_%s_Ord1_Psi1_UNBINNED",pt,method.Data()) );
    PV2BINNED[pt] = (TH1D*) listR->FindObject( Form("PV2_%d_%s_Ord1_Psi1_BINNED",pt,method.Data()) );
    //PV2UNBINNED[pt]->SetTitle( Form("[%.1f - %.1f ]",minPt,maxPt) );
    //PV2BINNED[pt]->SetTitle( Form("[%.1f - %.1f ]",minPt,maxPt) );
    PV2BINNED[pt]->SetLineColor(kRed-3);
    mixPT[pt]->SetLineColor(kRed-3);
    mixPTL[pt]->SetLineColor(kYellow-3);
    mixPTR[pt]->SetLineColor(kYellow-3);
    main->cd(pt-2);
    binPT[pt]->Draw();
    mixPTL[pt]->Draw("SAME");
    mixPTR[pt]->Draw("SAME");
    mixPT[pt]->Draw("SAME");
    main2->cd(pt-2);
    sgnPT[pt]->Draw();
    main3->cd(pt-2);
    PV2UNBINNED[pt]->Draw();
    PV2BINNED[pt]->Draw("SAME");
    // recording
    mainPUB->cd(1);
    binPT[pt]->Draw();
    mixPT[pt]->Draw("SAME");
    mainPUB->cd(2);
    PV2BINNED[pt]->Draw();
    mainPUB->SaveAs( Form("fit/%s_%s_PB%02d.root",cut.Data(),method.Data(),pt), "root" );
  }

}
