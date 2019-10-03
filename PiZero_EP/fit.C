#include "DoubleFitter.h"
#include "DoubleFitter.cxx"
#include "Model3.h"
#include "Model3.cxx"

Model3 *fitter1;

void FCN(Int_t& a, Double_t& b, Double_t &val, Double_t *par, Int_t e) {
  fitter1->FCN( a, b, val, par, e );
}

int fit(int ptbin=18, TString method="EP", TString cut="NOM") {
  Double_t minX = 0.055;
  Double_t maxX = 0.240;
  gStyle->SetOptStat(0);
  TFile *file = new TFile( Form("fit/%s_%s_PB%02d.root",cut.Data(),method.Data(),ptbin) );
  TCanvas *main = (TCanvas*) file->Get("pub");
  main->Draw();
  TList *list1 = ((TPad*) main->GetListOfPrimitives()->At(0))->GetListOfPrimitives();
  TList *list2 = ((TPad*) main->GetListOfPrimitives()->At(1))->GetListOfPrimitives();
  list1->ls();
  list2->ls();
  TH1D *yield = (TH1D*) list1->At(0);
  TString title = yield->GetTitle();
  TH1D *bgr = (TH1D*) list1->At(1);
  TH1D *flow = (TH1D*) list2->At(0);
  minX = flow->GetXaxis()->GetBinLowEdge(2);
  maxX = flow->GetXaxis()->GetBinLowEdge( flow->GetXaxis()->GetNbins() );
  TLatex *tex = new TLatex();
  tex->SetTextColor(kGreen-3);
  //cout << yield->GetTitle() << endl; 
  //cout << bgr->GetTitle() << endl; 
  //cout << flow->GetTitle() << endl; 
  //========================================
  //========================================
  //========================================
  TCanvas *mainBGR = new TCanvas();
  bgr->Draw("E");
  double bgrcounts = bgr->Integral();
  cout << "******* PT " << ptbin << " ENTRIES BGR "  << bgr->Integral() << endl;
  TString pol3 = Form("%s+%s*x+%s*(2*x*x-1)+%s*(4*x*x*x-3*x)","[0]","[1]","[2]","[3]");
  TString pol5 = pol3 + Form("+%s*(8*pow(x,4)-8*x*x+1)+%s*(16*pow(x,5)-20*pow(x,3)+5*x)","[4]","[5]");
  TString pol7 = pol5 + Form("+%s*(32*pow(x,6)-48*pow(x,4)+18*x*x-1)+%s*(64*pow(x,7)-112*pow(x,5)+56*x*x*x-7*x)","[6]","[7]");
  cout << pol7.Data() << endl;
  TF1 *bgrF = new TF1("bgrF",pol7.Data(),minX,maxX);
  if(bgrcounts<1000) {
    for(int ii=1; ii!=8; ++ii) {
      bgrF->SetParameter(ii,0);
      bgrF->SetParLimits(ii,+1,-1);
    }
  }
  bgr->Fit(bgrF,"LQR","",minX,maxX);
  bgr->Fit(bgrF,"EMIRW","",minX,maxX);
  cout << "Reduced chi squared " << bgrF->GetChisquare() / bgrF->GetNDF() << endl;
  //========================================
  //========================================
  //========================================
  TCanvas *mainR = new TCanvas();
  mainR->Divide(2,1);
  yield->SetMarkerStyle(20);
  yield->SetMarkerColor(kBlack);
  yield->SetLineColor(kBlack);
  yield->GetXaxis()->SetTitle( "Inv.Mass  [GeV/c^{2}]" );
  yield->SetTitle( Form( "#pi^{0}  candidates   %s GeV/c",yield->GetTitle()) );
  flow->SetMarkerStyle(20);
  flow->SetMarkerColor(kBlack);
  flow->SetLineColor(kBlack);
  flow->GetXaxis()->SetTitle( "Inv.Mass  [GeV/c^{2}]" );
  flow->SetTitle( "Anisotropy" );
  //////////////////////////
  TString pol3fix = Form("%f+%f*x+%f*(2*x*x-1)+%f*(4*x*x*x-3*x)",
  			 bgrF->GetParameter(0),
			 bgrF->GetParameter(1),
			 bgrF->GetParameter(2),
			 bgrF->GetParameter(3));
  TString pol5fix = pol3fix + Form("+%f*(8*pow(x,4)-8*x*x+1)+%f*(16*pow(x,5)-20*pow(x,3)+5*x)",
				   bgrF->GetParameter(4),
				   bgrF->GetParameter(5));
  TString pol7fix = pol5fix + Form("+%f*(32*pow(x,6)-48*pow(x,4)+18*x*x-1)+%f*(64*pow(x,7)-112*pow(x,5)+56*x*x*x-7*x)",
				   bgrF->GetParameter(6),
				   bgrF->GetParameter(7));
  cout << pol7fix.Data() << endl;  
  //////////////////////////
  //FIT ALTERNATIVO
  TF1 *yieldF = new TF1("yieldF",Form("[0]*TMath::Gaus(x,[1],[2])+[3]*(%s)",pol7fix.Data()));
  yieldF->SetLineColor(kGreen-3);
  yieldF->SetLineWidth(3);
  yieldF->SetParameter(1,0.138); yieldF->SetParLimits(1,0.132,0.139);
  yieldF->SetParameter(2,0.010); yieldF->SetParLimits(2,0.001,0.020);
  ///////////////////////////
  mainR->cd(1);
  yield->Draw("E");
  cout << "First attempt of fitting yield <L>" << endl;
  yield->Fit(yieldF,"LQR","",minX,maxX);
  cout << "Second attempt of fitting yield <EMI>" << endl;
  yield->Fit(yieldF,"EMIR","",minX,maxX);
  TString Spol7fix = Form("%f*(%s)",yieldF->GetParameter(3),pol7fix.Data());
  TF1 *SbgrF = new TF1("SbgrF",Spol7fix.Data(),minX,maxX);
  SbgrF->Draw("same");
  tex->DrawLatexNDC(0.15,0.25, Form("muSGN  %.3f GeV/c^{2}",yieldF->GetParameter(1)) );
  tex->DrawLatexNDC(0.15,0.20, Form("sgSGN  %.3f GeV/c^{2}",yieldF->GetParameter(2)) );
  //////////////////////////
  /*
  fitter1 = new Model3();
  fitter1->SetMassHistogram(yield);
  fitter1->SetFlowHistogram(flow);
  fitter1->Init();
  // sgn
  fitter1->SetParameterValue( 5, yieldF->GetParameter(0),true,0.1);
  fitter1->SetParameterValue( 6, yieldF->GetParameter(1));
  fitter1->SetParameterValue( 7, yieldF->GetParameter(2));
  // bgr
  fitter1->SetParameterValue( 8, bgrF->GetParameter(0),true);
  fitter1->SetParameterValue( 9, bgrF->GetParameter(1),true);
  fitter1->SetParameterValue(10, bgrF->GetParameter(2),true);
  fitter1->SetParameterValue(11, bgrF->GetParameter(3),true);
  fitter1->SetParameterValue(12, bgrF->GetParameter(4),true);
  fitter1->SetParameterValue(13, bgrF->GetParameter(5),true);
  fitter1->SetParameterValue(14, bgrF->GetParameter(6),true);
  fitter1->SetParameterValue(15, bgrF->GetParameter(7),true);
  fitter1->LoadIni( Form("fit/ini/EP_PB%02d_DF.ini",ptbin) );
  fitter1->Print( Form("fit/ini/EP_PB%02d_DF.ini",ptbin) );
  TVirtualFitter::SetDefaultFitter("Minuit");
  TVirtualFitter *minuit1 = TVirtualFitter::Fitter( NULL, fitter1->GetNumberOfParameters() );
  minuit1->SetFCN( FCN );
  fitter1->Fit( minuit1 );
  fitter1->SaveAs( Form("fit/res/EP_PB%02d_DF.res",ptbin)  );
  TF1 *exported_mass, *exported_flow, *exported_sgn, *exported_bgr;
  fitter1->ExportTF1(&exported_mass,&exported_flow,&exported_sgn,&exported_bgr);
  */
  mainR->cd(2);
  flow->Draw("E");
  TString gaus = Form("[0]*TMath::Gaus(x,[1],[2])");
  TF1 *flowF = new TF1("flowF", Form("([3]*%s+([4]+[5]*(x-[1])+[6]*pow(x-[1],2) )*(%s))/(%s+%s)",
				     gaus.Data(),Spol7fix.Data(),
				     gaus.Data(),Spol7fix.Data()) );
  flowF->SetLineColor(kGreen-3);
  flowF->SetParameter( 0, yieldF->GetParameter(0)); flowF->SetParLimits(0,+1,-1); flowF->SetParError( 0, yieldF->GetParError(0));
  flowF->SetParameter( 1, yieldF->GetParameter(1)); flowF->SetParLimits(1,+1,-1); flowF->SetParError( 1, yieldF->GetParError(1));
  flowF->SetParameter( 2, yieldF->GetParameter(2)); flowF->SetParLimits(2,+1,-1); flowF->SetParError( 2, yieldF->GetParError(2));
  flowF->SetParameter( 3, 0.1); flowF->SetParLimits(3,-1,+1);
  flowF->SetParameter( 4, 0.1); flowF->SetParLimits(4,-1,+1);
  flowF->SetParameter( 6, 0.0); flowF->SetParLimits(6,+1,-1);  // <=== -1, +1
  cout << endl;
  cout << "First attempt of fitting flowsgn <LR>" << endl;
  flow->Fit(flowF,"WLR","",minX,maxX);
  flow->Fit(flowF,"WLR","",minX,maxX);
  flow->Fit(flowF,"WLR","",minX,maxX);
  cout << "Second attempt of fitting flowsgn <MIR>" << endl;
  flow->Fit(flowF,"EMIR","",minX,maxX);
  TF1 *flowBgrF = new TF1("flowBgrF", Form("([3]*%s+([4]+[5]*(x-[1])+[6]*pow(x-[1],2) )*(%s))/(%s+%s)",
					   gaus.Data(),Spol7fix.Data(),
					   gaus.Data(),Spol7fix.Data()),minX,maxX );
  flowBgrF->SetLineColor(kRed-3);
  flowBgrF->SetParameter( 0, 0);
  flowBgrF->SetParameter( 1, flowF->GetParameter(1));
  flowBgrF->SetParameter( 2, flowF->GetParameter(2));
  flowBgrF->SetParameter( 3, flowF->GetParameter(3));
  flowBgrF->SetParameter( 4, flowF->GetParameter(4));
  flowBgrF->SetParameter( 5, flowF->GetParameter(5));
  flowBgrF->SetParameter( 6, flowF->GetParameter(6));
  flowBgrF->Draw("SAME");

  // bgr
  //new TCanvas();
  //exported_mass->Draw();
  //flow->Fit(v2all);
  Double_t nchi2 = flowF->GetChisquare()/flowF->GetNDF();
  tex->DrawLatexNDC(0.15,0.85, Form("FlowSGN  %.3f",flowF->GetParameter(3)) );
  tex->DrawLatexNDC(0.15,0.80, Form("Chi2/NDF  %.2f",nchi2) );
  mainR->SaveAs( Form("fit/res/%s_%s_PB%02d_DF.png",cut.Data(),method.Data(),ptbin),"png" );

  ofstream fout( Form("fit/res/%s_%s_PB%02d_DF.res",cut.Data(),method.Data(),ptbin)  );
  TString nocortitle = title(1,title.Length()-2);
  TString sptmin = nocortitle(0,nocortitle.Length()/2-2);
  TString sptmax = nocortitle(nocortitle.Length()/2+1,nocortitle.Length()-2);
  float ptmin = sptmin.Atof();
  float ptmax = sptmax.Atof();
  cout << "*" << ptmin << "*" << ptmax << "*" << endl;
  fout << ptbin << " " << ptmin << " " << ptmax << endl;
  for(int i=0; i!=3; ++i) {
    fout << Form( "%f %f", yieldF->GetParameter(i), yieldF->GetParError(i)) << endl;
  }
  fout << Form( "%f %f", flowF->GetParameter(3), flowF->GetParError(3)) << endl;
  fout.close();
  cout << "File  " << Form("fit/res/%s_%s_PB%02d_DF.res",cut.Data(),method.Data(),ptbin) << "  was saved." << endl;

}
