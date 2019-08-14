#include <TFile.h>
#include <TList.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TString.h>
#include "QC_EventPlane.h"
#include "QC_ScalarProduct.h"
#include "QC_Cumulants.h"

int main(int nv, char** vars) {
  TFile *fin;
  TString set="NOM";
  if(nv==2) {
    set = vars[1];
  }
  fin = new TFile( Form("PiZero_EP/allfiles/all%s.root",set.Data()) );
  fin->ls();
  TList *listEP = (TList*) fin->Get("AT_EP");
  //return 0;
  TList *list;
  {
    QC_EventPlane *epc = new QC_EventPlane();
    list = (TList*) fin->Get("EP_Ord1_Psi1");
    list->ls();
    epc->SetList( list );
    epc->Results();
    TList *mass = new TList();
    mass->SetName("yields");
    mass->SetOwner();
    list->Add( mass );
    //decoding yields
    TH2D *map = (TH2D*) list->FindObject( "binmap_EP_Ord1_Psi1" );
    TH1D *binPT[100];
    TH1D *mixPT[100];
    for(int pt=0; pt!=map->GetXaxis()->GetNbins(); ++pt) {
      Double_t minPt = map->GetXaxis()->GetBinLowEdge( pt+1 );
      Double_t maxPt = map->GetXaxis()->GetBinLowEdge( pt+2 );
      binPT[pt] = map->ProjectionY( Form("Yield_PB%d",pt), pt+1, pt+1 );
      binPT[pt]->SetTitle( Form("[%.1f - %.1f ]",minPt,maxPt) );
      mixPT[pt] = (TH1D*) listEP->FindObject( Form("hMass2_PB%d",pt) );
      mass->Add( binPT[pt] );
      mass->Add( mixPT[pt] );
    }
    epc->Save( Form("PiZero_EP/allfiles/all%s_EP_Ord1.root",set.Data()) );
  }

  if(0){
    QC_ScalarProduct *spc = new QC_ScalarProduct();
    list = (TList*) fin->Get("SP_Ord1_Psi1");
    list->ls();
    spc->SetList( list );
    spc->Results();
    TList *mass = new TList();
    mass->SetName("yields");
    mass->SetOwner();
    list->Add( mass );
    //decoding yields
    TH2D *map = (TH2D*) list->FindObject( "binmap_SP_Ord1_Psi1" );
    TH1D *binPT[100];
    TH1D *mixPT[100];
    for(int pt=0; pt!=map->GetXaxis()->GetNbins(); ++pt) {
      Double_t minPt = map->GetXaxis()->GetBinLowEdge( pt+1 );
      Double_t maxPt = map->GetXaxis()->GetBinLowEdge( pt+2 );
      binPT[pt] = map->ProjectionY( Form("Yield_PB%d",pt), pt+1, pt+1 );
      binPT[pt]->SetTitle( Form("[%.1f - %.1f ]",minPt,maxPt) );
      mixPT[pt] = (TH1D*) listEP->FindObject( Form("hMass2_PB%d",pt) );
      mass->Add( binPT[pt] );
      mass->Add( mixPT[pt] );
    }
    spc->Save( Form("PiZero_EP/allfiles/all%s_SP_Ord1.root",set.Data()) );
  }

  if(0){
    QC_Cumulants *qcc = new QC_Cumulants();
    list = (TList*) fin->Get("QC_Ord1_Psi1");
    list->ls();
    qcc->SetList( list );
    qcc->Results();
    TList *mass = new TList();
    mass->SetName("yields");
    mass->SetOwner();
    list->Add( mass );
    //decoding yields
    TH2D *map = (TH2D*) list->FindObject( "binmap_QC_Ord1_Psi1" );
    TH1D *binPT[100];
    TH1D *mixPT[100];
    for(int pt=0; pt!=map->GetXaxis()->GetNbins(); ++pt) {
      Double_t minPt = map->GetXaxis()->GetBinLowEdge( pt+1 );
      Double_t maxPt = map->GetXaxis()->GetBinLowEdge( pt+2 );
      binPT[pt] = map->ProjectionY( Form("Yield_PB%d",pt), pt+1, pt+1 );
      binPT[pt]->SetTitle( Form("[%.1f - %.1f ]",minPt,maxPt) );
      mixPT[pt] = (TH1D*) listEP->FindObject( Form("hMass2_PB%d",pt) );
      mass->Add( binPT[pt] );
      mass->Add( mixPT[pt] );
    }
    qcc->Save( Form("PiZero_EP/allfiles/all%s_QC_Ord1.root",set.Data()) );
  }

  fin->Close();


  /*  
  TFile *file;
  TList *listA;

  file = new TFile("aaa.root");
  listA = (TList*) file->Get("Results");
  TH1D *aa = (TH1D*) listA->FindObject("Results_PT8");
  TProfile *bb = (TProfile*) listA->FindObject("_pc_PT8");

  file = new TFile("bbb.root");
  listA = (TList*) file->Get("Results");
  TH1D *cc = (TH1D*) listA->FindObject("Results_PT8");

  aa->SetLineColor(kRed-3);
  bb->SetLineColor(kBlue-3);
  cc->SetLineColor(kGreen-3);
  TCanvas *main = new TCanvas();
  aa->Draw("H");
  bb->Draw("same");
  cc->Draw("same");

  main->SaveAs("aaa.jpg","JPG");
  */
  return 0;
}
