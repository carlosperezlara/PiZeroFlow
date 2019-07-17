#include <iostream>
#include <fstream>
#include <vector>
#include <TString.h>
#include <TTree.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include "Analysis.h"
#include "AT_EventChecker.h"
#include "EmcIndexer.h"
#include "EmcIndexer.C"
#include "PbGlIndexer.C"
#include "PbScIndexer.C"

AT_EventChecker::AT_EventChecker() : AT_ReadTree() {
  for(int i=0; i!=2; ++i) {
    hTrackMultiplicity[i] = NULL;
    hBBCMultiplicity[i] = NULL;
    hBBCTmean[i] = NULL;
    hBBCTrmsS[i] = NULL;
    hBBCTrmsN[i] = NULL;
    hBBCmeanT[i] = NULL;
    hBBCrmsT[i] = NULL;
    hBBCsgn[i] = NULL;
    hBBCsgnD[i] = NULL;
    hBBCvtx[i] = NULL;
    hFrac[i] = NULL;
  }

}
AT_EventChecker::~AT_EventChecker() {
}
void AT_EventChecker::MyInit() {
  hEvents->GetXaxis()->SetBinLabel(10,"AT_EventChecker init");
  hEvents->GetXaxis()->SetBinLabel(11,"AT_EventChecker no_errors");
  hEvents->GetXaxis()->SetBinLabel(12,"AT_EventChecker vtx_time");
  hEvents->GetXaxis()->SetBinLabel(13,"AT_EventChecker rms");
  hEvents->GetXaxis()->SetBinLabel(14,"AT_EventChecker energy");
  for(int i=0; i!=2; ++i) {
    hTrackMultiplicity[i] = new TH1F( Form("hTrackMultiplicity%d",i),
				      Form("hTrack multiplicity %d",i),
				      100, -0.5, 99.5 );
    hBBCMultiplicity[i] = new TH1F( Form("hBBCMultiplicity%d",i),
				    Form("hBBC multiplicity %d",i),
				    100, 0, 200 );
    hBBCTmean[i] = new TH1F( Form("hBBCTimeMean%d",i),
			     Form("hBBBC mean Time %d; TimeS  -  TimeN  (a.u.)",i),
                             1000, -5, +5 );
    hBBCTrms[i] = new TH1F( Form("hBBCTimeRMS%d",i),
			     Form("hBBBC rms Time %d; TRMS_S  +  TRMS_N  (a.u.)",i),
                             1000, 0, 3 );
    hBBCTrmsS[i] = new TH1F( Form("hBBCTimeRMS_S%d",i),
			     Form("hBBBC rms Time South %d;rmsS  (a.u.)",i),
                             1000, +0.0, +0.5 );
    hBBCTrmsN[i] = new TH1F( Form("hBBCTimeRMS_N%d",i),
			     Form("hBBBC rms Time North %d;rmsN  (a.u.)",i),
                             1000, +0.0, +0.5 );
    hBBCTnrmsS[i] = new TH1F( Form("hBBCTimeNRMS_S%d",i),
			      Form("hBBBC rms Time South %d;norm rmsS  (a.u.)",i),
			      1000, +0.0, +5 );
    hBBCTnrmsN[i] = new TH1F( Form("hBBCTimeNRMS_N%d",i),
			      Form("hBBBC rms Time North %d;norm rmsN  (a.u.)",i),
			      1000, +0.0, +5 );
    hBBCmeanT[i] = new TH2F( Form("hBBCTimeMean2D%d",i),
			     Form("hBBBC mean Time %d;time mean S  (a.u.);time mean N  (a.u)",i),
                             100, -10, +30, 100, -10, +30 );
    hBBCrmsT[i] = new TH2F( Form("hBBCTimeRMS2D%d",i),
			    Form("hBBC rms Time %d;time rms S  (a.u.);time rms N  (a.u.)",i),
                            1000, 0.0, 3.0,
                            1000, 0.0, 3.0 );
    hBBCsgn[i] = new TH2F( Form("hBBCSgn2D%d",i),
			   Form("hBBC signals %d;signal S  (a.u.);signal N  (a.u.)",i),
			   100, 0.0, 170,
			   100, 0.0, 120 );
    hBBCsgnD[i] = new TH1F( Form("hBBCSgn%d",i),
			    Form("hBBC signals %d;signal S-N  (a.u.)",i),
                            100, -30., +150. );
    hBBCvtx[i] = new TH1F( Form("hBBCvtx%d",i),
			   Form("hBBC Vertex %d;vtxZ  (cm)",i),
			   100, -30., 30. );
    hFrac[i] = new TH1F( Form("hFrac%d",i),
			 Form("hFrac %d",i),
			 100, 0., 1. );
  }


}
void AT_EventChecker::MyFinish() {
  /*
  TF1 *fitBBC[2], *fitTRK[2];
  for(int i=0; i!=2; ++i) {
    fitBBC[i] = new TF1( Form("fitBBC%d",i),"[0]*TMath::Gaus(x,[1],[2],1)");
    fitTRK[i] = new TF1( Form("fitTRK%d",i),"[0]*TMath::Gaus(x,[1],[2],1)");
    fitTRK[i]->SetParameter(1,12); fitTRK[i]->SetParLimits(1,6,16);
    fitTRK[i]->SetParameter(2,5);  fitTRK[i]->SetParLimits(2,3,8);
    fitBBC[i]->SetParameter(1,80); fitBBC[i]->SetParLimits(1,65,95);
    fitBBC[i]->SetParameter(2,15); fitBBC[i]->SetParLimits(2,10,20);
    hBBCMultiplicity[i]->Fit( fitBBC[i], "R", "", 60, 100 );
    hTrackMultiplicity[i]->Fit( fitTRK[i], "R", "", 5, 20 );
  }
  Analysis *ana = Analysis::Instance();
  int run = ana->RunNumber();
  int seg = ana->SegmentNumber();
  std::ofstream fout( Form("EventChecker/dat/%d_%02d.dat",run,seg) );
  fout << Form("%d_%02d",run,seg);
  if(0) {
    fout << " " << fitBBC[0]->GetParameter(1);
    fout << " " << fitBBC[0]->GetParameter(2);
    fout << " " << fitBBC[1]->GetParameter(1);
    fout << " " << fitBBC[1]->GetParameter(2);
    fout << " " << fitTRK[0]->GetParameter(1);
    fout << " " << fitTRK[0]->GetParameter(2);
    fout << " " << fitTRK[1]->GetParameter(1);
    fout << " " << fitTRK[1]->GetParameter(2);
  } else {
    fout << " " << hBBCMultiplicity[0]->GetMean();
    fout << " " << hBBCMultiplicity[0]->GetMeanError();
    fout << " " << hBBCMultiplicity[1]->GetMean();
    fout << " " << hBBCMultiplicity[1]->GetMeanError();
    fout << " " << hTrackMultiplicity[0]->GetMean();
    fout << " " << hTrackMultiplicity[0]->GetMeanError();
    fout << " " << hTrackMultiplicity[1]->GetMean();
    fout << " " << hTrackMultiplicity[1]->GetMeanError();
  }
  fout << std::endl;
  fout.close();
  */  
  for(int i=0; i!=2; ++i) {
    hBBCMultiplicity[i]->Write();
    hTrackMultiplicity[i]->Write();
    hBBCTmean[i]->Write();
    hBBCTrms[i]->Write();
    hBBCTrmsS[i]->Write();
    hBBCTrmsN[i]->Write();
    hBBCTnrmsS[i]->Write();
    hBBCTnrmsN[i]->Write();
    hBBCmeanT[i]->Write();
    hBBCrmsT[i]->Write();
    hBBCsgn[i]->Write();
    hBBCsgnD[i]->Write();
    hBBCvtx[i]->Write();
    hFrac[i]->Write();
  }

}
void AT_EventChecker::MyExec() {
  //std::cout << "AAAA" << std::endl;
  float vtx = fGLB.vtxZ;
  float cen = fGLB.cent;
  int bvtx = BinVertex( vtx );
  int bcen = BinCentrality( cen );
  float frac = fGLB.frac;
  double rmsS = fGLB2.bbcsTrms;
  double rmsN = fGLB2.bbcnTrms;
  float meanS = fGLB2.bbcsTmean;
  float meanN = fGLB2.bbcnTmean;
  float sgnS = fGLB.bbcs;
  float sgnN = fGLB2.bbcn;
  //if(sgnS>93||sgnS<71) return;
  //if(sgnN>50||sgnN<5) return;
  if(sgnS>kBBCADCMax_S||sgnS<kBBCADCMin_S) return;
  if(sgnN>kBBCADCMax_N||sgnN<kBBCADCMin_N) return;
  hEvents->Fill(9); // init
  if( TMath::IsNaN(meanS) || TMath::IsNaN(meanN) ) return;
  if( TMath::IsNaN(rmsS) || TMath::IsNaN(rmsN) ) {
    std::cout << "FUNNYYYY!!!" << std::endl;
    return;
  }
  hEvents->Fill(10); // noerrors

  hBBCvtx[0]->Fill( vtx );
  hBBCTmean[0]->Fill( meanS - meanN );
  hBBCTrms[0]->Fill( rmsS/fTime[1] + rmsN/fTime[2] );
  hBBCTrmsS[0]->Fill( rmsS );
  hBBCTrmsN[0]->Fill( rmsN );
  hBBCTnrmsS[0]->Fill( rmsS/fTime[1] );
  hBBCTnrmsN[0]->Fill( rmsN/fTime[2] );
  hBBCmeanT[0]->Fill( meanS, meanN );
  hBBCrmsT[0]->Fill( rmsS, rmsN );
  hBBCsgn[0]->Fill( sgnS, sgnN );
  hBBCsgnD[0]->Fill( sgnS - sgnN );
  hTrackMultiplicity[0]->Fill( fGLB2.alltrks);
  hBBCMultiplicity[0]->Fill( sgnS + sgnN );
  hFrac[0]->Fill( frac );

  if( TMath::Abs(meanS-meanN+fTime[0]) > 5*0.60 ) return;
  hEvents->Fill(11); // vtx time
  if( rmsS/fTime[1] + rmsN/fTime[2] > 1 ) return;
  hEvents->Fill(12); // rms time

  hBBCvtx[1]->Fill( vtx );
  hBBCTmean[1]->Fill( meanS - meanN );
  hBBCTrms[1]->Fill( rmsS/fTime[1] + rmsN/fTime[2] );
  hBBCTrmsS[1]->Fill( rmsS );
  hBBCTrmsN[1]->Fill( rmsN );
  hBBCTnrmsS[1]->Fill( rmsS/fTime[1] );
  hBBCTnrmsN[1]->Fill( rmsN/fTime[2] );
  hBBCmeanT[1]->Fill( meanS, meanN );
  hBBCrmsT[1]->Fill( rmsS, rmsN );
  hBBCsgn[1]->Fill( sgnS, sgnN );
  hBBCsgnD[1]->Fill( sgnS - sgnN );
  hTrackMultiplicity[1]->Fill( fGLB2.alltrks);
  hBBCMultiplicity[1]->Fill( sgnS + sgnN );
  hFrac[1]->Fill( frac );

  hEvents->Fill(13); // energy
}
