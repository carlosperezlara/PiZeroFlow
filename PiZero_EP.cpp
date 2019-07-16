#include <iostream>
#include "Analysis.h"
#include "AT_PiZero.h"
#include "AT_EP.h"

int main(int argc, char *argv[]){
  if(argc<3) {
    return 1;
  }
  TString run = argv[1];
  TString snev = argv[2];
  int nev = snev.Atoi();
  TString spar3 = argv[3];

  AT_PiZero *tsk = new AT_PiZero();
  tsk->CentralitySelection(0,5);
  tsk->SkipTracks();
  TString ssys = "";
  if(spar3.Contains("FD0")) {
    tsk->SetDist(7.0); ssys="FD0";
  } else if(spar3.Contains("D0")) {
    tsk->SetDist(7.5); ssys="D0";
  } else if(spar3.Contains("FD1")) {
    tsk->SetDist(9.0); ssys="FD1";
  } else if(spar3.Contains("D1")) {
    tsk->SetDist(8.5); ssys="D1";
  } else if(spar3.Contains("FA0")) {
    tsk->SetAlpha(0.65); ssys="FA0";
  } else if(spar3.Contains("A0")) {
    tsk->SetAlpha(0.75); ssys="A0";
  } else if(spar3.Contains("FA1")) {
    tsk->SetAlpha(0.90); ssys="FA1";
  } else if(spar3.Contains("A1")) {
    tsk->SetAlpha(0.85); ssys="A1";
  } else if(spar3.Contains("FT0")) {
    tsk->SetTime(4.0); ssys="FT0";
  } else if(spar3.Contains("T0")) {
    tsk->SetTime(4.5); ssys="T0";
  } else if(spar3.Contains("FT1")) {
    tsk->SetTime(6.0); ssys="FT1";
  } else if(spar3.Contains("T1")) {
    tsk->SetTime(5.5); ssys="T1";
  }

  Analysis *ana = Analysis::Instance();
  ana->InputFileName( Form("trees/%s.root",run.Data()) );
  ana->OutputFileName( Form("PiZero_EP/out%s/out_%s.root",ssys.Data(),run.Data()) );
  ana->DataSetTag( run );
  ana->NumberOfEventsToAnalyze( nev );
  ana->AddTask( tsk );

  AT_EP *tsk2 = new AT_EP();
  ana->AddTask( tsk2 );

  ana->Run();

}
