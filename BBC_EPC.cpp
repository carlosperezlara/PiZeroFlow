#include "Analysis.h"
#include "AT_ReadTree.h"
#include "AT_BBC_EPC.h"

int main(int argc, char *argv[]){
  if(argc<3) {
    return 1;
  }
  TString run = argv[1];
  TString snev = argv[2];
  int nev = snev.Atoi();

  Analysis *ana = Analysis::Instance();
  ana->InputFileName( Form("trees/%s.root",run.Data()) );
  ana->OutputFileName( Form("BBC_EPC/out/out_%s.root",run.Data()) );
  ana->DataSetTag( run );
  ana->NumberOfEventsToAnalyze( nev );

  AT_BBC_EPC *tsk = new AT_BBC_EPC();
  tsk->SkipBBCQCal();
  tsk->SkipTracks();
  ana->AddTask( tsk );

  ana->Run();
}
