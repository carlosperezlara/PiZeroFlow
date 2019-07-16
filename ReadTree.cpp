#include "Analysis.h"
#include "AT_ReadTree.h"

int main(int argc, char *argv[]){
  if(argc<3) {
    return 1;
  }
  TString run = argv[1];
  TString snev = argv[2];
  int nev = snev.Atoi();

  Analysis *ana = Analysis::Instance();
  ana->InputFileName( Form("trees/%s.root",run.Data()) );
  ana->OutputFileName( Form("Basic/out/out_%s.root",run.Data()) );
  ana->DataSetTag( run );
  ana->NumberOfEventsToAnalyze( nev );

  AT_ReadTree *tsk = new AT_ReadTree();
  ana->AddTask( tsk );

  ana->Run();
}
