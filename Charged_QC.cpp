#include "Analysis.h"
#include "AT_Charged.h"
#include "AT_QC.h"

int main(int argc, char *argv[]){
  if(argc<3) {
    return 1;
  }
  TString run = argv[1];
  TString snev = argv[2];
  int nev = snev.Atoi();

  Analysis *ana = Analysis::Instance();
  ana->InputFileName( Form("trees/%s.root",run.Data()) );
  ana->OutputFileName( Form("Charged/out/out_%s.root",run.Data()) );
  ana->DataSetTag( run );
  ana->NumberOfEventsToAnalyze( nev );

  AT_Charged *sel = new AT_Charged();
  ana->AddTask( sel );
  AT_QC *tsk = new AT_QC();
  ana->AddTask( tsk );

  ana->Run();
}
