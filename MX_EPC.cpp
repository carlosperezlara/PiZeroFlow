#include "Analysis.h"
#include "AT_ReadTree.h"
#include "AT_MX_EPC.h"

int main(int argc, char *argv[]){
  if(argc<3) {
    return 1;
  }
  TString run = argv[1];
  TString snev = argv[2];
  int nev = snev.Atoi();

  Analysis *ana = Analysis::Instance();
  ana->InputFileName( Form("trees/%s.root",run.Data()) );
  ana->OutputFileName( Form("MX_EPC/out/out_%s.root",run.Data()) );
  ana->DataSetTag( run );
  ana->NumberOfEventsToAnalyze( nev );

  AT_MX_EPC *tsk = new AT_MX_EPC();
  tsk->SkipBBCQCal();
  ana->AddTask( tsk );

  ana->Run();

  //AT_ReadTree *tskchk = new AT_ReadTree();
  //tskchk->CheckEP1();
  //tskchk->CheckEP2();
}
