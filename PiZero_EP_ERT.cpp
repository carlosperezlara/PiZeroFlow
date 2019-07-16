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

  Analysis *ana = Analysis::Instance();
  ana->InputFileName( Form("treesERT/%s.root",run.Data()) );
  ana->OutputFileName( Form("PiZero_EP/outERT060/out_%s.root",run.Data()) );
  ana->DataSetTag( run );
  ana->NumberOfEventsToAnalyze( nev );

  AT_PiZero *tsk = new AT_PiZero();
  unsigned int trigger_BBCLL1narrowcent  = 0x00000008;
  unsigned int trigger_BBCLL1narrow      = 0x00000010;
  unsigned int trigger_FVTXNSBBCScentral = 0x00100000;
  unsigned int trigger_FVTXNSBBCS        = 0x00400000;
  unsigned int trigger_ERT4x4B           = 0x00000040;
  unsigned int trigger_MPC_N_A           = 0x00010000;
  unsigned int trigger_MPC_N_B           = 0x00020000;
  unsigned int msk = trigger_BBCLL1narrowcent | trigger_BBCLL1narrow | trigger_ERT4x4B;
  tsk->TriggerMask( msk );
  tsk->CentralitySelection(0,60);
  ana->AddTask( tsk );

  AT_EP *tsk2 = new AT_EP();
  ana->AddTask( tsk2 );

  ana->Run();

}
