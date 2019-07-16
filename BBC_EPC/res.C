void res(int run=454777) {
  TFile *file = new TFile( Form("out/run%d.root",run) );
  TList *listA = (TList*) file->Get("BBCEPC");
  TList *list = (TList*) listA->FindObjects("Resolution");
  TF1 *fit = new TF1( "fit", "[0]+[1]*(x-20)", 10, 30 );
  ofstream fout( Form("tables/BBC_R_%d.dat",run) );
  for(int ord=0; ord!=6; ++ord) {
    for(int i=0; i!=60; ++i) {
      TProfile *QH = (TProfile*) list->FindObject( Form("Res_ORD%d_CB%02d",ord,i) );
      fit->SetParameter(0,0);
      fit->SetParameter(1,0);
      QH->Fit( fit, "RL", "", 10.5, 29.5 );
      fout << fit->GetParameter( 0 ) << " " << fit->GetParError( 0 ) << " ";
      fout << fit->GetParameter( 1 ) << " " << fit->GetParError( 1 ) << endl;
    }
    fout << endl;
  }
  fout.close();
  file->Close();
}
