int coef(int run=454777) {
  TFile *file = new TFile( Form("out/run%d.root",run) );
  TProfile2D *QH;
  Double_t coe;
  ofstream fout[3];
  fout[0].open( Form("tables/BBC_A_%d.dat",run) );
  fout[1].open( Form("tables/BBC_B_%d.dat",run) );
  fout[2].open( Form("tables/BBC_C_%d.dat",run) );
  //TList *list = (TList*) file->Get("BBCEPC");
  TList *listA = (TList*) file->Get("BBCEPC");
  TList *list = (TList*) listA->FindObject("Calib");
  //list->ls();
  for(int se=0; se!=3; ++se) {
    for(int ord=0; ord!=6; ++ord) {
      for(int i=0; i!=60; ++i) {
	QH = (TProfile2D*) list->FindObject( Form("PsiC_SE%d_ORD%d_CB%02d",se,ord,i) );
	//yes, they are reversed!
	int nbins = QH->GetXaxis()->GetNbins(); // vtx
	//cout << " VTX BINS " << nbins << endl;
	// ttable of 32rows x 40columns
	for(int in=0; in!=32; ++in) {
	  for(int j=0; j!=nbins; ++j) {
	    //if( QHH->GetBinEntries( j+1 ) < 30 ) coe=0.0;
	    //else
	    coe = QH->GetBinContent( j+1, in+1 );
	    if( TMath::IsNaN( coe ) ) {
	      cout << "ERROR IN " << QH->GetName() << " bins " << j+1 << " " << in+1 << endl;
	    }
	    fout[se] << Form(" %e", coe);
	  }
	  fout[se] << endl;
	}
	fout[se] << endl;
	//==
	QH = (TProfile2D*) list->FindObject( Form("PsiS_SE%d_ORD%d_CB%02d",se,ord,i) );
	//yes, they are reversed!
	int nbins = QH->GetXaxis()->GetNbins(); // vtx
	// ttable of 32rows x 40columns
	for(int in=0; in!=32; ++in) {
	  for(int j=0; j!=nbins; ++j) {
	    //if( QHH->GetBinEntries( j+1 ) < 30 ) coe=0.0;
	    //else
	    coe = QH->GetBinContent( j+1, in+1 );
	    if( TMath::IsNaN( coe ) ) {
	      cout << "ERROR IN " << QH->GetName() << " bins " << j+1 << " " << in+1 << endl;
	    }
	    fout[se] << Form(" %e", coe);
	  }
	  fout[se] << endl;
	}
	fout[se] << endl;
      }
    }
    fout[se].close();
  }
  return 0;
}
