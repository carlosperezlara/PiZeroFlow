int qcent(int run=454777) {
  TFile *file = new TFile( Form("out/run%d.root",run) );
  TList *list2 = (TList*) file->Get("BBCEPC");
  TList *list = (TList*) list2->FindObject("Calib");
  //list->ls();
  TH2D *QH;
  TH1D *h;
  ofstream fout;
  char xy[2] = {'x','y'};

  //========= RAW
  fout.open( Form("tables/BBC_%d.dat",run) );
  for(int ord=0; ord!=9; ++ord) {
    for(int ix=0; ix!=2; ++ix) {
      for(int se=0; se!=3; ++se) {
	//printing 60rows x 40 columns table
	for(int i=0; i!=60; ++i) {
	  QH = NULL;
	  QH = (TH2D*) list->FindObject( Form("Q%cVtx_ORD%d_SE%d_CB%02d_ST0",xy[ix],ord,se,i) );
	  int nbins = QH->GetXaxis()->GetNbins(); // NVtx
	  for(int j=0; j!=nbins; ++j) {
	    h = QH->ProjectionY( Form("%s_P%d",QH->GetName(),j), j+1, j+1 );
	    Double_t mean = h->GetMean();
	    if(h->GetEntries()<100) mean = 0;
	    fout << Form(" %e", mean);
	  }
	  fout << endl;
	}
	fout << endl;
      }
    }
  }
  fout.close();

  //========= PRIMED
  fout.open( Form("tables/BBC_P_%d.dat",run) );
  for(int ord=0; ord!=9; ++ord) {
    for(int ix=0; ix!=2; ++ix) {
      for(int se=0; se!=3; ++se) {
	//printing 60rows x 40 columns table
	for(int i=0; i!=60; ++i) {
	  QH = NULL;
	  QH = (TH2D*) list->FindObject( Form("Q%cVtx_ORD%d_SE%d_CB%02d_ST1",xy[ix],ord,se,i) );
	  int nbins = QH->GetXaxis()->GetNbins(); // NVtx
	  for(int j=0; j!=nbins; ++j) {
	    h = QH->ProjectionY( Form("%s_P%d",QH->GetName(),j), j+1, j+1 );
	    Double_t mean = h->GetMean();
	    if(h->GetEntries()<100) mean = 0;
	    fout << Form(" %e", mean);
	  }
	  fout << endl;
	}
	fout << endl;
      }
    }
  }
  fout.close();

  //========= PRI-PRIMED
  fout.open( Form("tables/BBC_Q_%d.dat",run) );
  for(int ord=0; ord!=9; ++ord) {
    for(int ix=0; ix!=2; ++ix) {
      for(int se=0; se!=3; ++se) {
	//printing 60rows x 40 columns table
	for(int i=0; i!=60; ++i) {
	  QH = NULL;
	  QH = (TH2D*) list->FindObject( Form("Q%cVtx_ORD%d_SE%d_CB%02d_ST2",xy[ix],ord,se,i) );
	  int nbins = QH->GetXaxis()->GetNbins(); // NVtx
	  for(int j=0; j!=nbins; ++j) {
	    h = QH->ProjectionY( Form("%s_P%d",QH->GetName(),j), j+1, j+1 );
	    Double_t mean = h->GetMean();
	    if(h->GetEntries()<100) mean = 0;
	    fout << Form(" %e", mean);
	  }
	  fout << endl;
	}
	fout << endl;
      }
    }
  }
  fout.close();


  return 0;
}
