fTree->SetBranchAddress("Event",&fGLB);
//=
fTree->SetBranchAddress("Q1ex",&pQ1ex);
fTree->SetBranchAddress("Q2ex",&pQ2ex);
fTree->SetBranchAddress("Q3ex",&pQ3ex);
fTree->SetBranchAddress("Q4ex",&pQ4ex);
fTree->SetBranchAddress("Q6ex",&pQ6ex);
fTree->SetBranchAddress("Q1fv",&pQ1fv);
fTree->SetBranchAddress("Q2fv",&pQ2fv);
fTree->SetBranchAddress("Q1bb",&pQ1bb);
fTree->SetBranchAddress("Q2bb",&pQ2bb);
//=
fTree->SetBranchAddress("EMCid",   &pEMCid);
fTree->SetBranchAddress("EMCtwrid",&pEMCtwrid);
fTree->SetBranchAddress("EMCx",    &pEMCx);
fTree->SetBranchAddress("EMCy",    &pEMCy);
fTree->SetBranchAddress("EMCz",    &pEMCz);
fTree->SetBranchAddress("EMCecore",&pEMCecore);
fTree->SetBranchAddress("EMCecent",&pEMCecent);
fTree->SetBranchAddress("EMCchisq",&pEMCchisq);
fTree->SetBranchAddress("EMCtimef",&pEMCtimef);
//=
fTree->SetBranchAddress("TRKqua",  &pTRKqua);
fTree->SetBranchAddress("TRKpt",   &pTRKpt);
fTree->SetBranchAddress("TRKphi",  &pTRKphi);
fTree->SetBranchAddress("TRKpz",   &pTRKpz);
fTree->SetBranchAddress("TRKecore",&pTRKecore);
fTree->SetBranchAddress("TRKetof", &pTRKetof);
fTree->SetBranchAddress("TRKplemc",&pTRKplemc);
fTree->SetBranchAddress("TRKtwrid",&pTRKtwrid);
fTree->SetBranchAddress("TRKchisq",&pTRKchisq);
fTree->SetBranchAddress("TRKdphi", &pTRKdphi);
fTree->SetBranchAddress("TRKdz",   &pTRKdz);
fTree->SetBranchAddress("TRKpc3sdphi",&pTRKpc3sdphi);
fTree->SetBranchAddress("TRKpc3sdz",  &pTRKpc3sdz);
fTree->SetBranchAddress("TRKzed",  &pTRKzed);
fTree->SetBranchAddress("TRKdisp", &pTRKdisp);
fTree->SetBranchAddress("TRKprob", &pTRKprob);
fTree->SetBranchAddress("TRKcid",  &pTRKcid);
//=
fTree->SetBranchAddress("MXSpt",  &pMXSpt);
fTree->SetBranchAddress("MXSpz",  &pMXSpz);
fTree->SetBranchAddress("MXSphi", &pMXSphi);
fTree->SetBranchAddress("MXSflyr",&pMXSflyr);
fTree->SetBranchAddress("MXSsingleD", &pMXSsingleD);
fTree->SetBranchAddress("MXSsingleP", &pMXSsingleP);
fTree->SetBranchAddress("MXSempccent",&pMXSempccent);
fTree->SetBranchAddress("MXSempc3x3", &pMXSempc3x3);

