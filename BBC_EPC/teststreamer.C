int teststreamer(int run=454777) {
  int n=0;
  double tmp;
  ifstream fin( Form("tables/BBC_A_%d.dat",run) );
  for( ; ; ++n ) {
    fin >> tmp;
    if(!fin.good()) break;
    cout << tmp << " ";
  }
  cout << "I could read " << n << " floats." << endl;
  return 0;
}
