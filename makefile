all:
	rootcint -f Dict.cpp -c qcQ.h LinkDef.h
	#g++ -o Run_ReadTree ReadTree.cpp AT_ReadTree.cxx Analysis.cxx qcQ.cxx Dict.cpp `root-config --cflags --glibs`
	#g++ -o Run_BBC_EPC BBC_EPC.cpp AT_BBC_EPC.cxx AT_ReadTree.cxx Analysis.cxx qcQ.cxx Dict.cpp `root-config --cflags --glibs`
	#g++ -o Run_MX_EPC MX_EPC.cpp AT_MX_EPC.cxx AT_ReadTree.cxx Analysis.cxx qcQ.cxx Dict.cpp `root-config --cflags --glibs`
	#g++ -o Run_PiZero PiZero.cpp AT_PiZero.cxx AT_ReadTree.cxx Analysis.cxx qcQ.cxx Dict.cpp `root-config --cflags --glibs`
	g++ -o Run_PiZero_EP PiZero_EP.cpp AT_EP.cxx AT_PiZero.cxx AT_ReadTree.cxx Analysis.cxx qcQ.cxx QC_EventPlane.cxx QC_ScalarProduct.cxx QC_Cumulants.cxx QCorrelation.cxx Dict.cpp `root-config --cflags --glibs`
	g++ -o Run_EventChecker EventChecker.cpp AT_EventChecker.cxx AT_ReadTree.cxx Analysis.cxx qcQ.cxx Dict.cpp `root-config --cflags --glibs`
	g++ -o redofinish redofinish.cpp qcQ.cxx QC_EventPlane.cxx QC_ScalarProduct.cxx QC_Cumulants.cxx QCorrelation.cxx Dict.cpp `root-config --cflags --glibs`
	#g++ -o toymodel toymodel.cpp qcQ.cxx QC_EventPlane.cxx QC_ScalarProduct.cxx QC_Cumulants.cxx QCorrelation.cxx Dict.cpp `root-config --cflags --glibs`
	rm Dict.*
