// This test runs a single root file over a version of babymaker that is tailored
// to run on condor. The major change is that the target file directory is flat
// author: alec aivazis

doTest() {
    // load the cms stuff
    gSystem->Load("libMiniFWLite.so");
    cout << "loaded miniFWLite" << endl;
    gROOT->ProcessLine(".L CORE/libCMS2NtupleMacrosCORE.so");
    cout << "loaded CMS2NtupleMacros" << endl;
    // load the babymaker library
    gROOT->ProcessLine(".L ScanChain.C++");  // compile ScanChain
    cout << "loaded ScanChain" << endl;
    // add the file to a TChain
    TChain *chain = new TChain("Events"); 
    chain->Add("source/*.root");
    // create the looper
    babyMaker* looper = new babyMaker();
    // loop over the file
    looper->ScanChain(chain, "test"); 

}
