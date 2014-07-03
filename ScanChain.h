#ifndef ScanChain_h
#define ScanChain_h

// C++ includes
#include <string>
#include <vector>

// ROOT includes
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "Math/VectorUtil.h"
#include "Math/Vector4D.h"

#include "CORE/ssSelections.h"
#include "CORE/muonSelections.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

class babyMaker {

 public:

  babyMaker() {};
  ~babyMaker() {
    delete BabyFile_;
    delete BabyTree_;
  };

  void ScanChain(TChain* chain, std::string baby_name = "testSample", int numEvents = -1, float customScale = -1, bool isData = false);

  void MakeBabyNtuple(const char *);
  void InitBabyNtuple();
  void FillBabyNtuple();
  void CloseBabyNtuple();

 private:

  TFile *BabyFile_;
  TTree *BabyTree_;

  //baby ntuple variables
  float met;
  float metPhi;

  //std::vector<LorentzVector> jets_p4;
  std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > jets_p4;

  std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > generated_p4;

  std::vector<int> generated_id;
  std::vector<int> generated_mother_id;


  int type;

  LorentzVector ll_p4;
  LorentzVector lt_p4;
  LorentzVector total_p4;

  int ll_id;
  int lt_id;
  int ll_charge;
  int lt_charge;
  int ll_index;
  int lt_index;
  

  float scale_1fb;
  std::vector<float> btagDiscriminant;
  
  int eventNumber;
  int runNumber;
  int lumiBlock;

  string file;

};

#endif

void babyMaker::MakeBabyNtuple(const char *BabyFilename){

  //
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();
  BabyFile_ = new TFile(Form("%s", BabyFilename), "RECREATE");
  BabyFile_->cd();
  BabyTree_ = new TTree("tree", "A Baby Ntuple");

  BabyTree_->Branch("met", &met );
  BabyTree_->Branch("metPhi", &metPhi );
  BabyTree_->Branch("jets_p4", &jets_p4 );
  BabyTree_->Branch("type", &type);

  BabyTree_->Branch("ll_p4", &ll_p4);
  BabyTree_->Branch("lt_p4", &lt_p4);
  BabyTree_->Branch("total_p4", &total_p4);

  BabyTree_->Branch("ll_id", &ll_id);
  BabyTree_->Branch("lt_id", &lt_id);
  BabyTree_->Branch("ll_charge" , &ll_charge);
  BabyTree_->Branch("lt_charge", &lt_charge);
  BabyTree_->Branch("ll_index", &ll_index);
  BabyTree_->Branch("lt_index", &lt_index);
  
  BabyTree_->Branch("eventNumber", &eventNumber);
  BabyTree_->Branch("runNumber", &runNumber);
  BabyTree_->Branch("lumiBlock", &lumiBlock);
  
  BabyTree_->Branch("scale_1fb", &scale_1fb);
  BabyTree_->Branch("btagDiscriminant", &btagDiscriminant);

  BabyTree_->Branch("generated_p4", &generated_p4);
  BabyTree_->Branch("generated_id", &generated_id);
  BabyTree_->Branch("generated_mother_id", &generated_mother_id);
  


  BabyTree_->Branch("file", &file);

  return;
}

void babyMaker::InitBabyNtuple () {

  met = -999.0;
  metPhi = 0;

  type = -1;

  ll_id = -1;
  lt_id = -1;
  ll_charge = -999;
  lt_charge = -999;
  ll_index = -1;
  lt_index = -1;

  scale_1fb = 1;

  eventNumber = -1;
  runNumber = -1;
  lumiBlock = 1;

  return;
}

void babyMaker::FillBabyNtuple(){
  BabyTree_->Fill();
  return;
}

void babyMaker::CloseBabyNtuple(){
  BabyFile_->cd();
  BabyTree_->Write();
  BabyFile_->Close();
  return;
}

