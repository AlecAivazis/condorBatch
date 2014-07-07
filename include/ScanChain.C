// C++
#include <iostream>
#include <vector>
#include <set>
#include <fstream>

// ROOT
#include "TDirectory.h"
#include "TTreeCache.h"
#include "Math/VectorUtil.h"
#include "Math/LorentzVector.h"

// CMS2
#include "CORE/CMS2.h"
#include "CORE/electronSelections.h"
#include "CORE/muonSelections.h"
#include "CORE/ssSelections.h"
#include "CORE/trackSelections.h"
#include "CORE/eventSelections.h"
#include "CORE/susySelections.h"


#include "goodrun.cc"
#include "Include.C"

// header
#include "ScanChain.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

using namespace std;
using namespace tas;
using namespace ROOT::Math::VectorUtil;

bool isValidPair(int hypIndex, int jetIndex){

    // require the generated pair to be either mu- b or mu+ bbar
    // 
 
    // if(genps_id().at(hypIndex) * genps_id().at(jetIndex) ) return false;
    
    // muons have id = 13
    // b's have id = 5
    // i need bbar and muon (or opposite)
    return genps_id().at(hypIndex) * genps_id().at(jetIndex) == -65;

}

void babyMaker::ScanChain(TChain* chain, std::string baby_name, int numEvents, float customScale, bool isData){

  if (numEvents != -1 ){
      cout << "Processing the first " << numEvents << " event(s)" << endl;
  }

  MakeBabyNtuple( Form("%s.root", baby_name.c_str()) );

  // File Loop
  int nDuplicates = 0;
  int nEvents = chain->GetEntries();
  unsigned int nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;

  double _hypPt20Counter = 0;
  double _osCounter = 0;
  double _typeCounter = 0;
  double _etaCounter = 0;
  double _muonIdCounter = 0;
  double _muonIsoCounter = 0;
  double _electronIdCounter = 0;
  double _electronIsoCounter = 0;
  double _eventsCounter = 0;
  double _muonCounter = 0;

  double hypCounter = 0;
  double eventHypCounter = 0;

  bool cont = true;

  while ( (currentFile = (TFile*)fileIter.Next()) && cont ) {

    // Get File Content
    TFile f( currentFile->GetTitle() );
    TTree *tree = (TTree*)f.Get("Events");
    TTreeCache::SetLearnEntries(10);
    tree->SetCacheSize(128*1024*1024);
    cms2.Init(tree);
    
    // set good run list
    if (evt_isRealData()) 
        set_goodrun_file("/home/users/jgran/analysis/sswh/fakes/json/final_19p49fb_cms2.txt");

    // event Loop
    unsigned int nEventsTree = tree->GetEntriesFast();

    for(unsigned int event = 0; event < nEventsTree; ++event) {

      // get event content
      tree->LoadTree(event);
      cms2.GetEntry(event);

      // count the number of events
      ++nEventsTotal;

      // if its less than the max
      if (numEvents != -1 && nEventsTotal > numEvents) {
          cont = false;
          break;
      }

       // select good vagina
      if(evt_isRealData() && !goodrun(evt_run(), evt_lumiBlock())) continue;

      if(evt_isRealData()){
        DorkyEventIdentifier id = {evt_run(), evt_event(), evt_lumiBlock()};
        if (is_duplicate(id)){
          continue;
        }
      }
      

      int index = -1;
      float _maxPt = 0.0;

      // count number of hypotheses
      if (hyp_p4().size() != 0) eventHypCounter++;

      hypCounter = hypCounter + hyp_p4().size();
      
      // apply cuts to hypotheses
      for (unsigned int i = 0; i< hyp_p4().size(); i++){

          if (hyp_ll_p4().at(i).pt() < 20) continue;
          if (hyp_lt_p4().at(i).pt() < 20) continue;
          _hypPt20Counter++;
          
          // oppositely charged leptons - NOTE: NOT WHAT ALEX DOES (uses id voodoo hoodoo)
          if (hyp_ll_charge().at(i)*hyp_lt_charge().at(i) > 0) continue;
          _osCounter++;

          // eta < 2.4
          if (fabs(hyp_ll_p4().at(i).eta()) > 2.4) continue;
          if (fabs(hyp_lt_p4().at(i).eta()) > 2.4) continue;
          _etaCounter++;
            
          // invariant mass > 20
          if ((hyp_ll_p4().at(i)+hyp_lt_p4().at(i)).M() < 20) continue;

          // count number of muons before ID/iso
          if (abs(hyp_ll_id().at(i)) == 13 && abs(hyp_lt_id().at(i)) == 13) {
	      _muonCounter++;
          }

          /* This iso is too tight. Saving in case they change their mind....
             if(!isGoodLepton(hyp_lt_id().at(i), hyp_lt_index().at(i))) continue;
             if(!isGoodLepton(hyp_ll_id().at(i), hyp_ll_index().at(i))) continue;
          */

          // electron id and iso
          if (abs(hyp_ll_id().at(i)) == 11 && 
              !passElectronSelection_ZMet2012_v3_Iso(hyp_ll_index().at(i))) continue;
          if (abs(hyp_lt_id().at(i)) == 11 && 
              !passElectronSelection_ZMet2012_v3_Iso(hyp_lt_index().at(i))) continue;
          _electronIdCounter++;

          // muon id and iso
          if (abs(hyp_ll_id().at(i)) == 13 && 
              !muonId(hyp_ll_index().at(i), ZMet2012_v1)) continue;
          if (abs(hyp_lt_id().at(i)) == 13 && 
              !muonId(hyp_lt_index().at(i), ZMet2012_v1)) continue;
          _muonIdCounter++;
          /*          
          // ll muon isolation
          if (abs(hyp_ll_id().at(i)) == 13){
              double chiso_ll = mus_isoR03_pf_ChargedHadronPt().at(hyp_ll_index().at(i));
              double nhiso_ll = mus_isoR03_pf_NeutralHadronEt().at(hyp_ll_index().at(i));
              double emiso_ll = mus_isoR03_pf_PhotonEt().at(hyp_ll_index().at(i));
              double dbeta_ll = mus_isoR03_pf_PUPt().at(hyp_ll_index().at(i));
              double iso_ll = (chiso_ll + max(0.0, nhiso_ll + emiso_ll - 0.5 * dbeta_ll))
              / hyp_ll_p4().at(i).pt();

              if (iso_ll > 0.2) continue; 
          }
          
          // lt muon isolation
          if (abs(hyp_lt_id().at(i)) == 13){
              double chiso_lt = mus_isoR03_pf_ChargedHadronPt().at(hyp_lt_index().at(i));
              double nhiso_lt = mus_isoR03_pf_NeutralHadronEt().at(hyp_lt_index().at(i));
              double emiso_lt = mus_isoR03_pf_PhotonEt().at(hyp_lt_index().at(i));
              double dbeta_lt = mus_isoR03_pf_PUPt().at(hyp_lt_index().at(i));
              double iso_lt = (chiso_lt + max(0.0, nhiso_lt + emiso_lt - 0.5 * dbeta_lt)) 
              / hyp_lt_p4().at(i).pt();

              if (iso_lt > 0.2) continue;
          }
          
          // electron isolation 
          if (abs(hyp_ll_id().at(i)) == 11){
              if(!samesign2011::isGoodLepton(hyp_ll_id().at(i), hyp_ll_index().at(i))) continue;
          };

          if (abs(hyp_lt_id().at(i)) == 11){
              if(!samesign2011::isGoodLepton(hyp_lt_id().at(i), hyp_lt_index().at(i))) continue;
          };

          // count the number of electrons and muons that pass
          // if it's not a mu/mu then it's an e/mu
          if (abs(hyp_lt_id().at(i)) == 13 && abs(hyp_ll_id().at(i)) == 13){
              _muonIsoCounter++;
          } else {
              _electronIsoCounter++;
          }
          */
          // select the highest pt hypothesis
          float sumPt = hyp_lt_p4().at(i).pt() + hyp_ll_p4().at(i).pt();
         
          if (sumPt > _maxPt){
              _maxPt = sumPt;
              index = i;
          }
          
      }

      if (index == -1) continue;
      _eventsCounter++;

      // create the ntuple
      InitBabyNtuple();

      // set the branch values
      met = evt_pfmet_type1cor();
      metPhi = evt_pfmetPhi();

      if (! isData){
          generated_p4 = genps_p4();
          generated_id = genps_id();
          generated_mother_id = genps_id_mother();
      }

      //  correct the jet pt at baby level
      std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > jets;
      for (unsigned int i = 0; i<pfjets_p4().size(); i++) {
          jets.push_back(pfjets_p4().at(i) * pfjets_corL1FastL2L3().at(i));
      }

      jets_p4 = jets;

      type = hyp_type().at(index);

      ll_charge = hyp_ll_charge().at(index);
      ll_id = hyp_ll_id().at(index);
      ll_index = hyp_ll_index().at(index);
      ll_p4 = hyp_ll_p4().at(index);

      lt_charge = hyp_lt_charge().at(index);
      lt_id = hyp_lt_id().at(index);
      lt_index = hyp_lt_index().at(index);
      lt_p4 = hyp_lt_p4().at(index);
      
      total_p4 = hyp_p4().at(index);

      btagDiscriminant = pfjets_combinedSecondaryVertexBJetTag();

      file = Form("%s", currentFile->GetTitle());

      eventNumber = evt_event();
      runNumber = evt_run();
      lumiBlock = evt_lumiBlock();

      // grab the appropriate scale
      if (customScale == -1) {
          scale_1fb = evt_scale1fb();
      } else {
          scale_1fb = customScale; 
      }

      // if they didn't specify a total number of events to run over
      if (numEvents == -1)
          // use the total number in the chain for the scale
          numEvents = nEventsChain;

      // scale the scale for minis
      scale_1fb *= (nEventsChain/numEvents);

      // fill it
      FillBabyNtuple();

      // track progress
      CMS2::progress( nEventsTotal, numEvents );

    } // end of loop over events in file

    delete tree;
    f.Close();

  } // end of loop over files

  if (nEventsChain != nEventsTotal){
      cout << "WARNING: The number of events added is not equal to the total number of events!" << endl;
  }

  cout << nDuplicates << " duplicate events were skipped." << endl;

  CloseBabyNtuple();

  ofstream stream;
  stream.open("cutflow.txt", ios::app);

  stream << baby_name << ": " << endl;
  stream << Form("Source (Events): %.0d", numEvents) << endl;
  stream << Form("Events with Hypothesis: %.0f (%.2f)", eventHypCounter, eventHypCounter/numEvents * 100) << endl;
  stream << Form("Source (Hypotheses): %.0f", hypCounter) << endl;
  stream << Form("Hypothesis Pt > 20: %.0f (%.2f, %.2f)",_hypPt20Counter, _hypPt20Counter/hypCounter * 100, 1-(_hypPt20Counter/hypCounter)) << endl;
  stream << Form("Oppositely Charged: %.0f (%.2f, %.2f)",_osCounter, _osCounter/hypCounter * 100, 1-(_osCounter/_hypPt20Counter)) << endl;
  stream << Form("Eta < 2.4: %.0f (%.2f, %.2f)",_etaCounter, _etaCounter/hypCounter * 100, 1-(_etaCounter/_typeCounter)) << endl;
  stream << Form("e/mu passing ID: %.0f (%.2f, %.2f)", _electronIdCounter, _electronIdCounter/hypCounter * 100, 1-(_electronIdCounter/_etaCounter)) << endl;
  stream << Form("# of muons : %.0f (%.2f)",_muonCounter, _muonCounter/hypCounter * 100) << endl;
  stream << Form("# muons passing ID: %.0f (%.2f, %.2f)",_muonIdCounter, _muonIdCounter/hypCounter * 100, 1-(_muonIdCounter/_etaCounter)) << endl;
  stream << "-" << endl;
  stream << Form("# of hypothesis passing ID/ISO: %.0f (%.2f)",(_muonIsoCounter + _electronIsoCounter), (_muonIsoCounter + _electronIsoCounter)/hypCounter * 100) << endl;
  stream << Form("# of events passing ID/ISO: %.0f (%.2f)",(_eventsCounter), (_eventsCounter)/eventHypCounter * 100) << endl;
//stream << Form("Muon Events: %.0f (%.2f)", muonCounter, muonCounter/numEvents * 100) << endl;
  stream << "--------------------------------" << endl;

  /*
    if (showControlRegions){v

    stream << Form("Control Region 1: %.1f", (CR1counter/numEvents) * 100 ) <<  " " << CR1counter/9 << endl; 
    stream << Form("Control Region 2: %.1f", (CR2counter/numEvents) * 100 ) << " " << CR2counter/4 << endl; 
    stream << Form("Control Region 3: %.1f", (CR3counter/numEvents) * 100 ) << " " << CR3counter/1 << endl; 
    stream << "--------------------------------" << endl;
    }
  */
  stream.close();
  

  return;
}

