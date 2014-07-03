#ifndef JETSMEARINGTOOLS_H
#define JETSMEARINGTOOLS_H

#include "CMS2.h"
#include <vector>
#include <utility>

// function to smear jet energy to account for differences in data-MC JER
class JetSmearer;
JetSmearer* makeJetSmearer(const char* ptFileName="$CMSSW_BASE/src/CondFormats/JetMETObjects/data/Spring10_PtResolution_AK5PF.txt",
                           const char* phiFileName="$CMSSW_BASE/src/CondFormats/JetMETObjects/data/Spring10_PhiResolution_AK5PF.txt",
                           const char* resFileName="$CMSSW_BASE/src/RecoMET/METProducers/python/METSigParams_cfi.py");
JetSmearer* makeJetSmearer(std::vector<std::string> &vector_of_file_names);
LorentzVector smearJet(const LorentzVector& p4, JetSmearer* jetSmearer, bool recoOnly = false);
double getJetResolution(const LorentzVector& p4, JetSmearer* jetSmearer);
std::vector<double> getJetResolutions(const std::vector<LorentzVector>& vp4s, JetSmearer* jetSmearer);
std::vector<LorentzVector> smearJets(const std::vector<LorentzVector>& vp4s, JetSmearer* jetSmearer, bool recoOnly = false);

//-----------------------------------------------------
// function to smear the MET for differences in jet
// resolution between data and MC
//-----------------------------------------------------
std::pair<double, double> smearMETforJER(std::pair<double, double> in_met, JetSmearer *jetSmearer, const std::vector<LorentzVector>& vp4s, bool recoOnly = false);

#endif
