// ----------------------------------------------------------------------------
// 'LambdaJetHunterOptions.h'
// Derek Anderson
// 01.25.2024
//
// Options for the SLambdaJetHunter module.
// ----------------------------------------------------------------------------

#ifndef LAMBDAJETHUNTEROPTIONS_H
#define LAMBDAJETHUNTEROPTIONS_H 

// c++ utilities
#include <string>
#include <utility>
// analysis utilities
#include "/sphenix/user/danderson/install/include/scorrelatorutilities/EvtTools.h"
#include "/sphenix/user/danderson/install/include/scorrelatorutilities/GenTools.h"
#include "/sphenix/user/danderson/install/include/scorrelatorutilities/JetTools.h"
#include "/sphenix/user/danderson/install/include/slambdajethunter/SLambdaJetHunter.h"
#include "/sphenix/user/danderson/install/include/slambdajethunter/SLambdaJetHunterConfig.h"

// make common namespacs implicit
using namespace std;
using namespace SColdQcdCorrelatorAnalysis;
using namespace SColdQcdCorrelatorAnalysis::SCorrelatorUtilities;



namespace LambdaJetHunterOptions {

  // acceptance cuts ----------------------------------------------------------

  // event acceptance
  const pair<float, float> vzEvtRange = {-10., 10.};
  const pair<float, float> vrEvtRange = {0.0,  0.418};

  // particle acceptance
  const pair<float, float> ptParRange  = {0.,   100.};
  const pair<float, float> etaParRange = {-1.1, 1.1};

  // jet acceptance
  const pair<float, float> ptJetRange  = {0.1, 100.};
  const pair<float, float> etaJetRange = {-0.7, 0.7};



  // bundle acceptances into pairs --------------------------------------------

  pair<ParInfo, ParInfo> GetParAccept() {

    pair<ParInfo, ParInfo> parAccept;
    parAccept.first.pt   = ptParRange.first;
    parAccept.first.eta  = etaParRange.first;
    parAccept.second.pt  = ptParRange.second;
    parAccept.second.eta = etaParRange.second;
    return parAccept;

  }  // end 'GetParAccept()'



  pair<JetInfo, JetInfo> GetJetAccept() {

    pair<JetInfo, JetInfo> cfg_jetAccept;
    cfg_jetAccept.first.pt   = ptJetRange.first;
    cfg_jetAccept.first.eta  = etaJetRange.first;
    cfg_jetAccept.second.pt  = ptJetRange.second;
    cfg_jetAccept.second.eta = etaJetRange.second;

  }  // end 'GetJetAccept()'



  // set up configuration -----------------------------------------------------

  SLambdaJetHunterConfig GetConfig(const int verbosity, const string outFile) {

    SLambdaJetHunterConfig cfg {
      .verbosity   = verbosity,
      .isDebugOn   = true,
      .isEmbed     = false,
      .moduleName  = "SLambdaJetHunter",
      .outTreeName = "LambdaJetTree",
      .outFileName = outFile,
      .isCharged   = true,
      .rJet        = 0.4,
      .jetAlgo     = "antikt_algorithm",
      .jetRecomb   = "pt_scheme",
      .vzAccept    = vzEvtRange,
      .vrAccept    = vrEvtRange,
      .parAccept   = GetParAccept(),
      .jetAccept   = GetJetAccept()
    };
    return cfg;

  }  // end 'GetConfig(int, string&)'

}  // end LambdaJetHunterOptions namespace

#endif

// end ------------------------------------------------------------------------
