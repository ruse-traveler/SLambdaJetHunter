// ----------------------------------------------------------------------------
// 'LambdaJetHunterOptions.h'
// Derek Anderson
// 01.25.2024
//
// Options for the SLambdaJetHunter module.
// ----------------------------------------------------------------------------

#ifndef LAMBDAJETHUNTEROPTIONS_H
#define LAMBDAJETHUNTEROPTIONS_H 

// analysis utilities
#include "src/SLambdaJetHunterConfig.h"
#include "/sphenix/user/danderson/install/include/slambdajethunter/SLambdaJetHunter.h"



namespace LambdaJetHunterOptions {

  /* TODO change to config + analysis types */

  // track & particle flow parameters
  const bool   runTracking(false);
  const bool   doTruthTableReco(false);
  const double nSigma(1.5);

  // jet tree general parameters
  const bool isEmbed(true);
  const bool doDebug(false);
  const bool doVtxCut(false);

  // jet tree jet parameters
  const double jetRes  = 0.4;
  const auto   jetAlgo = SLambdaJetHunter::Algo::AntiKt;
  const auto   jetReco = SLambdaJetHunter::Recomb::Pt;

  // event acceptance
  const pair<double, double> vzEvtRange = {-10., 10.};
  const pair<double, double> vrEvtRange = {0.0,  0.418};

  // particle acceptance
  const pair<double, double> ptParRange  = {0.,   9999.};
  const pair<double, double> etaParRange = {-1.1, 1.1};

  // module configuration
  SLambdaJetHunterConfig Config;

}  // end LambdaJetHunterOptions namespace

#endif

// end ------------------------------------------------------------------------
