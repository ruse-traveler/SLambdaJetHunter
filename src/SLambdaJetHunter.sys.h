// ----------------------------------------------------------------------------
// 'SLambdaJetHunter.sys.h'
// Derek Anderson
// 01.25.2024
//
// A minimal analysis module to find lambda-tagged
// jets in pythia events.
// ----------------------------------------------------------------------------

#ifndef SLAMBDAJETHUNTER_SYS_H
#define SLAMBDAJETHUNTER_SYS_H

// c++ utilities
#include <cassert>

// make common namespaces implicit
using namespace std;



namespace SColdQcdCorrelatorAnalysis {

  // system methods -----------------------------------------------------------

  void SLambdaJetHunter::InitTree() {

    // print debug statement
    if (m_config.isDebugOn) {
      cout << "SLambdaJetHunter::InitTree() Initializing output tree" << endl;
    }

    // set up output tree
    //  - FIXME use utility structs instead when i/o is ready
    m_outTree = new TTree(m_config.outTreeName.data(), "A tree of lambda-tagged jets");
    m_outTree -> Branch("EvtNJets",      &m_evtNJets,         "EvtNJets/I");
    m_outTree -> Branch("EvtNLambdas",   &m_evtNLambdas,      "EvtNLambdas/I");
    m_outTree -> Branch("EvtNChrgPars",  &m_evtNChrgPars,     "EvtNChrgPars/I");
    m_outTree -> Branch("EvtNNeuPars",   &m_evtNNeuPars,      "EvtNNeuPars/I");
    m_outTree -> Branch("EvtSumEPar",    &m_evtSumEPar,       "EvtSumEPar/D");
    m_outTree -> Branch("EvtVtxX",       &m_evtVtxX,          "EvtVtxX/D");
    m_outTree -> Branch("EvtVtxY",       &m_evtVtxY,          "EvtVtxY/D");
    m_outTree -> Branch("EvtVtxZ",       &m_evtVtxZ,          "EvtVtxZ/D");
    m_outTree -> Branch("PartonA_ID",    &m_evtPartID.first,  "PartonA_ID/I");
    m_outTree -> Branch("PartonB_ID",    &m_evtPartID.second, "PartonB_ID/I");
    m_outTree -> Branch("PartonA_Px",    &m_evtPartPx.first,  "PartonA_Px/D");
    m_outTree -> Branch("PartonA_Py",    &m_evtPartPy.first,  "PartonA_Py/D");
    m_outTree -> Branch("PartonA_Pz",    &m_evtPartPz.first,  "PartonA_Pz/D");
    m_outTree -> Branch("PartonA_E",     &m_evtPartE.first,   "PartonA_E/D");
    m_outTree -> Branch("PartonB_Px",    &m_evtPartPx.second, "PartonB_Px/D");
    m_outTree -> Branch("PartonB_Py",    &m_evtPartPy.second, "PartonB_Py/D");
    m_outTree -> Branch("PartonB_Pz",    &m_evtPartPz.second, "PartonB_Pz/D");
    m_outTree -> Branch("PartonB_E",     &m_evtPartE.second,  "PartonB_E/D");
    m_outTree -> Branch("JetHasLambda",  &m_jetHasLambda);
    m_outTree -> Branch("JetNCst",       &m_jetNCst);
    m_outTree -> Branch("JetID",         &m_jetID);
    m_outTree -> Branch("JetE",          &m_jetE);
    m_outTree -> Branch("JetPt",         &m_jetPt);
    m_outTree -> Branch("JetEta",        &m_jetEta);
    m_outTree -> Branch("JetPhi",        &m_jetPhi);
    m_outTree -> Branch("CstID",         &m_cstID);
    m_outTree -> Branch("CstJetID",      &m_cstJetID);
    m_outTree -> Branch("CstEmbedID",    &m_cstEmbedID);
    m_outTree -> Branch("CstZ",          &m_cstZ);
    m_outTree -> Branch("CstDr",         &m_cstDr);
    m_outTree -> Branch("CstEnergy",     &m_cstE);
    m_outTree -> Branch("CstPt",         &m_cstPt);
    m_outTree -> Branch("CstEta",        &m_cstEta);
    m_outTree -> Branch("CstPhi",        &m_cstPhi);
    m_outTree -> Branch("LambdaID",      &m_lambdaID);
    m_outTree -> Branch("LambdaJetID",   &m_lambdaJetID);
    m_outTree -> Branch("LambdaEmbedID", &m_lambdaEmbedID);
    m_outTree -> Branch("LambdaZ",       &m_lambdaZ);
    m_outTree -> Branch("LambdaDr",      &m_lambdaDr);
    m_outTree -> Branch("LambdaEnergy",  &m_lambdaE);
    m_outTree -> Branch("LambdaPt",      &m_lambdaPt);
    m_outTree -> Branch("LambdaEta",     &m_lambdaEta);
    m_outTree -> Branch("LambdaPhi",     &m_lambdaPhi);
    return;

  }  // end 'InitTree()'



  void SLambdaJetHunter::InitOutput() {

    // print debug statement
    if (m_config.isDebugOn) {
      cout << "SLambdaJetHunter::InitOuput() Initializing output file" << endl;
    }

    m_outFile = new TFile(m_config.outFileName.data(), "recreate");
    if (!m_outFile) {
      cerr << "PANIC: couldn't open SLambdaJetHunter output file!" << endl;
      assert(m_outFile);
    }
    return;

  }  // end 'InitOutput()'



  void SLambdaJetHunter::SaveAndCloseOutput() {

    // print debug statement
    if (m_config.isDebugOn) {
      cout << "SLambdaJetHunter::SaveAndCloseOuput() Saving to output file and closing" << endl;
    }

    m_outFile -> cd();
    m_outTree -> Write();
    m_outFile -> Close();
    return;

  }  // end 'SaveAndCloseOutput()'



  void SLambdaJetHunter::ResetOutput() {

    // print debug statement
    if (m_config.isDebugOn) {
      cout << "SLambdaJetHunter::ResetOuput() Resetting output containers" << endl;
    }

    // reset output variables
    m_genEvtInfo.Reset();
    m_jetInfo.clear();
    m_cstInfo.clear();
    m_lambdaInfo.clear();

    // reset lists/associations
    m_vecSubEvts.clear();
    m_mapCstJetAssoc.clear();
    m_mapLambdaJetAssoc.clear();

    // FIXME remove when i/o of utility structs is ready
    m_evtNJets = 999.;
    m_evtNLambdas = 999.;
    m_evtNTaggedJets = 999.;
    m_evtNChrgPars = 999.;
    m_evtNNeuPars  = 999.;;
    m_evtSumEPar = -999.;
    m_evtVtxX = -999.;
    m_evtVtxY = -999.;
    m_evtVtxZ = -999.;
    m_evtPartID = {-999., -999.};
    m_evtPartPx = {-999., -999.};
    m_evtPartPy = {-999., -999.};
    m_evtPartPz = {-999., -999.};
    m_evtPartE = {-999., -999.};
    m_jetHasLambda.clear();
    m_jetNCst.clear();
    m_jetID.clear();
    m_jetE.clear();
    m_jetPt.clear();
    m_jetEta.clear();
    m_jetPhi.clear();
    m_cstID.clear();
    m_cstJetID.clear();
    m_cstEmbedID.clear();
    m_cstZ.clear();
    m_cstDr.clear();
    m_cstE.clear();
    m_cstPt.clear();
    m_cstEta.clear();
    m_cstPhi.clear();
    m_lambdaID.clear();
    m_lambdaJetID.clear();
    m_lambdaEmbedID.clear();
    m_lambdaZ.clear();
    m_lambdaDr.clear();
    m_lambdaE.clear();
    m_lambdaPt.clear();
    m_lambdaEta.clear();
    m_lambdaPhi.clear();
    return;

  }  // end 'ResetOutput()'

}  // end SColdQcdCorrelatorAnalysis namespace

#endif

// end ------------------------------------------------------------------------

