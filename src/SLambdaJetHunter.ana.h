// ----------------------------------------------------------------------------
// 'SLambdaJetHunter.ana.h'
// Derek Anderson
// 01.25.2024
//
// A minimal analysis module to find lambda-tagged
// jets in pythia events.
// ----------------------------------------------------------------------------

#ifndef SLAMBDAJETHUNTER_ANA_H
#define SLAMBDAJETHUNTER_ANA_H

// fastjet libraries
#include <fastjet/PseudoJet.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/FunctionOfPseudoJet.hh>
// analysis utilities
#include "/sphenix/user/danderson/install/include/scorrelatorutilities/EvtTools.h"
#include "/sphenix/user/danderson/install/include/scorrelatorutilities/JetTools.h"
#include "/sphenix/user/danderson/install/include/scorrelatorutilities/CstTools.h"
#include "/sphenix/user/danderson/install/include/scorrelatorutilities/GenTools.h"

// make common namespaces implicit
using namespace std;
using namespace fastjet;
using namespace findNode;
using namespace SColdQcdCorrelatorAnalysis::SCorrelatorUtilities;



namespace SColdQcdCorrelatorAnalysis {

  // analysis methods ---------------------------------------------------------

  void SLambdaJetHunter::MakeJets(PHCompositeNode* topNode) {

    // print debug statement
    if (m_config.isDebugOn) {
      cout << "SLambdaJetHunter::MakeJets() Making jets" << endl;
    }

    // for fastjet
    vector<PseudoJet> vecPseudoJets;

    // loop over subevents
    vector<int> vecSubEvents = GrabSubevents(topNode);
    for (const int subEvt : vecSubEvents) {

      // loop over particles
      HepMC::GenEvent* genEvt = GetGenEvent(topNode, subEvt);
      for (
        HepMC::GenEvent::particle_const_iterator hepPar = genEvt -> particles_begin();
        hepPar != genEvt -> particles_end();
        ++hepPar
      ) {

        // grab particle info
        ParInfo particle(*hepPar, subEvt);

        // check if particle is good & final state
        const bool isParGood    = IsGoodParticle(particle);
        const bool isFinalState = IsFinalState(particle.status);
        if (!isParGood || !isFinalState) continue;

        // make pseudojet
        PseudoJet pseudo(particle.px, particle.py, particle.pz, particle.ene);
        pseudo.set_user_index(particle.barcode);
        vecPseudoJets.push_back(pseudo);

      }  // end particle loop
    }  // end subevent loop

    // set jet definition
    JetDefinition definition(
      m_mapJetAlgo[m_config.jetAlgo],
      m_config.rJet,
      m_mapRecomb[m_config.jetRecomb],
      fastjet::Best
    );

    // run clustering
    ClusterSequence* sequence = new ClusterSequence(vecPseudoJets, definition);

    // grab jets and return
    m_jets = sequence -> inclusive_jets();
    return;

  }  // end 'MakeJets(PHCompositeNode*)'



  bool SLambdaJetHunter::IsGoodParticle(ParInfo& particle) {

    // print debug statement
    if (m_config.isDebugOn) {
      cout << "SLambdaJetHunter::IsGoodParticle(ParInfo) checking if particle is good" << endl;
    }

    // run checks and return
    const bool isInAccept = IsInAcceptance(particle, m_config.parAccept.first, m_config.parAccept.second);
    return isInAccept;

  }  // end 'IsGoodParticle(ParInfo&)'



  bool SLambdaJetHunter::IsLambda() {

    /* TODO lambda checking goes here */
    return true;

  }  // end 'IsLambda()'



  ParInfo SLambdaJetHunter::FindLambda(const int barcode) {

    /* TODO lambda hunting goes here */
    ParInfo lambda;
    return lambda;

  }  // end 'FindLambda(int)'

}  // end SColdQcdCorrelatorAnalysis namespace

#endif

// end ------------------------------------------------------------------------
