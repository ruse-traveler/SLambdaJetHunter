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

  void SLambdaJetHunter::GrabEventInfo(PHCompositeNode* topNode) {

    // print debug statement
    if (m_config.isDebugOn) {
      cout << "SLambdaJetHunter::GrabEventInfo() Grabbing event info" << endl;
    }

    m_vecSubEvts = GrabSubevents(topNode);
    m_genEvtInfo.SetInfo(topNode, m_config.isEmbed, m_vecSubEvts);
    return;

  }  // end 'GrabEventInfo(PHCompositeNode*)'



  void SLambdaJetHunter::MakeJets(PHCompositeNode* topNode) {

    // print debug statement
    if (m_config.isDebugOn) {
      cout << "SLambdaJetHunter::MakeJets() Making jets" << endl;
    }

    // for fastjet
    vector<PseudoJet> vecPseudoJets;

    // loop over subevents
    for (const int subEvt : m_vecSubEvts) {

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



  void SLambdaJetHunter::HuntLambdas(PHCompositeNode* topNode) {

    // print debug statement
    if (m_config.isDebugOn) {
      cout << "SLambdaJetHunter::HuntLambdas() hunting down lambda-taggged jets" << endl;
    }

    // loop over fastjet output
    for (size_t iJet = 0; iJet < m_jets.size(); iJet++) {

      // grab jet info
      m_jetInfo.emplace_back( m_jets[iJet] );
      m_jetInfo.back().jetID = iJet;

      // loop over constituents
      m_csts.clear();
      for (size_t iCst = 0; iCst < m_jets[iJet].constituents().size(); iCst++) {

        // grab cst info
        m_csts.emplace_back( m_jets[iJet].constituents()[iCst] );
        m_csts.back().cstID = iCst;

        /* TODO analysis steps
         *   (1) retrieve particle based on user_index
         *   (2) loop through parents until a lambda is found
         *   (3) if lambda hasn't already been added, shove into list
         *   (4) continue on
         */ 

      }  // end cst loop

      // add constituents to output
      m_cstInfo.push_back( m_csts );

    }  // end jet loop
    return;

  }  // end 'HuntLambdas(PHCompositeNode*)'



  void SLambdaJetHunter::FillOutputTree() {

    // print debug statement
    if (m_config.isDebugOn) {
      cout << "SLambdaJetHunter::FillOutputTree() collecting output information and filling tree" << endl;
    }

    // collect event info
    m_evtNJets       = m_jetInfo.size();
    m_evtNLambdas    = 0;  // TODO fill in when ready
    m_evtNTaggedJets = 0;  // TODO fill in when ready
    m_evtNChrgPars   = m_genEvtInfo.nChrgPar;
    m_evtNNeuPars    = m_genEvtInfo.nNeuPar;
    m_evtSumEPar     = m_genEvtInfo.eSumChrg + m_genEvtInfo.eSumNeu;
    m_evtVtxX        = m_genEvtInfo.partons.first.vx;
    m_evtVtxY        = m_genEvtInfo.partons.first.vy;
    m_evtVtxZ        = m_genEvtInfo.partons.first.vz;

    // collect parton info
    m_evtPartID = make_pair( m_genEvtInfo.partons.first.pid, m_genEvtInfo.partons.second.pid );
    m_evtPartPx = make_pair( m_genEvtInfo.partons.first.px,  m_genEvtInfo.partons.second.px  );
    m_evtPartPy = make_pair( m_genEvtInfo.partons.first.py,  m_genEvtInfo.partons.second.py  );
    m_evtPartPz = make_pair( m_genEvtInfo.partons.first.pz,  m_genEvtInfo.partons.second.pz  );
    m_evtPartE  = make_pair( m_genEvtInfo.partons.first.ene, m_genEvtInfo.partons.second.ene );

    // collect jet information
    for (JetInfo& jet : m_jetInfo) {
      m_jetNCst.push_back( jet.nCsts );
      m_jetID.push_back (  jet.jetID );
      m_jetE.push_back(    jet.ene   );
      m_jetPt.push_back(   jet.pt    );
      m_jetEta.push_back(  jet.eta   );
      m_jetPhi.push_back(  jet.phi   );
    }

    m_outTree -> Fill();
    return;

  }  // end 'GrabOutputInfo()'



  bool SLambdaJetHunter::IsGoodParticle(ParInfo& particle) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 5)) {
      cout << "SLambdaJetHunter::IsGoodParticle(ParInfo) checking if particle is good" << endl;
    }

    // run checks and return
    const bool isInAccept = IsInAcceptance(particle, m_config.parAccept.first, m_config.parAccept.second);
    return isInAccept;

  }  // end 'IsGoodParticle(ParInfo&)'



  bool SLambdaJetHunter::IsLambda(ParInfo& particle) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 5)) {
      cout << "SLambdaJetHunter::IsLambda() checking if particle is a lambda" << endl;
    }

    const bool isLambda = (particle.pid == m_const.pidLambda);
    return isLambda;

  }  // end 'IsLambda(ParInfo&)'



  ParInfo SLambdaJetHunter::FindLambda(const int barcode) {

    /* TODO lambda hunting goes here */
    ParInfo lambda;
    return lambda;

  }  // end 'FindLambda(int)'

}  // end SColdQcdCorrelatorAnalysis namespace

#endif

// end ------------------------------------------------------------------------
