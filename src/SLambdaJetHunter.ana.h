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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

// c++ utilities
#include <map>
#include <cmath>
#include <limits>
// root libraries
#include <TMath.h>
// fastjet libraries
#include <fastjet/PseudoJet.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/FunctionOfPseudoJet.hh>
// hepmc libraries
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#include <HepMC/GenParticle.h>
// PHG4 libraries
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Particlev2.h>
#include <g4main/PHG4TruthInfoContainer.h>
// analysis utilities
#include "/sphenix/user/danderson/install/include/scorrelatorutilities/EvtTools.h"
#include "/sphenix/user/danderson/install/include/scorrelatorutilities/JetTools.h"
#include "/sphenix/user/danderson/install/include/scorrelatorutilities/CstTools.h"
#include "/sphenix/user/danderson/install/include/scorrelatorutilities/GenTools.h"

#pragma GCC diagnostic pop

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
      cout << "SLambdaJetHunter::GrabEventInfo(PHCompositeNode*) Grabbing event info" << endl;
    }

    // FIXME turn on subevents after figuring out why
    //   IsGoodParticle and IsFinalState aren't removing
    //   non FS particles
    //m_vecSubEvts = GrabSubevents(topNode);
    m_vecSubEvts = {SubEvt::NotEmbedSignal};
    m_genEvtInfo.SetInfo(topNode, m_config.isEmbed, m_vecSubEvts);
    return;

  }  // end 'GrabEventInfo(PHCompositeNode*)'



  void SLambdaJetHunter::FindLambdas(PHCompositeNode* topNode) {

    // print debug statement
    if (m_config.isDebugOn) {
      cout << "SLambdaJetHunter::FindLambdas(PHCompositeNode*) Finding all lambdas in event" << endl;
    }

    // loop over subevents
    for (const int subEvt : m_vecSubEvts) {

      // loop over particles
      HepMC::GenEvent* genEvt = GetGenEvent(topNode, subEvt);
      for (
        HepMC::GenEvent::particle_const_iterator hepPar = genEvt -> particles_begin();
        hepPar != genEvt -> particles_end();
        ++hepPar
      ) {

        // check if lambda and if found yet
        const bool isLambda = IsLambda( (*hepPar) -> pdg_id() );
        const bool isNew    = IsNewLambda( (*hepPar) -> barcode() );
        if (!isLambda || !isNew) continue;

        // grab info
        ParInfo lambda( (*hepPar), subEvt );

        // check if good
        const bool isGood = IsGoodLambda(lambda);
        if (!isGood) continue;

        // if a new good lambda, add to output vector and lambda-jet map
        m_lambdaInfo.push_back(lambda);
        m_mapLambdaJetAssoc.insert(
          make_pair( (*hepPar) -> barcode(), -1 )
        );

      }  // end particle loop
    }  // end subevent loop
    return;

  }  // end 'FindLambdas(PHCompositeNode*)'



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
    m_vecFastJets = sequence -> inclusive_jets();
    return;

  }  // end 'MakeJets(PHCompositeNode*)'



  void SLambdaJetHunter::CollectJetOutput(PHCompositeNode* topNode) {

    // print debug statement
    if (m_config.isDebugOn) {
      cout << "SLambdaJetHunter::CollectOutput() hunting down lambda-taggged jets" << endl;
    }

    // reserve space for fastjet output
    m_jetInfo.resize( m_vecFastJets.size() );
    m_cstInfo.resize( m_vecFastJets.size() );

    // loop over fastjet output
    for (size_t iJet = 0; iJet < m_vecFastJets.size(); iJet++) {

      // grab jet info
      m_jetInfo[iJet].SetInfo( m_vecFastJets[iJet] );
      m_jetInfo[iJet].jetID = iJet;

      // loop over constituents
      m_cstInfo[iJet].resize( m_vecFastJets[iJet].constituents().size() );
      for (size_t iCst = 0; iCst < m_vecFastJets[iJet].constituents().size(); iCst++) {

        // grab cst info
        m_cstInfo[iJet][iCst].SetInfo( m_vecFastJets[iJet].constituents()[iCst] );

        // get particle PDG
        HepMC::GenParticle* hepCst = GetHepMCGenParticleFromBarcode(m_cstInfo[iJet][iCst].cstID, topNode);
        m_cstInfo[iJet][iCst].pid  = hepCst -> pdg_id();

        // run calculations
        const float pCst  = hypot( m_cstInfo[iJet][iCst].px, m_cstInfo[iJet][iCst].py, m_cstInfo[iJet][iCst].pz );
        const float pJet  = hypot( m_jetInfo[iJet].px, m_jetInfo[iJet].py, m_jetInfo[iJet].pz );
        const float dfCst = m_cstInfo[iJet][iCst].phi - m_jetInfo[iJet].phi;
        const float dhCst = m_cstInfo[iJet][iCst].eta - m_jetInfo[iJet].eta;

        // grab remaining cst info
        m_cstInfo[iJet][iCst].jetID = m_jetInfo[iJet].jetID;
        m_cstInfo[iJet][iCst].z     = pCst / pJet;
        m_cstInfo[iJet][iCst].dr    = hypot( dfCst, dhCst );

      }  // end cst loop
    }  // end jet loop
    return;

  }  // end 'CollectJetOutput(PHCompositeNode*)'



  void SLambdaJetHunter::AssociateLambdasToJets(PHCompositeNode* topNode) {

    // print debug statement
    if (m_config.isDebugOn) {
      cout << "SLambdaJetHunter::AssociateLambdasToJets() associating found lambdas to jets" << endl;
    }

    // loop over lambdas
    for (ParInfo& lambda : m_lambdaInfo) {

      // if needed, first try associating by decay chain
      int iAssocJet = -1;
      if (!m_config.useOnlyDistHunt) {
        iAssocJet = HuntLambdasByDecayChain(lambda, topNode);
      }

      // otherwise, associate by distance
      if (iAssocJet == -1) {
        iAssocJet = HuntLambdasByDistance(lambda);
      }

      // assign associated jet
      m_mapLambdaJetAssoc[lambda.barcode] = iAssocJet;

    }  // end lambda loop
    return;

  }   // end 'AssociateLambdasToJets(PHCompositeNode*)'



  void SLambdaJetHunter::FillOutputTree() {

    // print debug statement
    if (m_config.isDebugOn) {
      cout << "SLambdaJetHunter::FillOutputTree() collecting output information and filling tree" << endl;
    }

    // collect event info
    //   - FIXME remove when i/o of utility structs is ready
    m_evtNJets       = m_jetInfo.size();
    m_evtNLambdas    = m_mapLambdaJetAssoc.size();
    m_evtNTaggedJets = 1.;   // TODO put in
    m_evtNChrgPars   = m_genEvtInfo.nChrgPar;
    m_evtNNeuPars    = m_genEvtInfo.nNeuPar;
    m_evtSumEPar     = m_genEvtInfo.eSumChrg + m_genEvtInfo.eSumNeu;
    m_evtVtxX        = m_genEvtInfo.partons.first.vx;
    m_evtVtxY        = m_genEvtInfo.partons.first.vy;
    m_evtVtxZ        = m_genEvtInfo.partons.first.vz;

    // collect parton info
    //   - FIXME remove when i/o of utility structs is ready
    m_evtPartID = make_pair( m_genEvtInfo.partons.first.pid, m_genEvtInfo.partons.second.pid );
    m_evtPartPx = make_pair( m_genEvtInfo.partons.first.px,  m_genEvtInfo.partons.second.px  );
    m_evtPartPy = make_pair( m_genEvtInfo.partons.first.py,  m_genEvtInfo.partons.second.py  );
    m_evtPartPz = make_pair( m_genEvtInfo.partons.first.pz,  m_genEvtInfo.partons.second.pz  );
    m_evtPartE  = make_pair( m_genEvtInfo.partons.first.ene, m_genEvtInfo.partons.second.ene );

    // collect lambda information
    //   - FIXME remove when i/o of utility structs is ready
    for (ParInfo& lambda : m_lambdaInfo) {

      // collect general information
      m_lambdaID.push_back( lambda.barcode );
      m_lambdaEmbedID.push_back( lambda.embedID );
      m_lambdaE.push_back( lambda.ene );
      m_lambdaPt.push_back( lambda.pt );
      m_lambdaEta.push_back( lambda.eta );
      m_lambdaPhi.push_back( lambda.phi );

      // collect associate jet information
      m_lambdaJetID.push_back( m_mapLambdaJetAssoc[lambda.barcode] );
      if (m_mapLambdaJetAssoc[lambda.barcode] > -1) {
        m_lambdaZ.push_back( GetLambdaAssocZ(lambda) );
        m_lambdaDr.push_back( GetLambdaAssocDr(lambda) );
      } else {
        m_lambdaZ.push_back( -1. );
        m_lambdaDr.push_back( -1. );
      }
    }

    // collect jet information
    //   - FIXME remove when i/o of utility structs is ready
    for (JetInfo& jet : m_jetInfo) {
      m_jetHasLambda.push_back( false );  // TODO fill in when ready
      m_jetNCst.push_back( jet.nCsts );
      m_jetID.push_back( jet.jetID );
      m_jetE.push_back( jet.ene );
      m_jetPt.push_back( jet.pt );
      m_jetEta.push_back( jet.eta );
      m_jetPhi.push_back( jet.phi );
    }  // end jet loop

    // collect cst information
    //   - FIXME remove when i/o of utility structs is ready
    m_cstID.resize( m_cstInfo.size() );
    m_cstPID.resize( m_cstInfo.size() );
    m_cstJetID.resize( m_cstInfo.size() );
    m_cstEmbedID.resize( m_cstInfo.size() );
    m_cstZ.resize( m_cstInfo.size() );
    m_cstDr.resize( m_cstInfo.size() );
    m_cstE.resize( m_cstInfo.size() );
    m_cstPt.resize( m_cstInfo.size() );
    m_cstEta.resize( m_cstInfo.size() );
    m_cstPhi.resize( m_cstInfo.size() );
    for (size_t iJet = 0; iJet < m_cstInfo.size(); iJet++) {
      m_cstID[iJet].resize( m_cstInfo[iJet].size() );
      m_cstPID[iJet].resize( m_cstInfo[iJet].size() );
      m_cstJetID[iJet].resize( m_cstInfo[iJet].size() );
      m_cstEmbedID[iJet].resize( m_cstInfo[iJet].size() );
      m_cstZ[iJet].resize( m_cstInfo[iJet].size() );
      m_cstDr[iJet].resize( m_cstInfo[iJet].size() );
      m_cstE[iJet].resize( m_cstInfo[iJet].size() );
      m_cstPt[iJet].resize( m_cstInfo[iJet].size() );
      m_cstEta[iJet].resize( m_cstInfo[iJet].size() );
      m_cstPhi[iJet].resize( m_cstInfo[iJet].size() );
      for (size_t iCst = 0; iCst < m_cstInfo[iJet].size(); iCst++) {
        m_cstID[iJet][iCst]      = m_cstInfo[iJet][iCst].cstID;
        m_cstPID[iJet][iCst]     = m_cstInfo[iJet][iCst].pid;
        m_cstJetID[iJet][iCst]   = m_cstInfo[iJet][iCst].jetID;
        m_cstEmbedID[iJet][iCst] = m_cstInfo[iJet][iCst].embedID;
        m_cstZ[iJet][iCst]       = m_cstInfo[iJet][iCst].z;
        m_cstDr[iJet][iCst]      = m_cstInfo[iJet][iCst].dr;
        m_cstE[iJet][iCst]       = m_cstInfo[iJet][iCst].ene;
        m_cstPt[iJet][iCst]      = m_cstInfo[iJet][iCst].pt;
        m_cstEta[iJet][iCst]     = m_cstInfo[iJet][iCst].eta;
        m_cstPhi[iJet][iCst]     = m_cstInfo[iJet][iCst].phi;
      }
    }  // end cst loop

    m_outTree -> Fill();
    return;

  }  // end 'FillOutputTree()'
 


  int SLambdaJetHunter::HuntLambdasByDecayChain(ParInfo& lambda, PHCompositeNode* topNode) {

    // print debug statement
    if (m_config.isDebugOn) {
      cout << "SLambdaJetHunter::HuntLambdasByDecayChain(ParInfo&, PHCompositeNode*) hunting for lambdas by inspecting decay chains" << endl;
    }

    // grab PHG4 and HepMC information
    HepMC::GenParticle* hepLambda = GetHepMCGenParticleFromBarcode(lambda.barcode, topNode);
    HepMC::GenVertex*   hepEndVtx = hepLambda -> end_vertex();

    // loop over constituents and check which cst.s the lambda decayed into
    int  iAssocJet     = -1;
    bool foundAssocJet = false;
    for (size_t iJet = 0; iJet < m_cstInfo.size(); iJet++) {
      for (size_t iCst = 0; iCst < m_cstInfo[iJet].size(); iCst++) {

        // check if constituent is in decay chain,
        bool isInDecayChain = false;
        if (hepEndVtx) {
          isInDecayChain = IsInHepMCDecayChain(m_cstInfo[iJet][iCst].cstID, hepEndVtx);
        } else {
          isInDecayChain = IsInPHG4DecayChain(m_cstInfo[iJet][iCst].cstID, lambda.barcode, topNode);
        }

        // if so, add lambda to lists
        if (isInDecayChain) {
          foundAssocJet = true;
          iAssocJet     = iJet;
          break;
        }
      } // end cst loop

      // if found association, exit loop
      if (foundAssocJet)  break;

    }  // end jet loop
    return iAssocJet;

  }  // end 'HuntLambdasByDecayChain(PHCompositeNode*)'



  int SLambdaJetHunter::HuntLambdasByDistance(ParInfo& lambda) {

    // print debug statement
    if (m_config.isDebugOn) {
      cout << "SLambdaJetHunter::HuntLambdasByDistance(ParInfo&) hunting for lambds by distance" << endl;
    }

    // loop over jets
    int    iAssocJet   = -1;
    double drAssocBest = numeric_limits<double>::max();
    for (size_t iJet = 0; iJet < m_jetInfo.size(); iJet++) {

      // calculate delta-phi
      double dfAssoc = lambda.phi - m_jetInfo[iJet].phi;
      if (dfAssoc < 0.)             dfAssoc += TMath::TwoPi();
      if (dfAssoc > TMath::TwoPi()) dfAssoc -= TMath::TwoPi();

      // calculate dr
      const double dhAssoc = lambda.eta - m_jetInfo[iJet].eta;
      const double drAssoc = hypot(dfAssoc, dhAssoc);

      // check if lambda is within rJet
      const bool isDistGood = (drAssoc < m_config.rJet);
      const bool isDistBest = (drAssoc < drAssocBest);
      if (isDistGood && isDistBest) {
        drAssocBest = drAssoc;
        iAssocJet   = iJet;
      }
    }
    return iAssocJet;

  }  // end 'HuntLambdasByDistance(ParInfo& lambda)'



  bool SLambdaJetHunter::IsInHepMCDecayChain(const int idToFind, HepMC::GenVertex* vtxStart) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 5)) {
      cout << "SLambdaJetHunter::IsInHepMCDecayChain(int, HepMC::GenVertex*) checking if particle is in HepMC decay chain" << endl;
    }

    // make sure vectors are clear and add initial vertex
    m_vecVtxChecking.clear();
    m_vecVtxToCheck.clear();
    m_vecVtxToCheck.push_back(vtxStart);

    // breadth-first search of all connected vertices for particle w/ barcode 'idToFind'
    int  nVtxChecked     = 0;
    bool isParticleFound = false;
    while (!isParticleFound && (nVtxChecked < m_const.maxVtxToCheck)) {

      // transfer vertices to-check to being-checked list
      m_vecVtxChecking.clear();
      for (auto vertex : m_vecVtxToCheck) {
        m_vecVtxChecking.push_back(vertex);
      }
      m_vecVtxToCheck.clear();

      // iterate over vertices being checked
      for (auto vertex : m_vecVtxChecking) {

        // iterate over particles in vertex
        HepMC::GenVertex::particles_out_const_iterator outPar;
        for (
          outPar = vertex -> particles_out_const_begin();
          outPar != vertex -> particles_out_const_end();
          ++outPar
        ) {
          if ((*outPar) -> barcode() == idToFind) {
            isParticleFound = true;
            break;
          } else {
            m_vecVtxToCheck.push_back((*outPar) -> end_vertex());
          }
        }  // end particle loop

        // if found particler, exit
        //   - otherwise increment no. of vertices checked
        if (isParticleFound) {
          break;
        } else {
          ++nVtxChecked;
        }
      }  // end vertex loop
    }  // end while (!isParticleFound)
    return isParticleFound;

  }  // end 'IsInHepMCDecayChain(int, HepMC::GenVertex*)'



  bool SLambdaJetHunter::IsInPHG4DecayChain(const int idToFind, const int idLambda, PHCompositeNode* topNode) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 5)) {
      cout << "SLambdaJetHunter::IsInPHG4DecayChain(int) checking if particle is in Geant decay chain" << endl;
    }

    // grab truth information
    PHG4TruthInfoContainer*            container = GetTruthContainer(topNode);
    PHG4TruthInfoContainer::ConstRange particles = container -> GetParticleRange();

    // loop over particles
    m_vecIDToCheck.clear();
    for (
      PHG4TruthInfoContainer::ConstIterator itPar = particles.first;
      itPar != particles.second;
      ++itPar
    ) {

      // grab current particle
      PHG4Particle* phPar = itPar -> second;

      // grab parent if information is available
      PHG4Particle* phParent = NULL;
      if (phPar -> get_parent_id() != 0) {
        phParent = container -> GetParticle(phPar -> get_parent_id());
      } else {
        continue;
      }

      // check if parent matches lambda barcode/pid
      const bool hasLambdaBarcode = (phParent -> get_barcode() == idLambda);
      const bool hasLambdaPID     = (phParent -> get_pid()     == m_const.pidLambda);
      if (hasLambdaBarcode && hasLambdaPID) {
            m_vecIDToCheck.push_back(phPar -> get_barcode());
      }
    }  // end particle loop

    // check if idToFind is among barcodes
    bool isParticleFound = false;
    for (const int idToCheck : m_vecIDToCheck) {
      if (idToCheck == idToFind) {
        isParticleFound = true;
        break;
      }
    }
    return isParticleFound;

  }  // end 'IsInPHG4DecayChain(int, PHG4Particle*, PHCompositeNode*)'



  bool SLambdaJetHunter::IsGoodParticle(ParInfo& particle) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 5)) {
      cout << "SLambdaJetHunter::IsGoodParticle(ParInfo&) checking if particle is good" << endl;
    }

    // check charge if needed
    bool isGoodCharge = true;
    if (m_config.isCharged) {
      isGoodCharge = (particle.charge != 0.);
    }

    // FIXME overloaded <, etc. operators aren't behaving as expected
    const bool isInPtAccept  = ((particle.pt > m_config.parAccept.first.pt) && (particle.pt < m_config.parAccept.second.pt));
    const bool isInEtaAccept = ((particle.eta > m_config.parAccept.first.eta) && (particle.eta < m_config.parAccept.second.eta));

    // run other checks and return
    const bool isInAccept = (isInPtAccept && isInEtaAccept);
    return (isGoodCharge && isInAccept);

  }  // end 'IsGoodParticle(ParInfo&)'



  bool SLambdaJetHunter::IsGoodLambda(ParInfo& lambda) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 5)) {
      cout << "SLambdaJetHunter::IsGoodLambda(ParInfo&) checking if lambda is good" << endl;
    }

    // FIXME overloaded <, etc. operators aren't behaving as expected
    const bool isInPtAccept  = ((lambda.pt > m_config.parAccept.first.pt) && (lambda.pt < m_config.parAccept.second.pt));
    const bool isInEtaAccept = ((lambda.eta > m_config.parAccept.first.eta) && (lambda.eta < m_config.parAccept.second.eta));

    // make sure lambda is in acceptance
    const bool isInAccept = (isInPtAccept && isInEtaAccept);
    return isInAccept;

  }  // end 'IsGoodLambda(ParInfo&)'



  bool SLambdaJetHunter::IsLambda(const int pid) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 5)) {
      cout << "SLambdaJetHunter::IsLambda(int) checking if particle is a lambda" << endl;
    }

    return (pid == m_const.pidLambda);

  }  // end 'IsLambda(int)'



  bool SLambdaJetHunter::IsNewLambda(const int id) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 5)) {
      cout << "SLambdaJetHunter::IsNewLambda(int) checking if lambda is already found" << endl;
    }

    return (m_mapLambdaJetAssoc.count(id) == 0);

  }  // end 'IsNewLambda(int)'



  bool SLambdaJetHunter::HasParentInfo(const int parent) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 5)) {
      cout << "SLambdaJetHunter::HasParentInfo(int) checking if PHG4particle has parent info" << endl;
    }

    return (parent != 0);

  }  // end 'HasParentInfo(int)'



  double SLambdaJetHunter::GetLambdaAssocZ(ParInfo& lambda) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 5)) {
      cout << "SLambdaJetHunter::GetLambdaAssocZ(ParInfo&) getting z of lambda relative to associated jet" << endl;
    }

    // get index of associated jet
    const int iAssoc = m_mapLambdaJetAssoc[lambda.barcode];

    // get total momenta
    const double pJet    = hypot(m_jetInfo.at(iAssoc).px, m_jetInfo.at(iAssoc).py, m_jetInfo.at(iAssoc).pz);
    const double pLambda = hypot(lambda.px, lambda.py, lambda.pz);

    // calculate z and return
    const double zLambda = pLambda / pJet;
    return zLambda;

  }  // end 'GetLambdaAssocZ(ParInfo& lambda)'



  double SLambdaJetHunter::GetLambdaAssocDr(ParInfo& lambda) {

    // print debug statement
    if (m_config.isDebugOn && (m_config.verbosity > 5)) {
      cout << "SLambdaJetHunter::GetLambdaAssocDr(ParInfo&) getting delta-R of lambda relative to associated jet" << endl;
    }

    // get index of associated jet
    const int iAssoc = m_mapLambdaJetAssoc[lambda.barcode];

    // calculate delta-phi
    double dfLambda = lambda.phi - m_jetInfo.at(iAssoc).phi;
    if (dfLambda < 0.)             dfLambda += TMath::TwoPi();
    if (dfLambda > TMath::TwoPi()) dfLambda -= TMath::TwoPi();

    // calculate dr and return
    const double dhLambda = lambda.eta - m_jetInfo.at(iAssoc).eta;
    const double drLambda = hypot(dfLambda, dhLambda);
    return drLambda;

  }  // end 'GetLambdaAssocZ(ParInfo& lambda)'

}  // end SColdQcdCorrelatorAnalysis namespace

#endif

// end ------------------------------------------------------------------------
