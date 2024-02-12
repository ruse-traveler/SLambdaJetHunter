// ----------------------------------------------------------------------------
// 'SLambdaJetHunter.h'
// Derek Anderson
// 01.25.2024
//
// A minimal analysis module to find lambda-tagged
// jets in pythia events.
// ----------------------------------------------------------------------------

#ifndef SLAMBDAJETHUNTER_H
#define SLAMBDAJETHUNTER_H

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

// c++ utilities
#include <map>
#include <string>
#include <vector>
// root libraries
#include <TFile.h>
#include <TTree.h>
// fastjet libraries
#include <fastjet/PseudoJet.hh>
#include <fastjet/JetDefinition.hh>
// hepmc libraries
#include <HepMC/GenVertex.h>
#include <HepMC/GenParticle.h>
// f4a utilities
#include <fun4all/SubsysReco.h>
// analysis utilities
#include "SLambdaJetHunterConfig.h"
#include "/sphenix/user/danderson/install/include/scorrelatorutilities/EvtTools.h"
#include "/sphenix/user/danderson/install/include/scorrelatorutilities/JetTools.h"
#include "/sphenix/user/danderson/install/include/scorrelatorutilities/CstTools.h"
#include "/sphenix/user/danderson/install/include/scorrelatorutilities/GenTools.h"

#pragma GCC diagnostic pop

// make common namespaces implicit
using namespace std;
using namespace fastjet;
using namespace SColdQcdCorrelatorAnalysis::SCorrelatorUtilities;

// forward declarations
class PHCompositeNode;



namespace SColdQcdCorrelatorAnalysis {

  // SLambdaJetHunter definition ----------------------------------------------

  class SLambdaJetHunter : public SubsysReco {

     public:

       // ctor/dtor
       SLambdaJetHunter(const string &name = "SLambdaJetHunter", const bool debug = false);
       SLambdaJetHunter(SLambdaJetHunterConfig& config);
       ~SLambdaJetHunter() override;

       // f4a methods
      int Init(PHCompositeNode *topNode)          override;
      int process_event(PHCompositeNode *topNode) override;
      int End(PHCompositeNode *topNode)           override;

      // setters
      void SetConfig(SLambdaJetHunterConfig& config) {m_config = config;}

      // getters
      SLambdaJetHunterConfig GetConfig() {return m_config;}

    private:

      // i/o members
      TFile* m_outFile = NULL;
      TTree* m_outTree = NULL;

      // module configuration
      SLambdaJetHunterConfig m_config;

      // output variables
      GenInfo                 m_genEvtInfo;
      vector<JetInfo>         m_jetInfo;
      vector<vector<CstInfo>> m_cstInfo;
      vector<vector<ParInfo>> m_lambdaInfo;

      // vectors for internal calculations
      vector<int>               m_vecSubEvts;
      vector<PseudoJet>         m_jets;
      vector<HepMC::GenVertex*> m_vecVtxToCheck;
      vector<HepMC::GenVertex*> m_vecVtxChecking;

      // jet-lambda associations
      map<int, int> m_mapLambdaJetAssoc;

      // output event variables
      //   - FIXME remove when i/o of utility structs is ready
      uint64_t m_evtNJets;
      uint64_t m_evtNLambdas;
      uint64_t m_evtNTaggedJets;
      uint64_t m_evtNChrgPars;
      uint64_t m_evtNNeuPars;
      double   m_evtSumEPar;
      double   m_evtVtxX;
      double   m_evtVtxY;
      double   m_evtVtxZ;
      // output parton variables
      pair<int,   int>     m_evtPartID;
      pair<double, double> m_evtPartPx;
      pair<double, double> m_evtPartPy;
      pair<double, double> m_evtPartPz;
      pair<double, double> m_evtPartE;
      // output jet variables
      vector<bool>     m_jetHasLambda;
      vector<uint64_t> m_jetNCst;
      vector<uint64_t> m_jetID;
      vector<double>   m_jetE;
      vector<double>   m_jetPt;
      vector<double>   m_jetEta;
      vector<double>   m_jetPhi;
      // output constituent variables
      vector<vector<int>>    m_cstID;
      vector<vector<int>>    m_cstJetID;
      vector<vector<int>>    m_cstEmbedID;
      vector<vector<double>> m_cstZ;
      vector<vector<double>> m_cstDr;
      vector<vector<double>> m_cstE;
      vector<vector<double>> m_cstPt;
      vector<vector<double>> m_cstEta;
      vector<vector<double>> m_cstPhi;
      // output lambda variables
      vector<vector<int>>    m_lambdaID;
      vector<vector<int>>    m_lambdaJetID;
      vector<vector<int>>    m_lambdaEmbedID;
      vector<vector<double>> m_lambdaZ;
      vector<vector<double>> m_lambdaDr;
      vector<vector<double>> m_lambdaE;
      vector<vector<double>> m_lambdaPt;
      vector<vector<double>> m_lambdaEta;
      vector<vector<double>> m_lambdaPhi;

      // analysis methods (*.ana.h)
      void GrabEventInfo(PHCompositeNode* topNode);
      void MakeJets(PHCompositeNode* topNode);
      void HuntLambdas(PHCompositeNode* topNode);
      void FillOutputTree();
      bool IsGoodParticle(ParInfo& particle);
      bool IsLambda(const int pid);
      bool IsNewLambda(const int id);
      bool IsInHepMCDecayChain(const int idToFind, HepMC::GenVertex* vtxToStart);
      bool IsInPHG4DecayChain(const int idToFind);

      // system methods (*.sys.h)
      void InitTree();
      void InitOutput();
      void SaveAndCloseOutput();
      void ResetOutput();

      // class-wide constants
      struct Const {
        int pidLambda;
        int maxVtxToCheck;
      }  m_const = {3122, 500};

      // map of user input onto fastjet options
      map<string, fastjet::JetAlgorithm> m_mapJetAlgo = {
        {"kt_algorithm",                    fastjet::JetAlgorithm::kt_algorithm},
        {"cambridge_algorithm",             fastjet::JetAlgorithm::cambridge_algorithm},
        {"antikt_algorithm",                fastjet::JetAlgorithm::antikt_algorithm},
        {"genkt_algorithm",                 fastjet::JetAlgorithm::genkt_algorithm},
        {"cambridge_for_passive_algorithm", fastjet::JetAlgorithm::cambridge_for_passive_algorithm},
        {"genkt_for_passive_algorithm",     fastjet::JetAlgorithm::genkt_for_passive_algorithm},
        {"ee_kt_algorithm",                 fastjet::JetAlgorithm::ee_kt_algorithm},
        {"ee_genkt_algorithm",              fastjet::JetAlgorithm::ee_genkt_algorithm},
        {"plugin_algorithm",                fastjet::JetAlgorithm::plugin_algorithm}
      };
      map<string, fastjet::RecombinationScheme> m_mapRecomb = {
        {"E_scheme",        fastjet::RecombinationScheme::E_scheme},
        {"pt_scheme",       fastjet::RecombinationScheme::pt_scheme},
        {"pt2_scheme",      fastjet::RecombinationScheme::pt2_scheme},
        {"Et_scheme",       fastjet::RecombinationScheme::Et_scheme},
        {"Et2_scheme",      fastjet::RecombinationScheme::Et2_scheme},
        {"BIpt_scheme",     fastjet::RecombinationScheme::BIpt_scheme},
        {"BIpt2_scheme",    fastjet::RecombinationScheme::BIpt2_scheme},
        {"WTA_pt_scheme",   fastjet::RecombinationScheme::WTA_pt_scheme},
        {"WTA_modp_scheme", fastjet::RecombinationScheme::WTA_modp_scheme},
        {"external_scheme", fastjet::RecombinationScheme::external_scheme}
      };

  };  // end SLambdaJetHunter

}  // end SColdQcdCorrelatorAnalysis namespace

#endif

// end ------------------------------------------------------------------------
