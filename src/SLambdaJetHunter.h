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
// f4a utilities
#include <fun4all/SubsysReco.h>
// analysis utilities
#include "SLambdaJetHunterConfig.h"
#include "/sphenix/user/danderson/install/include/scorrelatorutilities/EvtTools.h"
#include "/sphenix/user/danderson/install/include/scorrelatorutilities/JetTools.h"
#include "/sphenix/user/danderson/install/include/scorrelatorutilities/CstTools.h"
#include "/sphenix/user/danderson/install/include/scorrelatorutilities/GenTools.h"

// make common namespaces implicit
using namespace std;
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

      // fastjet output
      vector<fastjet::PseudoJet> m_jets;

      // analysis methods (*.ana.h)
      void    MakeJets(PHCompositeNode* topNode);
      void    GrabSubEvents();
      void    GrabOutputInfo();
      bool    IsGoodParticle(ParInfo& particle);
      bool    IsLambda();
      ParInfo FindLambda(const int barcode);

      // system methods (*.sys.h)
      void InitTree();
      void InitOutput();
      void SaveAndCloseOutput();
      void ResetOutput();

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
