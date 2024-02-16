// ----------------------------------------------------------------------------
// 'Fun4All_RunLambdaJetHunter.C'
// Derek Anderson
// 01.25.2024
//
// F4A macro the SLambdaJetHunter module.
// ----------------------------------------------------------------------------

#define FUN4ALL_RUNLAMBDAJETHUNTER_C

// c++ utilities
#include <vector>
#include <string>
// f4a/sphenix libraries
#include <FROG.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <g4main/Fun4AllDstPileupInputManager.h>
// analysis specific utilities
#include "LambdaJetHunterOptions.h"
#include "/sphenix/user/danderson/install/include/slambdajethunter/SLambdaJetHunter.h"
#include "/sphenix/user/danderson/install/include/slambdajethunter/SLambdaJetHunterConfig.h"

// load libraries
R__LOAD_LIBRARY(/sphenix/user/danderson/install/lib/libslambdajethunter.so)
R__LOAD_LIBRARY(/sphenix/user/danderson/install/lib/libscorrelatorutilities.so)

// make common namespaces implicit
using namespace std;
using namespace SColdQcdCorrelatorAnalysis;

// global constants
static const int            NEvtDefault      = 10;
static const int            VerbDefault      = 0;
static const string         SOutDirDefault   = "/sphenix/user/danderson/eec/SLambdaJetHunter/submit/ppJet10GeV/full/output";
static const string         SRecoNameDefault = "LambdaJetTree";
static const vector<string> SInDefault       = {"test.list"};



// macro body -----------------------------------------------------------------

void Fun4All_RunLambdaJetHunterOnCondor(
  const vector<string>& sInputLists = SInDefault,
  const int             nEvents     = NEvtDefault,
  const int             verbosity   = VerbDefault
) {

  // load libraries and create f4a server
  gSystem -> Load("libg4dst.so");
  gSystem -> Load("libFROG.so");

  FROG*          frog      = new FROG();
  Fun4AllServer* ffaServer = Fun4AllServer::instance();
  ffaServer -> Verbosity(verbosity);

  // figure out folder revisions, file numbers, etc
  string outDir   = SOutDirDefault;
  string recoName = SRecoNameDefault;
  if (outDir.substr(outDir.size() - 1, 1) != "/") {
    outDir += "/";
  }
  outDir += recoName + "/";

  string fileNumber   = sInputLists[0];
  size_t findLastDash = fileNumber.find_last_of("-");
  if (findLastDash != string::npos) {
    fileNumber.erase(0, findLastDash + 1);
  }
  string remove_this = ".list";

  size_t pos = fileNumber.find(remove_this);
  if (pos != string::npos) {
    fileNumber.erase(pos, remove_this.length());
  }

  // create output name
  string outputFileName = "outputData_" + recoName + "_" + fileNumber + ".root";
  string outputRecoDir  = outDir + "/inReconstruction/";
  string makeDirectory  = "mkdir -p " + outputRecoDir;

  system(makeDirectory.c_str());
  string outputRecoFile = outputRecoDir + outputFileName;

  // add input files 
  for (size_t iInput = 0; iInput < sInputLists.size(); iInput++) {
    Fun4AllDstInputManager* inManager = new Fun4AllDstInputManager("InputDstManager" + to_string(iInput));
    inManager -> AddListFile(sInputLists.at(iInput));
    ffaServer -> registerInputManager(inManager);
  }

  // set module configuration
  SLambdaJetHunterConfig config = LambdaJetHunterOptions::GetConfig(verbosity, outputRecoFile);

  // instantiate & register lambda jet hunter
  SLambdaJetHunter* hunter = new SLambdaJetHunter(config);
  ffaServer -> registerSubsystem(hunter);

  // move output and clean up logs
  ifstream file(outputRecoFile.c_str());
  if (file.good()) {
    string moveOutput = "mv " + outputRecoFile + " " + outDir;
    system(moveOutput.c_str());
  } else {
    string rmOutput = "rm " + outDir + recoName + fileNumber + ".root";
    system(rmOutput.c_str());
  }

  // run reconstruction & close f4a
  ffaServer -> run(nEvents);
  ffaServer -> End();
  delete ffaServer;

  // announce end & exit
  gSystem -> Exit(0);
  return;

}

// end ------------------------------------------------------------------------
