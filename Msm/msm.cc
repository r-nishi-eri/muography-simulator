#include "MsmDetectorConstruction.hh"
#include "MsmActionInitialization.hh"

#include "G4SystemOfUnits.hh"

#ifdef G4MULTITHREADED
 #include "G4MTRunManager.hh"
#else
 #include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "FTFP_BERT.hh"
#include "Randomize.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif


namespace {
  void PrintUsage() {
    G4cerr << " Usage: " << G4endl;
    G4cerr << " exec_msm [-m macro] [-u UIsession] [-t demfile] [-d detectorfile] [-i injectionfile] [-log logfileprefix]" << G4endl;
  }
}


int main(int argc,char** argv)
{

  // Evaluate arguments
  if (argc < 2){
    PrintUsage();
    return 1;
  }
  
  G4String macro;
  G4String session;
  G4String str_numcore;
  G4String str_demfile = "demfile";
  G4String str_detectorfile = "detectorfile";
  G4String str_injectionfile = "injectionfile";
  G4String str_dest_injectionlog = "run";
  double cutlength = 1.0 * mm;
  int flg_numcore = false;
  int numcore;
  for ( G4int i=1; i<argc; i=i+2 ) {
    if      ( G4String(argv[i]) == "-m" ) macro = argv[i+1];
    else if ( G4String(argv[i]) == "-u" ) session = argv[i+1];
    else if ( G4String(argv[i]) == "-c" ) {
      str_numcore = argv[i+1];
      numcore = atoi(str_numcore.c_str());
      flg_numcore = true;
    }
    else if ( G4String(argv[i]) == "-t" ) {
      str_demfile = argv[i+1];
      G4cout << str_demfile;
    }
    else if ( G4String(argv[i]) == "-d") {
      str_detectorfile = argv[i+1];
      G4cout << str_detectorfile;
    }
    else if ( G4String(argv[i]) == "-i") {
      str_injectionfile = argv[i+1];
    }
    else if ( G4String(argv[i]) == "-log"){
      str_dest_injectionlog = argv[i+1];
    }
    else {
      PrintUsage();
      return 1;
    }
  }    

  // Random Engine Initialization
  CLHEP::MTwistEngine randomEngine;
  G4Random::setTheEngine(&randomEngine);
  G4int seed = time(NULL);
  G4Random::setTheSeed(seed);
  
  // Construct default run manager
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;
  if (flg_numcore) runManager->SetNumberOfThreads(numcore);
  else runManager->SetNumberOfThreads(1);
#else
  G4RunManager* runManager = new G4RunManager;
#endif

  
  // Set mandatory initialization classes
  // Detector Construction
  MsmDetectorConstruction* detConstruction = new MsmDetectorConstruction(str_demfile, str_detectorfile);
  runManager->SetUserInitialization(detConstruction);

  // Physics List
  G4VUserPhysicsList* physicsList = new FTFP_BERT();
  // Register Physics List
  runManager->SetUserInitialization(physicsList);

  // User Action Initialization
  runManager->SetUserInitialization(new MsmActionInitialization(str_injectionfile,str_dest_injectionlog));
    
  // Initialize G4 kernel
  runManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  
#ifdef G4VIS_USE
  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif


  if ( macro.size() ) {
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command+macro);
  }
  else  {  

#ifdef G4UI_USE
    G4UIExecutive* ui = new G4UIExecutive(argc, argv, session);
#ifdef G4VIS_USE
    UImanager->ApplyCommand("/control/execute init_vis.mac");
#else
    UImanager->ApplyCommand("/control/execute init.mac"); 
#endif
    if (ui->IsGUI())
      UImanager->ApplyCommand("/control/execute gui.mac");
    ui->SessionStart();
    delete ui;
#endif
  }

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}


