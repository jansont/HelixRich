/***** The most basic main file ****/

// Our required includes
#include "G4RunManager.hh"
#include "G4UImanager.hh"

// Classes we will be pulling from for our construction
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "StackingAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"

//Project/McGillPhysics/HelixRich/Code/Applications/HelixRich/HR-build
//source ../../../geant4-install/bin/geant4.sh

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#include "G4TrajectoryDrawByParticleID.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

// add "int argc, char** argv" to int main if you'd like to run arguments
int main(int argc,char** argv)
{
  // Construct the default run manager
  // App Dev Guide notes that this is the ONLY manager class in G4 kernel that should be explicitly constructed in main()
  G4RunManager* runManager = new G4RunManager;

  // Create visualization manager (Not NEEDED, but useful)
  //G4VisManager* visManager = VisManager;
  //visManager -> initialize();

  // Set MANDATORY initialization classes
  runManager -> SetUserInitialization(new DetectorConstruction);
  runManager -> SetUserInitialization(new PhysicsList);

  // Set MANDATORY user action class
  runManager -> SetUserAction(new PrimaryGeneratorAction);
  //Other action classes
  runManager->SetUserAction(new StackingAction);
  RunAction* run_action = new RunAction;
  runManager->SetUserAction(new RunAction);
  EventAction* event_action = new EventAction(run_action);
  runManager->SetUserAction(event_action);
  SteppingAction* stepping_action = new SteppingAction(event_action);
  runManager->SetUserAction(stepping_action);

  #ifdef G4VIS_USE
    // visualization manager
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();
    G4TrajectoryDrawByParticleID* model = new G4TrajectoryDrawByParticleID;
    //  G4TrajectoryDrawByParticleID* model2 = new G4TrajectoryDrawByParticleID("test");

    model->SetDefault("cyan");
    model->Set("opticalphoton", "cyan");
    model->Set("gamma", "green");
    model->Set("e+", "magenta");
    model->Set("e-", G4Colour(1.0, 0.5, 0.0)); //orange
    model->Set("alpha", "yellow");

    visManager->RegisterModel(model);
  #endif

  // Run the kernel
  runManager -> Initialize();
  // Get the pointer to the User Interface manager
  //
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if (argc==1)   // Define UI session for interactive mode
    {
#ifdef G4UI_USE
      G4UIExecutive * ui = new G4UIExecutive(argc,argv);
#ifdef G4VIS_USE
      UImanager->ApplyCommand("/control/execute vis.mac");
#endif
      ui->SessionStart();
      delete ui;
#endif
    }
  else         // Batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
    }


    
  // get pointer to UI manager and set verbosities
  G4UImanager* UI = G4UImanager::GetUIpointer();
  UI -> ApplyCommand("/run/verbose 1");
  UI -> ApplyCommand("/event/verbose 1");
  UI -> ApplyCommand("/tracking/verbose 1");

  // start a run
  int numberOfEvents = 1;
  runManager -> BeamOn(numberOfEvents);

  // TERMINATE job
  // Remember to always clean up. Any time we create something with new we should delete it
  #ifdef G4VIS_USE
  delete visManager;
#endif

  delete runManager;
  //delete visManager; 
  return 0;
}
/*****************
* DETECTOR CONSTRUCTION specifies GEOMETRY, MATERIALS, SENSITIVE REGIONS, and READOUT schemes of sensitive regions
* PARTICLE LIST (from G4VUserPhysicsList) requires user define PARTICLES being used, RANGE CUTS for particles, and ALL PHYSICS PROCESSES to be simulated
*
* GEANT does not check for mandatory classes till initialize() and BeamOn() are invoked 
*
* OPTIONAL USER CLASSES
* UserRunAction, EventAction, StackingAction, TrackingAction, SteppingAction
******************/
