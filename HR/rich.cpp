//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// Helix Rich
// Author: Theodore Janson (theodore.janson@mail.mcgill.ca)



// Our required includes
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "StackingAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "SimulationConstants.hh"
#include "G4Cerenkov.hh"
#include "G4OpticalPhysics.hh"
#include "G4Material.hh"
#include "SimulationConstants.hh"
#include <string> 


#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#include "G4TrajectoryDrawByParticleID.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

void aRun(int argc,char** argv)
{
		// Construct the default run manager
	// App Dev Guide notes that this is the ONLY manager class in G4 kernel that should be explicitly constructed in main()
	G4RunManager* runManager = new G4RunManager;

	// Create visualization manager (Not NEEDED, but useful)
	//G4VisManager* visManager = VisManager;
	//visManager -> initialize();


    //choose the Random engine
    CLHEP::HepRandom::setTheEngine(new CLHEP::HepJamesRandom());;
    
    //set (somewhat) random seed with system time
    G4long seed = std::abs(time(NULL)*(time_t)&runManager);
    //G4long seed = 3;
    CLHEP::HepRandom::setTheSeed(seed);

	// Set MANDATORY initialization classes

	DetectorConstruction* detector = new DetectorConstruction();
	runManager -> SetUserInitialization(detector);

	PhysicsList* plist = new PhysicsList();
	// G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
	// plist->RegisterPhysics(opticalPhysics);

	runManager -> SetUserInitialization(plist);

	// Set MANDATORY user action class
	PrimaryGeneratorAction* primary_generator = new PrimaryGeneratorAction;
	runManager -> SetUserAction(primary_generator);
	// Set MANDATORY initialization classes


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

	if (SimulationConstants::ShowTerminalUI)
	{
		//Get the pointer to the User Interface manager
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
	}
	// start a run
	int numberOfEvents = SimulationConstants::runCount;
	runManager -> BeamOn(numberOfEvents);

	// TERMINATE job
	// Remember to always clean up. Any time we create something with new we should delete it
	#ifdef G4VIS_USE
		delete visManager;
	#endif

 	delete runManager;
}


int main(int argc,char** argv)
{
	// G4double* AerogelRindex = [1.150, 1.151, 1.152, 1.153, 1.154, 1.155, 1.156, 1.157, 1.158, 1.159, 1.160];
	// G4int length = (sizeof(AerogelRindex)/sizeof(*AerogelRindex));
	aRun(argc,argv);
	return 0;
}



