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

#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

RunAction::RunAction(): G4UserRunAction(),fEdep(0.),fEdep2(0.)
{ 
	// Register accumulable to the accumulable manager
	G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
	accumulableManager->RegisterAccumulable(fEdep);
	accumulableManager->RegisterAccumulable(fEdep2); 
}

RunAction::~RunAction()
{}


void RunAction::BeginOfRunAction(const G4Run*)
{ 
	// inform the runManager to save random number seed
	G4RunManager::GetRunManager()->SetRandomNumberStore(false);

	// reset accumulables to their initial values
	G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
	accumulableManager->Reset();
}

void RunAction::EndOfRunAction(const G4Run* run)
{
	G4int nofEvents = run->GetNumberOfEvent();
	if (nofEvents == 0)
		return;

	// Merge accumulables 
	G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
	accumulableManager->Merge();

	// Compute dose = total energy deposit in a run and its variance
	G4double edep  = fEdep.GetValue();
	G4double edep2 = fEdep2.GetValue();

	G4double rms = edep2 - edep*edep/nofEvents;
	if (rms > 0.)
		rms = std::sqrt(rms);
	else
		rms = 0.;  
}

void RunAction::AddEdep(G4double edep)
{
	fEdep  += edep;
	fEdep2 += edep*edep;
}



