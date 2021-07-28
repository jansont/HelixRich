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

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "SimulationConstants.hh"
#include "G4ThreeVector.hh"
#include "G4OpticalPhoton.hh"
#include "G4Electron.hh"
#include <fstream>

SteppingAction::SteppingAction(EventAction* eventAction):
G4UserSteppingAction(),
fEventAction(eventAction),
fScoringVolume(0)
{
	eFile.open("ElectronTracking.out", std::ofstream::app);
	if (eFile)
	{
		eFile << "# Pre Step Z (mm)" << "\t\t\t";
		// eFile << "Post Step Z (mm)" << "\t\t\t";
		eFile << "Post Step Energy (mm)" << "\t\t\t"<< std::endl;
	}
	eFile.close();
}

SteppingAction::~SteppingAction()
{}

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
	G4StepPoint* pPreStepPoint  = aStep ->GetPreStepPoint();
	G4StepPoint* pPostStepPoint = aStep ->GetPostStepPoint();
	const G4ThreeVector prePos= pPreStepPoint->GetPosition();
	const G4ThreeVector postPos= pPostStepPoint->GetPosition();
	const G4double preEnergy  = pPreStepPoint->GetKineticEnergy();
	const G4double  postEnergy  =  pPostStepPoint->GetKineticEnergy();

	G4Track* eTrack = aStep -> GetTrack();
	const G4DynamicParticle* aParticle = eTrack->GetDynamicParticle();
	const G4double energy = aParticle->GetKineticEnergy();

	if (aParticle->GetDefinition() != G4OpticalPhoton::OpticalPhoton())
	{
		G4String preVol = pPreStepPoint->GetPhysicalVolume()->GetName();
		// G4String postVol = pPostStepPoint->GetPhysicalVolume()->GetName();

		G4String vol = "Tile";
		if (preVol == vol)
		{
			eFile.open("ElectronTracking.out", std::ofstream::app);
			if (eFile){
				eFile << (G4double) prePos.z() << "\t\t\t";
				// eFile << (G4double) postPos.z() << "\t\t\t";
				eFile << (G4double) postEnergy << "\t\t\t"<< std::endl;
			}
			eFile.close();
		}
	}
}

