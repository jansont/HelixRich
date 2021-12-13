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
	if (SimulationConstants::TRACKING)
	{
		eFile.open("ElectronTracking.out", std::ofstream::app);
		if (eFile)
		{
			eFile << "# Pre Step Z (mm)" << "\t\t\t";
			// eFile << "Post Step Z (mm)" << "\t\t\t";
			eFile << "Post Step Energy (mm)" << "\t\t\t"<< std::endl;
		}
		eFile.close();

		oFile.open("PhotonTracking.out", std::ofstream::app);
		if (oFile)
		{
			oFile << "# Track ID" << "\t\t\t";

			oFile << "Step Number" << "\t\t\t";
			oFile << "X (mm)" << "\t\t\t";
			oFile << "Y (mm)" << "\t\t\t";
			oFile << "Z (mm)" << "\t\t\t"<< std::endl;
		}
		oFile.close();
	}

	fCerenkovCounter   = 0;

}

SteppingAction::~SteppingAction()
{}

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
	if (SimulationConstants::TRACKING)
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

		if (aParticle->GetDefinition() == G4OpticalPhoton::OpticalPhoton())
		{
			G4int track_id = eTrack->GetTrackID();
			const G4ThreeVector vertex = eTrack->GetVertexPosition();
			G4double z_loc = (G4double) vertex.z();
			G4double x_loc = (G4double) vertex.x();
			G4double y_loc = (G4double) vertex.y();
			G4int stepNum = eTrack->GetCurrentStepNumber();


			eFile.open("PhotonTracking.out", std::ofstream::app);
			if (eFile)
			{
				if (stepNum == 1)
				{
					eFile << track_id << "\t\t\t";
					eFile << stepNum  << "\t\t\t";
					eFile << x_loc  << "\t\t\t";
					eFile << y_loc  << "\t\t\t";
					eFile << z_loc << "\t\t\t"<< std::endl;
				}
			}
			eFile.close();
		}

	}


	G4int eventNumber = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
     G4cout << " Number of Cerenkov Photons in previous event: " << fCerenkovCounter << G4endl;
       

       G4Track* track = aStep->GetTrack();

  G4String ParticleName = track->GetDynamicParticle()->
                                 GetParticleDefinition()->GetParticleName();




       if (ParticleName == "opticalphoton") return;



  const std::vector<const G4Track*>* secondaries =
                                            aStep->GetSecondaryInCurrentStep();

  if (secondaries->size()>0) {
     for(unsigned int i=0; i<secondaries->size(); ++i) {
        if (secondaries->at(i)->GetParentID()>0) {
           if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition()
               == G4OpticalPhoton::OpticalPhotonDefinition()){
              if (secondaries->at(i)->GetCreatorProcess()->GetProcessName()== "Cerenkov")
              	fCerenkovCounter++;
           }
        }
     }
  }
}

