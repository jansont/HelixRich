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

#include "StackingAction.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Track.hh"

StackingAction::StackingAction(): photonCounter(0) {}
StackingAction::~StackingAction(){}

G4ClassificationOfNewTrack StackingAction::ClassifyNewTrack(const G4Track * aTrack)
{
	if(aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition())
	{ // particle is optical photon
		if(aTrack->GetParentID()>0) // particle is secondary
			photonCounter++;
	}
	return fUrgent;
}

void StackingAction::NewStage()
{
	G4cout << "Number of optical photons produced in this event : " << photonCounter << G4endl;
	outFile.open("HelixRich.out", std::ofstream::app);
	outFile << "# " << photonCounter << " Photons produced." << std::endl;
	outFile.close();
}

void StackingAction::PrepareNewEvent()
{
	photonCounter = 0;
}

G4int StackingAction::GetPhotonCounter()
{
	return photonCounter;
}

