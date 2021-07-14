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

#include "TrackerSD.hh"
#include "TrackerHit.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "SimulationConstants.hh"


TrackerSD::TrackerSD(G4String detector_name, G4VPhysicalVolume *detector)
 : G4VSensitiveDetector(detector_name), savetime(true), saveposition(true),saveenergy(true), detector(detector)
{
	G4String HCname;
	collectionName.insert(HCname="HelixRichHitsCollection");
}

TrackerSD::~TrackerSD(){}

void TrackerSD::Initialize(G4HCofThisEvent* HCE)
{
  photonCollection = new TrackerHitsCollection(SensitiveDetectorName,collectionName[0]);
  static G4int HCID = -1;
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection( HCID, photonCollection );
  outFile.open("HelixRich.out", std::ofstream::app);
}

G4bool TrackerSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  G4ThreeVector position = aStep->GetPostStepPoint()->GetPosition();
  G4double x = position.X;
  G4double y = position.Y;
  G4double z = position.Z;
  // G4double energy = aStep->GetPostStepPoint()->GetTotalEnergy();
  G4double energy = aStep->GetPostStepPoint()->GetKineticEnergy();


  if (SimulationConstants::particle_to_detect == "primary"){
    // if (aStep->GetPostStepPoint()->GetPhysicalVolume() != detector)
    //   return false;
    if (aStep->GetTrack()->GetDefinition()->GetParticleType() == "opticalphoton")
      return false;
  } 
  else if (SimulationConstants::particle_to_detect == "optical"){
    if (aStep->GetPostStepPoint()->GetPhysicalVolume() != detector)
      return false;
     if (aStep->GetTrack()->GetDefinition()->GetParticleType() != "opticalphoton")
      return false; 
  }
  else if (SimulationConstants::particle_to_detect == "all"){
    if (aStep->GetPostStepPoint()->GetPhysicalVolume() != detector)
      return false;
  }
  else{
    return false;
  }

  TrackerHit* newHit = new TrackerHit();
  newHit->SetEnergy(energy);
  newHit->SetTime  (aStep->GetPostStepPoint()->GetGlobalTime());
  newHit->SetPos   (aStep->GetPostStepPoint()->GetPosition());
  photonCollection->insert( newHit );
  return true;
}

void  TrackerSD::EndOfEvent(G4HCofThisEvent* HCE){
  // save results
  G4int NbHits = photonCollection->entries();

  if (outFile) {
    G4cout<< "Here1 "<<G4endl;
    outFile << "# " << NbHits << " Hits detected." << std::endl;
    outFile << "# ";
    if (savetime)
      outFile << "Time (ns)";
    if (saveposition) {
      if (savetime)
      outFile << "\tPosition (x in mm)\tPosition (y in mm)\tPosition (z in mm)\t";
    }
    if (saveenergy) {
      if (savetime || saveposition)
      outFile << "Energy (eV)";
    }

    outFile << std::endl;
    for (G4int i=0;i<NbHits;i++) (*photonCollection)[i]->Print(outFile, savetime, saveposition, saveenergy);
    outFile << std::endl << std::endl;
    outFile.close();
  }

  else {

    G4cout << "=============" << G4endl << "Sensitive Detector : " << NbHits << " optical photons detected" << G4endl;
    
    outFile.open("HelixRich.out", std::ofstream::app);
    if (outFile.is_open()) G4cout<< " File Open" <<G4endl;
    else G4cout<< "File not open"<<G4endl;
    G4cout << "# " << NbHits << " Hits detected." <<G4endl;

    outFile << "# " << NbHits << " Hits detected." << std::endl;
    outFile << "# ";


    if (savetime)
      outFile << "Time (ns)";
    if (saveposition) {
      if (savetime)
      outFile << "\tPosition (x in mm)\tPosition (y in mm)\tPosition (z in mm)\t";
    }
    if (saveenergy) {
      if (savetime || saveposition)
      outFile << "Energy (eV)";
    }
    outFile << std::endl;
    for (G4int i=0;i<NbHits;i++) (*photonCollection)[i]->Print(outFile, savetime, saveposition, saveenergy);
    outFile << std::endl << std::endl;

    outFile.close();
  }


}

// void  TrackerSD::clear(){} 

// void  TrackerSD::DrawAll(){ } 

// void  TrackerSD::PrintAll(){ } 