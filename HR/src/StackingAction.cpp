#include "StackingAction.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Track.hh"


StackingAction::StackingAction(): photonCounter(0) {}
StackingAction::~StackingAction(){}

G4ClassificationOfNewTrack StackingAction::ClassifyNewTrack(const G4Track * aTrack){
  if(aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()){ // particle is optical photon
    if(aTrack->GetParentID()>0){ // particle is secondary
      photonCounter++;
    }
  }
  return fUrgent;
}

void StackingAction::NewStage(){
  G4cout << "Number of optical photons produced in this event : " << photonCounter << G4endl;
  outFile.open("HelixRich.out", std::ofstream::app);
  outFile << "# " << photonCounter << " Photons produced." << std::endl;
  outFile.close();
}

void StackingAction::PrepareNewEvent()
{ photonCounter = 0; }

G4int StackingAction::GetPhotonCounter()
{return photonCounter; }

