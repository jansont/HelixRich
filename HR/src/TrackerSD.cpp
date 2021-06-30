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


  if (SimulationConstants::electron){
    // if (aStep->GetTrack()->GetDefinition()->GetParticleType() != "e-")
    //   return false; 
    // if  (x > windowX/2 || y > windowY/2 )
    //   return false;
    // if(z > 200.00*mm || z < 200.05*mm )
    // return false;
  } else {
    if (aStep->GetPostStepPoint()->GetPhysicalVolume() != detector)
      return false;
     if (aStep->GetTrack()->GetDefinition()->GetParticleType() == "e-")
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
        outFile << "\t";
      outFile << "Position (x in mm)\tPosition (y in mm)\tPosition (z in mm)\t";
    }
    if (saveenergy) {
      if (savetime || saveposition)
        outFile << "\t";
      outFile << "Energy (MeV)";
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
        outFile << "\t";
      outFile << "Position (x in mm)\tPosition (y in mm)\tPosition (z in mm)\t";
    }
    if (saveenergy) {
      if (savetime || saveposition)
        outFile << "\t";
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