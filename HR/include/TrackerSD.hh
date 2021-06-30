#ifndef TrackerSD_h
#define TrackerSD_h 1

#include "G4VSensitiveDetector.hh"
#include "TrackerHit.hh"
#include "StackingAction.hh"

class G4Step;
class G4HCofThisEvent;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class TrackerSD : public G4VSensitiveDetector
{
  public:
      TrackerSD(G4String, G4VPhysicalVolume *detector);
     ~TrackerSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*, G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);

      enum SaveData {
        Time,
        Position,
        Energy
      };

      // void setSaveOption(SaveData type, bool save);

  private:
      TrackerHitsCollection* photonCollection;
      bool savetime;
      bool saveposition;
      bool saveenergy;
      const G4VPhysicalVolume *detector;
      std::ofstream outFile;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

