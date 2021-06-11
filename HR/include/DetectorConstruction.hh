
#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4MaterialTable.hh"


// These two classes are called within the construct function
class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Box;

class DetectorConstruction : public G4VUserDetectorConstruction
{
   public:
      DetectorConstruction();
      virtual ~DetectorConstruction();
      virtual G4VPhysicalVolume* Construct();

   private:
      void DefineMaterials();

      G4Box* solid_world;
      G4LogicalVolume* logical_world;
      G4VPhysicalVolume* physical_world;

      G4Box* solid_tile;
      G4LogicalVolume* logical_tile;
      G4VPhysicalVolume* physical_tile;

      G4Box* solid_detector;
      G4LogicalVolume* logical_detector;
      G4VPhysicalVolume* physical_detector;


};
#endif
