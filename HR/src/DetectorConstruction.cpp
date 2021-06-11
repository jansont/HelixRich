/******* We will create the most basics of a detector construction file  ********/
/*** 
 * To do this we will create a world box and place another box inside it and fill it. We need 
 * to create the logical volume, use the solid, and add other attributes (like composition)
 * ***/

#include "DetectorConstruction.hh"
#include "G4NistManager.hh"         // Contains the NIST data for particles
#include "G4Box.hh"                 // Needed to create box shapes
#include "G4LogicalVolume.hh"       // Will be required in every Detector construction file
#include "G4PVPlacement.hh"        
#include "G4RotationMatrix.hh"      
#include "G4Transform3D.hh"         
#include "G4PhysicalConstants.hh"   
#include "G4SystemOfUnits.hh"
#include "TrackerSD.hh"
#include "G4SDManager.hh"
#include "G4MaterialTable.hh"
#include "SimulationConstants.hh"


// constructor 
DetectorConstruction::DetectorConstruction()  : G4VUserDetectorConstruction(), physical_detector(0)
{
   DefineMaterials(); // Construct materials for detector
}

// destructor
DetectorConstruction::~DetectorConstruction()
{}  // Don't need to destruct anything, but we keep good habits

void DetectorConstruction::DefineMaterials()
{
  /*--------------------------------------Define Elements-----------------------------------*/
  G4double a, z;
  G4String name, symbol;
  // element 1: silicon
  a = 28.09*g/mole;
  G4Element* elemental_silicon = new G4Element(name="Silicon", symbol="Si", z=14., a);
  //element 2: oxygen
  a = 16.00*g/mole;
  G4Element* elemental_oxygen = new G4Element(name="Oxygen", symbol="O", z = 8., a);
  //element 3: hydrgoen
  a = 1.01*g/mole;
  G4Element* elemental_hydrogen = new G4Element(name="Hydrogen", symbol="H", z = 1., a);
  //element 4: carbon
  a = 12.01*g/mole;
  G4Element* elemental_carbon = new G4Element(name="Carbon", symbol="C", z = 6., a);

  G4NistManager* nist = G4NistManager::Instance();
  G4Element *N = nist->FindOrBuildElement("N");
  G4Element *C = nist->FindOrBuildElement("C");
  G4Element *H = nist->FindOrBuildElement("H");
  G4Element *Si = nist->FindOrBuildElement("Si");
  G4Element *O = nist->FindOrBuildElement("O");
  G4Element *B = nist->FindOrBuildElement("B");
  G4Element *Na = nist->FindOrBuildElement("Na");
  G4Element *Al = nist->FindOrBuildElement("Al");
  G4Element *K = nist->FindOrBuildElement("K");
  G4Element* F = nist->FindOrBuildElement("F");
  G4Element *Sb = nist->FindOrBuildElement("Sb");
  G4Element *Rb = nist->FindOrBuildElement("Rb");
  G4Element *Cs = nist->FindOrBuildElement("Cs");
  G4bool isotopes = false;

/*--------------------------------------Define Materials-----------------------------------*/
  // secondary materials: 
  G4int ncomponents;
  G4double density, fractionmass;
  G4String natoms;
  
  // water
  density = 1.000*g/cm3;
  G4Material* H2O = new G4Material(name="Water", density, ncomponents=2);
  H2O -> AddElement(H, 67*perCent);
  H2O -> AddElement(O, 33*perCent);
  // silica
  density = 2.200*g/cm3;
  G4Material* SiO2 = new G4Material(name="Silica", density, ncomponents=2);
  SiO2 -> AddElement(Si, 67*perCent);
  SiO2 -> AddElement(O , 33*perCent);

  // air 
  density = 1.290 * mg/cm3;
   G4Material* Air = new G4Material("Air", density, 2);
   Air -> AddElement(N, 70*perCent);
   Air -> AddElement(O, 30*perCent);

/*---------------------------------Aerogel: Material and properties-----------------------------------*/
  // Aerogel type
  G4Material* Aerogel;
  switch (Aerogel_type){
    case 1: 
      density = 0.200*g/cm3;
      Aerogel = new G4Material(name="Aerogel", density, ncomponents=3); 
      Aerogel -> AddMaterial(SiO2, fractionmass=97.*perCent);
      Aerogel -> AddMaterial(H2O , fractionmass=3.0*perCent);
      Aerogel -> AddElement (C, fractionmass= 0.1*perCent);
      break;

  }

  //define index of refraction for aerogel
  switch (Aerogel_properties){
    case 1:
      //understand everything below
        //PhotonEnergy
      G4double PhotonMinEnergy=1.3*CLHEP::eV;
      G4double PhotonMaxEnergy=7.3*CLHEP::eV;
      G4int NumPhotWaveLengthBins = 1000;
      G4int ibin=0;
      G4double PhotonEnergyStep = (PhotonMaxEnergy-PhotonMinEnergy)/ NumPhotWaveLengthBins;
      G4double* PhotonMomentum = new G4double[NumPhotWaveLengthBins];
      for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
        PhotonMomentum[ibin] = PhotonMinEnergy + PhotonEnergyStep * ibin;
      }
      G4MaterialPropertiesTable* AerogelMPT = new G4MaterialPropertiesTable(); //create a table for aerogel properties
      G4double* AerogelRindex = new G4double[NumPhotWaveLengthBins];        
      G4double* AerogelAbsorpLength=new G4double[NumPhotWaveLengthBins];
      for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
          AerogelAbsorpLength[ibin]=1.E32*mm;
          AerogelRindex[ibin]=1.15;
        }
      AerogelMPT->AddProperty("ABSLENGTH",PhotonMomentum,
      AerogelAbsorpLength,NumPhotWaveLengthBins);
      AerogelMPT->AddProperty("RINDEX", PhotonMomentum, AerogelRindex, NumPhotWaveLengthBins);
      Aerogel->SetMaterialPropertiesTable(AerogelMPT);
      break;
  }
}


// Now we will be creating the actual world and object inside the world here
G4VPhysicalVolume* DetectorConstruction::Construct()
{
   G4NistManager* nist = G4NistManager::Instance();
  //--------------------------------World Geometry------------------------------------------------------------
  G4Material* world_material;
  switch (world_material_type){
    case 1:
      world_material = nist -> FindOrBuildMaterial("Air"); //manager finds materials defined above 
      break;
  }

  solid_world = new G4Box("World", worldX, worldY, worldZ);   //simple volume (name, x_size, y_size, z_size)
  logical_world = new G4LogicalVolume(solid_world, world_material, "World"); //logical volume (solid volume, material, name)
  //physical volume
  physical_world = new G4PVPlacement(0, //No rotation
                                    G4ThreeVector(),//Null translation
                                    logical_world,  //logical volume
                                    "World",        //name
                                    0,              //null pointer to mother volume
                                    false,          //boolean operation
                                    0);              //copy number


  //--------------------------------Aerogel tile Geometry-------------------------------------------------------


   /* create a solid tile*/
   G4Material* tile_material = nist -> FindOrBuildMaterial("Aerogel"); 
   
   solid_tile = new G4Box("Tile", tileX, tileY, tileZ);  //simple tile volume (name, x_size, y_size, z_size)
   logical_tile = new G4LogicalVolume(solid_tile, tile_material, "Tile"); //logical tile volume (solid volume, material, name)
   
   // Physical volume is a placed instance of the logical volume in its mother logical volume
   physical_tile = new G4PVPlacement(0,  // Rotation
                                    G4ThreeVector(),   // its location
                                    logical_tile,      // the logical volume
                                    "Tile",            // its name
                                    logical_world,     // its mother volume 
                                    false,             // boolean operations
                                    0);                // its copy number

  //--------------------------------Detector Geometry-------------------------------------------------------
  //material: AIR
  G4Material* window_material = nist -> FindOrBuildMaterial("Air");
   
   solid_detector = new G4Box("Detector", windowX,windowY, windowZ); //simple tile volume (name, x_size, y_size, z_size)
   logical_detector = new G4LogicalVolume(solid_detector,window_material,"Detector"); //logical tile volume (solid volume, material, name)
   //physical volume
   physical_detector = new G4PVPlacement(0,                         // Rotation
                                        G4ThreeVector(0,0,20.*cm),  // its location
                                        logical_detector,           // the logical volume
                                        "Detector",                 // its name
                                        logical_world,              // its mother volume 
                                        false,                      // boolean operations
                                        0);                         // its copy number

  //set the sensitivity of the detector (uses TrackerSD.cpp)
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String sensitiveDetectorName = "/detector/sensitiveDetector";
  TrackerSD* theTrackerSD = new TrackerSD(sensitiveDetectorName, physical_detector); 
  SDman->AddNewDetector(theTrackerSD);
  logical_detector->SetSensitiveDetector(theTrackerSD);

   return physical_world; // Always return the world 
}
