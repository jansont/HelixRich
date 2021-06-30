/******* We will create the most basics of a detector construction file  ********/
/*** 
 * To do this we will create a world box and place another box inside it and fill it. We need 
 * to create the logical volume, use the solid, and add other attributes (like composition)
 * ***/

#include "DetectorConstruction.hh"
#include "G4NistManager.hh"         // Contains the NIST data for particles
#include "G4Box.hh"                 // Needed to create box shapes
#include "G4Tubs.hh"
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
#include <random>

// constructor 
DetectorConstruction::DetectorConstruction()  : G4VUserDetectorConstruction(), physical_detector(0)
{
   DefineMaterials(); // Construct materials for detector
}

// destructor
DetectorConstruction::~DetectorConstruction()
{}  // Don't need to destruct anything, but we keep good habits


std::vector<G4double> InitializePhotonMomentumVector() {
  G4int ibin;
  //PhotonEnergy
  G4int NumPhotWaveLengthBins = 1000;
  G4double PhotonMinEnergy=1.3*CLHEP::eV; //sim constant
  G4double PhotonMaxEnergy=7.3*CLHEP::eV; //sim constant
  G4double PhotonEnergyStep=(PhotonMaxEnergy-PhotonMinEnergy)/ NumPhotWaveLengthBins;
  
  std::vector<G4double>PhotMomVect(NumPhotWaveLengthBins);
  for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
    PhotMomVect[ibin]=PhotonMinEnergy+PhotonEnergyStep*ibin;
  }
  return PhotMomVect;
}

std::vector<G4double> InitAgelPhotW() {
  return InitializePhotonMomentumVector() ; 
}

void DetectorConstruction::DefineMaterials()
{
  /*--------------------------------------Define Elements-----------------------------------*/
  G4String name, symbol;
  G4NistManager* nist = G4NistManager::Instance();
  G4Element *N = nist->FindOrBuildElement("N");
  G4Element *C = nist->FindOrBuildElement("C");
  G4Element *H = nist->FindOrBuildElement("H");
  G4Element *Si = nist->FindOrBuildElement("Si");
  G4Element *O = nist->FindOrBuildElement("O");
  G4Element *Ar = nist->FindOrBuildElement("Ar");
  G4Element *Sb = nist->FindOrBuildElement("Sb");
  G4Element *Rb = nist->FindOrBuildElement("Rb");
  G4Element *Cs = nist->FindOrBuildElement("Cs");

/*--------------------------------------Define Materials-----------------------------------*/
  // secondary materials: 
  G4int ncomponents, natoms;
  G4double density, fractionmass, a, z;
  
  // water
  density = 1.0*g/cm3;
  G4Material* H2O = new G4Material(name="Water", density, ncomponents=2);
  H2O -> AddElement(H, natoms=2);
  H2O -> AddElement(O, natoms=1);
    //o2
  density = 1.43*g/L;
  G4Material* O2 = new G4Material(name="Oxygen2", density, ncomponents=1);
  O2 -> AddElement(O, natoms=2);
  // silica
  G4Material* SiO2 = new G4Material(name="Silica", SimulationConstants::silica_density, ncomponents=2);
  SiO2 -> AddElement(Si, natoms=1);
  SiO2 -> AddElement(O, natoms=2);

  //Atmospheric composition
  G4double pressure;
  G4double temperature;
  G4double atmospheric_density;
  switch(SimulationConstants::altitude){
    case 0: //0m 
      pressure = 101325.00*pascal;
      temperature = 288.15*kelvin;
      atmospheric_density = 1.225*kg/m3;
      break;
    case 1: //110000m
      pressure = 22632.10*pascal;
      temperature = 216.65*kelvin;
      atmospheric_density =0.4135*kg/m3;
      break;
    case 2: //20000m
      pressure = 5474.89*pascal;
      temperature = 216.65*kelvin;
      atmospheric_density = 0.08891*kg/m3;
      break;
    case 3: //32000;
      pressure = 868.02*pascal;
      temperature = 228.65*kelvin;
      atmospheric_density = 0.01841*kg/m3;
      break;
    case 4: //47000m
      pressure = 110.91*pascal;
      temperature = 270.65*kelvin;
      atmospheric_density = 0.001027*kg/m3;
      break;
    case 5: //51000m
      pressure = 66.94*pascal;
      temperature = 270.65*kelvin;
      atmospheric_density = 0.001027*kg/m3;
      break;
    case 6: //71000m
      pressure = 13.96*pascal;
      temperature = 214.65*kelvin;
      atmospheric_density = 0.00008283*kg/m3;
      break;
  }

  //n2
  density = 1.25*g/L;
  G4Material* N2 = new G4Material(name="NitrogenGas", density, ncomponents=1, kStateGas,temperature,pressure);
  N2 -> AddElement(N,natoms=2);
  //co2
  density = 1.96*g/L;
  G4Material* CO2 = new G4Material(name="CO2", density, ncomponents=2, kStateGas,temperature,pressure);
  CO2 -> AddElement(O, natoms=2);
  CO2 -> AddElement(C, natoms=1);
  // air 
  density = 1.290 * mg/cm3;
  G4Material* Air = new G4Material("Air", atmospheric_density, ncomponents = 4, kStateGas, temperature, pressure);
  Air -> AddMaterial(N2, 78.0840*perCent);
  Air -> AddMaterial(O2, 20.9460*perCent);
  Air -> AddMaterial(CO2, 0.0360*perCent);
  Air -> AddElement(Ar, 0.9340*perCent);

  G4int ibin;
  //PhotonEnergy
  G4double PhotonMinEnergy=1.3*CLHEP::eV; //sim constant
  G4double PhotonMaxEnergy=7.3*CLHEP::eV; //sim constant
  G4int NumPhotWaveLengthBins = 1000;
  G4double PhotonEnergyStep = (PhotonMaxEnergy-PhotonMinEnergy) / NumPhotWaveLengthBins;
  G4double* PhotonMomentum = new G4double[NumPhotWaveLengthBins];
  for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
    PhotonMomentum[ibin] = PhotonMinEnergy + PhotonEnergyStep * ibin;
  }

  //vacuum
  density = universe_mean_density; //from PhysicalConstants.h pressure = 1.e-19*pascal;
  temperature = 0.1*kelvin;
  G4Material* vacuum = new G4Material(name="Galactic", z=1., a=1.01*g/mole, density,kStateGas,temperature,pressure);

  //Materila properties of air
  G4double* AirAbsorpLength = new G4double[NumPhotWaveLengthBins];
  G4double* AirRindex = new G4double[NumPhotWaveLengthBins];

  for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
    AirAbsorpLength[ibin]=1.E32*mm;
    AirRindex[ibin]=1.000273;
  }

  G4MaterialPropertiesTable* AirMPT = new G4MaterialPropertiesTable();
  Air->SetMaterialPropertiesTable(AirMPT);
  AirMPT->AddProperty("RINDEX", PhotonMomentum, AirRindex,NumPhotWaveLengthBins);
  Air->SetMaterialPropertiesTable(AirMPT);


/*---------------------------------Aerogel: Material and properties-----------------------------------*/
//density, state, temperature, pressure,
  
  // Aerogel type
  G4Material* Aerogel = new G4Material(name="Aerogel", SimulationConstants::aerogel_density, ncomponents=2); 
  Aerogel -> AddMaterial(SiO2, SimulationConstants::silica_prop);
  Aerogel -> AddMaterial(H2O , SimulationConstants::water_prop);

  //define index of refraction for aerogel
  switch (SimulationConstants::Aerogel_properties){
    case 1:

      G4MaterialPropertiesTable* AerogelMPT = new G4MaterialPropertiesTable(); //create a table for aerogel properties

      G4double* AerogelRindex = new G4double[NumPhotWaveLengthBins];        
      G4double* AerogelAbsorpLength = new G4double[NumPhotWaveLengthBins];
      G4double* AerogelRScatLength = new G4double[NumPhotWaveLengthBins];
      G4double* AgelRaleighScatLength = new G4double[NumPhotWaveLengthBins];
    
      std::random_device rd{};
      std::mt19937 gen{rd()};
      float std_dev = 0.0007; //sim constant
      float mean = 1.15;       // sim constant
      std::normal_distribution<float> d{mean, std_dev}; 


      static const G4double PhotMomWaveConv=1243.125;
      static const G4double AerogelClarity =0.00719*micrometer*micrometer*micrometer*micrometer/cm;
      std::vector<G4double>AgelPhotW = InitAgelPhotW();
      G4double aClarity=AerogelClarity/(micrometer*micrometer*micrometer*micrometer);
      for(G4int ibinw=0; ibinw<NumPhotWaveLengthBins; ibinw++ ){
        G4double ephoton=AgelPhotW[ibinw]/eV;
        //In the following the 1000 is to convert form nm to micrometer
        G4double wphoton=(PhotMomWaveConv/ephoton)/1000.0;
        AgelRaleighScatLength[ibinw]=(std::pow(wphoton,4))/aClarity;
      }

      for (ibin=0; ibin<NumPhotWaveLengthBins; ibin++){
          AerogelAbsorpLength[ibin]=1.E32*mm;
          AerogelRindex[ibin]= d(gen);
          AerogelRScatLength[ibin] = AgelRaleighScatLength[ibin];
        }
      AerogelMPT->AddProperty("ABSLENGTH",PhotonMomentum, AerogelAbsorpLength,NumPhotWaveLengthBins);
      AerogelMPT->AddProperty("RAYLEIGH",PhotonMomentum,AerogelRScatLength,NumPhotWaveLengthBins);
      AerogelMPT->AddProperty("RINDEX", PhotonMomentum, AerogelRindex, NumPhotWaveLengthBins);
      Aerogel->SetMaterialPropertiesTable(AerogelMPT);
      break;
  }

  G4Material* BialkaliCathode = new G4Material("photoCathode", 3*g/cm3, 3, kStateSolid);
  BialkaliCathode->AddElement(Sb, 1);
  BialkaliCathode->AddElement(Rb, 1);
  BialkaliCathode->AddElement(Cs, 1);
  G4MaterialPropertiesTable* BialkaliCathodeProp = new G4MaterialPropertiesTable();
  const G4int numentriesbialkalicath = 2;
  G4double bialkalicathodeenergy[numentriesbialkalicath] = {1.2*eV, 6.5*eV};
  G4double bialkalicathodeabsorp[numentriesbialkalicath] = {1.e-6*mm, 1.e-6*mm}; // absorb all
  BialkaliCathodeProp->AddProperty("ABSLENGTH", bialkalicathodeenergy, bialkalicathodeabsorp, numentriesbialkalicath);
   const  G4int nGlass = 4;
  G4double BoroSiGlassEnergy[nGlass]= {3.1*eV, 3.87*eV, 4.13*eV, 4.96*eV} ;// {1.2*eV, 3.1*eV, 6.5*eV};
  G4double BoroSiGlassRindex[nGlass]= {1.15, 1.15, 1.15, 1.15};
  BialkaliCathodeProp->AddProperty("RINDEX", BoroSiGlassEnergy,BoroSiGlassRindex ,nGlass ); // use values from window to prevent refraction
  BialkaliCathode->SetMaterialPropertiesTable(BialkaliCathodeProp);
}


// Now we will be creating the actual world and object inside the world here
G4VPhysicalVolume* DetectorConstruction::Construct()
{
   G4NistManager* nist = G4NistManager::Instance();
  //--------------------------------World Geometry------------------------------------------------------------
  G4Material* world_material;
  switch (SimulationConstants::world_material_type){
    case 1:
      world_material = nist -> FindOrBuildMaterial("Air"); //manager finds materials defined above 
      break;
    case 2: 
      world_material = nist -> FindOrBuildMaterial("Galactic");
  }

  solid_world = new G4Box("World", SimulationConstants::worldX, SimulationConstants::worldY, SimulationConstants::worldZ);   //simple volume (name, x_size, y_size, z_size)
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
  if (SimulationConstants::curved){
    solid_tile_curved = new G4Tubs("Tile",
                                    SimulationConstants::pRMin, //Inner radius
                                    SimulationConstants::pRMax, //Outer radius
                                    SimulationConstants::pDz,  //Half length in z
                                    SimulationConstants::pSPhi, //Starting phi angle in radians
                                    SimulationConstants::pDPhi); //Angle of the segment in radians
    logical_tile_curved = new G4LogicalVolume(solid_tile_curved, tile_material, "Tile"); //logical tile volume (solid volume, material, name)
    physical_tile_curved = new G4PVPlacement(0,  // Rotation
                                      G4ThreeVector(),   // its location
                                      logical_tile_curved,      // the logical volume
                                      "Tile",            // its name
                                      logical_world,     // its mother volume 
                                      false,             // boolean operations
                                      0);                // its copy number
  }else{
    solid_tile_box = new G4Box("Tile", SimulationConstants::tileX, SimulationConstants::tileY, SimulationConstants::tileZ);  //simple tile volume (name, x_size, y_size, z_size)
    logical_tile_box = new G4LogicalVolume(solid_tile_box, tile_material, "Tile"); //logical tile volume (solid volume, material, name)
    // Physical volume is a placed instance of the logical volume in its mother logical volume
    physical_tile_box = new G4PVPlacement(0,  // Rotation
                                    G4ThreeVector(),   // its location
                                    logical_tile_box,      // the logical volume
                                    "Tile",            // its name
                                    logical_world,     // its mother volume 
                                    false,             // boolean operations
                                    0);                // its copy number
  }

    //--------------------------------Detector Geometry-------------------------------------------------------
  solid_detector = new G4Box("detector", SimulationConstants::windowX,SimulationConstants::windowY, SimulationConstants::windowZ); //simple tile volume (name, x_size, y_size, z_siz

  //material: AIR
  if (SimulationConstants::air){
  G4Material* window_material = nist -> FindOrBuildMaterial("Air");
  logical_detector = new G4LogicalVolume(solid_detector,window_material,"Detector"); //logical tile volume (solid volume, material, name)
  physical_detector = new G4PVPlacement(0,                         // Rotation
                                        G4ThreeVector(0,0,20.*cm),  // its location
                                        logical_detector,           // the logical volume
                                        "Detector",                 // its name
                                        logical_world,              // its mother volume 
                                        false,                      // boolean operations
                                        0);                         // its copy number
  } else {
  G4Material* window_material2 = nist -> FindOrBuildMaterial("photoCathode");
  logical_detector = new G4LogicalVolume(solid_detector,window_material2, "detector", 0,0,0);
  physical_detector = new G4PVPlacement(0, G4ThreeVector(0,0,20.*cm), logical_detector, "detector", logical_world, false, 0);
  }

  //set the sensitivity of the detector (uses TrackerSD.cpp)
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String sensitiveDetectorName = "/detector/sensitiveDetector";
  TrackerSD* theTrackerSD = new TrackerSD(sensitiveDetectorName, physical_detector); 
  SDman->AddNewDetector(theTrackerSD);
  logical_detector->SetSensitiveDetector(theTrackerSD);

   return physical_world; // Always return the world 
}


