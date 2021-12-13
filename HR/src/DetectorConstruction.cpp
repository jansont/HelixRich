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
#include "G4OpticalSurface.hh"
#include "PhysicsList.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4Isotope.hh"
#include <random>
#include "G4ParticleGun.hh"
#include "G4Cerenkov.hh"
#include <string> 

// constructor 
DetectorConstruction::DetectorConstruction(): G4VUserDetectorConstruction(), physical_detector(0)
{
  DefineMaterials(); // Construct materials for detector
}

// destructor
DetectorConstruction::~DetectorConstruction()
{}

//Photon momentum vector: photons created by cherenkov radiation have varying properties based on 
//their momentum. For example, constants for photon absorption will vary. Function below creates a
//bounded vector of photon momentums. 
std::vector<G4double> InitializePhotonMomentumVector()
{
  //Used to create rindex and abslength vectors 
  G4int numBins = 1000;
  G4double PhotonMinEnergy=SimulationConstants::PHOTON_MIN_ENERGY; 
  G4double PhotonMaxEnergy=SimulationConstants::PHOTON_MAX_ENERGY; 
  G4double step=(PhotonMaxEnergy-PhotonMinEnergy)/ numBins;
  std::vector<G4double>PhotonMomemtums(numBins);
  for (int i=0; i<numBins; i++)
  {
    PhotonMomemtums[i]=PhotonMinEnergy+step*i;
    }
  return PhotonMomemtums;
}


void DetectorConstruction::DefineMaterials()
{
  /*--------------------------------------Define Elements-----------------------------------*/
  G4String name, symbol;
  G4double density, a, z;
  G4int ncomponents, natoms;

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

  // water
  density = 1.00*g/cm3;
  G4Material* H2O = new G4Material(name="Water", density, ncomponents=2);
  H2O -> AddElement(H, natoms=2);
  H2O -> AddElement(O, natoms=1);
  //o2
  density = 1.43*g/L;
  G4Material* O2 = new G4Material(name="Oxygen2", density, ncomponents=1);
  O2 -> AddElement(O, natoms=2);
  // silica
  G4Material* SiO2 = new G4Material(name="Silica", SimulationConstants::SILICA_DENSITY, ncomponents=2);
  SiO2 -> AddElement(Si, natoms=1);
  SiO2 -> AddElement(O, natoms=2);

  //Atmospheric composition
  G4double pressure, temperature, atmospheric_density = 0.0;

  if (SimulationConstants::ALTITUDE ==  "sea")
  { //0m 
    pressure = 101325.00*pascal;
    temperature = 288.15*kelvin;
    atmospheric_density = 1.225*kg/m3;
  }
  else if (SimulationConstants::ALTITUDE == "troposphere")
  { //110000m
    pressure = 22632.10*pascal;
    temperature = 216.65*kelvin;
    atmospheric_density =0.4135*kg/m3;
  }
  else if (SimulationConstants::ALTITUDE == "stratosphere")
  { //20000m
    pressure = 5474.89*pascal;
    temperature = 216.65*kelvin;
    atmospheric_density = 0.08891*kg/m3;
  }
  else if (SimulationConstants::ALTITUDE == "mesosphere1")
  { //32000;
    pressure = 868.02*pascal;
    temperature = 228.65*kelvin;
    atmospheric_density = 0.01841*kg/m3;
  }
  else if (SimulationConstants::ALTITUDE == "mesosphere2")
  { //47000m
    pressure = 110.91*pascal;
    temperature = 270.65*kelvin;
    atmospheric_density = 0.001027*kg/m3;
  }
  else if (SimulationConstants::ALTITUDE == "thermosphere1")
  { //51000m
    pressure = 66.94*pascal;
    temperature = 270.65*kelvin;
    atmospheric_density = 0.001027*kg/m3;
  }
  else
  { //71000m
    pressure = 13.96*pascal;
    temperature = 214.65*kelvin;
    atmospheric_density = 0.00008283*kg/m3;
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

  // ----------------------------- Air ---------------------------------

  density = 1.290 * mg/cm3;
  G4Material* Air = new G4Material("Air", atmospheric_density, ncomponents = 4, kStateGas, temperature, pressure);
  Air -> AddMaterial(N2, 78.0840*perCent);
  Air -> AddMaterial(O2, 20.9460*perCent);
  Air -> AddMaterial(CO2, 0.0360*perCent);
  Air -> AddElement(Ar, 0.9340*perCent);


  //photon energy vector for world material
  G4double PhotonMinEnergy=SimulationConstants::PHOTON_MIN_ENERGY; 
  G4double PhotonMaxEnergy=SimulationConstants::PHOTON_MAX_ENERGY; 
  G4int numBins = 1000;
  G4double PhotonEnergyStep = (PhotonMaxEnergy-PhotonMinEnergy) / numBins;
  G4double* PhotonMomentum = new G4double[numBins];
  for (G4int i=0; i<numBins; i++)
  {
    PhotonMomentum[i] = PhotonMinEnergy + PhotonEnergyStep * i;
  }
  //Material properties of air
  G4double* AirAbsorpLength = new G4double[numBins];
  G4double* AirRindex = new G4double[numBins];
  for (G4int i=0; i<numBins; i++)
  {
    AirAbsorpLength[i]=1.E32*mm; //absoption 
    AirRindex[i]=1.000273;    // refractive index 
  }
  G4MaterialPropertiesTable* AirMPT = new G4MaterialPropertiesTable();
  Air->SetMaterialPropertiesTable(AirMPT);
  AirMPT->AddProperty("RINDEX", PhotonMomentum, AirRindex,numBins);
  Air->SetMaterialPropertiesTable(AirMPT);


  // ------------------------------------ Vaccuum ------------------------------------

  //Material properties of vaccuum
  G4double* VacRindex = new G4double[numBins];
  G4double* VacAbsorpLength = new G4double[numBins];
  for (G4int i=0; i<numBins; i++)
  {
    VacAbsorpLength[i]=1.E32*mm;
    VacRindex[i]=1.0;
  }
  //vacuum
  density = universe_mean_density; //from PhysicalConstants.h 
  pressure = 1.e-19*pascal;
  temperature = 0.1*kelvin;
  G4Material* vacuum = new G4Material(name="Galactic", z=1., a=1.e-8*g/mole, density,kStateGas,temperature,pressure);
  G4MaterialPropertiesTable* vacMPT = new G4MaterialPropertiesTable();
  vacMPT->AddProperty("RINDEX", PhotonMomentum, VacRindex,numBins);
  vacuum->SetMaterialPropertiesTable(vacMPT);

  

  /*---------------------------------Aerogel: Material and properties-----------------------------------*/  
  // Aerogel type
  Aerogel = new G4Material(name="Aerogel", SimulationConstants::AEROGEL_DENSITY, ncomponents=3, kStateGas, 293.*kelvin,  101325.*pascal); 
  Aerogel -> AddElement(Si, 1);
  Aerogel -> AddElement(Si, 4);
  Aerogel -> AddElement(Si, 4);


  //define index of refraction for aerogel
  G4MaterialPropertiesTable* AerogelMPT = new G4MaterialPropertiesTable(); //create a table for aerogel properties

  G4double* AerogelRindex = new G4double[numBins];        //refractive index 
  G4double* AerogelAbsorpLength = new G4double[numBins];  //Absorption
  G4double* AerogelRScatLength = new G4double[numBins]; //rayleigh scattering
  G4double* AerogelRaleighScatLength = new G4double[numBins];

  // ----------------------------Rayleigh scattering Length from RICH -----------------------------------
  // https://gitlab.cern.ch/geant4/geant4/tree/8a51e0bc408e73c60512e61fbc51c2be9c22904d/examples/advanced/Rich

  static const G4double PhotMomWaveConv=1243.125;
  static const G4double AerogelClarity =0.00719*micrometer*micrometer*micrometer*micrometer/cm;
  std::vector<G4double>AgelPhotW = InitializePhotonMomentumVector();
  G4double aClarity=AerogelClarity/(micrometer*micrometer*micrometer*micrometer);
  for(G4int i = 0; i<numBins; i ++ )
  {
    G4double ephoton=AgelPhotW[i]/eV;
    //In the following the 1000 is to convert form nm to micrometer
    G4double wphoton=(PhotMomWaveConv/ephoton)/1000.0;
    AerogelRaleighScatLength[i]=(std::pow(wphoton,4))/aClarity;
  }

  // --------------------------------Aerogel Refractive Index Configurations-------------------------------------------

  if (SimulationConstants::AEROGEL_PROPERTIES == "GaussianRindex")
  {
    //rindex is gaussian around 1.15
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<float> d{SimulationConstants::REFRACTIVE_INDEX, SimulationConstants::REFRACTIVE_INDEX_SDEV}; 
    for (G4int i=0; i<numBins; i++)
    {
      AerogelAbsorpLength[i]=1.E32*mm;
      AerogelRindex[i]= d(gen);
      AerogelRScatLength[i] = AerogelRaleighScatLength[i];
    }
  }
  else if (SimulationConstants::AEROGEL_PROPERTIES == "ExactRindex"){
  //rindex is constant 1.15 
    for (G4int i=0; i<numBins; i++)
    {
      AerogelAbsorpLength[i]=1.E32*mm;
      AerogelRindex[i]= SimulationConstants::REFRACTIVE_INDEX;
      AerogelRScatLength[i] = AerogelRaleighScatLength[i];
    }
  }

  double aero_n[]  = {1.155, 1.155};
  G4double energyRange[] = {1239.84193*eV/900, 1239.84193*eV/400};


  AerogelMPT->AddProperty("ABSLENGTH",PhotonMomentum, AerogelAbsorpLength,numBins)->SetSpline(true);
  AerogelMPT->AddProperty("RAYLEIGH",PhotonMomentum,AerogelRScatLength,numBins)->SetSpline(true);
  AerogelMPT->AddProperty("RINDEX", energyRange, aero_n, 2)->SetSpline(true);

  Aerogel->SetMaterialPropertiesTable(AerogelMPT);



  /*---------------------------------Detector Plane (photocathode): Material and properties-----------------------------------*/  


  G4Material* BialkaliCathode = new G4Material("photoCathode",
                         3*g/cm3,
                         3,
                         kStateSolid);
  BialkaliCathode->AddElement(Sb, 1);
  BialkaliCathode->AddElement(Rb, 1);
  BialkaliCathode->AddElement(Cs, 1);
  G4MaterialPropertiesTable* BialkaliCathodeProp = new G4MaterialPropertiesTable();
  const G4int num_bialkalicath = 2;
  G4double bialkalicathodeenergy[num_bialkalicath] = {1.2*eV, 6.5*eV};
  G4double bialkalicathodeabsorp[num_bialkalicath] = {1.e-6*mm, 1.e-6*mm}; // absorb all
  BialkaliCathodeProp->AddProperty("ABSLENGTH",
                  bialkalicathodeenergy,
                  bialkalicathodeabsorp,
                  num_bialkalicath);
  const  G4int nGlass = 4;
  G4double BoroSiGlassEnergy[nGlass]= {3.1*eV, 3.87*eV, 4.13*eV, 4.96*eV};
  G4double BoroSiGlassRindex[nGlass]= {1,1,1,1};
  BialkaliCathodeProp->AddProperty("RINDEX",
                   BoroSiGlassEnergy,
                   BoroSiGlassRindex,
                   nGlass ); // use values from window to prevent refraction
  BialkaliCathode->SetMaterialPropertiesTable(BialkaliCathodeProp);




 //--------------------------------Configuration of Aergoel Rindex Gradient-----------------------------------------

  if (SimulationConstants::ADD_TILE_GRADIENT && !(SimulationConstants::ADD_TILE))
  {
  
  //subdivide 100*100*10 mm^3 aerogel tile into 1mm^3 cubes
  // G4Material* Aerogel_cubes[100][100][10];

  G4int counter = 0;
  for (G4int x = 0; x < 100; x++){
    for (G4int t = 0; t < 100; t++){
      for (G4int u = 0; u < 10; u++){

          G4Material* Aerogel_cubes;
          std::string mat_name_str = std::to_string(counter);
          G4String mat_name = G4String(mat_name_str);
          // std::string mat_name = std::to_string(counter);


          Aerogel_cubes = new G4Material(name=mat_name, SimulationConstants::AEROGEL_DENSITY, ncomponents=3, kStateGas, 293.*kelvin,  101325.*pascal); 
          Aerogel_cubes -> AddElement(Si, 1);
          Aerogel_cubes -> AddElement(Si, 4);
          Aerogel_cubes -> AddElement(Si, 4);

          // rindex is gaussian around 1.155
          std::random_device rd{};
          std::mt19937 gen{rd()};
          float mean_index = SimulationConstants::REFRACTIVE_INDEX + (u / 1000);
          std::normal_distribution<float> d{mean_index, SimulationConstants::REFRACTIVE_INDEX_SDEV}; 
          
          for (G4int i=0; i<numBins; i++)
          {
            AerogelAbsorpLength[i]=1.E32*mm;
            AerogelRindex[i]= d(gen);
            AerogelRScatLength[i] = AerogelRaleighScatLength[i];
          }

          G4MaterialPropertiesTable* Aerogel_cubes_properties =  nullptr;

          Aerogel_cubes_properties->AddProperty("ABSLENGTH",PhotonMomentum, AerogelAbsorpLength,numBins)->SetSpline(true);
          Aerogel_cubes_properties->AddProperty("RAYLEIGH",PhotonMomentum,AerogelRScatLength,numBins)->SetSpline(true);
          Aerogel_cubes_properties->AddProperty("RINDEX", PhotonMomentum, AerogelRindex, 2)->SetSpline(true);

          Aerogel_cubes->SetMaterialPropertiesTable(Aerogel_cubes_properties);

          counter +=1 ;
        }
      }
    }
  }
}



// Now we will be creating the actual world and object inside the world here
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  G4NistManager* nist = G4NistManager::Instance();
  //--------------------------------World Geometry------------------------------------------------------------
  G4Material* world_material;
  if (SimulationConstants::WORLD_MATERIAL == "Air")
    world_material = nist -> FindOrBuildMaterial("Air"); //manager finds materials defined above 
  else if (SimulationConstants::WORLD_MATERIAL == "Galactic")
    world_material = nist -> FindOrBuildMaterial("Galactic");
  else{
    world_material = nist -> FindOrBuildMaterial("Galactic");
  }

  solid_world = new G4Box("World", 
              SimulationConstants::WORLDX,
              SimulationConstants::WORLDY,
              SimulationConstants::WORLDZ);   //simple volume (name, x_size, y_size, z_size)
  logical_world = new G4LogicalVolume(solid_world,
                    world_material,
                    "World"); //logical volume (solid volume, material, name)
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
  if (SimulationConstants::ADD_TILE)
  {
    if (SimulationConstants::CURVED)
    {
      solid_tile_curved = new G4Tubs("Tile",
      SimulationConstants::pRMin, //Inner radius
      SimulationConstants::pRMax, //Outer radius
      SimulationConstants::pDz,  //Half length in z
      SimulationConstants::pSPhi, //Starting phi angle in radians
      SimulationConstants::pDPhi); //Angle of the segment in radians

      G4double curved_edge = 100.*m;
      G4double curved_center = curved_edge - SimulationConstants::pRMin;
      logical_tile_curved = new G4LogicalVolume(solid_tile_curved,
                          tile_material,
                           "Tile"); //logical tile volume (solid volume, material, name)

      G4RotationMatrix* rotationMatrix_curved = new G4RotationMatrix();
        rotationMatrix_curved->rotateY(90.*deg);
      physical_tile_curved = new G4PVPlacement(rotationMatrix_curved,  // Rotation
                          G4ThreeVector(0,0,curved_center),   // its location
                          logical_tile_curved,      // the logical volume
                          "Tile",            // its name
                          logical_world,     // its mother volume 
                          false,             // boolean operations
                          0);                // its copy number
    }
    else
    {
      G4double tilex = SimulationConstants::TILEX/2;
      G4double tiley = SimulationConstants::TILEY/2;
      G4double tilez = SimulationConstants::TILEZ/2;
      solid_tile_box = new G4Box("Tile",
                    tilex,
                    tiley,
                    tilez);  //simple tile volume (name, x_size, y_size, z_size)

      logical_tile_box = new G4LogicalVolume(solid_tile_box,
                          tile_material,
                          "Tile"); //logical tile volume (solid volume, material, name)
      // Physical volume is a placed instance of the logical volume in its mother logical volume


      G4RotationMatrix* rotationMatrix = new G4RotationMatrix();
        rotationMatrix->rotateZ(90.*deg);

      physical_tile_box = new G4PVPlacement(rotationMatrix,  // Rotation
                        G4ThreeVector(0,0,100.*mm),   // its location
                        logical_tile_box,      // the logical volume
                        "Tile",            // its name
                        logical_world,     // its mother volume 
                        false,             // boolean operations
                        0);                // its copy number
    }

  }

  if (SimulationConstants::ADD_TILE_GRADIENT && !(SimulationConstants::ADD_TILE))
    {

    //--------------------------------Aergoel Rindex Gradient-----------------------------------------
      G4int x,t,u;

      G4double tilex = SimulationConstants::TILEX/(2*100);
      G4double tiley = SimulationConstants::TILEY/(2*100);
      G4double tilez = SimulationConstants::TILEZ/(2*10);
      solid_tile_box = new G4Box("Tile", tilex, tiley, tilez);  

    G4int count = 0;
    G4VPhysicalVolume* aerogel_cubes[100][100][10];

    for ( x = 0; x < 100; x++){
      for ( t = 0; t < 100; t++){
        for ( u = 0; u < 10; u++){


          std::string name_str = std::to_string(count);
          G4String name = G4String(name_str);          // std::string name = std::to_string(count);

          G4Material* cube_material = nist -> FindOrBuildMaterial(name); //mat_name

            G4RotationMatrix* rotationMatrix = new G4RotationMatrix();
            rotationMatrix->rotateZ(90.*deg);

            G4Box* solid_cube = new G4Box(name, tilex, tiley, tilez); 
            G4LogicalVolume* logical_cube = new G4LogicalVolume(solid_cube, cube_material, name); 

            G4double posZ = (95.5 + u);
            posZ = posZ*mm;
            G4double posXY = (-50 + s);
            posXY = posXY*mm;

            G4VPhysicalVolume* cube = new G4PVPlacement(rotationMatrix,  // Rotation
                              G4ThreeVector(posXY*mm,posXY*mm,posZ*mm),   // its location 
                              logical_cube,      // the logical volume
                              name,            // its name
                              logical_world,     // its mother volume 
                              false,             // boolean operations
                              0);                // its copy number

            aerogel_cubes[x][t][u] = cube;


          count+=1;

        }
      }
    }
  }

  //--------------------------------Detector Geometry-------------------------------------------------------
  solid_detector = new G4Box("detector",
               SimulationConstants::WINDOWX,
               SimulationConstants::WINDOWY,
               SimulationConstants::WINDOWZ); //simple tile volume (name, x_size, y_size, z_siz
  //material: AIR
  if (SimulationConstants::AIR)
  {
    G4Material* window_material = nist -> FindOrBuildMaterial("Air");
    logical_detector = new G4LogicalVolume(solid_detector,
                    window_material,
                    "detector"); //logical tile volume (solid volume, material, name)
    physical_detector = new G4PVPlacement(0,                         // Rotation
                      G4ThreeVector(0,0,383*mm),  // its location
                      logical_detector,           // the logical volume
                      "detector",                 // its name
                      logical_world,              // its mother volume 
                      false,                      // boolean operations
                      0);                         // its copy number
  } 
  else
  {
    G4Material* window_material2 = nist -> FindOrBuildMaterial("photoCathode");
    logical_detector = new G4LogicalVolume(solid_detector,
                        window_material2,
                        "detector", 
                        0,0,0);
    physical_detector = new G4PVPlacement(0,  
                      G4ThreeVector(0,0,383*mm),
                      logical_detector,
                      "detector",
                      logical_world,
                      false,
                       0);
  }

  //set the sensitivity of the detector (TrackerSD.cpp)
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String sensitiveDetectorName = "/detector/sensitiveDetector";
  TrackerSD* theTrackerSD = new TrackerSD(sensitiveDetectorName, physical_detector); 
  SDman->AddNewDetector(theTrackerSD);
  logical_detector->SetSensitiveDetector(theTrackerSD);

  
  // stop compiler freaking out about unsed variable
  (void)physical_detector;

  /*------------------------------------------------Aerogel Surface----------------------------------------------------------------------*/
  int numBins = 1000;
  int i;
  //photon energy vector for world material
  G4double PhotonMinEnergy=SimulationConstants::PHOTON_MIN_ENERGY; 
  G4double PhotonMaxEnergy=SimulationConstants::PHOTON_MAX_ENERGY; 
  G4double PhotonEnergyStep = (PhotonMaxEnergy-PhotonMinEnergy) / numBins;
  G4double* PhotonMomentum = new G4double[numBins];
  for (i=0; i<numBins; i++)
  {
    PhotonMomentum[i] = PhotonMinEnergy + PhotonEnergyStep * i;
  }

  G4double* AerogelRindex = new G4double[numBins];        
  G4double* AerogelAbsorpLength = new G4double[numBins];
  G4double* AerogelRScatLength = new G4double[numBins];
  G4double* AerogelRaleighScatLength = new G4double[numBins];

  static const G4double PhotMomWaveConv=1243.125;
  static const G4double AerogelClarity =0.00719*micrometer*micrometer*micrometer*micrometer/cm;
  std::vector<G4double>AgelPhotW = InitializePhotonMomentumVector();
  G4double aClarity=AerogelClarity/(micrometer*micrometer*micrometer*micrometer);

  for(i=0; i<numBins; i++ )
  {
    G4double ephoton=AgelPhotW[i]/eV;
    //In the following the 1000 is to convert form nm to micrometer
    G4double wphoton=(PhotMomWaveConv/ephoton)/1000.0;
    AerogelRaleighScatLength[i]=(std::pow(wphoton,4))/aClarity;
  }

  if (SimulationConstants::AEROGEL_PROPERTIES == "GaussianRindex")
  {
    //rindex is gaussian around 1.15
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<float> d{SimulationConstants::REFRACTIVE_INDEX, SimulationConstants::REFRACTIVE_INDEX_SDEV}; 
    for (i=0; i<numBins; i++)
    {
      AerogelAbsorpLength[i]=1.E32*mm;
      AerogelRindex[i]= d(gen);
      AerogelRScatLength[i] = AerogelRaleighScatLength[i];
    }
  }
  else if (SimulationConstants::AEROGEL_PROPERTIES == "ExactRindex")
  {
    //rindex is constant 1.15 
    for (G4int j=0; j<numBins; j++)
    {
      AerogelAbsorpLength[j]=1.E32*mm;
      AerogelRindex[j]= 1.155;
      AerogelRScatLength[j] = AerogelRaleighScatLength[j];
    }
  }

  // Surface roughness defined using unified model
  if (SimulationConstants::SURFACE_MODEL == "unified")
  {
    G4MaterialPropertiesTable* AerogelSurfaceMPT = new G4MaterialPropertiesTable();
    G4OpticalSurface* AerogelSurface = new G4OpticalSurface("AerogelSurface");
    G4LogicalBorderSurface* AerogelSurface_logical = new G4LogicalBorderSurface("AerogelSurface",
                                          physical_tile_box,
                                          physical_world,
                                          AerogelSurface);
    AerogelSurface->SetType(dielectric_dielectric);
    AerogelSurface->SetModel(unified);
    AerogelSurface->SetFinish(ground);
    AerogelSurface->SetSigmaAlpha(SimulationConstants::AEROGEL_ROUGHNESS);

    std::vector<G4double> pp = {1.*eV, 10.*eV};
    std::vector<G4double> specularlobe = {0.0, 0.0};
    std::vector<G4double> specularspike = {0.0, 0.0};
    std::vector<G4double> backscatter = {0.0, 0.0};
    std::vector<G4double> reflectivity = {0.0, 0.0};
    std::vector<G4double> efficiency = {0.0, 0.0};

    AerogelSurfaceMPT->AddProperty("RINDEX", PhotonMomentum, AerogelRindex, numBins);
    AerogelSurfaceMPT->AddProperty("SPECULARLOBECONSTANT", pp, specularlobe);
    AerogelSurfaceMPT->AddProperty("SPECULARSPIKECONSTANT", pp, specularspike);
    AerogelSurfaceMPT->AddProperty("BACKSCATTERCONSTANT", pp, backscatter);
    AerogelSurfaceMPT->AddProperty("EFFICIENCY", pp, efficiency);
    AerogelSurface->SetMaterialPropertiesTable(AerogelSurfaceMPT);  
  }

  // Surface roughness defined using glisur model
  if (SimulationConstants::SURFACE_MODEL == "glisur")
  {
    std::vector<G4double> pp = {1.*eV, 10.*eV};
    std::vector<G4double> specularlobe = {0.0, 0.0};
    std::vector<G4double> specularspike = {0.0, 0.0};
    std::vector<G4double> backscatter = {0.0, 0.0};
    std::vector<G4double> efficiency = {0.0, 0.0};                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  

    G4OpticalSurface* OpSurface = new G4OpticalSurface("AerogelSurface");
    G4LogicalSkinSurface* Surface = new G4LogicalSkinSurface("AerogelSurface",
                              logical_tile_box,
                              OpSurface);
    OpSurface->SetType(dielectric_dielectric);
    OpSurface->SetFinish(ground);
    OpSurface->SetModel(glisur);
    OpSurface->SetPolish(SimulationConstants::AEROGEL_ROUGHNESS);

    G4MaterialPropertiesTable* OpSurfaceProperty = new G4MaterialPropertiesTable();
    OpSurfaceProperty->AddProperty("RINDEX", PhotonMomentum, AerogelRindex, numBins);
    OpSurfaceProperty->AddProperty("EFFICIENCY", pp, efficiency);
    OpSurface->SetMaterialPropertiesTable(OpSurfaceProperty);
  }
  // Surface roughness defined using LUT Davis model
  else if (SimulationConstants::SURFACE_MODEL == "davis")
  {
    G4OpticalSurface* OpSurface = new G4OpticalSurface("AerogelSurface");
    OpSurface->SetType(dielectric_LUTDAVIS);
    if (SimulationConstants::DAVIS_ROUGHNESS == "rough")
      OpSurface->SetFinish(Rough_LUT);
    else if (SimulationConstants::DAVIS_ROUGHNESS == "polished")
      OpSurface->SetFinish(Polished_LUT);
    OpSurface->SetModel(DAVIS);
    G4LogicalBorderSurface* aerogelSurface = new G4LogicalBorderSurface("AerogelSurface",
                                    physical_tile_box,
                                  physical_world,
                                    OpSurface);
  }
  else //for surface model = none, do  nothing
  {}
  return physical_world; // Always return the world 
}