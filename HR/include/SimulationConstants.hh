#ifndef SimulationConstants_h
#define SimulationConstants_h 1

/*-------------------------------Detector geometry------------------------------*/
//World volume
static const G4double worldX = 2.*m;
static const G4double worldY = 2.*m;
static const G4double worldZ = 2.*m;
//Tile Volume
static const G4double tileX = 15*cm; 
static const G4double tileY = 15*cm;
static const G4double tileZ = 1.*cm;  
//Detector Volume
static const G4double windowX = 400*mm; 
static const G4double windowY = 400*mm;
static const G4double windowZ = 1.*mm; 

/*-------------------------------Detector Materials------------------------------*/
static const int Aerogel_type = 1;
static const int Aerogel_properties = 1;
static const int world_material_type = 1;

/*-------------------------------Physics Materials------------------------------*/
//Electromagnetic processes
static const bool multiple_scattering = true;
static const bool ionization = true;
static const bool bremsstrahlung = true;
//Optical processes 
static bool bool cerenkov = true;
static bool bool rayleigh_scattering = true;
G4int MaxNumPhotons = 100;

#endif
