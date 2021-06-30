#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace SimulationConstants{
 /*------------------------------Detector geometry------------------------------*/
    //World volume
    extern const G4double worldX;
    extern const G4double worldY;
    extern const G4double worldZ;
    extern const bool curved;
    //Tile Volume: box
    extern const G4double tileX; 
    extern const G4double tileY;
    extern const G4double tileZ; 
    //Tile Volume: curved 
    extern const G4double pRMin; //Inner radius
    extern const G4double pRMax; //Outer radius
    extern const G4double pDz; //Half length in z
    extern const G4double pSPhi; //Starting phi angle in radians
    extern const G4double pDPhi; //Angle of the segment in radians 
    //Detector Volume
    extern const G4double windowX; 
    extern const G4double windowY;
    extern const G4double windowZ; 

    /*-------------------------------Detector Materials------------------------------*/
    extern const G4double aerogel_density;
    extern const G4double silica_prop;
    extern const G4double water_prop;
    extern const int Aerogel_properties;
    extern const int world_material_type;
    extern const G4double silica_density;
    extern const int altitude;

    /*-------------------------------Sensitive Detector------------------------------*/
    extern const bool air;
    extern const bool electron;

    /*-------------------------------Particle Generator------------------------------*/
    extern const G4int n_particle; 
    extern const G4double electronEnergy;
    extern const int distribution;
    extern const double mean;
    extern const double std_dev;

    /*-------------------------------Physics Processes------------------------------*/
    //Electromagnetic processes
    extern const bool multiple_scattering;
    extern const bool ionization;
    extern const bool bremsstrahlung;
    //Optical processes 
    extern const bool cerenkov;
    extern const bool rayleigh_scattering;
    extern const bool absorption;
    extern const G4int MaxNumPhotons;
}
