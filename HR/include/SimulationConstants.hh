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
    extern const G4double silica_density;
    extern const G4double silica_prop;
    extern const G4double water_prop;
    extern const float refractive_index;
    extern const std::string Aerogel_properties;
    extern const float rindex_sdev;
    extern const G4double PhotonMinEnergy_;
    extern const G4double PhotonMaxEnergy_;
    extern const std::string world_material_type;
    extern const std::string altitude;

    /*-------------------------------Sensitive Detector------------------------------*/
    extern const bool air;
    extern const bool electron;

    /*-------------------------------Particle Generator------------------------------*/
    extern const std::string sourceParticle; 
    extern const G4int n_particle; 
    extern const G4double particleEnergy;
    extern const std::string energyDistribution;
    extern const G4double uniformEnergyRadius;  
    extern const float E_sdev;
    extern const std::string momentumDistribution;
    extern const double mean_momentum_z;
    extern const double sdev_momentum_z;
    extern const double mean_momentum_x;
    extern const double sdev_momentum_x;
    extern const double mean_momentum_y;
    extern const double sdev_momentum_y;
    extern const std::string source;
    extern const G4double uniform_source_center;
    extern const G4double uniform_source_radius;
    extern const double mean_source_x;
    extern const double sdev_source_x;
    extern const double mean_source_y;
    extern const double sdev_source_y;

    /*-------------------------------Physics Processes------------------------------*/
    //Electromagnetic processes
    extern const bool multiple_scattering;
    extern const bool ionization;
    extern const bool bremsstrahlung;
    extern const bool rayleigh;
    extern const bool photoElectric;
    extern const bool compton;
    //Optical processes 
    extern const bool cerenkov;
    extern const bool rayleigh_scattering;
    extern const bool absorption;
    extern const bool mie;
    extern const G4int MaxNumPhotons;
    extern const G4int maxBetaChange;
}
