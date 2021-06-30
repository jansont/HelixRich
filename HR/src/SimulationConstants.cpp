#include "SimulationConstants.hh"
#include <math.h>

namespace SimulationConstants {
    /*------------------------------Detector geometry------------------------------*/
    //World volume
    const G4double worldX = 2.*m;
    const G4double worldY = 2.*m;
    const G4double worldZ = 2.*m;
    const bool curved = false;
    //Tile Volume: box
    const G4double tileX = 15*cm; 
    const G4double tileY = 15*cm;
    const G4double tileZ = 1.*cm; 
    //Tile Volume: curved 
    const G4double pRMin = 20*cm; //Inner radius
    const G4double pRMax = 16.*cm; //Outer radius
    const G4double pDz = 15.*cm; //Half length in z
    const G4double pSPhi = 0.; //Starting phi angle in radians
    const G4double pDPhi = M_PI/2; //Angle of the segment in radians 
    //Detector Volume
    const G4double windowX = 500*mm; 
    const G4double windowY = 500*mm;
    const G4double windowZ = 1.*mm; 

    /*-------------------------------Detector Materials------------------------------*/
    const G4double aerogel_density = 0.200*g/cm3;
    const G4double silica_prop = 97*perCent;
    const G4double water_prop =  3*perCent;
    const int Aerogel_properties = 1;
    const int world_material_type = 1;
    const G4double silica_density = 2.200*g/cm3;
    const int altitude = 0;


    /*-------------------------------Sensitive Detector------------------------------*/
    const bool air = false;
    const bool electron = false;

    /*-------------------------------Particle Generator------------------------------*/
    const G4int n_particle = 1; 
    const G4double electronEnergy = 35.0 * MeV;
    const int distribution = 1;
    const double mean = 0.0;
    const double std_dev = 30; //degrees (opening angle)

    /*-------------------------------Physics Processes------------------------------*/
    //Electromagnetic processes
    const bool multiple_scattering = true;
    const bool ionization = true;
    const bool bremsstrahlung = false;
    //Optical processes 
    const bool cerenkov = true;
    const bool rayleigh_scattering = true;
    const bool absorption = true;
    const G4int MaxNumPhotons = 1000;
};