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

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace SimulationConstants{
 /*------------------------------Detector geometry------------------------------*/
    //Global
    extern const G4int runCount;
    extern const bool ShowTerminalUI;
    extern const bool TRACKING;
    //World volume
    extern const G4double WORLDX;
    extern const G4double WORLDY;
    extern const G4double WORLDZ;
    extern const bool CURVED;
    //Tile Volume: box
    extern const G4double TILEX; 
    extern const G4double TILEY;
    extern const G4double TILEZ; 
    //Tile Volume: curved 
    extern const G4double pRMin; //Inner radius
    extern const G4double pRMax; //Outer radius
    extern const G4double pDz; //Half length in z
    extern const G4double pSPhi; //Starting phi angle in radians
    extern const G4double pDPhi; //Angle of the segment in radians 
    //Detector Volume
    extern const G4double WINDOWX; 
    extern const G4double WINDOWY;
    extern const G4double WINDOWZ; 

    /*-------------------------------Detector Materials------------------------------*/
    extern const bool ADD_TILE;
    extern const bool ADD_TILE_GRADIENT;
    extern const G4double AEROGEL_DENSITY;
    extern const G4double SILICA_DENSITY;
    extern const G4double SILICA_PROP;
    extern const std::string SURFACE_MODEL;
    extern const G4double AEROGEL_ROUGHNESS;
    extern const std::string DAVIS_ROUGHNESS;
    extern const float REFRACTIVE_INDEX;
    extern const std::string AEROGEL_PROPERTIES;
    extern const float REFRACTIVE_INDEX_SDEV;
    extern const G4double PHOTON_MIN_ENERGY;
    extern const G4double PHOTON_MAX_ENERGY;
    extern const std::string WORLD_MATERIAL;
    extern const std::string ALTITUDE;

    /*-------------------------------Sensitive Detector------------------------------*/
    extern const bool AIR;
    extern const std::string PARTICLE_TO_DETECT;

    /*-------------------------------Particle Generator------------------------------*/
    extern const std::string SOURCE_PARTICLE; 
    extern const G4int N_PARTICLE; 
    extern const G4double PARTICLE_ENERGY;
    extern const std::string ENERGY_DISTRIBUTION;
    extern const G4double UNIFORM_ENERGY_RADIUS;  
    extern const float ENERGY_SDEV;
    extern const std::string MOMENTUM_DISTRIBUTION;
    extern const double MEAN_MOMENTUM_X;
    extern const double SDEV_MOMENTUM_X;
    extern const double MEAN_MOMENTUM_Y;
    extern const double SDEV_MOMENTUM_Y;
    extern const double MEAN_THETA;
    extern const double SDEV_THETA;
    extern const std::string SOURCE_LOCATION;
    extern const G4double UNIFORM_SOURCE_CENTRE;
    extern const G4double UNIFORM_SOURCE_RADIUS;
    extern const double MEAN_SOURCE_X;
    extern const double SDEV_SOURCE_X;
    extern const double MEAN_SOURCE_Y;
    extern const double SDEV_SOURCE_Y;

    /*-------------------------------Physics Processes------------------------------*/
    //Electromagnetic processes
    extern const bool MULTIPLE_SCATTERING;
    extern const bool IONIZATION_PROCESS;
    extern const bool BREMSSTRAHLUNG_PROCESS;

    //Optical processes 
    extern const bool RAYLEIGH_SCATTERING;
    extern const bool ABSORPTION_PROCESS;
    extern const bool BOUNDARY_PROCESS;
    extern const bool COMPTON_SCATTERING;
    extern const bool SCINTILLATION_PROCESS;
    extern const G4int MAX_NUM_PHOTONS;
    extern const G4int MAX_BETA_CHANGE;
}
