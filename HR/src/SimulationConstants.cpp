#include "SimulationConstants.hh"
#include <math.h>

namespace SimulationConstants {
	/*------------------------------Detector geometry------------------------------*/
	//Global 
	const G4int runCount = 10000;                            //$ Run Count>
	const bool ShowTerminalUI = false;						//$ Show in terminal>
	const bool TRACKING = false;							//$ Electron and Photon Tracking>
	//World volume
	const G4double WORLDX = 2.*m;                         //$ World size x>
	const G4double WORLDY = 2.*m;                         //$ World size y>
	const G4double WORLDZ = 2.*m;                         //$ World size z>
	const bool CURVED = false;                            //$ Curved tile geometry>
	//Tile Volume: box
	const G4double TILEX = 100*mm;                         //$ Tile size X>
	const G4double TILEY = 100*mm;                         //$ Tile size Y>
	const G4double TILEZ = 10*mm;                        //$ Tile size Z>
	//Tile Volume: curved 
	const G4double pRMin = 4*cm;                         //$ Curved tile inner radius>
	const G4double pRMax = 5*cm;                        //$ Curved tile outer radius>
	const G4double pDz = 15.*cm;                          //$ Curved tile half length in z>
	const G4double pSPhi = 0;                             //$ Curved tile starting phi angle (rad)>
	const G4double pDPhi = M_PI;                          //$ Curved tile angle of the segment (rad)>
	//Detector Volume
	const G4double WINDOWX = 500*mm;                      //$ Sensitive detector size X>
	const G4double WINDOWY = 500*mm;                      //$ Sensitive detector size Y>
	const G4double WINDOWZ = 1.*mm;                       //$ Sensitive detector size Z>

	/*-------------------------------Detector Materials------------------------------*/
	//Aerogel
	const bool ADD_TILE = true;							//$ Add aerogel> 0.200*g/cm3
	const bool ADD_TILE_GRADIENT = false;							//$ Add aerogel gradient> 0.200*g/cm3
	const G4double AEROGEL_DENSITY = 531*kg/m3;         //$ Aerogel density> 0.200*g/cm3
	const std::string SURFACE_MODEL = "none";           //$ Aerogel surface roughness model> "unified" , "glisur" , "davis" 
	const G4double AEROGEL_ROUGHNESS = 0.1;               //$ Aerogel surface roughness parameter> 
	const std::string DAVIS_ROUGHNESS = "none";          //$ LUT Davis roughness parameter>  "rough" , "polished", "none"
	const G4double SILICA_DENSITY = 0.200*g/cm3;          //$ Silica density> 2.200*g/cm3 
	const G4double SILICA_PROP = 97*perCent;              //$ Aerogel silica proportion> 97
	const float REFRACTIVE_INDEX = 1.155;                 //$ Refractive index of Aerogel>
	const std::string AEROGEL_PROPERTIES = "ExactRindex"; //$ Refractive index distribution> "GaussianRindex", "ExactRindex"  
	const float REFRACTIVE_INDEX_SDEV = 0.001;                     //$ Refractive index standard deviation>
	const G4double PHOTON_MIN_ENERGY = 1.3*CLHEP::eV;        //$ Minimum threshold for photon generation> 1.3 = lower end of visible
	const G4double PHOTON_MAX_ENERGY = 7.3*CLHEP::eV;         //$ Maximum threshold for photon generation> 3.18 = upper end of visible
	//Environment
	const std::string WORLD_MATERIAL = "Galactic";        //$ Environment type> "Air", "Galactic" 
	const std::string ALTITUDE = "sea";                   //$ Detector altitude> "sea", "troposhere", "stratosphere", "mesosphere1", "mesosphere2", "thermosphere1", "thermosphere2"

	/*-------------------------------Sensitive Detector------------------------------*/
	const bool AIR = false;                               //$ Detector is cross section of air>
	const std::string PARTICLE_TO_DETECT = "optical";     //$ Particle to detect> "all", "primary" , 'optical'

	/*-------------------------------Particle Generator------------------------------*/
	const std::string SOURCE_PARTICLE = "e-";              //$ Primary particle> "e-" ,  "Be9"  , "Be10", "photon"
	const G4int N_PARTICLE = 1;                           //$ Primaries generated>
	//energy
	const G4double PARTICLE_ENERGY = 35*MeV;              //$ Primary particle energy>
	const std::string ENERGY_DISTRIBUTION = "exact";        //$ Primary energy distribution> "exact" , "uniform", "gaussian"
	const float ENERGY_SDEV = 1;                              //$ Standard deviation of Gaussian energy>
	const G4double UNIFORM_ENERGY_RADIUS = 50.0;              //$ Radius of uniform energy>
	//momentum
	const std::string MOMENTUM_DISTRIBUTION = "divergence";      //$ Momentum distribution> "exact" , "uniform", "gaussian", "divergence"
	const double MEAN_MOMENTUM_X = 0;                      //$ Mean x for Gaussian momentum>                           
	const double SDEV_MOMENTUM_X = 0.05;                  //$ Standard deviation of x for Gaussian momentum>
	const double MEAN_MOMENTUM_Y = 0;                      //$ Mean y for Gaussian momentum>
	const double SDEV_MOMENTUM_Y = 0.05;                  //$ Standard deviation of y for Gaussian momentum>
	const double MEAN_THETA = 0;                            //$ Mean theta for beam divergence> 0, 0.0872665, 0.174533, 0.261799
	const double SDEV_THETA = 0.0191986;                     //$ Standard deviation of theta for beam divergence>
	//location
	const std::string SOURCE_LOCATION = "gaussian_plane";                    //$ Particle source location> "exact", "gaussian_plane", "uniform_plane"
	const G4double UNIFORM_SOURCE_CENTRE = 0*cm;           //$ Uniform source location center>
	const G4double UNIFORM_SOURCE_RADIUS = 1*mm;           //$ Uniform source location radius>
	const double MEAN_SOURCE_X = 0;                        //$ Gaussian source mean location x>
	const double SDEV_SOURCE_X = 0.02;                     //$ Gaussian source standard deviation location x>
	const double MEAN_SOURCE_Y = 0;                        //$ Gaussian source mean location y>
	const double SDEV_SOURCE_Y = 0.02;                     //$ Gaussian source standard deviation location y>
	/*-------------------------------Physics Processes------------------------------*/
	//Electromagnetic processes
	const bool MULTIPLE_SCATTERING = true;                  //$ Multiple scattering>
	const bool IONIZATION_PROCESS = true;                   //$ Ionization>
	const bool BREMSSTRAHLUNG_PROCESS = true;				//$ Brehmsstrahlung Process>
	//Optical processes 
	const bool RAYLEIGH_SCATTERING = true;                 //$ Rayleigh scattering of photons> this causes an issue in which every tiny step gets printed to the terminal
	const bool ABSORPTION_PROCESS = true;                           //$ Absorption> no photons detected if false
	const bool BOUNDARY_PROCESS = true;							//$ Boundary Process>
	const bool SCINTILLATION_PROCESS = true;					//$ Scintillation Process>
	const bool COMPTON_SCATTERING = true; 						//$ Compton Scattering>
	const G4int MAX_NUM_PHOTONS = 1000;                         //$ Maximum number of photons produced per primary particle>
	const G4int MAX_BETA_CHANGE = 100;                       //$ Max beta change>

};