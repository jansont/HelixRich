#include "SimulationConstants.hh"
#include <math.h>

namespace SimulationConstants {
	/*------------------------------Detector geometry------------------------------*/
	//Global 
	const G4int runCount = 1;                            //$ Run Count>
	//World volume
	const G4double worldX = 2.*m;                         //$ World size x>
	const G4double worldY = 2.*m;                         //$ World size y>
	const G4double worldZ = 2.*m;                         //$ World size z>
	const bool curved = false;                            //$ Curved tile geometry>
	//Tile Volume: box
	const G4double tileX = 15*cm;                         //$ Tile size X>
	const G4double tileY = 15*cm;                         //$ Tile size Y>
	const G4double tileZ = 1.0*cm;                        //$ Tile size Z>
	//Tile Volume: curved 
	const G4double pRMin = 19*cm;                         //$ Curved tile inner radius>
	const G4double pRMax = 20.*cm;                        //$ Curved tile outer radius>
	const G4double pDz = 15.*cm;                          //$ Curved tile half length in z>
	const G4double pSPhi = 0;                             //$ Curved tile starting phi angle (rad)>
	const G4double pDPhi = M_PI;                          //$ Curved tile angle of the segment (rad)>
	//Detector Volume
	const G4double windowX = 500*mm;                      //$ Sensitive detector size X>
	const G4double windowY = 500*mm;                      //$ Sensitive detector size Y>
	const G4double windowZ = 1.*mm;                       //$ Sensitive detector size Z>

	/*-------------------------------Detector Materials------------------------------*/
	//Aerogel
	const G4double aerogel_density = 0.200*g/cm3;         //$ Aerogel density> 0.200*g/cm3
	const std::string surfaceModel = "none";           //$ Aerogel surface roughness model> "unified" , "glisur" , "davis" - unified currently causing some issues in terminal
	const G4double aerogel_roughness = 0.3;               //$ Aerogel surface roughness parameter> 
	const std::string davis_roughness = "polished";          //$ LUT Davis roughness parameter>  "rough" , "polished", "none"
	const G4double silica_density = 2.200*g/cm3;          //$ Silica density> 2.200*g/cm3 
	const G4double silica_prop = 97*perCent;              //$ Aerogel silica proportion> 97
	const G4double water_prop =  3*perCent;               //$ Aerogel water proportion> 3
	const float refractive_index = 1.15;                 //$ Refractive index of Aerogel>
	const std::string Aerogel_properties = "ExactRindex"; //$ Refractive index distribution> "GaussianRindex", "ExactRindex"  
	const float rindex_sdev = 0.0007;                     //$ Refractive index standard deviation>
	const G4double PhotonMinEnergy_= 1.3*CLHEP::eV;        //$ Minimum threshold for photon generation> 1.3 = lower end of visible
	const G4double PhotonMaxEnergy_= 7.3*CLHEP::eV;         //$ Maximum threshold for photon generation> 3.18 = upper end of visible
	//Environment
	const std::string world_material_type = "Galactic";        //$ Environment type> "Air", "Galactic" (nothing detected for galactic??)
	const std::string altitude = "sea";                   //$ Detector altitude> "sea", "troposhere", "stratosphere", "mesosphere1", "mesosphere2", "thermosphere1", "thermosphere2"

	/*-------------------------------Sensitive Detector------------------------------*/
	const bool air = false;                               //$ Detector is cross section of air>
	const std::string particle_to_detect = "optical";     //$ Particle to detect> "all", "primary" , 'optical'

	/*-------------------------------Particle Generator------------------------------*/
	const std::string sourceParticle = "e-";              //$ Primary particle> "e-" ,  "Be9"  , "Be10"
	const G4int n_particle = 1;                           //$ Primaries generated>
	//energy
	const G4double particleEnergy = 35*MeV;              //$ Primary particle energy>
	const std::string energyDistribution = "exact";        //$ Primary energy distribution> "exact" , "uniform", "gaussian"
	const float E_sdev = 20.0;                              //$ Standard deviation of Gaussian energy>
	const G4double uniformEnergyRadius = 30.0;              //$ Radius of uniform energy>
	//momentum
	const std::string momentumDistribution = "exact";      //$ Momentum distribution> "exact" , "uniform", "gaussian", "divergence"
	const double mean_momentum_z = 1;                     //$ Mean z for Gaussian momentum> 
	const double sdev_momentum_z = 0;                   //$ Standard deviation of z for Gaussian momentum>
	const double mean_momentum_x = 0;                      //$ Mean x for Gaussian momentum>                           
	const double sdev_momentum_x = 0.05;                  //$ Standard deviation of x for Gaussian momentum>
	const double mean_momentum_y = 0;                      //$ Mean y for Gaussian momentum>
	const double sdev_momentum_y = 0.05;                  //$ Standard deviation of y for Gaussian momentum>
	const double mean_theta = 0;                            //$ Mean theta for beam divergence> 0, 0.0872665, 0.174533, 0.261799
	const double sdev_theta = 0.261799;                     //$ Standard deviation of theta for beam divergence>
	//location
	const std::string source = "exact";                    //$ Particle source location> "exact", "gaussian_plane", "uniform_plane"
	const G4double uniform_source_center = 0*cm;           //$ Uniform source location center>
	const G4double uniform_source_radius = 1*cm;           //$ Uniform source location radius>
	const double mean_source_x = 0;                        //$ Gaussian source mean location x>
	const double sdev_source_x = 0.05;                     //$ Gaussian source standard deviation location x>
	const double mean_source_y = 0;                        //$ Gaussian source mean location y>
	const double sdev_source_y = 0.05;                     //$ Gaussian source standard deviation location y>
	/*-------------------------------Physics Processes------------------------------*/
	//Electromagnetic processes
	const bool multiple_scattering = false;                  //$ Multiple scattering>
	const bool ionization = true;                           //$ Ionization>
	const bool bremsstrahlung = false;                       //$ Bremsstrahlung> no photon generated if false
	const bool rayleigh = false;                             //$ Rayleigh scattering>
	const bool photoElectric = false;                        //$ Photoelectric effect>
	const bool compton = false;                              //$ Compton scattering>
	//Optical processes 
	const bool cerenkov = true;                             //$ Cerenkov process> no photons generated if false
	const bool rayleigh_scattering = false;                 //$ Rayleigh scattering of photons> this causes an issue in which every tiny step gets printed to the terminal
	const bool absorption = true;                           //$ Absorption> no photons detected if false
	const bool mie = false;                                 //$ Mie scattering>
	const G4int MaxNumPhotons = 1000;                         //$ Maximum number of photons produced per primary particle>
	const G4int maxBetaChange = 10.0;                       //$ Max beta change>

};