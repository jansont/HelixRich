/**
 * Remember to use the header file in the /include folder
 * This is the third REQUIRED file (Detector Const, and Phys list)
 */

#include "PrimaryGeneratorAction.hh"

// Now to create some events in our world
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ThreeVector.hh"
#include "G4Geantino.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"
#include "G4GeneralParticleSource.hh"
#include "SimulationConstants.hh"
#include "G4ChargedGeantino.hh"
#include <random>


/* We'll use the Geantino (non-interacting particle) for the gun, can be changed. 
 * This particle is generally used for testing. Think of it as similar to a neutrino. 
 *
 * In our gun design all particles will fire from a stationary spot
 * */
PrimaryGeneratorAction::PrimaryGeneratorAction()
{
	G4ParticleDefinition* particleDefinition;
	particleGun = new G4ParticleGun(SimulationConstants::n_particle);  // creation of particle gu

}


// Create destructor
PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
   delete particleGun;
}

// Primary event
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	G4ParticleDefinition* particle = particleGun->GetParticleDefinition();
	G4int z;
	G4double A, excitEnergy, ionCharge;
	//define particle type
	if (SimulationConstants::sourceParticle == "Be9"){
		z = 4; //atomic number
		A = 9.012182; //atomic mass
		excitEnergy = 0.*keV;
		ionCharge = 4*eplus;
		G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(z, A, excitEnergy);

		particleGun->SetParticleCharge(ionCharge);
		particleGun->SetParticleDefinition(ion);
	}
	else if (SimulationConstants::sourceParticle == "Be10"){
		z = 4; //atomic number
		A = 10.0135347; //atomic mass
		excitEnergy = 0.*keV;
		ionCharge = 4*eplus;
		G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(z, A, excitEnergy);
		particleGun->SetParticleCharge(ionCharge);
		particleGun->SetParticleDefinition(ion);
	}
	else{
		G4ParticleDefinition* particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("e-");
		particleGun->SetParticleDefinition(particleDefinition);
	}

	//define particle energy
	if (SimulationConstants::energyDistribution == "uniform"){
  	G4double energy_center = SimulationConstants::particleEnergy;
  	G4double energy_radius = SimulationConstants::uniformEnergyRadius;
  	energy_center += energy_radius*(G4UniformRand());
	}
	else if (SimulationConstants::energyDistribution == "gaussian"){
		std::random_device rd{};
		std::mt19937 gen{rd()};
		std::normal_distribution<double> E_rand{SimulationConstants::particleEnergy, SimulationConstants::E_sdev}; 
		G4double energy = E_rand(gen);
    	particleGun->SetParticleEnergy(energy);    
    }
    
	else if (SimulationConstants::energyDistribution == "exact"){
    	particleGun->SetParticleEnergy(SimulationConstants::particleEnergy);    
    }

	//define particle momentum
	if (SimulationConstants::momentumDistribution == "exact")
		particleGun -> SetParticleMomentumDirection(G4ThreeVector(0,0,1));
	else if (SimulationConstants::momentumDistribution == "uniform"){
		G4double Phi = CLHEP::twopi * G4UniformRand();
  		G4double Theta = 0.5 * CLHEP::pi * G4UniformRand();	
		particleGun -> SetParticleMomentumDirection(G4ThreeVector(cos(Theta), sin(Theta)*cos(Phi), sin(Theta)*sin(Phi)));
	}
	else if (SimulationConstants::momentumDistribution == "gaussian"){
		std::random_device rd{};
		std::mt19937 gen{rd()};
		std::normal_distribution<double> z_rand{SimulationConstants::mean_momentum_z, SimulationConstants::sdev_momentum_z}; 
		G4double z_dir = z_rand(gen);

		std::normal_distribution<double> x_rand{SimulationConstants::mean_momentum_x, SimulationConstants::sdev_momentum_x}; 
		G4double x_dir = x_rand(gen);

		std::normal_distribution<double> y_rand{SimulationConstants::mean_momentum_y, SimulationConstants::sdev_momentum_y}; 
		G4double y_dir = y_rand(gen);

		particleGun -> SetParticleMomentumDirection(G4ThreeVector(x_dir,y_dir,z_dir));
      }

     //define particle source location
      if (SimulationConstants::source == "uniform_plane"){
		G4double x0 = SimulationConstants::uniform_source_center;
		G4double y0 = SimulationConstants::uniform_source_center;
		// G4double z0 = SimulationConstants::uniform_source_center;	  
		G4double source_radius = SimulationConstants::uniform_source_radius;
		x0 += source_radius*(G4UniformRand()-0.5);
		y0 += source_radius*(G4UniformRand()-0.5);
		// z0 += source_radius*(G4UniformRand()-0.5);
		x0 *= cm;
		y0 *= cm;
		particleGun->SetParticlePosition(G4ThreeVector(x0,y0,0.*cm)); //constant z-position
      }
      else if (SimulationConstants::source == "gaussian_plane"){
		std::random_device rd{};
		std::mt19937 gen{rd()};
		std::normal_distribution<double> x_rand{SimulationConstants::mean_source_x, SimulationConstants::sdev_source_x}; 
		G4double X = x_rand(gen);
		X*=cm;
		std::normal_distribution<double> y_rand{SimulationConstants::mean_source_y, SimulationConstants::sdev_source_y}; 
		G4double Y = y_rand(gen);
		Y *= cm;
		particleGun -> SetParticleMomentumDirection(G4ThreeVector(X,Y,0.));
      }
      else{
      	particleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
      }
  	particleGun -> GeneratePrimaryVertex(anEvent);  // creates the initial momentum
}










