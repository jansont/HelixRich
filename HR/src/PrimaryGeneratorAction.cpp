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


#include "PrimaryGeneratorAction.hh"
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

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
	G4ParticleDefinition* particleDefinition;
	particleGun = new G4ParticleGun(SimulationConstants::n_particle);  // creation of particle gun
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