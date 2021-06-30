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
#include "Randomize.hh"
#include "G4GeneralParticleSource.hh"
#include "SimulationConstants.hh"


/* We'll use the Geantino (non-interacting particle) for the gun, can be changed. 
 * This particle is generally used for testing. Think of it as similar to a neutrino. 
 *
 * In our gun design all particles will fire from a stationary spot
 * */
PrimaryGeneratorAction::PrimaryGeneratorAction()
{
	particleGun = new G4ParticleGun(SimulationConstants::n_particle);  // creation of particle gu
	G4ParticleDefinition* particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("e-");
	particleGun->SetParticleDefinition(particleDefinition);
	particleGun -> SetParticleEnergy(SimulationConstants::electronEnergy);  
}

// Create destructor
PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
   delete particleGun;
}

// Primary event
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	G4ThreeVector* v;
	if (SimulationConstants::distribution == 1)
		particleGun -> SetParticleMomentumDirection(G4ThreeVector(0,0,1));
	else if (SimulationConstants::distribution == 2){
		G4double Phi = CLHEP::twopi * G4UniformRand();
  		G4double Theta = 0.5 * CLHEP::pi * G4UniformRand();	
		particleGun -> SetParticleMomentumDirection(G4ThreeVector(cos(Theta), sin(Theta)*cos(Phi), sin(Theta)*sin(Phi)));
	}
	else if (SimulationConstants::distribution == 3){
		// // https://github.com/DavidSarria89/TGF-TEB-Propagation-Geant4/blob/master/src/src/PrimaryGeneratorAction.cc
		// R_max = 10000.;      // -> maximum angle is atan(10000) = 89.9943 degrees
  //       sigma_sample_R = std::tan(SimulationConstants::std_dev * degree);
  //       R_try = R_max + 10.; // just for initialization
  //       while (R_try > R_max) {
  //           X_try = CLHEP::RandGauss::shoot(SimulationConstants::mean, sigma_sample_R); // gaussian position sample
  //           Y_try = CLHEP::RandGauss::shoot(SimulationConstants::mean, sigma_sample_R); // gaussian position sample
  //           R_try = sqrt(X_try * X_try + Y_try * Y_try);
	}
  	// Pure z momentum
  	particleGun -> GeneratePrimaryVertex(anEvent);  // creates the initial momentum
}

