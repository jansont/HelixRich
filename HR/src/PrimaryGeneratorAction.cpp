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

/* We'll use the Geantino (non-interacting particle) for the gun, can be changed. 
 * This particle is generally used for testing. Think of it as similar to a neutrino. 
 *
 * In our gun design all particles will fire from a stationary spot
 * */
PrimaryGeneratorAction::PrimaryGeneratorAction()
{
	G4int n_particle = 100;   // Number of particles fired per beamOn run
	particleGun = new G4ParticleGun(n_particle);  // creation of particle gun

	G4ParticleDefinition* particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("e-");
	particleGun->SetParticleDefinition(particleDefinition);
	// 1GeV energy of gun, this can be changed from the command line
	particleGun -> SetParticleEnergy(35.0 * MeV);  
	particleGun -> SetParticlePosition(G4ThreeVector(0.*m, 0.*m, -25.*cm));   // Set gun to be at furthest z in world (far left to standard orientation)
}

// Create destructor
PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
   delete particleGun;
}

// Primary event
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
   particleGun -> SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));  // Pure z momentum
 
   particleGun -> GeneratePrimaryVertex(anEvent);  // creates the initial momentum
}
