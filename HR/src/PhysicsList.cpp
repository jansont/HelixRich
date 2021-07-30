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

/* Second REQUIRED file (Detector Construction, Primary Generator Action) */
#include "PhysicsList.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4Positron.hh"
#include "G4eMultipleScattering.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eIonisation.hh"
#include "G4ProcessManager.hh"
#include "G4Cerenkov.hh"
#include "G4OpAbsorption.hh"
#include "G4ParticleDefinition.hh"
#include "G4OpRayleigh.hh"
#include "G4OpBoundaryProcess.hh"
#include "SimulationConstants.hh"
#include "G4OpMieHG.hh"
#include "G4RayleighScattering.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4ComptonScattering.hh"
#include "G4Material.hh"
#include "G4IonConstructor.hh"


#include "DetectorConstruction.hh"
#include "G4NistManager.hh"         // Contains the NIST data for particles
#include "G4Box.hh"                 // Needed to create box shapes
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"       // Will be required in every Detector construction file
#include "G4PVPlacement.hh"        
#include "G4RotationMatrix.hh"      
#include "G4Transform3D.hh"         
#include "G4PhysicalConstants.hh"   
#include "G4SystemOfUnits.hh"
#include "TrackerSD.hh"
#include "G4SDManager.hh"
#include "G4MaterialTable.hh"
#include "SimulationConstants.hh"
#include "G4OpticalSurface.hh"
#include "PhysicsList.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4Isotope.hh"
#include <random>
#include "G4ParticleGun.hh"
#include "G4Cerenkov.hh"


PhysicsList::PhysicsList(G4Material* agel): G4VModularPhysicsList()
{
	Agel = agel;
	theParticleTable = G4ParticleTable::GetParticleTable();
	theParticleIterator = theParticleTable->GetIterator();
}

PhysicsList::~PhysicsList()
{}

void PhysicsList::SetCuts()
{
	G4VUserPhysicsList::SetCuts();
}

void PhysicsList::ConstructProcess()
{        
	AddTransportation();    //particle transportation
	ConstructEM();          //electromagnetic processes     
	ConstructOp();          //optical processes
}

void PhysicsList::ConstructParticle()
{
	ConstructLeptons();
	ConstructBosons();
	ConstructIons();
}

void PhysicsList::ConstructIons()
{

	//Construct light ions
	G4IonConstructor pConstructor;
	pConstructor.ConstructParticle();  
}


void PhysicsList::ConstructBosons()
{
 	G4OpticalPhoton::OpticalPhotonDefinition();
	G4Gamma::GammaDefinition();
}

void PhysicsList::ConstructLeptons()
{
	G4Electron::ElectronDefinition();
	G4Positron::PositronDefinition();
}

//Electromagnetic processes
void PhysicsList::ConstructEM()
{
	auto particleIterator=GetParticleIterator();
	particleIterator->reset();
	G4cout<<" Now creating EM processes"<<G4endl;
	while((*particleIterator)())
	{
		G4ParticleDefinition* particle = particleIterator->value();
		G4ProcessManager* pmanager = particle->GetProcessManager();
		G4String particleName = particle->GetParticleName();
		if (particleName == "e-")
		{
			if (SimulationConstants::multiple_scattering)
				pmanager->AddProcess(new G4eMultipleScattering());
				// pmanager->AddProcess(new G4eMultipleScattering(),-1,1,1);
			if (SimulationConstants::ionization)
				pmanager->AddProcess(new G4eIonisation());
				// pmanager->AddProcess(new G4eIonisation(),-1,2,2);
			// if (SimulationConstants::bremsstrahlung)
			// 	pmanager->AddProcess(new G4eBremsstrahlung());
				// pmanager->AddProcess(new G4eBremsstrahlung(),-1,3,3);
		}
	} 
}



//Optical processes (ie: Cerenkov radiation)
G4Cerenkov* PhysicsList::ConstructOp()
{
	theCerenkovProcess = new G4Cerenkov("Cerenkov");
	G4OpAbsorption* theAbsorptionProcess = new G4OpAbsorption();
	G4OpBoundaryProcess* theBoundaryProcess = new G4OpBoundaryProcess();
	G4OpticalSurfaceModel themodel = glisur;
	// theBoundaryProcess->SetModel(themodel);

	if (SimulationConstants::cerenkov)
	{
		theCerenkovProcess->SetVerboseLevel(0);
		theCerenkovProcess->SetMaxBetaChangePerStep(SimulationConstants::maxBetaChange);
		theCerenkovProcess->SetTrackSecondariesFirst(true);
		theCerenkovProcess->SetMaxNumPhotonsPerStep(SimulationConstants::MaxNumPhotons);
	}
	if (SimulationConstants::absorption)
	{
		// theAbsorptionProcess->SetVerboseLevel(0);
	}
	auto particleIterator=GetParticleIterator();
	particleIterator->reset();
	while((*particleIterator)())
	{
		G4ParticleDefinition* particle = particleIterator->value();
		G4ProcessManager* pmanager = particle->GetProcessManager();
		G4String particleName = particle->GetParticleName();
		if (SimulationConstants::cerenkov)
		{
			if (particle->GetParticleName () == "e-")
			{
				G4cout << "$$$$check$$$$$$$" << particle->GetParticleName ()  << G4endl;
				pmanager->AddDiscreteProcess(theCerenkovProcess); 
			    pmanager->SetProcessOrdering(theCerenkovProcess,idxPostStep);
			}
		}
		if (particleName == "opticalphoton")
		{
			// if (SimulationConstants::absorption)
				// pmanager->AddDiscreteProcess(theAbsorptionProcess);
			// if (SimulationConstants::rayleigh_scattering)
				// pmanager->AddDiscreteProcess(theRayleighScatteringProcess);
				// G4OpRayleigh* theRayleighScatteringProcess = new G4OpRayleigh();
			// if (SimulationConstants::mie)
			// 	pmanager->AddDiscreteProcess(theMieHGScatteringProcess);
			// 	// G4OpMieHG* theMieHGScatteringProcess = new G4OpMieHG();
      		// pmanager->AddDiscreteProcess(theBoundaryProcess);
		}
	}

	/*--------------------------------------Define Elements-----------------------------------*/
	G4String name, symbol;
	G4double density, fractionmass, a, z;
	G4int ncomponents, natoms, n, iz;

	G4NistManager* nist = G4NistManager::Instance();
	G4Element *N = nist->FindOrBuildElement("N");
	G4Element *C = nist->FindOrBuildElement("C");
	G4Element *H = nist->FindOrBuildElement("H");
	G4Element *Si = nist->FindOrBuildElement("Si");
	G4Element *O = nist->FindOrBuildElement("O");
	G4Element *Ar = nist->FindOrBuildElement("Ar");
	G4Element *Sb = nist->FindOrBuildElement("Sb");
	G4Element *Rb = nist->FindOrBuildElement("Rb");
	G4Element *Cs = nist->FindOrBuildElement("Cs");
	G4Element *Be = nist->FindOrBuildElement("Be");  

	/*--------------------------------------Define Isotopes-----------------------------------*/
	//Isotope = name, atomic number, nucleon count, mass per mole 
	G4Isotope* Be9 = new G4Isotope(name = "Be9", iz = 4, n = 9, a = 9.012182*g/mole);
	G4Isotope* Be10 = new G4Isotope(name = "Be10", iz = 4, n = 10, a = 10.0135347*g/mole);

	/*--------------------------------------Define Materials-----------------------------------*/
	// secondary materials: 

	// water
	density = 1.00*g/cm3;
	G4Material* H2O = new G4Material(name="aqua", density, ncomponents=2);
	H2O -> AddElement(H, natoms=2);
	H2O -> AddElement(O, natoms=1);
	// silica
	G4Material* SiO2 = new G4Material(name="sili", SimulationConstants::silica_density, ncomponents=2);
	SiO2 -> AddElement(Si, natoms=1);
	SiO2 -> AddElement(O, natoms=2);


	G4Material* Aerogel = new G4Material(name="agel", SimulationConstants::aerogel_density, ncomponents=2); 
	Aerogel -> AddMaterial(SiO2, SimulationConstants::silica_prop);
	Aerogel -> AddMaterial(H2O , SimulationConstants::water_prop);

		// int numBins = 1000;
		// int i;
		// //photon energy vector for world material
		// G4double PhotonMinEnergy=SimulationConstants::PhotonMinEnergy_; 
		// G4double PhotonMaxEnergy=SimulationConstants::PhotonMaxEnergy_; 
		// G4double PhotonEnergyStep = (PhotonMaxEnergy-PhotonMinEnergy) / numBins;
		// G4double* PhotonMomentum = new G4double[numBins];
		// for (i=0; i<numBins; i++)
		// 	PhotonMomentum[i] = PhotonMinEnergy + PhotonEnergyStep * i;


		// G4double* AerogelRindex = new G4double[numBins];        
		// for (i=0; i<numBins; i++)
		// 	AerogelRindex[i]= SimulationConstants::refractive_index;
	

		// G4PhysicsOrderedFreeVector* rindex_vector  = new G4PhysicsOrderedFreeVector(PhotonMomentum, AerogelRindex, numBins);
			

		// const G4double photon_count = theCerenkovProcess->GetAverageNumberOfPhotons(-1*eplus,
  //                             									0.99, 
  //                             									Aerogel,
		//                        									rindex_vector);
		// G4cout << "Number of optical photons expected in this event : " << photon_count << G4endl;
}



