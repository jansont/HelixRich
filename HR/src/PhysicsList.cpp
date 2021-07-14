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
// #include "G4MultipleScattering.hh"
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
#include "G4IonConstructor.hh"

PhysicsList::PhysicsList(): G4VModularPhysicsList()
{
  //setVerboseLevel(1);
  theParticleTable = G4ParticleTable::GetParticleTable();
  theParticleIterator = theParticleTable->GetIterator();
}

PhysicsList::~PhysicsList()
{}

/* Cuts are a threshold value at which a secondary will not be created. From the 
 * documentation they note that this is because some of the elctromagnetic processes
 * will create infrared divergence if this is not set. There is also a max threshold
 * energy set at 10 GeV and can be changed by running "/cuts/setMaxCutEnergy" in the
 * command line of the run.
 * We will just keep the defaults */
void PhysicsList::SetCuts()
{
  G4VUserPhysicsList::SetCuts();
  //   SetCutsWithDefault();   
}

void PhysicsList::ConstructProcess(){        
  AddTransportation();    //particle transportation
  ConstructEM();          //electromagnetic processes     
  ConstructOp();          //optical processes
}

void PhysicsList::ConstructParticle(){
  ConstructLeptons();
  ConstructBosons();
  ConstructIons();
}

void PhysicsList::ConstructIons(){
  {
    //  Construct light ions
    G4IonConstructor pConstructor;
    pConstructor.ConstructParticle();  
}

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
  theParticleIterator->reset();
  G4cout<<" Now creating EM processes"<<G4endl;
  while((*theParticleIterator)() )
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    if (SimulationConstants::multiple_scattering)
      pmanager->AddProcess(new G4eMultipleScattering(),-1,1,1);
    if (SimulationConstants::ionization)
      pmanager->AddProcess(new G4eIonisation(),-1,2,2);
    if (SimulationConstants::bremsstrahlung)
      pmanager->AddProcess(new G4eBremsstrahlung(),-1,-1,3);
    if (SimulationConstants::rayleigh)
      pmanager->AddProcess(new G4RayleighScattering());
    if (SimulationConstants::photoElectric)
  pmanager->AddProcess(new G4PhotoElectricEffect());
    if (SimulationConstants::compton)
  pmanager->AddProcess(new G4ComptonScattering());
  } 
}

//Optical processes (ie: Cerenkov radiation)
void PhysicsList::ConstructOp()
{
  G4Cerenkov* theCerenkovProcess = new G4Cerenkov("Cerenkov");
  G4OpAbsorption* theAbsorptionProcess = new G4OpAbsorption();
  G4OpRayleigh* theRayleighScatteringProcess = new G4OpRayleigh();
  G4OpMieHG* theMieHGScatteringProcess = new G4OpMieHG();
  G4OpBoundaryProcess* theBoundaryProcess = new G4OpBoundaryProcess();
  G4OpticalSurfaceModel themodel = unified;
  // theBoundaryProcess->SetModel(themodel);

  if (SimulationConstants::cerenkov){
    theCerenkovProcess->SetVerboseLevel(0);
    theCerenkovProcess->SetMaxBetaChangePerStep(SimulationConstants::maxBetaChange);
    theCerenkovProcess->SetTrackSecondariesFirst(true);
    theCerenkovProcess->SetMaxNumPhotonsPerStep(SimulationConstants::MaxNumPhotons);
  }
  if (SimulationConstants::absorption){
    theAbsorptionProcess->SetVerboseLevel(0);
  }
  if (SimulationConstants::rayleigh_scattering){
  theRayleighScatteringProcess->SetVerboseLevel(0);
  }

  // theBoundaryProcess->SetVerboseLevel(0);
  //G4OpticalSurfaceModel themodel = unified;
  // theBoundaryProcess->SetModel(themodel);

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    if (SimulationConstants::cerenkov)
      pmanager->AddDiscreteProcess(theCerenkovProcess); 
    if (particleName == "opticalphoton") {
      if (SimulationConstants::absorption)
        pmanager->AddDiscreteProcess(theAbsorptionProcess);
      if (SimulationConstants::rayleigh_scattering)
        pmanager->AddDiscreteProcess(theRayleighScatteringProcess);
      if (SimulationConstants::mie)
        pmanager->AddDiscreteProcess(theMieHGScatteringProcess);
    }
    if (SimulationConstants::cerenkov)
      pmanager->AddDiscreteProcess(theCerenkovProcess); 
      pmanager->AddDiscreteProcess(theBoundaryProcess);
   }
}







