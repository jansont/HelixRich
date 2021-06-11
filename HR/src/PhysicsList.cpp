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
#include "G4OpRayleigh.hh"
#include "G4OpBoundaryProcess.hh"
#include "SimulationConstants.hh"


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
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    if (multiple_scattering)
      pmanager->AddProcess(new G4eMultipleScattering(),-1,1,1);
    if (ionization)
      pmanager->AddProcess(new G4eIonisation(),-1,2,2);
    if (bremsstrahlung)
      pmanager->AddProcess(new G4eBremsstrahlung(),-1,-1,3);
    }
}

//Optical processes (ie: Cerenkov radiation)
void PhysicsList::ConstructOp()
{
  G4Cerenkov* theCerenkovProcess = new G4Cerenkov("Cerenkov");
  theCerenkovProcess->SetVerboseLevel(0);
  theCerenkovProcess->SetMaxBetaChangePerStep(10.0);
  theCerenkovProcess->SetTrackSecondariesFirst(true);
  theCerenkovProcess->SetMaxNumPhotonsPerStep(MaxNumPhotons);

  // G4OpBoundaryProcess* theBoundaryProcess = new G4OpBoundaryProcess();
  G4OpAbsorption* theAbsorptionProcess = new G4OpAbsorption();
  G4OpRayleigh*   theRayleighScatteringProcess = new G4OpRayleigh();
  theAbsorptionProcess->SetVerboseLevel(0);
  theRayleighScatteringProcess->SetVerboseLevel(0);
  // theBoundaryProcess->SetVerboseLevel(0);
  //G4OpticalSurfaceModel themodel = unified;
  // theBoundaryProcess->SetModel(themodel);

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    pmanager->AddDiscreteProcess(theCerenkovProcess); 

    if (particleName == "opticalphoton") {
      if (cerenkov)
        pmanager->AddDiscreteProcess(theAbsorptionProcess);
      if (rayleigh_scattering)
        pmanager->AddDiscreteProcess(theRayleighScatteringProcess);
    //  pmanager->AddDiscreteProcess(theBoundaryProcess);
    }
  }
}







