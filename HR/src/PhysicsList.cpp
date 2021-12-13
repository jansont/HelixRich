#include "PhysicsList.hh"
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4EmSaturation.hh"
#include "G4LossTableManager.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"
#include "SimulationConstants.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4Scintillation.hh"
//#include "G4eplusAnnihilation"

PhysicsList::PhysicsList() : G4VUserPhysicsList() {
}
PhysicsList::~PhysicsList() {
}

void PhysicsList::ConstructParticle() {
	G4BosonConstructor bConstructor;
	bConstructor.ConstructParticle();

	G4LeptonConstructor lConstructor;
	lConstructor.ConstructParticle();
	G4IonConstructor iConstructor;
	iConstructor.ConstructParticle();
	ConstructIons();
	ConstructBosons();
	ConstructLeptons();
	ConstructBaryons();
}

void PhysicsList::ConstructIons()
{
	G4GenericIon::GenericIonDefinition();
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

void PhysicsList::ConstructBaryons()
{
	G4Proton::ProtonDefinition();
	G4AntiProton::AntiProtonDefinition();
}

void PhysicsList::ConstructProcess() {
  AddTransportation();
  ConstructEM();
  ConstructOp();
}

void PhysicsList::ConstructOp()
{

	// create Cerenkov
	cerenkovProcess = new G4Cerenkov("Cerenkov");
	cerenkovProcess->SetMaxNumPhotonsPerStep(50);
	cerenkovProcess->SetMaxBetaChangePerStep(1);
	cerenkovProcess->SetTrackSecondariesFirst(true);

	// create boundary process to be applied to opticalphotons
	boundaryProcess = new G4OpBoundaryProcess();

	// create scinitillation
	if (SimulationConstants::SCINTILLATION_PROCESS)
	{
		G4Scintillation *scintProcess = new G4Scintillation("Scintilation");
		scintProcess->SetScintillationYieldFactor(1.0);
		scintProcess->SetTrackSecondariesFirst(true);
	}

	G4OpAbsorption *absProcess = new G4OpAbsorption();
	G4OpRayleigh *rayProcess = new G4OpRayleigh();

	GetParticleIterator()->reset();
	while( (*GetParticleIterator())())
	{
		G4ParticleDefinition* particle = GetParticleIterator()->value();
		G4ProcessManager* pmanager = particle->GetProcessManager();
		G4String particleName = particle->GetParticleName();

		if (cerenkovProcess->IsApplicable(*particle) || particleName == "GenericIon") {
		  pmanager->AddProcess(cerenkovProcess);
		  pmanager->SetProcessOrdering(cerenkovProcess, idxPostStep);
		}
		if (particleName == "opticalphoton")
		{
				pmanager->AddDiscreteProcess(boundaryProcess);

			if (SimulationConstants::ABSORPTION_PROCESS)
				pmanager->AddDiscreteProcess(absProcess);

			if (SimulationConstants::RAYLEIGH_SCATTERING)
				pmanager->AddDiscreteProcess(rayProcess);

			if (SimulationConstants::COMPTON_SCATTERING)
				pmanager->AddProcess(new G4ComptonScattering());

			if (SimulationConstants::MULTIPLE_SCATTERING)
				pmanager->AddProcess(new G4hMultipleScattering());
		}
	}
}

void PhysicsList::ConstructEM()
{
	GetParticleIterator()->reset();
	while( (*GetParticleIterator())() ){
		G4ParticleDefinition* particle = GetParticleIterator()->value();
		G4ProcessManager* pmanager = particle->GetProcessManager();
		G4String particleName = particle->GetParticleName();

		if (particleName == "gamma")
		{
			if (SimulationConstants::COMPTON_SCATTERING)
				pmanager->AddDiscreteProcess(new G4ComptonScattering());
		}
		else if (particleName == "e-")
		{
			if (SimulationConstants::MULTIPLE_SCATTERING)
				pmanager->AddProcess(new G4eMultipleScattering(),-1, 1, 1);
			if (SimulationConstants::IONIZATION_PROCESS)
				pmanager->AddProcess(new G4eIonisation(),-1, 2, 2);
			if (SimulationConstants::BREMSSTRAHLUNG_PROCESS)
				pmanager->AddProcess(new G4eBremsstrahlung(),   -1, -1, 3);
	    } 
	    else if (particleName == "e+")
	    {
			if (SimulationConstants::MULTIPLE_SCATTERING)
				pmanager->AddProcess(new G4eMultipleScattering(),-1, 1, 1);
			if (SimulationConstants::IONIZATION_PROCESS)
				pmanager->AddProcess(new G4eIonisation(),-1, 2, 2);
			if (SimulationConstants::BREMSSTRAHLUNG_PROCESS)
				pmanager->AddProcess(new G4eBremsstrahlung(),-1, -1, 3);
	    }
	    else if (particleName == "proton" || particleName == "anti_proton")
	    {
			if (SimulationConstants::MULTIPLE_SCATTERING)	
				pmanager->AddProcess(new G4hMultipleScattering(),-1, 1, 1);
			if (SimulationConstants::IONIZATION_PROCESS)
				pmanager->AddProcess(new G4hIonisation(),       -1, 2, 2);
		}
		else if ((particleName == "alpha") && particle->GetPDGCharge() != 0.0)
		{
			if (SimulationConstants::IONIZATION_PROCESS)
				pmanager->AddProcess(new G4ionIonisation(), -1, 2, 2);
			if (SimulationConstants::MULTIPLE_SCATTERING)
				pmanager->AddProcess(new G4hMultipleScattering(),-1, 1, 1);
		}
		else if ((particleName == "GenericIon") && particle->GetPDGCharge() != 0.0)
		{
			pmanager->AddProcess(new G4ionIonisation(), -1, 2, 2);
		}
		else
		{
			if ((particle->GetPDGCharge() != 0.0) && (particle->GetParticleName() != "chargedgeantino")
				&& !particle->IsShortLived())
			{
				if (SimulationConstants::MULTIPLE_SCATTERING)
					pmanager->AddProcess(new G4hMultipleScattering(),-1,1,1);
				if (SimulationConstants::IONIZATION_PROCESS)	
					pmanager->AddProcess(new G4hIonisation(),       -1,2,2);
			}
		}
	}
}

void PhysicsList::SetCuts() {

  SetCutsWithDefault();
}
