/* We have to have a .hh file for our primary generator action file */
#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh" //This is a Geant file, not created by user
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class G4GeneralParticleSource;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction();
    virtual ~PrimaryGeneratorAction();  // We must be able to destroy what we can create

  public:
    virtual void GeneratePrimaries(G4Event*); // See .cpp file
    const G4ParticleGun* GetParticleGun() const { return particleGun;}
    void SetOptPhotonPolar();
    void SetOptPhotonPolar(G4double);

  private:
    G4ParticleGun* particleGun;
    G4double timeConstant;

};

#endif
