/* Physics List header */
#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "G4ParticleTable.hh"


class PhysicsList : public G4VModularPhysicsList
{
  public:
    PhysicsList();
    virtual ~PhysicsList();
    void ConstructParticle();
    void ConstructProcess();

    //these methods Construct particles
    // void ConstructBosons();
    void ConstructLeptons();
    void ConstructBosons();
    void ConstructEM();
    void ConstructIons();
    void ConstructOp();
    virtual void SetCuts();

   private:
    G4ParticleTable* theParticleTable;
  	G4ParticleTable::G4PTblDicIterator* theParticleIterator;


};
#endif
