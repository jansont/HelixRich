#ifndef StackingAction_H
#define StackingAction_H 1

#include "globals.hh"
#include "G4UserStackingAction.hh"
#include <fstream>

class StackingAction : public G4UserStackingAction
{
  public:
    StackingAction();
   ~StackingAction();

  public:
    G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
    void NewStage();
    void PrepareNewEvent();
    G4int GetPhotonCounter();
  private:
    G4int photonCounter;
    std::ofstream outFile;
};
#endif
