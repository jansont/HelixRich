#ifndef TrackerHit_h
#define TrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class TrackerHit : public G4VHit
{
  public:

      TrackerHit();
     ~TrackerHit();
      TrackerHit(const TrackerHit&);
      const TrackerHit& operator=(const TrackerHit&);
      G4int operator==(const TrackerHit&) const;

      inline void* operator new(size_t);
      inline void  operator delete(void*);

      void Print(std::ostream &stream = G4cout, bool printtime = true, bool printposition = false, bool printenergy = false);

  public:

      void SetTime(G4double t)      { time = t; }
      void SetEnergy   (G4double en)      { energy = en; }
      void SetPos      (G4ThreeVector xyz){ position = xyz; }

      G4double GetTime()  { return time; }
      G4double GetEnergy()    { return energy; }
      G4ThreeVector GetPos(){ return position; }
      
  private:
  
      G4double      time;
      G4double      energy;
      G4ThreeVector position;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<TrackerHit> TrackerHitsCollection;

extern G4Allocator<TrackerHit> TrackerHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* TrackerHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) TrackerHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void TrackerHit::operator delete(void *aHit)
{
  TrackerHitAllocator.FreeSingle((TrackerHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
