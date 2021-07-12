#include "TrackerHit.hh"
#include "G4UnitsTable.hh"
#include <typeinfo>

G4Allocator<TrackerHit> TrackerHitAllocator;


TrackerHit::TrackerHit() {}

TrackerHit::~TrackerHit() {}

TrackerHit::TrackerHit(const TrackerHit& hit): G4VHit()
{
  time = hit.time;
  energy = hit.energy;
  position  = hit.position;
}

const TrackerHit& TrackerHit::operator=(const TrackerHit& hit)
{
  time = hit.time;
  energy = hit.energy;
  position = hit.position;
  return *this;
}


G4int TrackerHit::operator==(const TrackerHit& hit) const
{
  return (this==&hit) ? 1 : 0;
}

void TrackerHit::Print(std::ostream &stream, bool printtime, bool printposition, bool printenergy)
{
  if (printtime)
    stream << "\t";
    stream << time <<'ns'<< "\t";
    stream << std::ends;
  if (printposition) {
    if (printtime)
      // stream << "\t";
    stream  << position.x() <<'mm' << "\t" << position.y() <<'mm'  << "\t" << position.z() * 0.001 <<'mm'<< "\t" ;
    stream << std::ends;
  }
  if (printenergy) {
    // G4float e m= (G4float) energy;
    // stream << typeid(e).name() <<'eV';
    energy *= 1000000;
    std::string e = std::to_string (energy);
    stream << e << "\n" <<'eV';
    stream << std::flush;
  }
}


