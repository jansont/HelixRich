#include "TrackerHit.hh"
#include "G4UnitsTable.hh"

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
    stream << time <<'ns';
  if (printposition) {
    if (printtime)
      stream << "\t";
    stream  << position.x() <<'mm' << "\t" << position.y() <<'mm'  << "\t" << position.z() <<'mm';
  }
  if (printenergy) {
    if (printtime || printposition)
      stream << "\t";
    stream << energy<<'eV';
  }
  stream << std::endl;
}


