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
    stream << time << "\t\t";
  if (printposition) {
    if (printtime)
    stream  << position.x() << "\t\t" << position.y() << "\t\t" << position.z();
  }
  if (printenergy) {
    // G4float e m= (G4float) energy;
    // stream << typeid(e).name() <<'eV';
    stream  << "\t";
    energy *= 1e6; //convert   to eV
    std::string e = std::to_string (energy);
    stream << e << "\n";
  }
  stream << std::endl;
}


