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

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4MaterialTable.hh"
#include "G4Tubs.hh"
#include "PrimaryGeneratorAction.hh"


// These two classes are called within the construct function
class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Box;

class DetectorConstruction : public G4VUserDetectorConstruction
{
   public:
      DetectorConstruction();
      virtual ~DetectorConstruction();
      virtual G4VPhysicalVolume* Construct();

      // G4LogicalVolume GetLogicalCurved()  { returnlogical_tile_curved;}
      // G4LogicalVolume GetLogicalBox()  { return logical_tile_box;}

   private:
      void DefineMaterials();

      G4Box* solid_world;
      G4LogicalVolume* logical_world;
      G4VPhysicalVolume* physical_world;

      G4Box* solid_tile_box;
      G4LogicalVolume* logical_tile_box;
      G4VPhysicalVolume* physical_tile_box;
      G4Tubs* solid_tile_curved;
      G4LogicalVolume* logical_tile_curved;
      G4VPhysicalVolume* physical_tile_curved;

      G4Box* solid_detector;
      G4LogicalVolume* logical_detector;
      G4VPhysicalVolume* physical_detector;
};
#endif
