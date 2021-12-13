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

/* Physics List header */
#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "G4Cerenkov.hh"
#include "G4ParticleTable.hh"
#include "G4OpBoundaryProcess.hh"



class PhysicsList : public G4VUserPhysicsList
{
  public:
    PhysicsList();
    virtual ~PhysicsList();
    void ConstructParticle();
    void ConstructProcess();

    //these methods Construct particles
    void ConstructLeptons();
    void ConstructBosons();
    void ConstructEM();
    void ConstructBaryons();
    void ConstructIons();
    void ConstructOp();
    virtual void SetCuts();

   private:
    G4ParticleTable* theParticleTable;
  	G4ParticleTable::G4PTblDicIterator* particleIterator;
    G4Cerenkov* cerenkovProcess;
    G4OpBoundaryProcess* boundaryProcess;


};
#endif
